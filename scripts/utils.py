"""Shared helpers for the antigen triage pipeline.

Author: Akhil Venkatesh <akhilvenkat197@gmail.com>
Part of the Gram-Negative Antigen Triage Pipeline. MIT licensed.
"""
import csv, gzip, json, os, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Multi-organism support. Set these to run more than one organism side by side:
#   PIPELINE_CONFIG=config/config_acinetobacter.yaml \
#   PIPELINE_ORG=acinetobacter python3 03_reannotate.py
# Resolved at import time so `from utils import *` picks up the right paths.
CONFIG_PATH = Path(os.environ.get("PIPELINE_CONFIG", ROOT / "config" / "config.yaml"))
if not CONFIG_PATH.is_absolute():
    CONFIG_PATH = ROOT / CONFIG_PATH
_ORG = os.environ.get("PIPELINE_ORG", "").strip()
RESULTS = (ROOT / "results" / _ORG) if _ORG else (ROOT / "results")
RESULTS.mkdir(parents=True, exist_ok=True)


def load_config(path=None):
    import yaml
    with open(path or CONFIG_PATH) as f:
        return yaml.safe_load(f)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def done(name):
    """True if a step's output already exists (pipeline is resumable)."""
    p = RESULTS / name
    if p.exists() and p.stat().st_size > 0:
        log(f"SKIP {name} (already present)")
        return True
    return False


# ---------------------------------------------------------------- FASTA I/O

def read_fasta(path):
    seqs, hdr, cur = {}, None, []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if hdr:
                    seqs[hdr] = "".join(cur)
                hdr, cur = line[1:].split()[0], []
            else:
                cur.append(line)
    if hdr:
        seqs[hdr] = "".join(cur)
    return seqs


def write_fasta(path, seqs, width=60):
    with open(path, "w") as f:
        for name, s in seqs.items():
            f.write(f">{name}\n")
            for i in range(0, len(s), width):
                f.write(s[i:i + width] + "\n")


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fieldnames=None):
    if not rows:
        Path(path).write_text("")
        return
    fieldnames = fieldnames or list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------- HTTP

def curl(args, timeout=120, retries=3):
    """Thin curl wrapper with retries. Returns stdout as text."""
    for attempt in range(retries):
        try:
            r = subprocess.run(["curl", "-sS", "-m", str(timeout)] + args,
                               capture_output=True, text=True, timeout=timeout + 30)
            if r.returncode == 0 and r.stdout:
                return r.stdout
        except subprocess.TimeoutExpired:
            pass
        if attempt < retries - 1:
            time.sleep(2 ** attempt)
    return ""


def uniprot_batch(accessions, fields="accession,reviewed,protein_name,gene_names", chunk=40):
    """Fetch UniProt annotation for many accessions. Returns {accession: {...}}."""
    out = {}
    accessions = sorted(set(a for a in accessions if a))
    for i in range(0, len(accessions), chunk):
        block = accessions[i:i + chunk]
        query = " OR ".join(f"accession:{a}" for a in block)
        txt = curl(["-G",
                    "--data-urlencode", f"query={query}",
                    "--data-urlencode", "format=tsv",
                    "--data-urlencode", f"fields={fields}",
                    "--data-urlencode", "size=500",
                    "https://rest.uniprot.org/uniprotkb/search"])
        lines = txt.strip().split("\n")
        if len(lines) < 2:
            continue
        hdr = lines[0].split("\t")
        for ln in lines[1:]:
            parts = ln.split("\t")
            if parts and parts[0]:
                out[parts[0]] = dict(zip(hdr, parts))
        log(f"  UniProt {min(i + chunk, len(accessions))}/{len(accessions)}")
        time.sleep(0.5)
    return out


def alphafold_pdb_url(accession):
    """Resolve the CURRENT AlphaFold model URL. Do not build the filename by
    hand -- the model version changes (v4 -> v6 during development)."""
    txt = curl([f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"])
    try:
        meta = json.loads(txt)
        return meta[0]["pdbUrl"] if meta else None
    except Exception:
        return None


# ---------------------------------------------------------------- structures

def parse_plddt(pdb_path):
    """Per-residue pLDDT from the B-factor column of CA atoms."""
    out = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    out[int(line[22:26])] = float(line[60:66])
                except ValueError:
                    pass
    return out


def require(binary):
    from shutil import which
    if which(binary) is None:
        sys.exit(f"ERROR: '{binary}' not found on PATH. See README prerequisites.")


def contiguous(indices, min_len=1):
    """Group sorted integers into consecutive runs of at least min_len."""
    runs, cur = [], []
    for i in sorted(indices):
        if cur and i == cur[-1] + 1:
            cur.append(i)
        else:
            if cur:
                runs.append(cur)
            cur = [i]
    if cur:
        runs.append(cur)
    return [r for r in runs if len(r) >= min_len]
