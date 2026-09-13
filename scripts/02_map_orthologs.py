#!/usr/bin/env python3
"""Step 02 -- map each protein to a UniProt accession via BLAST against a
reference proteome, so downstream steps have something to look structures up by.

Skip entirely (orthologs.enabled: false) when your organism already has its own
populated UniProtKB entries. D23580 does not, which is why this step exists.

Author: Akhil <akhilvenkat197@gmail.com>
Part of the Gram-Negative Antigen Triage Pipeline. MIT licensed.
"""
from utils import *

def fetch_proteome(proteome_id, taxon_id, dest):
    if dest.exists() and dest.stat().st_size > 1000:
        log(f"using cached {dest.name}")
        return
    if proteome_id:
        q = f"proteome:{proteome_id}"
    else:
        q = f"taxonomy_id:{taxon_id}"
    log(f"downloading proteome ({q}) ...")
    subprocess.run(["curl", "-sS", "-G",
                    "--data-urlencode", f"query={q}",
                    "--data-urlencode", "format=fasta",
                    "--data-urlencode", "compressed=true",
                    "-o", str(dest) + ".gz",
                    "https://rest.uniprot.org/uniprotkb/stream"], check=True, timeout=900)
    subprocess.run(["gunzip", "-f", str(dest) + ".gz"], check=True)


def blast(query, db, out, evalue, max_targets=1):
    subprocess.run(["blastp", "-query", str(query), "-db", str(db), "-out", str(out),
                    "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
                    "-evalue", str(evalue), "-max_target_seqs", str(max_targets)],
                   check=True)


def main():
    cfg = load_config()
    oc = cfg["orthologs"]
    if not oc.get("enabled", True):
        log("orthologs.enabled is false -- skipping step 02")
        return
    require("blastp"); require("makeblastdb")

    query = RESULTS / "01_hypothetical.fasta"
    ref = RESULTS / "02_reference_proteome.fasta"
    fetch_proteome(oc["reference_proteome"], oc["species_taxon_id"], ref)

    db = RESULTS / "02_ref_db"
    subprocess.run(["makeblastdb", "-in", str(ref), "-dbtype", "prot", "-out", str(db)],
                   check=True, capture_output=True)

    def collect(path, into):
        for line in open(path):
            q, s, pid, ln, qlen, slen, ev, bs = line.rstrip("\n").split("\t")
            pid, ln, qlen = float(pid), int(ln), int(qlen)
            cov = 100.0 * ln / qlen
            if pid < oc["min_identity"] or cov < oc["min_coverage"]:
                continue
            acc = s.split("|")[1] if "|" in s else s
            if q not in into or float(bs) > into[q]["bitscore"]:
                into[q] = {"entry": q, "accession": acc, "identity": round(pid, 1),
                           "coverage": round(cov, 1), "evalue": ev, "bitscore": float(bs),
                           "source": path.stem}

    raw = RESULTS / "02_blast_hits.tsv"
    log("pass 1: blastp against the reference proteome ...")
    blast(query, db, raw, oc["evalue"], max_targets=1)
    best = {}
    collect(raw, best)
    log(f"  pass 1 mapped {len(best)}")

    # ---- pass 2: species-wide fallback for whatever pass 1 missed ----------
    # Without this the audit undercounts: unmapped proteins get filed as
    # NO_ORTHOLOG in step 03 rather than being properly classified.
    all_entries = set(read_fasta(query).keys())
    missing = sorted(all_entries - set(best))
    if missing and oc.get("species_taxon_id"):
        log(f"pass 2: {len(missing)} unmapped -- falling back to species-wide search ...")
        sp_fasta = RESULTS / "02_species_proteome.fasta"
        fetch_proteome(None, oc["species_taxon_id"], sp_fasta)
        sp_db = RESULTS / "02_species_db"
        subprocess.run(["makeblastdb", "-in", str(sp_fasta), "-dbtype", "prot",
                        "-out", str(sp_db)], check=True, capture_output=True)
        miss_fa = RESULTS / "02_unmapped.fasta"
        seqs = read_fasta(query)
        write_fasta(miss_fa, {m: seqs[m] for m in missing})
        raw2 = RESULTS / "02_blast_hits_species.tsv"
        blast(miss_fa, sp_db, raw2, oc["evalue"], max_targets=5)
        before = len(best)
        collect(raw2, best)
        log(f"  pass 2 recovered {len(best) - before}")

    src = {r["entry"]: r for r in read_tsv(RESULTS / "01_hypothetical.tsv")}
    rows = []
    for entry in src:
        b = best.get(entry)
        rows.append({
            "entry": entry,
            "product": src[entry]["product"],
            "length": src[entry]["length"],
            "accession": b["accession"] if b else "",
            "identity": b["identity"] if b else "",
            "coverage": b["coverage"] if b else "",
            "ortholog_source": b.get("source", "") if b else "",
            "tier": ("strong_ortholog" if b and b["identity"] >= 90
                     else "weak_ortholog" if b else "no_ortholog"),
        })
    write_tsv(RESULTS / "02_orthologs.tsv", rows)
    mapped = sum(1 for r in rows if r["accession"])
    log(f"mapped {mapped}/{len(rows)} proteins to a UniProt accession")
    log("NOTE: check the identity column. Substituting a distant homolog's "
        "structure for your protein is only defensible at high identity.")


if __name__ == "__main__":
    main()
