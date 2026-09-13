#!/usr/bin/env python3
"""Step 08 -- how widely is each candidate conserved?

Downloads complete RefSeq proteomes for the target taxon and BLASTs the
candidates against all of them. An antigen absent from part of the population
is not useful, so this step is a hard filter, not a nice-to-have.

Author: Akhil Venkatesh <akhilvenkat197@gmail.com>
Part of the Gram-Negative Antigen Triage Pipeline. MIT licensed.
"""
from utils import *


def download_genomes(taxon, level, workdir):
    datasets = workdir / "datasets"
    if not datasets.exists():
        log("fetching the NCBI datasets CLI ...")
        subprocess.run(["curl", "-sS", "-o", str(datasets),
                        "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets"],
                       check=True, timeout=600)
        datasets.chmod(0o755)
    zip_path = workdir / "genomes.zip"
    if not zip_path.exists():
        log(f"downloading complete RefSeq proteomes for taxon {taxon} (this can take a while) ...")
        subprocess.run([str(datasets), "download", "genome", "taxon", str(taxon),
                        "--assembly-source", "RefSeq", "--assembly-level", level,
                        "--include", "protein", "--filename", str(zip_path)],
                       check=True, timeout=3600)
    extracted = workdir / "genomes"
    if not extracted.exists():
        extracted.mkdir()
        subprocess.run(["unzip", "-q", str(zip_path), "-d", str(extracted)], check=True)
    return extracted


def main():
    cfg = load_config()
    cc = cfg["conservation"]
    require("blastp"); require("makeblastdb")
    workdir = RESULTS / "08_conservation"
    workdir.mkdir(exist_ok=True)

    extracted = download_genomes(cc["taxon_id"], cc["assembly_level"], workdir)

    combined = workdir / "all_proteins.faa"
    if not combined.exists():
        import glob
        n_g = 0
        with open(combined, "w") as out:
            for f in glob.glob(str(extracted / "ncbi_dataset" / "data" / "*" / "protein.faa")):
                acc = os.path.basename(os.path.dirname(f))
                n_g += 1
                for line in open(f):
                    out.write(f">{acc}|{line[1:]}" if line.startswith(">") else line)
        log(f"combined {n_g} proteomes")

    n_genomes = len(set(l.split("|")[0][1:] for l in open(combined) if l.startswith(">")))
    log(f"reference set: {n_genomes} genomes")

    db = workdir / "db"
    if not (workdir / "db.pin").exists():
        subprocess.run(["makeblastdb", "-in", str(combined), "-dbtype", "prot",
                        "-out", str(db)], check=True, capture_output=True)

    query = RESULTS / "04_surface_candidates.fasta"
    hits_path = workdir / "hits.tsv"
    if not hits_path.exists():
        log("running blastp (this is the slow part) ...")
        subprocess.run(["blastp", "-query", str(query), "-db", str(db), "-out", str(hits_path),
                        "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
                        "-evalue", "1e-5", "-max_target_seqs", "5000"], check=True)

    from collections import defaultdict
    best = defaultdict(dict)
    for line in open(hits_path):
        q, s, pid, ln, qlen, slen, ev, bs = line.rstrip("\n").split("\t")
        if 100.0 * int(ln) / int(qlen) >= cc["min_coverage"]:
            g = s.split("|")[0]
            best[q][g] = max(best[q].get(g, 0), float(pid))

    rows = []
    for r in read_tsv(RESULTS / "04_surface_candidates.tsv"):
        g = best.get(r["entry"], {})
        near = sum(1 for v in g.values() if v >= cc["conserved_identity"])
        rows.append({**r, "genomes_searched": n_genomes,
                     "genomes_present": len(g),
                     "pct_present": round(100 * len(g) / n_genomes, 1),
                     "genomes_conserved": near,
                     "pct_conserved": round(100 * near / n_genomes, 1)})
    rows.sort(key=lambda r: -r["pct_conserved"])
    write_tsv(RESULTS / "08_conservation.tsv", rows)
    for r in rows:
        log(f"  {r['entry']:<16}{r['pct_conserved']:>6.1f}% conserved")


if __name__ == "__main__":
    main()
