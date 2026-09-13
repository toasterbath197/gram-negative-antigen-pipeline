#!/usr/bin/env python3
"""Step 03 -- THE MOST IMPORTANT STEP. Re-check every protein against current
UniProt annotation.

A genome's /product field freezes at deposition. In the worked example, 96% of
proteins labelled 'hypothetical' in a 2009 reference had acquired a name or a
full functional assignment by the time of analysis. Filtering on the deposited
label alone means screening proteins that are no longer uncharacterised.

Author: Akhil <akhilvenkat197@gmail.com>
Part of the Gram-Negative Antigen Triage Pipeline. MIT licensed.
"""
import re
from utils import *

GENERIC = {
    "Cytoplasmic protein", "Inner membrane protein", "Periplasmic protein",
    "Outer membrane protein", "Outer membrane lipoprotein", "Membrane protein",
    "Periplasmic or exported protein", "Inner membrane lipoprotein",
    "Exported protein", "Secreted protein", "Putative inner membrane protein",
}
UNCHAR = re.compile(r"^(Putative )?[Uu]ncharacterized protein\s*[A-Za-z0-9]*$")


def classify(name, reviewed):
    if not name:
        return "NO_ANNOTATION_FOUND"
    if UNCHAR.match(name.strip()):
        return "TRULY_UNCHARACTERIZED"
    if reviewed == "reviewed":
        return "CHARACTERIZED_reviewed"
    if name.strip() in GENERIC:
        return "GENERIC_localization_only"
    return "NAMED_function_unreviewed"


def main():
    cfg = load_config()
    src = read_tsv(RESULTS / "02_orthologs.tsv")
    accs = [r["accession"] for r in src if r["accession"]]
    log(f"querying UniProt for {len(set(accs))} accessions ...")
    ann = uniprot_batch(accs)

    rows = []
    for r in src:
        a = ann.get(r["accession"], {})
        name = a.get("Protein names", "")
        rev = a.get("Reviewed", "")
        rows.append({**r,
                     "uniprot_name": name,
                     "uniprot_gene": a.get("Gene Names", ""),
                     "uniprot_reviewed": rev,
                     "annotation_class": classify(name, rev) if r["accession"] else "NO_ORTHOLOG"})
    write_tsv(RESULTS / "03_reannotated.tsv", rows)

    from collections import Counter
    counts = Counter(r["annotation_class"] for r in rows)
    total = len(rows)
    still = counts.get("TRULY_UNCHARACTERIZED", 0)
    log("")
    log("=== ANNOTATION AUDIT ===")
    for k, v in counts.most_common():
        log(f"  {v:5d}  {k}")
    log("")
    log(f"  Still genuinely uncharacterised: {still}/{total} ({100*still/total:.1f}%)")
    log(f"  Annotation now STALE for:        {total-still}/{total} ({100*(total-still)/total:.1f}%)")
    log("")

    # Working pool: no assigned function, and not cytoplasmic.
    pool = [r for r in rows
            if r["annotation_class"] in ("GENERIC_localization_only", "TRULY_UNCHARACTERIZED")
            and r["uniprot_name"] != "Cytoplasmic protein"]
    write_tsv(RESULTS / "03_candidate_pool.tsv", pool)
    seqs = read_fasta(RESULTS / "01_hypothetical.fasta")
    write_fasta(RESULTS / "03_candidate_pool.fasta",
                {r["entry"]: seqs[r["entry"]] for r in pool if r["entry"] in seqs})
    log(f"candidate pool (unassigned function, non-cytoplasmic): {len(pool)}")


if __name__ == "__main__":
    main()
