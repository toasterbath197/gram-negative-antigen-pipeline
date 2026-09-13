#!/usr/bin/env python3
"""Step 04 -- predicted subcellular localization.

DeepTMHMM runs on BioLib's cloud (no local install, no API key needed).
The lipoprotein +2 sorting rule is Enterobacteriaceae-specific; disable it in
config for anything outside that family.

Author: Akhil Venkatesh <akhilvenkat197@gmail.com>
Part of the Gram-Negative Antigen Triage Pipeline. MIT licensed.
"""
import re
from utils import *


def run_deeptmhmm(fasta, outdir):
    import biolib
    app = biolib.load("DTU/DeepTMHMM")
    # NOTE: machine='local' is no longer supported by BioLib.
    job = app.cli(args=f"--fasta {fasta}")
    job.wait()
    job.save_files(str(outdir))
    return outdir


def parse_3line(path):
    txt = open(path).read().strip().split("\n")
    recs = {}
    for i in range(0, len(txt), 3):
        name, cls = [x.strip() for x in txt[i][1:].split("|")]
        seq, topo = txt[i + 1], txt[i + 2]
        m = re.match(r"^S+", topo)
        recs[name] = {
            "class": cls, "seq": seq, "topology": topo,
            "signal_len": len(m.group(0)) if m else 0,
            "n_tm": len(re.findall(r"M+", topo)),
            "exposed_frac": (topo.count("O") + topo.count("P")) / len(topo),
        }
    return recs


def plus2_rule(seq, signal_len):
    """Asp at +2 after the lipidated Cys => retained in the inner membrane."""
    cys = seq.find("C", max(0, signal_len - 3))
    if cys == -1:
        return {"is_lipoprotein": False, "plus2": "", "sorting": "not_a_lipoprotein"}
    p2 = seq[cys + 2] if cys + 2 < len(seq) else ""
    return {"is_lipoprotein": True, "plus2": p2,
            "sorting": "inner_membrane_retained" if p2 == "D" else "outer_membrane_sorted"}


def main():
    cfg = load_config()
    lc = cfg["localization"]
    fasta = RESULTS / "03_candidate_pool.fasta"
    outdir = RESULTS / "04_deeptmhmm"

    if not (outdir / "predicted_topologies.3line").exists():
        log("submitting to DeepTMHMM (BioLib cloud) ...")
        run_deeptmhmm(fasta, outdir)
    recs = parse_3line(outdir / "predicted_topologies.3line")

    rows = []
    for r in read_tsv(RESULTS / "03_candidate_pool.tsv"):
        d = recs.get(r["entry"])
        if not d:
            continue
        surface = ("SP" in d["class"] and d["n_tm"] <= lc["max_tm_helices"]) or \
                  (d["n_tm"] <= lc["max_tm_helices"] and d["exposed_frac"] >= 0.5)
        row = {**r, "tmhmm_class": d["class"], "signal_len": d["signal_len"],
               "n_tm_helices": d["n_tm"], "exposed_fraction": round(d["exposed_frac"], 2),
               "surface_candidate": "yes" if surface else "no"}
        if lc.get("apply_plus2_rule") and "SP" in d["class"]:
            row.update(plus2_rule(d["seq"], d["signal_len"]))
        rows.append(row)

    write_tsv(RESULTS / "04_localization.tsv", rows)
    keep = [r for r in rows if r["surface_candidate"] == "yes"]
    seqs = read_fasta(fasta)
    write_fasta(RESULTS / "04_surface_candidates.fasta",
                {r["entry"]: seqs[r["entry"]] for r in keep if r["entry"] in seqs})
    write_tsv(RESULTS / "04_surface_candidates.tsv", keep)
    log(f"surface candidates: {len(keep)}/{len(rows)}")
    log("NOTE: automatic 'membrane protein' names often reflect a hydrophobic "
        "signal peptide misread as a TM helix. Trust the topology, not the name.")


if __name__ == "__main__":
    main()
