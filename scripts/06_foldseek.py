#!/usr/bin/env python3
"""Step 06 -- structural similarity search against PDB100 via the Foldseek
web API (undocumented but stable endpoints).

A hit is only called 'known fold' when probability AND E-value agree. A
probability of 1.0 alongside an E-value of 0.02 is not a confident match.

Author: Akhil Venkatesh <akhilvenkat197@gmail.com>
Part of the Gram-Negative Antigen Triage Pipeline. MIT licensed.
"""
from utils import *


def submit(pdb_path, mode, databases):
    args = ["-X", "POST", "https://search.foldseek.com/api/ticket",
            "-F", f"q=@{pdb_path}", "-F", f"mode={mode}"]
    for db in databases:
        args += ["-F", f"database[]={db}"]
    try:
        return json.loads(curl(args)).get("id")
    except Exception:
        return None


def wait(ticket, tries=24, delay=5):
    for _ in range(tries):
        time.sleep(delay)
        try:
            st = json.loads(curl([f"https://search.foldseek.com/api/ticket/{ticket}"])).get("status")
        except Exception:
            st = None
        if st in ("COMPLETE", "ERROR"):
            return st
    return "TIMEOUT"


def best_hit(ticket):
    try:
        res = json.loads(curl([f"https://search.foldseek.com/api/result/{ticket}/0"], timeout=180))
    except Exception:
        return None
    best = None
    for block in res.get("results", []):
        for group in (block.get("alignments") or []):
            for hit in group:
                if best is None or hit.get("prob", -1) > best.get("prob", -1):
                    best = hit
    return best


def verdict(hit, fc):
    if hit is None:
        return "no_structural_hit"
    p, e = float(hit.get("prob", 0)), float(hit.get("eval", 1))
    if p >= fc["strong_prob"] and e <= fc["strong_evalue"]:
        return "known_fold_strong"
    if p >= fc["strong_prob"] and e <= fc["weak_evalue"]:
        return "known_fold_weak_evalue"
    if p >= fc["twilight_prob"]:
        return "twilight_inconclusive"
    return "no_significant_hit"


def main():
    cfg = load_config()
    fc = cfg["foldseek"]
    rows = read_tsv(RESULTS / "05_structures.tsv")
    out_path = RESULTS / "06_foldseek.tsv"

    prior = {r["entry"]: r for r in read_tsv(out_path)} if out_path.exists() else {}
    results = []
    for r in rows:
        if r["entry"] in prior and prior[r["entry"]].get("foldseek_verdict"):
            results.append(prior[r["entry"]]); continue
        if not r.get("structure_file"):
            results.append({**r, "foldseek_verdict": "no_structure"}); continue

        pdb = RESULTS / "05_structures" / r["structure_file"]
        log(f"searching {r['entry']} ...")
        ticket = submit(pdb, fc["mode"], fc["databases"])
        if not ticket:
            results.append({**r, "foldseek_verdict": "submit_failed"}); continue
        if wait(ticket) != "COMPLETE":
            results.append({**r, "foldseek_verdict": "job_failed"}); continue
        hit = best_hit(ticket)
        results.append({**r,
            "foldseek_top_hit": (hit or {}).get("target", "")[:100],
            "foldseek_prob": (hit or {}).get("prob", ""),
            "foldseek_evalue": (hit or {}).get("eval", ""),
            "foldseek_seqid": (hit or {}).get("seqId", ""),
            "foldseek_verdict": verdict(hit, fc)})
        write_tsv(out_path, results)   # checkpoint after every query
        time.sleep(2)

    write_tsv(out_path, results)
    from collections import Counter
    for k, v in Counter(r.get("foldseek_verdict") for r in results).most_common():
        log(f"  {v:4d}  {k}")


if __name__ == "__main__":
    main()
