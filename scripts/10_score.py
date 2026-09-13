#!/usr/bin/env python3
"""Step 10 -- assemble the final scorecard.

The weights here are a transparent, deliberately simple tally, not a validated
model. Treat the ranking as a way to order follow-up work, never as evidence
that a protein is an antigen.

Epitope terms are OFF by default. Turn them on only if step 07's negative
controls showed real separation for YOUR data.

Author: Akhil <akhilvenkat197@gmail.com>
Part of the Gram-Negative Antigen Triage Pipeline. MIT licensed.
"""
from collections import Counter
from utils import *


def low_complexity(seq, window, ):
    worst = 0.0
    for i in range(max(1, len(seq) - window + 1)):
        w = seq[i:i + window]
        if w:
            worst = max(worst, max(Counter(w).values()) / len(w))
    return worst


def main():
    cfg = load_config()
    sc = cfg["scoring"]
    base = {r["entry"]: r for r in read_tsv(RESULTS / "09_expression.tsv")}
    for extra in ("06_foldseek.tsv", "07_epitopes.tsv"):
        p = RESULTS / extra
        if p.exists():
            for r in read_tsv(p):
                if r.get("entry") in base:
                    base[r["entry"]].update({k: v for k, v in r.items() if k != "entry"})

    seqs = read_fasta(RESULTS / "04_surface_candidates.fasta")

    def f(row, key, default=0.0):
        try:
            return float(row.get(key) or default)
        except (TypeError, ValueError):
            return default

    rows = []
    for e, r in base.items():
        seq = seqs.get(e, "")
        lc = low_complexity(seq, sc["low_complexity_window"]) if seq else 0.0

        score, reasons = 0, []
        if r.get("foldseek_verdict") not in ("known_fold_strong", "known_fold_weak_evalue"):
            score += 1; reasons.append("no_known_fold")
        if f(r, "pct_conserved") >= 95:
            score += 1; reasons.append("conserved")
        if f(r, "mean_plddt") >= 70:
            score += 1; reasons.append("confident_structure")
        pct = r.get("infection_percentile_length_matched") or r.get("infection_percentile_all")
        if pct and float(pct) >= 50:
            score += 1; reasons.append("expressed_in_infection")
        if r.get("annotation_class") == "TRULY_UNCHARACTERIZED":
            score += 1; reasons.append("no_annotation_at_all")
        if sc.get("include_epitope_terms"):
            if f(r, "discotope_patches_ge4") >= 2:
                score += 1; reasons.append("conformational_epitopes")
            if f(r, "bepipred_epitopes_ge6") >= 3:
                score += 1; reasons.append("linear_epitopes")
        if lc >= sc["low_complexity_max_fraction"]:
            score -= 1; reasons.append("PENALTY_low_complexity")
        if r.get("sorting") == "inner_membrane_retained":
            score -= 1; reasons.append("PENALTY_periplasm_facing_lipoprotein")

        rows.append({**r, "low_complexity": round(lc, 2),
                     "score": score, "score_reasons": ";".join(reasons)})

    rows.sort(key=lambda r: (-r["score"], -float(r.get("pct_conserved") or 0)))
    cols = ["entry", "uniprot_name", "annotation_class", "length",
            "tmhmm_class", "sorting", "mean_plddt", "foldseek_verdict",
            "foldseek_prob", "foldseek_evalue", "pct_conserved",
            "infection_tpm", "infection_percentile_length_matched",
            "discotope_patches_ge4", "bepipred_epitopes_ge6", "kt_mean",
            "low_complexity", "score", "score_reasons"]
    write_tsv(RESULTS / "10_FINAL_SCORECARD.tsv", rows, cols)

    log("")
    log("=== FINAL RANKING ===")
    log(f"{'entry':<16}{'score':<7}{'cons%':<8}{'TPM':<10}{'fold':<24}{'name'}")
    for r in rows:
        log(f"{r['entry']:<16}{r['score']:<7}{str(r.get('pct_conserved','')):<8}"
            f"{str(r.get('infection_tpm','')):<10}{str(r.get('foldseek_verdict',''))[:23]:<24}"
            f"{r.get('uniprot_name','')}")
    log("")
    log("Reminder: every column here is a prediction. None of it establishes")
    log("that a protein is surface-exposed in vivo, translated, or immunogenic.")


if __name__ == "__main__":
    main()
