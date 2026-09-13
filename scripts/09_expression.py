#!/usr/bin/env python3
"""Step 09 -- expression in infection-relevant conditions. OPTIONAL.

Skipped when expression.geo_url is null. Most organisms have no comparable
dataset; that absence is a real limit on the pipeline, not a bug.

TPM normalises by feature length, so short genes inflate. Percentiles are
therefore reported both genome-wide and against a length-matched background.
"""
import bisect, statistics
from utils import *


def main():
    cfg = load_config()
    xc = cfg["expression"]
    if not xc.get("geo_url"):
        log("expression.geo_url is null -- skipping step 09")
        rows = read_tsv(RESULTS / "08_conservation.tsv")
        write_tsv(RESULTS / "09_expression.tsv", rows)
        return

    raw = RESULTS / "09_expression_raw.tsv"
    if not raw.exists():
        gz = str(raw) + ".gz"
        log("downloading expression table ...")
        subprocess.run(["curl", "-sS", "-o", gz, xc["geo_url"]], check=True, timeout=900)
        subprocess.run(["gunzip", "-f", gz], check=True)

    import csv as _csv
    table = list(_csv.reader(open(raw, encoding=xc.get("encoding", "utf-8")), delimiter="\t"))
    hdr = table[0]
    first = xc["first_data_column"]
    conds = [h.strip() for h in hdr[first:]]
    if xc["infection_condition"] not in conds:
        log(f"WARNING: '{xc['infection_condition']}' not in {conds}")
        return
    ci = conds.index(xc["infection_condition"])

    data, lengths = {}, {}
    for r in table[1:]:
        try:
            data[r[xc["locus_tag_column"]]] = [float(x) for x in r[first:]]
            lengths[r[xc["locus_tag_column"]]] = abs(int(r[2]) - int(r[1])) + 1
        except (ValueError, IndexError):
            continue

    allv = sorted(v[ci] for v in data.values())
    med = statistics.median(allv)
    log(f"genome-wide {xc['infection_condition']} median = {med:.1f} (n={len(allv)})")

    rows = []
    for r in read_tsv(RESULTS / "08_conservation.tsv"):
        v = data.get(r["entry"])
        if v is None:
            rows.append({**r, "expression_note": "not_in_table"}); continue
        val = v[ci]
        pct_all = 100 * bisect.bisect_left(allv, val) / len(allv)
        out = {**r, "infection_tpm": round(val, 2),
               "infection_percentile_all": round(pct_all, 1),
               "max_tpm": round(max(v), 2), "min_tpm": round(min(v), 2),
               "above_genome_median": "yes" if val >= med else "no"}
        if xc.get("length_matched_percentile"):
            L = lengths.get(r["entry"], 0)
            lo, hi = L * 0.5, L * 2.0
            matched = sorted(data[g][ci] for g in data if lo <= lengths.get(g, 0) <= hi)
            if len(matched) >= 30:
                out["infection_percentile_length_matched"] = round(
                    100 * bisect.bisect_left(matched, val) / len(matched), 1)
        rows.append(out)

    write_tsv(RESULTS / "09_expression.tsv", rows)
    for r in sorted(rows, key=lambda x: -float(x.get("infection_tpm") or 0)):
        log(f"  {r['entry']:<16}TPM={r.get('infection_tpm','-'):<10}"
            f"pct(all)={r.get('infection_percentile_all','-')}  "
            f"pct(len-matched)={r.get('infection_percentile_length_matched','-')}")


if __name__ == "__main__":
    main()
