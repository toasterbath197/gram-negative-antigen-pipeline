#!/usr/bin/env python3
"""Step 01 -- pull CDS features annotated as 'hypothetical' from an ENA/EMBL flat file.

WARNING: the /product field reflects the annotation as deposited, which may be
many years old. Step 03 exists precisely because this label goes stale.

Author: Akhil <akhilvenkat197@gmail.com>
Part of the Gram-Negative Antigen Triage Pipeline. MIT licensed.
"""
from utils import *

def download_embl(accession, dest):
    if dest.exists() and dest.stat().st_size > 1000:
        log(f"using cached {dest.name}")
        return
    log(f"downloading {accession} from ENA ...")
    url = f"https://www.ebi.ac.uk/ena/browser/api/embl/{accession}?download=true"
    subprocess.run(["curl", "-sS", "-L", "-o", str(dest), url], check=True, timeout=900)


def parse_embl(path, keyword):
    """State machine over `FT   CDS` blocks. Handles multi-line quoted values."""
    records, cur, in_cds, qual, buf = [], None, False, None, []

    def flush():
        nonlocal qual, buf
        if qual and cur is not None:
            cur[qual] = "".join(buf).strip('"')
        qual, buf = None, []

    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith("FT   ") and not line.startswith("FT      "):
                key = line[5:21].strip()
                flush()
                if cur:
                    records.append(cur)
                cur = {"feature": key} if key else None
                in_cds = (key == "CDS")
                continue
            if not in_cds or cur is None or not line.startswith("FT      "):
                continue
            body = line[21:].rstrip("\n")
            if body.startswith("/"):
                flush()
                if "=" in body:
                    qual, first = body[1:].split("=", 1)
                    buf = [first]
                else:
                    cur[body[1:]] = True
            elif qual:
                buf.append(body)
    flush()
    if cur:
        records.append(cur)

    hits = []
    for r in records:
        if r.get("feature") != "CDS":
            continue
        product = str(r.get("product", ""))
        if keyword.lower() not in product.lower():
            continue
        if r.get("pseudo") or not r.get("translation"):
            continue
        hits.append({
            "entry": r.get("locus_tag", "").strip(),
            "gene": r.get("gene", "").strip(),
            "product": product,
            "length": len(r["translation"]),
            "sequence": r["translation"],
        })
    return hits


def main():
    cfg = load_config()
    acc = cfg["annotation"]["ena_accession"]
    kw = cfg["annotation"]["hypothetical_keyword"]

    embl = RESULTS / f"{acc}.embl"
    download_embl(acc, embl)

    log(f"parsing {embl.name} for /product containing '{kw}' ...")
    hits = parse_embl(embl, kw)
    hits = [h for h in hits if h["entry"]]
    log(f"found {len(hits)} candidate CDS")

    write_fasta(RESULTS / "01_hypothetical.fasta", {h["entry"]: h["sequence"] for h in hits})
    write_tsv(RESULTS / "01_hypothetical.tsv", hits,
              ["entry", "gene", "product", "length"])
    log(f"wrote 01_hypothetical.fasta and 01_hypothetical.tsv ({len(hits)} proteins)")


if __name__ == "__main__":
    main()
