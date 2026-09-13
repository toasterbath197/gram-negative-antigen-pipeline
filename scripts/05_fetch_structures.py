#!/usr/bin/env python3
"""Step 05 -- download AlphaFold models for the surface candidates.

Always resolve the model URL from the API. The filename version is not stable
(it moved v4 -> v6 during development, and hand-built v4 URLs now 404).

Author: Akhil <akhilvenkat197@gmail.com>
Part of the Gram-Negative Antigen Triage Pipeline. MIT licensed.
"""
from utils import *


def main():
    cfg = load_config()
    rows = read_tsv(RESULTS / "04_surface_candidates.tsv")
    outdir = RESULTS / "05_structures"
    outdir.mkdir(exist_ok=True)

    manifest = []
    for r in rows:
        acc = r.get("accession", "")
        if not acc:
            manifest.append({**r, "structure_file": "", "note": "no_accession"})
            continue
        dest = outdir / f"{r['entry']}_{acc}.pdb"
        if not dest.exists() or dest.stat().st_size < 500:
            url = alphafold_pdb_url(acc)
            if not url:
                manifest.append({**r, "structure_file": "", "note": "no_alphafold_entry"})
                continue
            subprocess.run(["curl", "-sS", "-o", str(dest), url], timeout=120)
            time.sleep(0.5)
        if dest.exists() and dest.stat().st_size > 500:
            pl = parse_plddt(dest)
            mean_plddt = round(sum(pl.values()) / len(pl), 1) if pl else ""
            manifest.append({**r, "structure_file": dest.name,
                             "mean_plddt": mean_plddt, "note": ""})
            log(f"  {r['entry']} ({acc}) mean pLDDT={mean_plddt}")
        else:
            manifest.append({**r, "structure_file": "", "note": "download_failed"})

    write_tsv(RESULTS / "05_structures.tsv", manifest)
    ok = sum(1 for m in manifest if m.get("structure_file"))
    log(f"structures retrieved: {ok}/{len(manifest)}")


if __name__ == "__main__":
    main()
