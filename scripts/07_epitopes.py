#!/usr/bin/env python3
"""Step 07 -- epitope prediction, WITH NEGATIVE CONTROLS.

Read this before trusting any number this step produces.

1. The BioLib DiscoTope-3 wrapper overwrites every B-factor with 100.00 during
   chain extraction, in all input modes. DiscoTope uses pLDDT as a feature for
   AlphaFold models, so the run silently treats every residue as maximally
   confident. This script re-filters epitope residues against the real pLDDT
   parsed from the downloaded models. See docs/KNOWN_ISSUES.md.

2. In the worked example, BepiPred-2.0 and Kolaskar-Tongaonkar scored known
   cytoplasmic proteins as highly as surface candidates -- i.e. they carried no
   discriminating signal. The controls run by default. If your candidates do not
   separate from the controls, do not use these metrics.

Author: Akhil Venkatesh <akhilvenkat197@gmail.com>
Part of the Gram-Negative Antigen Triage Pipeline. MIT licensed.
"""
from utils import *

IEDB = "https://tools-cluster-interface.iedb.org/tools_api/bcell/"


def iedb_bcell(method, seq, extra=""):
    return curl(["--data", f"method={method}&sequence_text={seq}{extra}", IEDB], timeout=180)


def linear_scores(seq):
    out = {}
    txt = iedb_bcell("Bepipred-2.0", seq)
    rows = [l.split("\t") for l in txt.strip().split("\n")[1:] if l.strip()]
    try:
        sc = [float(r[2]) for r in rows]
        above = [i for i, v in enumerate(sc) if v >= 0.5]
        out["bepipred_residues"] = len(above)
        out["bepipred_epitopes_ge6"] = len(contiguous(above, 6))
        out["bepipred_max"] = round(max(sc), 3)
    except Exception:
        out["bepipred_error"] = txt[:80]
    time.sleep(2)
    txt = iedb_bcell("Kolaskar-Tongaonkar", seq, "&window_size=7")
    rows = [l.split("\t") for l in txt.strip().split("\n")[1:] if l.strip()]
    try:
        sc = [float(r[5]) for r in rows]
        out["kt_mean"] = round(sum(sc) / len(sc), 3)
        out["kt_max"] = round(max(sc), 3)
    except Exception:
        out["kt_error"] = txt[:80]
    time.sleep(2)
    return out


def run_discotope(zip_path, outdir):
    import biolib
    app = biolib.load("DTU/DiscoTope-3")
    # --struc_type is REQUIRED. --pdb_dir is not honoured by the BioLib build;
    # a zip of PDBs is.
    job = app.cli(args=f"--list_or_pdb_or_zip_file {zip_path} "
                       f"--struc_type alphafold --out_dir output", blocking=True)
    if job.get_status() != "completed":
        log("DiscoTope failed:"); log(job.get_stdout().decode()[-800:]); return None
    job.save_files(str(outdir))
    return outdir


def summarize_discotope(csv_dir, struct_dir, min_plddt):
    import glob, csv as _csv
    rows = []
    for f in sorted(glob.glob(str(csv_dir / "output" / "*.csv"))):
        base = os.path.basename(f).split("_A_discotope3")[0]
        entry = "_".join(base.split("_")[:2])
        src = struct_dir / f"{base}.pdb"
        real = parse_plddt(src) if src.exists() else {}
        d = list(_csv.DictReader(open(f)))
        ep = [r for r in d if r["epitope"].strip().lower() == "true"]
        conf = [r for r in ep if real.get(int(r["res_id"]), 0) >= min_plddt]
        patches = contiguous([int(r["res_id"]) for r in conf], 4)
        rows.append({"entry": entry, "length": len(d),
                     "discotope_epitope_residues_raw": len(ep),
                     "discotope_epitope_residues_plddt_filtered": len(conf),
                     "discotope_patches_ge4": len(patches),
                     "discotope_biggest_patch": max((len(p) for p in patches), default=0),
                     "mean_plddt": round(sum(real.values()) / len(real), 1) if real else ""})
    return rows


def main():
    cfg = load_config()
    ec = cfg["epitopes"]
    seqs = read_fasta(RESULTS / "04_surface_candidates.fasta")
    merged = {e: {"entry": e} for e in seqs}

    if ec.get("run_linear", True):
        log(f"IEDB linear predictions for {len(seqs)} candidates ...")
        for e, s in seqs.items():
            merged[e].update(linear_scores(s))
            log(f"  {e} done")

        log("running NEGATIVE CONTROLS (known cytoplasmic proteins) ...")
        controls = []
        for acc in ec.get("negative_controls", []):
            try:
                seq = json.loads(curl([f"https://rest.uniprot.org/uniprotkb/{acc}.json"]))["sequence"]["value"]
            except Exception:
                continue
            row = {"entry": f"CONTROL_{acc}", "length": len(seq)}
            row.update(linear_scores(seq))
            controls.append(row)
            log(f"  CONTROL {acc}: {row.get('bepipred_epitopes_ge6')} epitopes, KT {row.get('kt_mean')}")
        write_tsv(RESULTS / "07_negative_controls.tsv", controls)

        if controls:
            cand_kt = [v.get("kt_mean", 0) for v in merged.values() if v.get("kt_mean")]
            ctrl_kt = [c.get("kt_mean", 0) for c in controls if c.get("kt_mean")]
            if cand_kt and ctrl_kt:
                log("")
                log(f"  candidate KT mean : {sum(cand_kt)/len(cand_kt):.3f}")
                log(f"  control   KT mean : {sum(ctrl_kt)/len(ctrl_kt):.3f}")
                log("  If these are the same, the metric is not discriminating. "
                    "Leave scoring.include_epitope_terms = false.")
                log("")

    if ec.get("run_conformational", True):
        struct_dir = RESULTS / "05_structures"
        zip_path = RESULTS / "07_structures.zip"
        subprocess.run(["zip", "-qrj", str(zip_path), str(struct_dir)], check=True)
        outdir = RESULTS / "07_discotope"
        if not (outdir / "output").exists():
            run_discotope(zip_path, outdir)
        if (outdir / "output").exists():
            for row in summarize_discotope(outdir, struct_dir, ec["min_plddt"]):
                merged.setdefault(row["entry"], {"entry": row["entry"]}).update(row)

    write_tsv(RESULTS / "07_epitopes.tsv", list(merged.values()))
    log("wrote 07_epitopes.tsv and 07_negative_controls.tsv")


if __name__ == "__main__":
    main()
