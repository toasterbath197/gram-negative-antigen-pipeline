# Known issues and gotchas

Recorded from an actual run. Each of these cost real debugging time.

## 1. DiscoTope-3 via BioLib destroys pLDDT

**Severity: high — silently wrong results.**

The BioLib `DTU/DiscoTope-3` wrapper overwrites every B-factor with `100.00`
during chain extraction. This happens in **all** input modes (single PDB, zip,
`--pdb_dir`). DiscoTope-3 uses pLDDT as a feature for AlphaFold models, so the
run treats every residue as maximally confident and inflates epitope scores in
disordered regions.

Verify it yourself:

```bash
# original AlphaFold model
grep "^ATOM" model.pdb | awk '{print substr($0,61,6)}' | sort -u | wc -l   # ~130
# what DiscoTope actually scored
grep "^ATOM" discotope_out/input_chains/model_A.pdb | awk '{print substr($0,61,6)}' | sort -u | wc -l   # 1
```

**Workaround (implemented in step 07):** parse pLDDT from the original models
and discard epitope residues below `epitopes.min_plddt`. For publication-grade
work, prefer the DTU web server over BioLib.

## 2. DiscoTope requires `--struc_type`, and `--pdb_dir` doesn't work

`--struc_type alphafold|solved` is mandatory; omitting it fails with a usage
error. `--pdb_dir` is advertised but the BioLib build ignores it — pass a zip
via `--list_or_pdb_or_zip_file`.

## 3. AlphaFold model version is not stable

`AF-{ACC}-F1-model_v4.pdb` returned HTTP 404 during development; the current
version is v6. Always resolve `pdbUrl` from
`https://alphafold.ebi.ac.uk/api/prediction/{accession}` instead of building the
filename.

## 4. BioLib rejects `machine="local"`

`app.cli(args=..., machine='local')` raises `BioLibError`. Cloud execution works
without an API key.

## 5. Foldseek: probability saturates

The web API's `prob` reaches 1.0 for hits whose E-value is only 0.02. Requiring
probability alone produces false confidence. Step 06 requires
`prob >= 0.95 AND evalue <= 1e-3` for a strong call.

Endpoints (undocumented but stable):
```
POST /api/ticket          -F q=@model.pdb -F mode=3diaa -F "database[]=pdb100"
GET  /api/ticket/{id}     poll until status COMPLETE
GET  /api/result/{id}/0   JSON
```

## 6. Epitope predictors may not discriminate at all

On the D23580 data, BepiPred-2.0 and Kolaskar-Tongaonkar scored known
cytoplasmic proteins as highly as surface candidates. Kolaskar-Tongaonkar
returned ~1.00–1.03 for everything, and BepiPred counts tracked sequence length
rather than antigenicity. Always run the controls in step 07 before believing
these numbers.

## 7. TPM inflates short genes

TPM normalises by feature length. In `GSE119724`, median macrophage TPM was 59.9
for features ≤250 nt versus 19.9 for those >1000 nt, and the genome-wide top
hits were almost all ncRNAs. A 165-nt gene sat at the 99.7th percentile
genome-wide but the 97.6th against length-matched genes. Step 09 reports both.

## 8. Some GEO tables are not UTF-8

`GSE119724_D23.17cond.tpm.tsv` is latin-1 (a `25°C` column header). Reading it
as UTF-8 raises `UnicodeDecodeError`. Set `expression.encoding`.

## 9. UniProt proteome records can be empty

Proteome `UP000002622` (D23580) reports a protein count but has zero
individually accessioned UniProtKB entries — `proteomeStatistics` shows
`reviewedProteinCount: 0, unreviewedProteinCount: 0`. Step 02 exists to work
around exactly this. Check before assuming your organism needs it.

## 10. Statistical discipline

An apparently striking pattern in the D23580 run — zero overlap between
"still uncharacterised" and "has a signal peptide" — evaporated under a test:
expected 1.3 by chance, Fisher exact **p = 0.63**. Test before you narrate.
