# Worked example: *S.* Typhimurium D23580

Outputs from a full run with the shipped config.

| File | Contents |
|---|---|
| `D23580_03_reannotated.tsv` | All 763 proteins with current UniProt annotation and audit class |
| `D23580_final_scorecard.tsv` | Ranked candidates with structure, conservation and expression evidence |
| `D23580_negative_controls.json` | Epitope scores for six known cytoplasmic proteins |

Headline numbers reproduced by `03_reannotate.py`:

```
    273  GENERIC_localization_only
    227  CHARACTERIZED_reviewed
    222  NAMED_function_unreviewed
     33  TRULY_UNCHARACTERIZED
      8  NO_ORTHOLOG

  Still genuinely uncharacterised: 33/763 (4.3%)
  Annotation now STALE for:        730/763 (95.7%)
```

Compare `D23580_negative_controls.json` against the candidate rows in the
scorecard. GroEL — a cytoplasmic chaperonin that cannot be a surface antigen —
returns 13 BepiPred epitopes and a Kolaskar-Tongaonkar mean of 1.022. The best
real candidate returns 7 and 1.025. That is why the epitope terms are switched
off in scoring by default.
