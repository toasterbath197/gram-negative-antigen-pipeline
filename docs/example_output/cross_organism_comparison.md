# Cross-organism validation

The audit in step 03 was run on two Gram-negative organisms with unrelated
genome annotations, to check that the pipeline is not tuned to one dataset.

| | *S.* Typhimurium D23580 | *A. baumannii* ATCC 17978 |
|---|---|---|
| Reference accession | `FN424405` (2009) | `CP000521` (2007) |
| CDS labelled "hypothetical" | 763 | 1141 |
| Mapped to a UniProt accession | 755 | 879 |
| No ortholog found | 8 | 262 |
| Swiss-Prot reviewed | 227 | 22 |
| Named function (unreviewed) | 222 | 653 |
| Generic localization only | 273 | 0 |
| **Still uncharacterized** | **33** | **204** |
| % of all "hypothetical" | 4.3% | 17.9% |
| % of mapped proteins | 4.4% | 23.2% |

## What this shows

**The effect reproduces.** In both organisms the large majority of proteins
labelled "hypothetical" in the reference genome are no longer uncharacterized:
95.7% stale in D23580, 82.1% in *A. baumannii*. Filtering on the deposited
`/product` field is unsafe in both.

**Its size depends on how well studied the organism is.** Comparing like with
like — uncharacterized as a fraction of proteins that mapped — D23580 gives
4.4% against 23.2% for *A. baumannii*, a
5.3-fold difference. The
Swiss-Prot counts tell the same story from the other side: 227 curated
entries for Salmonella against 22 for *A. baumannii*. Salmonella has
been a laboratory workhorse since the 1950s and its proteome has been picked over;
*A. baumannii* has not. This matches the UniProt "Uncharacterized protein"
fractions quoted in the README (8.7% vs 11.4%) and supports the advice there on
choosing an organism.

## Two caveats, stated plainly

**Ortholog coverage differs and the numbers are not perfectly comparable.**
D23580 has a near-identical, well-annotated relative (LT2) so 755/763
mapped. *A. baumannii* ATCC 17978 has no populated UniProt proteome record and
only a patchy species-level set, so 262 proteins (23%)
found no ortholog at all. Those are counted as *not* uncharacterized, so
17.9% is a **floor** — the true figure for *A. baumannii* is
higher.

**The `GENERIC_localization_only` class is annotation-style dependent.** It
captured 273 Salmonella proteins and **0** in *A. baumannii*.
Names like "Inner membrane protein" are a convention in the enterics UniProt
entries that *A. baumannii* entries do not use. The classifier's `GENERIC` bucket
therefore does not transfer between annotation styles, and the candidate pool for
*A. baumannii* comes entirely from the truly-uncharacterized class. Anyone porting
this to a third organism should inspect the actual protein-name distribution
before trusting that bucket.

## Reproducing

```bash
PIPELINE_CONFIG=config/config_acinetobacter.yaml PIPELINE_ORG=AB17978 \
  python3 scripts/01_extract_hypothetical.py
PIPELINE_CONFIG=config/config_acinetobacter.yaml PIPELINE_ORG=AB17978 \
  python3 scripts/02_map_orthologs.py
PIPELINE_CONFIG=config/config_acinetobacter.yaml PIPELINE_ORG=AB17978 \
  python3 scripts/03_reannotate.py
```
