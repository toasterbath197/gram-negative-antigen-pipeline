# Gram-Negative Antigen Triage Pipeline

A sequential, resumable pipeline that takes a Gram-negative bacterial genome
annotation and narrows its unannotated proteins down to a ranked shortlist of
plausible vaccine-antigen candidates — using structure, localization,
conservation and infection-condition expression.

It was built and validated on *Salmonella enterica* serovar Typhimurium **D23580**
(ST313, African invasive non-typhoidal lineage), and is parameterised by config
so it can be pointed at another organism.

**What this pipeline does not do:** prove anything. Every output is a prediction.
Nothing here establishes that a protein is surface-exposed in vivo, translated,
or immunogenic. It is a triage tool for deciding what to test in a lab.

---

## The finding that motivated it

Running step 03 on D23580 produced the result that reshaped the whole project:

> Of 763 proteins annotated `hypothetical protein` in the 2009 reference genome
> (`FN424405`), only **33 (4.3%)** are still uncharacterised in current UniProt.
> 227 are Swiss-Prot reviewed with full functional assignments. One shortlisted
> "hypothetical" protein matched a solved crystal structure of a *named,
> published Salmonella protein* at 100% sequence identity.

A genome's `/product` field freezes at deposition. Any pipeline that filters on
that field alone is screening proteins that stopped being unknown years ago.
**Step 03 exists to catch this, and it should be the first thing you run on a
new organism**, before investing in anything downstream.

---

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 01 | `01_extract_hypothetical.py` | Parse ENA/EMBL flat file, extract CDS annotated as hypothetical |
| 02 | `02_map_orthologs.py` | BLAST to a reference proteome to obtain UniProt accessions |
| 03 | `03_reannotate.py` | **Re-check every protein against current UniProt.** Audit + candidate pool |
| 04 | `04_localization.py` | DeepTMHMM topology; Enterobacteriaceae lipoprotein +2 sorting rule |
| 05 | `05_fetch_structures.py` | Download AlphaFold models, extract per-residue pLDDT |
| 06 | `06_foldseek.py` | Foldseek search vs PDB100; probability **and** E-value both required |
| 07 | `07_epitopes.py` | DiscoTope-3 + IEDB linear predictors, **with negative controls** |
| 08 | `08_conservation.py` | Presence across all complete RefSeq genomes for the taxon |
| 09 | `09_expression.py` | Expression in infection conditions (optional) |
| 10 | `10_score.py` | Composite scorecard |

```bash
./run_all.sh                      # everything, in order
python3 scripts/03_reannotate.py  # or one step at a time
```

Every step writes to `results/` and skips work already done, so the pipeline is
safe to interrupt and resume.

---

## Install

```bash
git clone <this repo> && cd gram-negative-antigen-pipeline
pip install -r requirements.txt
```

Also required on `PATH`:

- **NCBI BLAST+** (`blastp`, `makeblastdb`) — `apt install ncbi-blast+` or conda
- `curl`, `unzip`, `zip`

DeepTMHMM and DiscoTope-3 run on BioLib's cloud via `pybiolib` — no local
install and no API key. Foldseek and IEDB are used through their public web
APIs. Nothing here needs a GPU.

---

## Running it on a different organism

Edit `config/config.yaml`. The parts that actually matter:

| Setting | Notes |
|---|---|
| `annotation.ena_accession` | Your strain's chromosome |
| `orthologs.enabled` | Set **false** if your organism has its own populated UniProtKB entries. D23580 does not, which is the only reason step 02 exists |
| `localization.gram` / `apply_plus2_rule` | The +2 rule is **Enterobacteriaceae-specific**. Disable it elsewhere. None of the localization logic transfers to Gram-positives or eukaryotes |
| `conservation.taxon_id` | Needs a few hundred complete genomes to mean anything |
| `expression.geo_url` | Set to `null` to skip. This is the least portable step — most organisms have no infection-condition RNA-seq |

**Choosing an organism.** The pipeline pays off where many genomes exist but
molecular characterisation is thin. Fraction of proteins still called
"Uncharacterized protein" in UniProt: *P. aeruginosa* 5.3%, *S.* Typhimurium
8.7%, *K. pneumoniae* 9.9%, *A. baumannii* 11.4%, *N. gonorrhoeae* 15.1%.
Salmonella sits near the bottom — it has been a lab workhorse since the 1950s,
so its proteome was always going to be picked over. *A. baumannii* or
*N. gonorrhoeae* are better targets.

---

## Two things this pipeline will tell you that most won't

**It runs negative controls.** Step 07 pushes known cytoplasmic proteins —
GroEL, DnaK, EF-Tu, enolase, RpoA, AhpC — through the same epitope predictors as
your candidates. On the D23580 data they scored the same: GroEL returned 13
BepiPred epitopes and a Kolaskar-Tongaonkar mean of 1.022, against 7 and 1.025
for the best real candidate. Those metrics carried no discriminating signal, so
they are excluded from scoring by default. Check your own controls before
turning them back on.

**It corrects a silent bug in a dependency.** The BioLib DiscoTope-3 wrapper
overwrites every B-factor with `100.00`, destroying the pLDDT values DiscoTope
uses as a feature for AlphaFold models. Step 07 re-filters epitope residues
against pLDDT parsed from the original models. See `docs/KNOWN_ISSUES.md`.

---

## Output

`results/10_FINAL_SCORECARD.tsv` — one row per candidate with annotation class,
topology, lipoprotein sorting, mean pLDDT, Foldseek verdict, conservation,
expression percentile, epitope counts, and a transparent score with the reason
for every point awarded.

The score is a simple tally, not a validated model. Use it to order follow-up
work, nothing more. As a calibration: reverse vaccinology's founding study on
meningococcus B went from ~600 in silico candidates to 28 that raised
bactericidal antibodies to a handful in the licensed vaccine.

---

## Citation of the tools used

AlphaFold DB · Foldseek · DeepTMHMM · DiscoTope-3 · IEDB Analysis Resource ·
SignalP · NCBI BLAST+ · NCBI Datasets · UniProt.
Worked-example expression data: GEO `GSE119724` (Canals et al., *PLOS Biology*
2019).

## License

MIT
