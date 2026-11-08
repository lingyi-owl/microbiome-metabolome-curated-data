# BenchmarkData

Harmonised inputs for benchmarking metagenome-to-metabolome abundance prediction tools
(e.g. MelonnPan, SparseNED, MiMeNet, mNODE, LOCATE, DMoVGPE, SIMBA-GNN).

## Contents

One folder per dataset of `data/processed_data/`, preprocessed following the DMoVGPE Methods
(Weng et al. 2025, *BMC Bioinformatics* 26:93, doi:10.1186/s12859-025-06110-7):

| File | Content |
|---|---|
| `microbe_filtered.tsv` | taxon relative abundance after the >50%-zero filter, re-closed to sum 1 per sample |
| `microbe_clr.tsv` | after filter, zero replacement and CLR |
| `metabolite_filtered.tsv` | metabolite abundance (original scale) after the >50%-zero filter |
| `metabolite_clr.tsv` | after filter, zero replacement and CLR |
| `metadata.tsv` | sample metadata (incl. `Study.Group`) |
| `microbe_proteom.csv` | one row per `microbe_clr.tsv` column: `gtdb_lineage` (exact column header); `proteome_match` = `exact` (shotgun: GTDB species representative; 16S: representative of the genus's GTDB type species), `fallback` (16S genus with no type-species genome in GTDB; path left empty) or `none` (unresolved genus); `gtdb_release` (r207 / r202) and `protein_faa_in_archive` (path inside `gtdb_proteins_aa_reps_<release>.tar.gz`), filled for `exact` rows |

All tables are samples × features, first column `Sample`, and share the same sample order within a dataset.
Files larger than 25 MB are gzip-compressed (`.tsv.gz`; `pandas.read_csv` reads them directly).
The two `iHMP_IBDMDB_2019` metabolite tables (50,481 metabolites) are split by columns into five parts
(`metabolite_*.partNof5.tsv.gz`, each with the `Sample` column) to stay under GitHub's file-size limits. Rejoin with:

```python
import pandas as pd
mtb = pd.concat([pd.read_csv(f"metabolite_clr.part{i}of5.tsv.gz", sep="\t", index_col=0) for i in range(1, 6)], axis=1)
```
`preprocessing_log.tsv` records per-dataset counts and replacement values; `preprocess_dmovgpe.py` is the code.

| DMoVGPE name | Folder | Samples | Microbe input | Taxa (in → kept) | Metabolites (in → kept) | Metabolite input scale |
|---|---|---|---|---|---|---|
| MFGM | `HE_INFANTS_MFGM_2019` | 277 | `genera.tsv` (genus) | 295 → 13 | 120 → 81 | linear |
| ASD | `KANG_AUTISM_2017` | 44 | `genera.tsv` (genus) | 321 → 55 | 61 → 53 | linear |
| INFANTS_DIABETES | `KOSTIC_INFANTS_DIABETES_2015` | 103 | `genera.tsv` (genus) | 584 → 33 | 184 → 175 | linear |
| PRETERMS | `WANDRO_PRETERMS_2018` | 75 | `genera.tsv` (genus) | 82 → 4 | 453 → 453 | linear |
| SINHA_CRC | `SINHA_CRC_2016` | 131 | `genera.tsv` (genus) | 86 → 31 | 530 → 530 | log |
| ADENOMAS | `KIM_ADENOMAS_2020` | 240 | `genera.tsv` (genus) | 499 → 54 | 462 → 462 | log |
| BIO_ML | `POYET_BIO_ML_2019` | 164 | `genera.tsv` (genus) | 353 → 59 | 489 → 471 | linear |
| IBD_RELATIVES | `JACOBS_IBD_FAMILIES_2016` | 90 | `genera.tsv` (genus) | 508 → 117 | 4626 → 1052 | linear |
| IBS | `MARS_IBS_2020` | 444 | `species.tsv` (species) | 11513 → 1186 | 43 → 39 | linear |
| CRC | `YACHIDA_CRC_2019` | 347 | `species.tsv.zip` (species) | 57702 → 6119 | 450 → 164 | linear |
| GC | `ERAWIJANTARI_GASTRIC_CANCER_2020` | 96 | `species.tsv` (species) | 48243 → 6441 | 524 → 144 | linear |
| ESRD | `WANG_ESRD_2020` | 287 | `species.tsv` (species) | 56961 → 7491 | 276 → 276 | linear |
| IBD | `FRANZOSA_IBD_2019` | 220 | `species.tsv` (species) | 55882 → 5604 | 8848 → 7156 | linear |
| IBDMDB | `iHMP_IBDMDB_2019` | 382 | `species.tsv` (species) | 42871 → 2875 | 81867 → 50481 | linear |

## Preprocessing

1. **NaN → 0**, as in the released `DMoVGPE.py` (`np.nan_to_num(X, nan=0.0)`).
2. **Feature filter:** keep a taxon or metabolite only if ≤ 50% of samples are zero.
3. **Zero replacement:** multiplicative replacement (Martín-Fernández et al. 2003, *Math Geol* 35:253) on the closed
   composition, with a per-feature value δ<sub>j</sub> = 0.65 × the smallest non-zero value of feature *j*
   (65% of the empirical detection limit). Non-zero parts are shrunk so each sample still sums to 1.
4. **CLR:** clr(x)<sub>j</sub> = log x<sub>j</sub> − mean<sub>k</sub>(log x<sub>k</sub>) (natural log).

## Deviations and assumptions

- **Taxonomic level.** DMoVGPE does not state it. Matching its Table 1 feature counts indicates genus (`genera.tsv`)
  for the 16S datasets and species (`species.tsv`) for the shotgun datasets; this mapping is used here.
- **Zero-replacement value.** DMoVGPE reports only "a zero-replacement function"; its value is neither in the paper
  nor in the released code. scikit-bio's default (δ = 1/D², D = number of features) was rejected, because for the
  low-dimensional 16S tables it exceeds the smallest observed abundance in 42–99% of samples and so inverts ranks.
- **Log-scaled metabolite tables.** `SINHA_CRC_2016` and `KIM_ADENOMAS_2020` metabolites are supplied already
  log-transformed (about half the values are negative). They get no zero replacement, and CLR is applied in log
  space as x − rowmean(x), which equals clr(exp(x)).
- Feature counts can differ slightly from DMoVGPE Table 1, which mixes filtered and unfiltered counts
  (see the dataset map in the analysis notes).

## Source

Muller E., Algavi Y.M., Borenstein E. (2022) The gut microbiome-metabolome dataset collection: a curated resource for
integrative meta-analysis. *npj Biofilms and Microbiomes* 8:79. doi:10.1038/s41522-022-00345-5
