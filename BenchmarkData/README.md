# BenchmarkData

Harmonised inputs for benchmarking metagenome-to-metabolome abundance prediction tools
(e.g. MelonnPan, SparseNED, MiMeNet, mNODE, MMINP, LOCATE, MICOM, AGORA2/mgPipe, SIMBA-GNN).

## Planned contents

For each dataset in `data/processed_data/<DATASET>/`:

- Paired microbiome (species / genus relative abundance) and metabolome tables,
  restricted to shared samples and aligned in the same sample order
- Matching sample metadata
- Train/test split definitions, including cross-dataset (train on one cohort, test on another) splits

## Source

Derived from the curated gut microbiome-metabolome collection
(Muller et al., npj Biofilms and Microbiomes 8:79, 2022; doi:10.1038/s41522-022-00345-5).
