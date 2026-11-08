"""
Preprocess the 14 curated gut microbiome-metabolome datasets following the
DMoVGPE (Weng et al. 2025, BMC Bioinformatics 26:93) Methods:

  1. Filter out features (taxa / metabolites) with > 50% zero abundance
  2. Replace remaining zeros with a small value
  3. Centered log-ratio (CLR) transform

Input : data/processed_data/<DATASET>/{genera.tsv|species.tsv, mtb.tsv, metadata.tsv}
Output: BenchmarkData/<DATASET>/
          microbe_filtered.tsv     relative abundance after step 1 (re-closed to sum 1)
          microbe_clr.tsv          after steps 1-3
          metabolite_filtered.tsv  original scale after step 1
          metabolite_clr.tsv       after steps 1-3
          metadata.tsv             samples in the same order

Choices / assumptions (DMoVGPE does not publish its preprocessing code):
  * Taxonomic level: genera.tsv for 16S datasets, species.tsv for shotgun
    datasets (inferred by matching DMoVGPE Table 1 feature counts).
  * NaN -> 0 before filtering (as in DMoVGPE.py: np.nan_to_num(X, nan=0.0)).
  * Zero test: a feature is kept if its fraction of zeros is <= 0.5.
  * Zero replacement: DMoVGPE only says "a zero-replacement function"; the value
    is not reported and is absent from the released code. Here each zero cell of
    the closed (relative-abundance, rows sum to 1) table is set to a fixed
    pseudocount of 1e-5 (PSEUDOCOUNT); non-zero cells are left unchanged.
    Metabolite tables on a linear scale are closed to relative abundance first,
    so the pseudocount is on the same scale for microbes and metabolites.
  * Metabolite tables that are already log-scaled (contain negative values:
    SINHA_CRC_2016, KIM_ADENOMAS_2020) are not zero-replaced; CLR is applied in
    log space as x - rowmean(x), which equals clr(exp(x)).
  * Samples whose filtered microbe or metabolite row is all zero are dropped.
"""
import numpy as np
import pandas as pd

ZERO_FRAC_MAX = 0.5


def zero_filter(df: pd.DataFrame) -> pd.DataFrame:
    keep = (df == 0).mean(axis=0) <= ZERO_FRAC_MAX
    return df.loc[:, keep]


def closure(df: pd.DataFrame) -> pd.DataFrame:
    return df.div(df.sum(axis=1), axis=0)


PSEUDOCOUNT = 1e-5


def pseudocount_replacement(df: pd.DataFrame, pseudo: float = PSEUDOCOUNT) -> pd.DataFrame:
    """Close rows to relative abundance, then set zero cells to `pseudo`."""
    xc = closure(df)
    return xc.mask(xc == 0, pseudo)


def clr(df: pd.DataFrame) -> pd.DataFrame:
    lg = np.log(df.values)
    return pd.DataFrame(lg - lg.mean(axis=1, keepdims=True), index=df.index, columns=df.columns)


def clr_logspace(df: pd.DataFrame) -> pd.DataFrame:
    v = df.values
    return pd.DataFrame(v - v.mean(axis=1, keepdims=True), index=df.index, columns=df.columns)


def process_pair(mic: pd.DataFrame, mtb: pd.DataFrame):
    """mic, mtb: samples x features, numeric, same index. Returns dict of outputs + log."""
    log = {}
    mic = mic.fillna(0.0)
    mtb_is_log = bool((mtb < 0).to_numpy().any())
    mtb = mtb.fillna(0.0)
    log.update(n_samples_in=len(mic), microbe_features_in=mic.shape[1], metabolite_features_in=mtb.shape[1],
               metabolite_input_scale="log" if mtb_is_log else "linear")

    mic_f = zero_filter(mic)
    mtb_f = zero_filter(mtb) if not mtb_is_log else mtb.loc[:, (mtb == 0).mean(axis=0) <= ZERO_FRAC_MAX]

    # drop samples with an all-zero filtered row
    bad = (mic_f.sum(axis=1) == 0)
    if not mtb_is_log:
        bad |= (mtb_f.sum(axis=1) == 0)
    mic_f, mtb_f = mic_f.loc[~bad], mtb_f.loc[~bad]
    log.update(samples_dropped_all_zero=int(bad.sum()), n_samples_out=len(mic_f),
               microbe_features_out=mic_f.shape[1], metabolite_features_out=mtb_f.shape[1],
               pseudocount=PSEUDOCOUNT,
               microbe_frac_zero_cells_replaced=float((mic_f == 0).to_numpy().mean()),
               metabolite_frac_zero_cells_replaced=None if mtb_is_log else float((mtb_f == 0).to_numpy().mean()),
               microbe_zero_frac_after_filter=float((mic_f == 0).to_numpy().mean()),
               metabolite_zero_frac_after_filter=float((mtb_f == 0).to_numpy().mean()))

    mic_clr = clr(pseudocount_replacement(mic_f))
    mtb_clr = clr_logspace(mtb_f) if mtb_is_log else clr(pseudocount_replacement(mtb_f))
    return dict(microbe_filtered=closure(mic_f), microbe_clr=mic_clr,
                metabolite_filtered=mtb_f, metabolite_clr=mtb_clr), log
