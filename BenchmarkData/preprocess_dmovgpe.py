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
    is not reported and is absent from the released code. We use multiplicative
    replacement (Martin-Fernandez et al. 2003, Math Geol 35:253) on the closed
    composition with a per-feature delta_j = 0.65 x (smallest non-zero value of
    feature j), i.e. 65% of the empirical detection limit. Non-zero parts in a
    row are shrunk by (1 - sum of deltas replaced in that row) so rows sum to 1.
    (scikit-bio's default delta = 1/D^2 was rejected: for the low-dimensional 16S
    tables it exceeds the smallest observed abundance in 42-99% of samples.)
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


DL_FRACTION = 0.65


def replacement_deltas(closed: pd.DataFrame) -> np.ndarray:
    """Per-feature delta_j = DL_FRACTION * min non-zero value of feature j."""
    return DL_FRACTION * closed.where(closed > 0).min(axis=0).values


def multiplicative_replacement(df: pd.DataFrame) -> pd.DataFrame:
    xc = closure(df)
    delta = replacement_deltas(xc)
    x = xc.values
    z = x == 0
    shrink = 1 - (z * delta[None, :]).sum(axis=1, keepdims=True)
    out = np.where(z, delta[None, :], x * shrink)
    return pd.DataFrame(out, index=df.index, columns=df.columns)


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
               microbe_delta_median=float(np.median(replacement_deltas(closure(mic_f)))),
               metabolite_delta_median=None if mtb_is_log else float(np.median(replacement_deltas(closure(mtb_f)))),
               microbe_zero_frac_after_filter=float((mic_f == 0).to_numpy().mean()),
               metabolite_zero_frac_after_filter=float((mtb_f == 0).to_numpy().mean()))

    mic_clr = clr(multiplicative_replacement(mic_f))
    mtb_clr = clr_logspace(mtb_f) if mtb_is_log else clr(multiplicative_replacement(mtb_f))
    return dict(microbe_filtered=closure(mic_f), microbe_clr=mic_clr,
                metabolite_filtered=mtb_f, metabolite_clr=mtb_clr), log
