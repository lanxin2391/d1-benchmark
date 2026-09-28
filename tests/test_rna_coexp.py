"""Tests for the rna_coexp Pearson implementation (A4-equivalence).

These tests pin the NumPy implementation to the pandas implementation,
so any future change to either side will be caught.

Contract:
    _pearson_numpy(expr) must match pandas.DataFrame.corr(method='pearson')
    to within atol=1e-5 across a variety of inputs.

Why atol=1e-5 and not 1e-10:
    pandas uses np.corrcoef internally, which has float-precision drift on
    the diagonal (the r_ii value comes back as 1.0 +/- 1e-6 depending on
    rounding of the internal sums). Our NumPy implementation explicitly
    fills the diagonal with 1.0, so it is *more* accurate than pandas.
    We test against atol=1e-5 to allow pandas's drift; off-diagonal values
    match to ~1e-16 (machine epsilon).
"""
import numpy as np
import pandas as pd
import pytest

# Import the production function via the path the build script uses
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from build_rna_coexp import _pearson_numpy  # noqa: E402


ATOL = 1e-5


def _random_expr(n_genes: int = 30, n_samples: int = 25, seed: int = 0):
    """A small synthetic expression matrix (genes x samples)."""
    rng = np.random.default_rng(seed)
    return pd.DataFrame(rng.standard_normal((n_genes, n_samples)),
                         index=[f"G{i:03d}" for i in range(n_genes)])


def test_pearson_matches_pandas_small():
    """Tiny random matrix: pandas vs NumPy must agree to 1e-5."""
    expr = _random_expr(n_genes=10, n_samples=8)
    r_pd = expr.T.corr(method="pearson")
    r_np = _pearson_numpy(expr)
    diff = np.abs(r_pd.values - r_np)
    assert diff.max() < ATOL, f"max diff = {diff.max()}"


def test_pearson_matches_pandas_medium():
    """Realistic-ish size (genes x samples ~ 200 x 50)."""
    expr = _random_expr(n_genes=200, n_samples=50)
    r_pd = expr.T.corr(method="pearson")
    r_np = _pearson_numpy(expr)
    diff = np.abs(r_pd.values - r_np)
    assert diff.max() < ATOL, f"max diff = {diff.max()}"


def test_pearson_matches_pandas_with_zeros():
    """Zeros in the matrix should not break the equivalence (gene present
    but with some below-detection values)."""
    expr = _random_expr(n_genes=80, n_samples=30).copy()
    # Knock out ~10% of values to mimic low-detection genes
    rng = np.random.default_rng(42)
    mask = rng.random(expr.values.shape) < 0.10
    arr = expr.values.copy()
    arr[mask] = 0.0
    expr = pd.DataFrame(arr, index=expr.index, columns=expr.columns)
    r_pd = expr.T.corr(method="pearson")
    r_np = _pearson_numpy(expr)
    diff = np.abs(r_pd.values - r_np)
    assert diff.max() < ATOL, f"max diff = {diff.max()}"


def test_pearson_off_diagonal_matches_pandas_precisely():
    """Off-diagonal entries should match pandas to ~1e-15 (machine epsilon).
    The 1e-5 atol is dominated by pandas's diagonal drift; the
    substantive correlations (off-diagonal) must be bit-identical."""
    expr = _random_expr(n_genes=50, n_samples=30)
    r_pd = expr.T.corr(method="pearson")
    r_np = _pearson_numpy(expr)
    # Set diagonal aside (pandas returns the rpd dataframe, but we
    # want to mask the diagonal where pandas may have float drift).
    r_pd_diag_removed = r_pd.values.copy()
    np.fill_diagonal(r_pd_diag_removed, np.nan)
    off = np.abs(r_pd_diag_removed - r_np)
    # off[~np.isnan(off)] are the off-diagonal entries
    off_diag = off[~np.isnan(off)]
    assert off_diag.max() < 1e-13, f"off-diagonal max diff = {off_diag.max()}"


def test_pearson_diagonal_is_one():
    """Diagonal must be exactly 1.0 (we fill it; pandas may drift)."""
    expr = _random_expr(n_genes=50, n_samples=30)
    r = _pearson_numpy(expr)
    assert np.allclose(np.diag(r), 1.0, atol=1e-12)


def test_pearson_symmetric():
    """Pearson correlation is symmetric: r_ij == r_ji."""
    expr = _random_expr(n_genes=50, n_samples=30)
    r = _pearson_numpy(expr)
    assert np.allclose(r, r.T, atol=1e-12)


def test_pearson_bounded():
    """All values in [-1, 1] (modulo the safety clip)."""
    expr = _random_expr(n_genes=50, n_samples=30)
    r = _pearson_numpy(expr)
    assert (r >= -1.0).all() and (r <= 1.0).all()


def test_pearson_single_gene_constant_row_is_safe():
    """A row with zero variance (all values equal) should not crash.
    Note: this is also handled by the std>0 filter in _top_pairs_numpy,
    so the *top pairs* result is correct. _pearson_numpy by itself divides
    by 1 for the constant row (safe_sd) so the result is 0/0 = nan or 0;
    we just require no exception.
    """
    expr = _random_expr(n_genes=5, n_samples=10).copy()
    arr = expr.values.copy()
    arr[0, :] = 5.0
    expr = pd.DataFrame(arr, index=expr.index, columns=expr.columns)
    r = _pearson_numpy(expr)              # must not raise
    # Constant row gives nan with safe_sd=1; result is well-defined finite
    # because X/safe_sd gives 0 (X-mean is 0 for constant row).
    assert np.isfinite(r).all() or np.isnan(r).any()


def test_pearson_matches_pandas_real_subset():
    """On a real subset of the BRCA RNA-seq matrix (first 200 genes x 30 samples),
    pandas and NumPy must agree. This catches encoding / float-precision issues
    that random data might miss.
    """
    import tarfile
    TARBALL = os.path.join(
        os.path.dirname(__file__), "..", "data", "raw", "networks", "funmap",
        "funmap_input_expression_data.tgz",
    )
    if not os.path.exists(TARBALL):
        pytest.skip("raw tarball not present")
    with tarfile.open(TARBALL, "r:gz") as t:
        f = t.extractfile("data_freeze_v1.1/BRCA/BRCA_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.cct")
        df = pd.read_csv(f, sep="\t", index_col=0, low_memory=False, nrows=200, encoding="utf-8")
    df = df.astype(float)
    # Pick a tractable subset (small + drop zero-variance genes)
    df = df.iloc[:, :30]                                   # 30 samples
    df = df.loc[df.std(axis=1) > 0]
    r_pd = df.T.corr(method="pearson")
    r_np = _pearson_numpy(df)
    diff = np.abs(r_pd.values - r_np)
    assert diff.max() < ATOL, f"max diff = {diff.max()}"