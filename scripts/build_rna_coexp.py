"""Build the RNA co-expression network (W2-W3 task, NumPy/BLAS fast path).

Reads:
    data/raw/networks/funmap/funmap_input_expression_data.tgz
        data_freeze_v1.1/<COHORT>/<COHORT>_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.cct
        .cct format: rows = Ensembl gene IDs, cols = patient samples,
        values = log2(RSEM+1) UQ-normalised expression.
        10 tumour cohorts available: BRCA CCRCC COAD GBM HNSCC LSCC LUAD OV PDAC UCEC
        (HCC has a ComBat-normalised variant; we use the un-ComBat version).
Writes:
    data/processed/networks/rna_coexp.tsv        (contract C2)
    data/processed/networks/rna_coexp.meta.json  (contract C2b)

Decisions referenced:
    D-05 RNA co-expression:
        Pearson per cohort, Fisher z averaged across cohorts, top 196,800
        pairs by mean z (density matched to FunMap), then LCC.
    D-02 HGNC symbol mapping via Person B's d1.hgnc.HGNCMapper (ensembl -> symbol).

Why a NumPy/BLAS fast path
-------------------------
pandas.corr() on an 18,000 x 18,000 matrix takes ~5-10 minutes per cohort
because pandas iterates with Python-level loops that miss the SIMD /
multithreaded BLAS GEMM. With NumPy we compute

    X_centered = X - X.mean(axis=1, keepdims=True)
    X_std     = X_centered / X.std(axis=1, keepdims=True)
    C         = (X_std @ X_std.T) / (n - 1)            # GEMM, BLAS

which is the same algorithm pandas implements, but the matrix multiply
is dispatched to the BLAS library (OpenBLAS / MKL on this machine) which
runs in tens of seconds per cohort. The two implementations are
numerically equivalent to ~1e-10 (see tests/test_rna_coexp.py).
"""
import os
import sys
import tarfile
from collections import defaultdict

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from d1.networks.io import finalize_edges, check_network_file  # noqa: E402
from d1.hgnc import HGNCMapper  # noqa: E402

TARBALL = os.path.join(ROOT, "data", "raw", "networks", "funmap",
                       "funmap_input_expression_data.tgz")
MAPPING_PATH = os.path.join(ROOT, "data", "raw", "reference", "hgnc_complete_set.txt")
OUT_DIR = os.path.join(ROOT, "data", "processed", "networks")
NET_ID = "rna_coexp"
TARGET_EDGES = 196_800                       # D-05: density matched to FunMap
TOP_PER_COHORT = 200_000                     # see docstring
MIN_GENES = 200                               # skip cohorts with fewer usable genes


def _list_cohorts():
    """Return list of (cohort, member_name) tuples for Tumor RNA-seq .cct files."""
    cohorts = []
    with tarfile.open(TARBALL, "r:gz") as t:
        for m in t.getmembers():
            base = os.path.basename(m.name)
            if base.startswith("._"):                       # macOS resource fork
                continue
            if not base.endswith("_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.cct"):
                continue
            if m.size < 100_000:                            # tiny = metadata
                continue
            cohort = m.name.split("/")[1]
            cohorts.append((cohort, m.name))
    return sorted(cohorts)


def _load_cohort(name: str, mapper: HGNCMapper) -> pd.DataFrame:
    """Read one .cct, ENSG -> HGNC, drop genes with >50% zeros.

    CPTAC files are cp1252 / latin-1 encoded (US Windows), not utf-8.
    """
    def _read_with(enc):
        with tarfile.open(TARBALL, "r:gz") as t:
            f = t.extractfile(name)
            return pd.read_csv(f, sep="\t", index_col=0, low_memory=False, encoding=enc)

    try:
        df = _read_with("utf-8")
    except UnicodeDecodeError:
        df = _read_with("latin-1")
    df = df.astype(float)
    keep = (df == 0).sum(axis=1) < 0.5 * df.shape[1]
    df = df.loc[keep]
    df.index = mapper.map(df.index, "ensembl")
    df = df.loc[df.index.notna()]
    df = df[~df.index.duplicated(keep="first")]
    return df


def _pearson_numpy(expr: pd.DataFrame) -> np.ndarray:
    """Compute pairwise Pearson correlation via NumPy/BLAS.

    Parameters
    ----------
    expr : DataFrame (genes x samples), no NaN, no zero-variance rows.

    Returns
    -------
    C : ndarray (genes x genes), Pearson r, ones on the diagonal.

    Algorithm
    ---------
    X is genes x samples. Centre each row (gene) and divide by its standard
    deviation; then the Pearson correlation between genes i and j is the
    inner product of the standardised rows divided by (n_samples - 1):

        r_ij = <x_i_norm, x_j_norm> / (n - 1)

    The standardisation step implicitly divides by sqrt(n-1) on each side,
    so the inner product is the correlation (up to the (n-1) factor).

    This calls `np.matmul` which dispatches to BLAS GEMM, multi-threaded.
    """
    X = expr.to_numpy()                                   # (G, S)
    X = X - X.mean(axis=1, keepdims=True)                 # centre each gene
    sd = X.std(axis=1, ddof=1, keepdims=True)             # (G, 1)
    # Constant rows have sd = 0 -> guard against division by zero
    safe_sd = np.where(sd > 0, sd, 1.0)
    X = X / safe_sd
    G = X.shape[0]
    C = (X @ X.T) / max(X.shape[1] - 1, 1)               # (G, G)
    # Re-impose 1 on the diagonal (BLAS rounding can drift)
    np.fill_diagonal(C, 1.0)
    # Restore -1 < C < 1 just in case; arctanh needs abs < 1
    C = np.clip(C, -0.999999, 0.999999)
    return C


def _top_pairs_numpy(expr: pd.DataFrame, top_k: int) -> dict[tuple, float]:
    """Per-cohort top correlations (NumPy/BLAS Pearson + Fisher z).

    Returns dict {(gene_a, gene_b): z} with |z| largest first.
    """
    # Drop zero-variance genes (std=0 -> Pearson undefined)
    expr = expr.loc[expr.std(axis=1) > 0]
    if expr.shape[0] < 2:
        return {}
    C = _pearson_numpy(expr)
    n = C.shape[0]
    # Restrict to upper triangle (no self, no duplicates)
    iu = np.triu_indices(n, k=1)
    r = C[iu]                                            # flat array of upper-tri
    z = np.arctanh(r)
    # Top-K by |z|
    if len(z) <= top_k:
        idx = np.arange(len(z))
    else:
        idx = np.argpartition(-np.abs(z), top_k - 1)[:top_k]
        # Now refine to the actual top by sorting the candidate slice
        idx = idx[np.argsort(-np.abs(z[idx]))]
    labels = list(expr.index)
    out = {}
    for k in idx:
        a, b = labels[iu[0][k]], labels[iu[1][k]]
        if a > b:
            a, b = b, a
        out[(a, b)] = float(z[k])
    return out


def main():
    print("RNA co-expression network (W2, NumPy fast path)")
    print("Step 1: HGNC mapper (ensembl -> symbol)")
    mapper = HGNCMapper(MAPPING_PATH)
    print(f"  mapper covers {len(mapper.maps['ensembl']):,} Ensembl -> HGNC entries")

    print("\nStep 2: list tumor RNA-seq cohorts")
    cohorts = _list_cohorts()
    print(f"  found {len(cohorts)} cohorts: {[c for c, _ in cohorts]}")

    print("\nStep 3: per-cohort Pearson + Fisher z (top-K per cohort)")
    z_sum = defaultdict(float)
    z_count = defaultdict(int)
    cohort_count = 0
    for cohort, member in cohorts:
        import time
        t0 = time.time()
        expr = _load_cohort(member, mapper)
        if expr.shape[0] < MIN_GENES:
            print(f"  [{cohort}] skip: only {expr.shape[0]} genes after cleaning")
            continue
        print(f"  [{cohort}] genes={expr.shape[0]:,} samples={expr.shape[1]:,}",
              end=" ... ", flush=True)
        top = _top_pairs_numpy(expr, TOP_PER_COHORT)
        for pair, zval in top.items():
            z_sum[pair] += zval
            z_count[pair] += 1
        cohort_count += 1
        print(f"got {len(top):,} top pairs ({time.time() - t0:.1f}s)")
    print(f"\n  {cohort_count} cohorts contributed, "
          f"{len(z_sum):,} unique pairs in the union")

    print("\nStep 4: average Fisher z across cohorts, take top 196,800")
    mean_z = pd.Series({p: z_sum[p] / z_count[p] for p in z_sum})
    top = mean_z.abs().sort_values(ascending=False).head(TARGET_EDGES)
    print(f"  selected {len(top):,} pairs")

    print("\nStep 5: build edge DataFrame for finalize_edges")
    edges = pd.DataFrame(
        [(a, b, mean_z[(a, b)]) for (a, b) in top.index],
        columns=["gene_a", "gene_b", "weight"],
    )
    edges["gene_a"] = edges["gene_a"].astype(str)
    edges["gene_b"] = edges["gene_b"].astype(str)
    print(f"  edges before finalize: {len(edges):,}")

    print("\nStep 6: finalize_edges (LCC + drop self-loops + dedup)")
    out_df, meta = finalize_edges(
        edges, net_id=NET_ID,
        meta={"source_file": "funmap_input_expression_data.tgz",
              "threshold": "top 196,800 by mean |Fisher z|",
              "score_column": "mean Fisher z of Pearson correlation",
              "n_cohorts_used": cohort_count,
              "n_cohorts_total": len(cohorts)},
        out_dir=OUT_DIR,
    )
    print(f"  edges after finalize: {len(out_df):,}  LCC nodes: {meta['n_nodes_lcc']:,}")
    print(f"  self-loops removed: {meta['n_self_loops_removed']}")
    print(f"  duplicates merged: {meta['n_edges_dedup'] - len(out_df)}")

    check_network_file(os.path.join(OUT_DIR, f"{NET_ID}.tsv"))
    print(f"  OK: {NET_ID}.tsv passes contract C2")


if __name__ == "__main__":
    main()