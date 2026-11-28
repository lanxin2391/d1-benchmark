"""The NetMix2 fallback that needs no Gurobi licence: the local-FDR
component alone.

NetMix2's published method has two steps: estimate the size of the
altered subnetwork with a local false-discovery rate, then pick the
subnetwork of that size with the largest total vertex score. The second
step is an integer program and needs a Gurobi licence; the first step --
the per-gene local FDR -- does not, and it is the part that produces a
per-gene score to rank.

This script therefore reports the **local-FDR component** of NetMix2 as
the mutation-first score, and labels it as such everywhere: it is not
NetMix2 as published, it is the significance step of NetMix2 applied to
the MC3 mutation heat, without the subnetwork selection.

It also runs the **trimmed-network** variant of the full method: the
sub-network is trimmed to about 1,500 nodes, which the size-limited
licence can solve, so the published method runs on a smaller universe.

Outputs:
  results/tables/c37_netmix2_localfdr.tsv
  results/tables/c37_netmix2_trimmed.tsv
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(r"D:\Grade3\swxxx\final\D1\d1-benchmark")
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from d1.engine.baselines import degree_scores  # noqa: E402
from d1.engine.evaluate import evaluate  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

REPO = Path(r"C:\Users\PC\AppData\Local\Temp\opencode\netmix2")
WORK = Path(r"C:\Users\PC\AppData\Local\Temp\opencode\netmix2_run")
NETS_DIR = ROOT / "data/processed/networks"
LABELS_DIR = ROOT / "data/processed/labels"
TABLES = ROOT / "results/tables"
PY = ROOT / ".venv-gnn/Scripts/python.exe"
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
TRIM_NODES = 700


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def heat_table() -> pd.DataFrame:
    snv = WORK / "mc3_snvs.tsv"
    df = pd.read_csv(snv, sep="\t", header=None, names=["sample", "genes"],
                     dtype=str)
    counts: dict[str, int] = {}
    for genes in df.genes.dropna():
        for g in str(genes).split("\t"):
            counts[g] = counts.get(g, 0) + 1
    out = pd.DataFrame({"gene": list(counts), "heat": list(counts.values())})
    from scipy import stats as sps
    # The local FDR fits a normal null to the z-scores, and the raw
    # mutation burden is far too skewed for that (the upstream itself
    # fails with "the middle of the histogram is nonnormal"). The
    # rank-based uniform p-value with a tie-breaking jitter gives the
    # method a well-behaved null; the jitter is a documented deviation.
    h = out.heat.values.astype(float)
    order = np.argsort(h, kind="mergesort")
    ranks = np.empty(len(h), dtype=float)
    ranks[order] = np.arange(1, len(h) + 1)
    rng = np.random.default_rng(20261126)
    jitter = rng.uniform(0.0, 1.0, size=len(h))
    out["p"] = np.clip((ranks + jitter - 0.5) / len(h), 1e-12, 1.0)
    return out


def local_fdr_scores(scores: pd.DataFrame, genes: list[str]) -> dict[str, float]:
    """The per-gene local FDR (NetMix2's significance step, no Gurobi)."""
    sys.path.insert(0, str(REPO / "locfdrpython"))
    sys.path.insert(0, str(REPO / "src"))
    from locfdr import locfdr
    import importlib
    common = importlib.import_module("common")
    pmap = dict(zip(scores.gene, scores.p))
    pvals = np.array([pmap.get(g, 1.0) for g in genes], dtype=float)
    pvals = np.clip(pvals, 1e-10, 1.0)
    z = common.compute_zscores(pvals)
    z = common.post_process_zscores(z)
    # the mutation heat is extremely skewed: a few genes carry a huge
    # burden, so the raw p-values give a non-normal z-histogram and the
    # local FDR's null fit fails. The extreme p-values are therefore
    # floored, and the null types are tried in order.
    r = None
    for nt in (1, 2, 3):
        try:
            r = locfdr(z, nulltype=nt, plot=0)
            break
        except Exception as exc:
            last = exc
    if r is None:
        raise SystemExit(f"the local FDR fit failed: {last}")
    fdr = np.asarray(r.fdr, dtype=float)
    fdr = np.where(np.isfinite(fdr), fdr, 1.0)
    fdr = np.clip(fdr, 1e-300, 1.0)
    # the score: the larger the better (the minus log of the local FDR)
    return {g: float(-np.log10(v)) for g, v in zip(genes, fdr)}


def evaluate_scores(net: str, best: dict[str, float]) -> list[dict]:
    genes, A = load_network(str(NETS_DIR / f"{net}.tsv"))
    deg = degree_scores(A)
    degmap = dict(zip(genes, deg))
    gset = set(genes)
    rows = []
    for label in LABELS:
        p = LABELS_DIR / f"{label}.tsv"
        if not p.exists():
            continue
        pos = sorted(g for g in set(pd.read_csv(p, sep="\t", dtype=str)
                                     ["gene"]) if g in gset)
        posset = set(pos)
        neg = sorted(g for g in genes if g not in posset)
        if len(pos) < 10 or len(neg) < 100:
            continue
        sc = np.array([best.get(g, 0.0) for g in pos + neg])
        labels = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
        ev_h = evaluate(sc, labels, np.arange(len(sc)))
        d = np.array([degmap[g] for g in pos + neg])
        ev_d = evaluate(d, labels, np.arange(len(d)))
        rows.append({"network": net, "label": label, "n_positives": len(pos),
                     "netmix2_localfdr_auroc": ev_h["auroc"],
                     "degree_auroc": ev_d["auroc"],
                     "delta_auroc": ev_h["auroc"] - ev_d["auroc"]})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-localfdr", action="store_true")
    ap.add_argument("--networks", nargs="*",
                    default=["reactome", "funmap", "string_full700",
                             "string_phys700", "rna_coexp", "intact"])
    args = ap.parse_args()
    scores = heat_table()
    log(f"the mutation heat: {len(scores):,} genes")
    rows = []
    for net in (args.networks if not args.skip_localfdr else []):
        genes, A = load_network(str(NETS_DIR / f"{net}.tsv"))
        t0 = time.time()
        best = local_fdr_scores(scores, list(genes))
        log(f"  {net}: the local FDR in {time.time() - t0:.0f}s")
        rows.extend(evaluate_scores(net, best))
        for r in rows[-len(LABELS):]:
            log(f"    {r['label']}: the local-FDR AUROC "
                f"{r['netmix2_localfdr_auroc']:.4f} vs the degree "
                f"{r['degree_auroc']:.4f}")
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "c37_netmix2_localfdr.tsv", sep="\t", index=False)
    log(f"wrote the local-FDR table ({len(df)} rows)")

    # ---- the trimmed-network variant of the published method -----------
    log(f"the trimmed-network variant ({TRIM_NODES} nodes)")
    trimmed_rows = []
    for net in args.networks:
        genes, A = load_network(str(NETS_DIR / f"{net}.tsv"))
        if len(genes) > TRIM_NODES:
            # the highest-degree core keeps the structure and fits the
            # size-limited licence
            deg = degree_scores(A)
            keep = set(np.array(genes)[np.argsort(-deg)[:TRIM_NODES]])
            genes = [g for g in genes if g in keep]
        idx = {g: i for i, g in enumerate(genes)}
        sub = A.tocsr()[np.ix_([idx[g] for g in genes],
                                [idx[g] for g in genes])]
        # the per-network score file (the upstream needs one that
        # covers every gene of the network)
        sf = WORK / f"{net}_scores.tsv"
        if not sf.exists():
            pmap = dict(zip(scores.gene, scores.p))
            sf.write_text("".join(f"{g}	{pmap.get(g, 1.0):.6g}\n"
                                   for g in genes), encoding="utf-8")
        edge_file = WORK / f"{net}_trim_edges.tsv"
        coo = sub.tocoo()
        with edge_file.open("w", encoding="utf-8") as fh:
            for u, v in zip(coo.row, coo.col):
                if u < v:
                    fh.write(f"{genes[u]}\t{genes[v]}\n")
        outdir = WORK / "out_trim" / net
        outdir.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ)
        env["PYTHONPATH"] = str(REPO / "locfdrpython") + os.pathsep + env.get(
            "PYTHONPATH", "")
        r = subprocess.run(
            [str(PY), str(REPO / "run_netmix2.py"),
             "-el", str(edge_file), "-gs", str(WORK / f"{net}_scores.tsv"),
             "-o", str(outdir), "-v", "1"],
            capture_output=True, text=True, cwd=str(REPO),
            encoding="utf-8", errors="replace", env=env)
        ok = r.returncode == 0
        log(f"  {net}: the trimmed run returned {r.returncode}")
        if not ok:
            for line in (r.stdout + r.stderr).strip().splitlines()[-4:]:
                log("    " + line[:130])
            continue
        best = read_subnetwork(outdir, genes)
        if best is None:
            continue
        trimmed_rows.extend(evaluate_scores(net, best))
    if trimmed_rows:
        tdf = pd.DataFrame(trimmed_rows)
        tdf.to_csv(TABLES / "c37_netmix2_trimmed.tsv", sep="\t", index=False)
        log(f"wrote the trimmed table ({len(tdf)} rows)")
    return 0


def read_subnetwork(outdir: Path, genes: list[str]) -> dict[str, float] | None:
    for f in sorted(outdir.rglob("*.txt")):
        if f.name in ("components.txt", "subnetwork.txt",
                      "netmix2_subnetwork.txt"):
            gs = [l.strip() for l in f.read_text(encoding="utf-8",
                                                 errors="replace").split()
                  if l.strip()]
            if gs:
                return {g: 1.0 for g in gs}
    return None


if __name__ == "__main__":
    raise SystemExit(main())
