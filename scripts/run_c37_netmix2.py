"""Control 37 -- the NetMix2 arm (the successor of HotNet2).

NetMix2 (Ruan et al., *Journal of Computational Biology* 2022, DOI
10.1089/cmb.2022.0336; the paper is in `data/papers/`) identifies the
altered subnetworks of an interaction network from gene scores, and its
"propagation family" is the subnetwork family that approximates the
network-propagation subnetworks. It is the right comparator for this
paper: unlike HotNet2 it does not need an O(n^2) matrix, so it runs on
the panel's networks.

The design is the one Control 37 already fixed:

* the **gene scores are the MC3 TCGA somatic mutation heat** -- the
  mutation evidence is independent of the labels, so this is a genuine
  mutation-first strategy, not the label re-read;
* the **network** is one of ours (the gene symbols as the edge list);
* the **evaluation** is the same endpoint as everywhere else: the AUROC
  of the per-gene score for the label's positives against the degree
  baseline, on the same genes.

NetMix2 takes the gene scores as p-values (its local-FDR step models the
null from them), so the mutation heat is converted into an upper-tail
empirical p-value per gene.
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
NETWORKS = ["reactome", "funmap", "string_full700", "string_phys700",
            "rna_coexp", "intact"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
# the repository's own (Python 3) port of locfdr has to be importable
ENV = dict(os.environ)
ENV["PYTHONPATH"] = str(REPO / "locfdrpython") + os.pathsep + ENV.get(
    "PYTHONPATH", "")


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def heat_scores() -> pd.DataFrame:
    """The MC3 mutation heat -> the per-gene upper-tail p-value."""
    snv = WORK / "mc3_snvs.tsv"
    if not snv.exists():
        log("  the MC3 SNV file is missing; run the HotNet2 pipeline first")
        return pd.DataFrame()
    df = pd.read_csv(snv, sep="\t", header=None, names=["sample", "genes"],
                     dtype=str)
    counts: dict[str, int] = {}
    for genes in df.genes.dropna():
        for g in str(genes).split("\t"):
            counts[g] = counts.get(g, 0) + 1
    out = pd.DataFrame({"gene": list(counts), "heat": list(counts.values())})
    # NetMix2 takes p-values (its local-FDR step computes the z-scores
    # itself), and the local FDR needs a smooth null: the empirical rank
    # has far too many ties, so the mutation burden gets the upper tail of
    # a Poisson with the panel's mean burden.
    from scipy import stats as sps
    out = out[np.isfinite(out.heat.values.astype(float))
              & (out.heat.values.astype(float) >= 0)]
    lam = float(out.heat.mean())
    p = sps.poisson.sf(out.heat.values.astype(float) - 1, lam)
    out["p"] = np.clip(p, 1e-300, 1.0)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--networks", nargs="*", default=NETWORKS)
    args = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    scores = heat_scores()
    if not len(scores):
        return 1
    log(f"the gene scores: {len(scores):,} genes with a mutation heat")
    score_file = WORK / "gene_scores.tsv"
    scores[["gene", "p"]].to_csv(score_file, sep="\t", index=False,
                                  header=False)

    rows = []
    for net in args.networks:
        net_file = NETS_DIR / f"{net}.tsv"
        if not net_file.exists():
            continue
        outdir = WORK / "out" / net
        outdir.mkdir(parents=True, exist_ok=True)
        edge_file = WORK / f"{net}_edges.tsv"
        if not edge_file.exists():
            genes, A = load_network(str(net_file))
            with edge_file.open("w", encoding="utf-8") as fh:
                for line in net_file.read_text(encoding="utf-8").splitlines()[1:]:
                    fh.write(line + "\n")
        # the score file must cover every gene of the network: the
        # upstream turns a missing gene into a NaN z-score, which breaks
        # the local-FDR step
        genes_net, _ = load_network(str(net_file))
        zmap = dict(zip(scores.gene, scores.p))
        per_net = WORK / f"{net}_scores.tsv"
        # a gene the score file does not cover gets the maximum p-value,
        # which is the neutral value
        per_net.write_text("".join(
            f"{g}\t{zmap.get(g, 1.0):.6g}\n" for g in genes_net),
            encoding="utf-8")
        log(f"  running NetMix2 on {net}")
        t0 = time.time()
        r = subprocess.run(
            [str(PY), str(REPO / "run_netmix2.py"),
             "-el", str(edge_file), "-gs", str(per_net),
             "-o", str(outdir), "-v", "1"],
            capture_output=True, text=True, cwd=str(REPO),
            encoding="utf-8", errors="replace", env=ENV)
        out = (r.stdout or "") + (r.stderr or "")
        log(f"    {time.time() - t0:.0f}s, the return code {r.returncode}")
        if r.returncode != 0:
            for line in out.strip().splitlines()[-8:]:
                log("    " + line[:140])
            continue
        # the per-gene score: the genes in the selected subnetwork, with
        # their z-scores if the output has them
        best = read_netmix2_scores(outdir)
        if best is None:
            log("    the output could not be read")
            continue
        genes, A = load_network(str(net_file))
        deg = degree_scores(A)
        degmap = dict(zip(genes, deg))
        gset = set(genes)
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
            rows.append({"network": net, "label": label,
                         "n_positives": len(pos),
                         "n_scored": sum(1 for g in pos + neg
                                         if best.get(g, 0.0) > 0),
                         "netmix2_auroc": ev_h["auroc"],
                         "degree_auroc": ev_d["auroc"],
                         "delta_auroc": ev_h["auroc"] - ev_d["auroc"]})
            log(f"    {label}: the NetMix2 AUROC {ev_h['auroc']:.4f} vs the "
                f"degree {ev_d['auroc']:.4f}")
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "c37_netmix2.tsv", sep="\t", index=False)
    log(f"wrote the Control 37 NetMix2 table ({len(df)} rows)")
    return 0


def read_netmix2_scores(outdir: Path) -> dict | None:
    """The per-gene score from the NetMix2 output."""
    import json
    for name in ("netmix2_results.json", "results.json"):
        f = outdir / name
        if f.exists():
            data = json.loads(f.read_text(encoding="utf-8"))
            genes = data.get("genes") or data.get("vertices") or []
            scores = data.get("zscores") or data.get("scores") or []
            if genes and scores and len(genes) == len(scores):
                return {str(g): float(s) for g, s in zip(genes, scores)}
    # the fallback: the selected subnetwork's genes score 1, the rest 0
    for f in outdir.rglob("*.txt"):
        if f.name in ("components.txt", "subnetwork.txt"):
            genes = [l.strip() for l in f.read_text(encoding="utf-8",
                                                   errors="replace").split()
                     if l.strip()]
            if genes:
                return {g: 1.0 for g in genes}
    return None


if __name__ == "__main__":
    raise SystemExit(main())
