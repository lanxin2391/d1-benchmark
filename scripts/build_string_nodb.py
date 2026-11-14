"""Control 18 -- STRING without the `database` channel.

Our network panel is not six independent tests: Reactome is the gold
standard of FunMap's edges, Reactome is an input to STRING's `database`
channel, and Reactome-derived and Cancer-Gene-Census-derived evidence is
also what the labels are built from. So a large part of what makes
STRING useful in this task may be the *same literature knowledge* the
labels encode.

Control 18 quantifies that. The STRING v12 detailed file carries the
per-channel evidence, so the `database` channel can be removed and the
network rebuilt from the remaining channels.

**This is an approximation, and it is labelled as one.** STRING's
combined_score is a Bayesian integration of the seven channels; we do
not re-implement it. Instead the edge weight is the **maximum over the
non-database channels** (neighbourhood, fusion, co-occurrence,
co-expression, experimental, text-mining) and the edges are cut to the
same edge count as the primary `string_full700` arm, so the two
networks are density-matched. The question is qualitative -- does the
margin survive when the database channel is gone? -- not a claim about
STRING's own scoring function.

Outputs:
  data/processed/networks/string_nodb700.tsv (+ meta.json)
  results/tables/c18_nodb_channel.tsv
  docs/nodb_channel.md
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.networks.io import check_network_file, finalize_edges  # noqa: E402

RAW = ROOT / "data/raw/networks/string"
DETAIL = RAW / "9606.protein.links.detailed.v12.0.txt.gz"
INFO = RAW / "9606.protein.info.v12.0.txt.gz"
OUT = ROOT / "data/processed/networks"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/nodb_channel.md"

CHANNELS = ["neighborhood", "fusion", "cooccurence", "coexpression",
             "experimental", "textmining"]
TARGET_NET = "string_full700"
NET_ID = "string_nodb700"
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load():
    info = pd.read_csv(INFO, sep="\t")
    info["ensp"] = info["#string_protein_id"].str.replace("9606.", "",
                                                          regex=False)
    mapping = dict(zip(info["ensp"], info["preferred_name"]))
    log("reading the detailed links (this takes a minute)")
    df = pd.read_csv(DETAIL, sep=" ")
    log(f"  {len(df):,} raw edges")
    df["gene_a"] = df["protein1"].str.replace("9606.", "", regex=False) \
        .map(mapping)
    df["gene_b"] = df["protein2"].str.replace("9606.", "", regex=False) \
        .map(mapping)
    df = df.dropna(subset=["gene_a", "gene_b"])
    log(f"  {len(df):,} edges with both ends mapped to a symbol")
    df["nodb_weight"] = df[CHANNELS].max(axis=1)
    df = df[df["nodb_weight"] > 0]
    log(f"  {len(df):,} edges with a non-database channel")
    return df


def build_matched(df, target_edges: int):
    """Cut the non-database network to the same edge count as the primary
    STRING arm (a density-matched comparison)."""
    df = df.sort_values("nodb_weight", ascending=False)
    lo, hi = 1, len(df)
    best = None
    for _ in range(18):
        k = (lo + hi) // 2
        sub = df.head(k)
        out, meta = finalize_edges(
            sub[["gene_a", "gene_b", "nodb_weight"]].rename(
                columns={"nodb_weight": "weight"}),
            net_id=NET_ID,
            meta={"source_file": DETAIL.name,
                  "channels_used": ",".join(CHANNELS),
                  "weight_rule": "max over the non-database channels "
                                 "(an approximation of STRING's combined "
                                 "score, not STRING's own formula)",
                  "edge_count_matched_to": TARGET_NET},
            out_dir=OUT)
        n = len(out)
        if best is None or abs(n - target_edges) < abs(best[1] - target_edges):
            best = (out, n, meta, k)
        if n < target_edges:
            lo = k + 1
        else:
            hi = k - 1
        if hi <= lo:
            break
    out, n, meta, k = best
    out.to_csv(OUT / f"{NET_ID}.tsv", sep="\t", index=False)
    import json
    (OUT / f"{NET_ID}.meta.json").write_text(json.dumps(meta, indent=2),
                                            encoding="utf-8")
    check_network_file(str(OUT / f"{NET_ID}.tsv"))
    log(f"  {NET_ID}: {n} edges, {meta['n_nodes_lcc']} nodes in the LCC "
        f"(target {target_edges})")
    return n, meta


def main() -> int:
    df = load()
    target = len(pd.read_csv(OUT / f"{TARGET_NET}.tsv", sep="\t"))
    log(f"the primary {TARGET_NET} has {target} edges")
    n, meta = build_matched(df, target)

    # the splits, with the same cap as the primary design
    import subprocess
    py = r"E:\anaconda\python.exe"
    for lab in LABELS:
        r = subprocess.run(
            [py, "scripts/build_splits.py", "--label", lab,
             "--universe", f"native_{NET_ID}", "--apply-d11-cap",
             "--n-rep", "40", "--out-dir", "data/processed/splits_v40"],
            capture_output=True, text=True, cwd=str(ROOT))
        if r.returncode != 0:
            print((r.stdout or r.stderr)[-1500:])
            raise SystemExit("the splits failed")
    log("the splits are built")

    runs = ROOT / "results/runs_v40_nodb"
    scores = ROOT / "results/scores_v40_nodb"
    runs.mkdir(parents=True, exist_ok=True)
    scores.mkdir(parents=True, exist_ok=True)
    r = subprocess.run(
        [py, "scripts/run_a4.py", "--only-net", NET_ID,
         "--splits-dir", "data/processed/splits_v40",
         "--runs-dir", str(runs), "--scores-dir", str(scores),
         "--table-suffix", "_nodb_v40", "--n-rep", "40"],
        capture_output=True, text=True, cwd=str(ROOT))
    if r.returncode != 0:
        print((r.stdout or r.stderr)[-2000:])
        raise SystemExit("the arm failed")
    log("the four arms are run")

    # the comparison
    nodb = pd.read_csv(TABLES / "a4_table1_delta_auroc_nodb_v40.tsv", sep="\t")
    base = pd.read_csv(TABLES / "a4_table1_delta_auroc_v40.tsv", sep="\t")
    base = base[base.network == TARGET_NET]
    rows = []
    for r in nodb.to_dict("records"):
        b = base[base.label_set == r["label_set"]]
        rows.append({
            "network": TARGET_NET, "label": r["label_set"],
            "with_the_database_channel": float(b["a0.5_mean"].iloc[0]),
            "without_the_database_channel": float(r["a0.5_mean"]),
            "change": float(r["a0.5_mean"]) - float(b["a0.5_mean"].iloc[0]),
        })
    comp = pd.DataFrame(rows)
    comp.to_csv(TABLES / "c18_nodb_channel.tsv", sep="\t", index=False)
    write_doc(comp, n, meta)
    log("wrote the comparison table and the doc")
    return 0


def write_doc(comp: pd.DataFrame, n: int, meta: dict):
    lines = [
        "# Control 18 -- STRING without the `database` channel",
        "",
        "**Date:** 2026-11-12  ",
        "**Script:** `scripts/build_string_nodb.py`  ",
        "**Table:** `results/tables/c18_nodb_channel.tsv`",
        "",
        "## Why",
        "",
        "Reactome is the gold standard of FunMap's edges, Reactome feeds "
        "STRING's `database` channel, and Reactome- and CGC-derived "
        "evidence is also what the labels are built from. A large part of "
        "what makes STRING useful here may therefore be the same "
        "literature knowledge the labels encode.",
        "",
        "## How (and the approximation)",
        "",
        "The edge weight is the **maximum over the non-database channels** "
        "(neighbourhood, fusion, co-occurrence, co-expression, "
        "experimental, text-mining) and the edges are cut to the same edge "
        f"count as the primary arm ({n} edges against "
        f"{len(pd.read_csv(OUT / f'{TARGET_NET}.tsv', sep='\\t'))}), so the "
        "two networks are density-matched. **This is not STRING's own "
        "combined score** (a Bayesian integration of the seven channels) "
        "and is labelled as an approximation; the question it answers is "
        "qualitative -- does the margin survive when the database channel "
        "is gone?",
        "",
        "| label | with the database channel | without it | change |",
        "|---|---|---|---|",
    ]
    for r in comp.itertuples():
        lines.append(f"| {r.label} | {r.with_the_database_channel:+.4f} | "
                     f"{r.without_the_database_channel:+.4f} | "
                     f"{r.change:+.4f} |")
    lines += [
        "",
        "## Caveats",
        "",
        "1. The approximation: the maximum over channels is not STRING's "
        "Bayesian combination, so the absolute weights differ; the "
        "comparison is about the *presence* of the database channel, not "
        "about STRING's scoring.",
        "2. The node universe changes with the edge set, so a margin "
        "difference mixes the channel effect with the universe effect.",
        "3. Exploratory family (amendment 4).",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
