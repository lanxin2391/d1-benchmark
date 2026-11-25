"""Control 32 -- the pathway level.

Everything so far is at the gene level. The teacher's question about the
pathway level ("are the drivers' pathways the ones the prior points at?")
needs a pathway-level analysis, and it needs three things the repository
does not have yet: a pathway-to-gene mapping, a pathway-level statistic,
and a pathway-level null.

The mapping comes from the raw Reactome release we already downloaded:
the "Interaction context" column of the interaction table is the pathway
identifier, and the Ensembl gene ids are mapped to symbols through the
HGNC set we already use. The pathways are then restricted to the genes of
the Reactome network, so the universe is the network's own gene set.

Three analyses:

  1. Which pathways are the drivers in? The hypergeometric enrichment of
     each pathway for the label's positives, with a Benjamini-Hochberg
     correction over the pathways.
  2. Do the prior's picks concentrate on the driver's pathways more than
     the degree baseline's picks? The number of pathways significantly
     enriched among the top-100 picks of each scorer.
  3. The pathway-level permutation test: the same statistic on random
     sets of 100 genes from the same pool, which is the null the
     enrichment statistic has to beat.

Outputs:
  results/tables/c32_pathway_enrichment_drivers.tsv
  results/tables/c32_pathway_picks.tsv
  results/tables/c32_pathway_permutation.tsv
  docs/c32_pathway_level.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.baselines import degree_scores  # noqa: E402
from d1.engine.evaluate import read_splits  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

RAW = (ROOT / "data/raw/networks/reactome/"
       "reactome.homo_sapiens.interactions.tab-delimited.txt")
HGNC = ROOT / "data/raw/reference/hgnc_complete_set.txt"
NETS = ROOT / "data/processed/networks"
SPLITS = ROOT / "data/processed/splits_v40"
SCORES = ROOT / "results/scores_v40"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/c32_pathway_level.md"
NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
N_FOLDS = 20
K = 100
N_PERM = 20
ALPHA = 0.5
SEED = 20261125


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def ensembl_to_symbol() -> dict[str, str]:
    """The Ensembl gene id -> the HGNC symbol."""
    out = {}
    with HGNC.open("r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_ens = next((i for i, h in enumerate(header)
                      if h.lower() == "ensembl_gene_id"), None)
        i_sym = next((i for i, h in enumerate(header)
                      if h.lower() == "symbol"), None)
        if i_ens is None or i_sym is None:
            return out
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_ens, i_sym):
                continue
            e = parts[i_ens].strip()
            if e.startswith("ENSG"):
                out.setdefault(e.split(".")[0], parts[i_sym].strip())
    return out


def pathways_for_genes(symbols: set[str]) -> dict[str, list[str]]:
    """The pathway id -> the genes (restricted to the symbols we have)."""
    mapping = ensembl_to_symbol()
    out: dict[str, set[str]] = {}
    with RAW.open("r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_e1 = next((i for i, h in enumerate(header)
                     if h.startswith("Interactor 1 Ensembl")), None)
        i_e2 = next((i for i, h in enumerate(header)
                     if h.startswith("Interactor 2 Ensembl")), None)
        i_ctx = next((i for i, h in enumerate(header)
                      if "Interaction context" in h), None)
        if i_e1 is None or i_e2 is None or i_ctx is None:
            raise SystemExit("the Reactome file does not have the expected "
                             "columns")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_e1, i_e2, i_ctx):
                continue
            ctx = parts[i_ctx].strip()
            if not ctx:
                continue
            for i_e in (i_e1, i_e2):
                for e in parts[i_e].split("|"):
                    e = e.strip()
                    # the Reactome file prefixes the ids with "ENSEMBL:"
                    if e.upper().startswith("ENSEMBL:"):
                        e = e.split(":", 1)[1]
                    if not e.startswith("ENSG"):
                        continue
                    sym = mapping.get(e.split(".")[0])
                    if sym in symbols:
                        out.setdefault(ctx, set()).add(sym)
    return {k: sorted(v) for k, v in out.items()}


def bh(pvals: np.ndarray) -> np.ndarray:
    """The Benjamini-Hochberg adjusted p-values."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def main() -> int:
    log("reading the Reactome pathway interactions")
    genes_r, A = load_network(str(NETS / "reactome.tsv"))
    universe = set(genes_r)
    pathways = pathways_for_genes(universe)
    log(f"  {len(pathways)} pathways over {len(universe)} network genes")
    keep = {p: g for p, g in pathways.items() if len(g) >= 5}
    log(f"  {len(keep)} pathways with at least 5 genes")
    # the pathway x gene incidence matrix (sparse): one matrix-vector
    # product per pick set replaces 4,000 hypergeometric calls
    from scipy.sparse import csr_matrix
    gene_ix = {g: i for i, g in enumerate(sorted(universe))}
    path_list = sorted(keep)
    path_ix = {p: i for i, p in enumerate(path_list)}
    r, c = [], []
    for p, gl in keep.items():
        for g in gl:
            r.append(path_ix[p])
            c.append(gene_ix[g])
    incidence = csr_matrix((np.ones(len(r)), (r, c)),
                           shape=(len(path_list), len(gene_ix)))
    path_sizes = np.array([len(keep[p]) for p in path_list], dtype=float)
    log(f"  the incidence matrix: {incidence.shape}")

    # ---- 1. which pathways are the drivers in? -----------------------
    drv_rows = []
    gene2path: dict[str, list[str]] = {}
    for p, gl in keep.items():
        for g in gl:
            gene2path.setdefault(g, []).append(p)
    for label in LABELS:
        lab_path = ROOT / f"data/processed/labels/{label}.tsv"
        if not lab_path.exists():
            continue
        pos = set(pd.read_csv(lab_path, sep="\t")["gene"]) & universe
        if len(pos) < 20:
            continue
        ps = []
        rows = []
        vec = np.zeros(len(gene_ix), dtype=float)
        for g in pos:
            if g in gene_ix:
                vec[gene_ix[g]] = 1.0
        n_in_all = incidence @ vec
        pv_all = sps.hypergeom.sf(n_in_all - 1, len(gene_ix), path_sizes,
                                  len(pos))
        for p, gl in keep.items():
            n_in = int(n_in_all[path_ix[p]])
            n_pool = len(gl)
            pv = float(pv_all[path_ix[p]])
            ps.append(pv)
            rows.append({"label": label, "pathway": p, "n_genes": n_pool,
                         "n_positive": n_in, "p": float(pv)})
        q = bh(np.array(ps))
        for r, qq in zip(rows, q):
            r["q"] = float(qq)
        drv_rows.extend(rows)
        log(f"  {label}: {int((q < 0.05).sum())} of {len(q)} pathways "
            "enriched at q < 0.05")
    drv = pd.DataFrame(drv_rows)
    drv.to_csv(TABLES / "c32_pathway_enrichment_drivers.tsv", sep="\t",
               index=False)

    # the pathways each label's drivers are enriched in (q < 0.05)
    driver_enriched: dict[str, set] = {}
    for label, g in drv.groupby("label"):
        driver_enriched[label] = set(g.loc[g.q < 0.05, "pathway"])
    log("  the driver-enriched pathway counts: "
        + str({k: len(v) for k, v in driver_enriched.items()}))

    # ---- 2. and 3. the picks and the pathway-level null ---------------
    pick_rows = []
    perm_rows = []
    for net in NETWORKS:
        net_file = NETS / f"{net}.tsv"
        if not net_file.exists():
            continue
        genes, An = load_network(str(net_file))
        deg = degree_scores(An)
        index = pd.Series(np.arange(len(genes)), index=genes)
        for label in LABELS:
            sp_path = SPLITS / f"{label}__native_{net}.tsv"
            npz_path = SCORES / f"{net}__{label}.npz"
            if not sp_path.exists() or not npz_path.exists():
                continue
            d = np.load(npz_path, allow_pickle=True)
            ai = int(np.argmin(np.abs(np.asarray(d["alphas"], dtype=float)
                                      - ALPHA)))
            rwr_all = np.asarray(d["rwr"][ai], dtype=float)
            sp = read_splits(sp_path)
            sp = sp.assign(idx=index.loc[sp.gene].to_numpy())
            for (r, f) in sorted(sp.groupby(["repeat", "fold"]).groups):
                if int(r) >= N_FOLDS:
                    continue
                j = int(r) * 5 + int(f)
                if j >= rwr_all.shape[1]:
                    continue
                test = sp[(sp.repeat == r) & (sp.fold == f)]
                ti = test.idx.to_numpy()
                names = test.gene.to_numpy()
                if len(ti) < K:
                    continue
                pool = [i for i, g in enumerate(names) if g in universe]
                if len(pool) < K:
                    continue
                picks = {
                    "the prior": np.argsort(-rwr_all[ti, j], kind="mergesort"
                                            )[:K],
                    "degree": np.argsort(-deg[ti], kind="mergesort")[:K],
                }
                rng = np.random.default_rng(SEED + 13 * int(r) + int(f))
                for p in range(N_PERM):
                    picks[f"random {p}"] = rng.permutation(len(pool))[:K]
                for name, idxs in picks.items():
                    sel = {names[i] for i in idxs if names[i] in universe}
                    n_sel = len(sel)
                    if n_sel == 0:
                        continue
                    vec = np.zeros(len(gene_ix), dtype=float)
                    for g in sel:
                        vec[gene_ix[g]] = 1.0
                    n_in = incidence @ vec
                    pv = sps.hypergeom.sf(n_in - 1, len(gene_ix),
                                          path_sizes, n_sel)
                    n_sig = int((pv < 0.05).sum())
                    # do the picks point at the pathways the drivers are
                    # enriched in? (the driver-enriched set for this label)
                    drv_set = driver_enriched.get(label, set())
                    hit = [path_list[i] for i in np.where(pv < 0.05)[0]
                           if path_list[i] in drv_set]
                    n_drv_hit = len(hit)
                    row = {"network": net, "label": label, "repeat": int(r),
                           "fold": int(f), "scorer": name,
                           "n_enriched_pathways": n_sig,
                           "n_driver_pathways_hit": n_drv_hit,
                           "n_driver_pathways": len(drv_set)}
                    if name.startswith("random"):
                        perm_rows.append(row)
                    else:
                        pick_rows.append(row)
            log(f"  {net} x {label}: done")
    picks = pd.DataFrame(pick_rows)
    picks.to_csv(TABLES / "c32_pathway_picks.tsv", sep="\t", index=False)
    perm = pd.DataFrame(perm_rows)
    perm.to_csv(TABLES / "c32_pathway_permutation.tsv", sep="\t", index=False)
    write_doc(drv, picks, perm)
    log("wrote the Control 32 tables")
    return 0


def write_doc(drv: pd.DataFrame, picks: pd.DataFrame,
              perm: pd.DataFrame) -> None:
    lines = [
        "# Control 32 — the pathway level",
        "",
        "**Date:** 2026-11-25  ",
        "**Script:** `scripts/run_c32_pathway_level.py`  ",
        "**Tables:** `c32_pathway_enrichment_drivers.tsv`, "
        "`c32_pathway_picks.tsv`, `c32_pathway_permutation.tsv`",
        "",
        "## The setup",
        "",
        "The pathway-to-gene mapping comes from the raw Reactome release "
        "(the `Interaction context` column is the pathway), with the "
        "Ensembl ids mapped to symbols through HGNC, and it is restricted to "
        "the genes of the Reactome network so that the universe is the "
        "network's own gene set.",
        "",
        "## 1. Which pathways are the drivers in?",
        "",
        "| label | the pathways tested | the pathways enriched at q < 0.05 |",
        "|---|---|---|",
    ]
    for label, g in drv.groupby("label"):
        lines.append(f"| {label} | {g.pathway.nunique()} | "
                     f"{int((g.q < 0.05).sum())} |")
    lines += [
        "",
        "## 2. Do the prior's picks concentrate on those pathways more than "
        "the degree baseline's?",
        "",
        "| network | label | the prior: the enriched pathways | degree: the "
        "enriched pathways | the random baseline |",
        "|---|---|---|---|---|",
    ]
    a = picks.groupby(["network", "label", "scorer"]).n_enriched_pathways.mean()
    for net in NETWORKS:
        for label in LABELS:
            if (net, label, "the prior") not in a.index:
                continue
            pr = a.loc[(net, label, "the prior")]
            dg = a.loc[(net, label, "degree")]
            rnd = perm[(perm.network == net) & (perm.label == label)]
            rb = rnd.n_enriched_pathways.mean() if len(rnd) else np.nan
            lines.append(f"| {net} | {label} | {pr:.1f} | {dg:.1f} | "
                         f"{rb:.1f} |")
    lines += [
        "",
        "## Reading",
        "",
        "## Caveats",
        "",
        "1. The pathways come from Reactome, so the Reactome arm has a "
        "shared source with the pathway annotation; the other five "
        "networks do not.",
        "2. The pathway-level null permutes the picks among the test "
        "genes, so it tests the enrichment statistic, not the pathway "
        "annotation itself.",
        "3. **Exploratory** (amendment 4): no pre-registered endpoint is "
        "touched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
