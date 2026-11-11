"""Direction 1 -- a consensus prior across networks, and its stability.

Why this model, and why now
----------------------------
Our own S1-S3 results say that the benchmark's one negative cell is a
COVERAGE failure, not a noise artefact: reactome's curated pathway
annotations lag behind the 2020-2024 driver discoveries, and the
failure survives three controls that equalise the positive set
(D-09/D-11/D-12). A single stale network cannot be repaired by a better
learner; it can only be out-voted.

So: build a CONSENSUS prior over the six networks on the shared
universe (1,605 genes), and ask three falsifiable questions.

  Q1  Does the consensus beat the degree baseline by at least the
      MEDIAN of its members? (If it only matches the median, the
      consensus is a robustness device, not a better method.)
  Q2  Is the consensus's per-gene ranking MORE STABLE than each
      member's? (This is the axis S1/S2 built, and the one that no
      other paper in this space has.)
  Q3  Does the consensus REPAIR the reactome temporal failure, and does
      dropping reactome from the consensus change that? (Leave-one-
      network-out consensus variants.)

Method
------
* Universe: the shared universe (D-12), so every network is scored on
  exactly the same 1,605 genes and the same 40-repeat x 5-fold split.
* Seeds: the positives of the training fold (D-06); the RWR always
  propagates on the network's own full LCC.
* Consensus: within one (repeat, fold) evaluation, each network's score
  is converted to a percentile rank over the shared genes, and the six
  ranks are averaged (rank aggregation is scale-free, so a network with
  a very different score scale cannot dominate).
* Reference: the degree baseline is the full-LCC degree of the shared
  genes, exactly as in D-12.
* Stability: the S1 metrics (score CV, rank SD) computed on the shared
  universe for the consensus and for every member.

Outputs:
  results/tables/consensus_table.tsv        (per label, per variant)
  results/tables/consensus_stability.tsv
  results/tables/consensus_member_table.tsv (per member, for reference)
  figures/figure7_consensus_prior.png
  docs/consensus_prior.md
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.evaluate import evaluate, run_arm  # noqa: E402

TABLES = ROOT / "results" / "tables"
FIG_DIR = ROOT / "figures"
SCORES_DIR = ROOT / "results" / "scores_shared"
NET_DIR = ROOT / "data" / "processed" / "networks"
SPLITS_DIR = ROOT / "data" / "processed" / "splits_v40"

NETWORKS = ["funmap", "rna_coexp", "string_phys700", "intact",
            "string_full700", "reactome"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
ALPHA = 0.5
SCORES_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def ensure_scores(net_id: str, label: str) -> Path:
    """Per-gene scores of one network on the shared universe."""
    dump = SCORES_DIR / f"{net_id}__{label}.npz"
    split = SPLITS_DIR / f"{label}_shared__shared.tsv"
    if dump.exists() and not split.exists():
        raise FileNotFoundError(split)
    if dump.exists():
        return dump
    log(f"  running {net_id} x {label} on the shared universe ...")
    run_arm(network=str(NET_DIR / f"{net_id}.tsv"), splits=str(split),
            network_id=net_id, label_id=f"{label}_shared",
            universe_id="native_" + net_id, alphas=(ALPHA,),
            dump_scores=str(dump))
    return dump


def load_shared(net_id: str, label: str, shared_genes: list[str]):
    """(scores per evaluation for the shared genes, degree, in_test)."""
    d = np.load(SCORES_DIR / f"{net_id}__{label}.npz", allow_pickle=True)
    genes = np.asarray(d["genes"]).astype(str)
    ai = int(np.argmin(np.abs(np.asarray(d["alphas"], dtype=float) - ALPHA)))
    rwr = np.asarray(d["rwr"][ai], dtype=float)          # (n_genes, n_combos)
    in_test = np.asarray(d["in_test"], dtype=bool)        # (n_combos, n_genes)
    degree = np.asarray(d["degree"], dtype=float)
    sel = {g: i for i, g in enumerate(genes)}
    idx = np.array([sel[g] for g in shared_genes])
    deg = degree[idx]
    # the test mask restricted to the shared genes
    it = in_test[:, idx]
    # ranks of the shared genes within each evaluation, over the genes
    # that are actually being scored in that evaluation
    # rwr is (n_genes, n_combos): gene indexing is axis 0
    scores = rwr[idx, :]
    ranks = np.empty_like(scores)
    for j in range(scores.shape[1]):
        m = it[j]
        if m.sum() == 0:
            ranks[m, j] = np.nan
            continue
        v = rwr[idx, j][m]
        order = np.argsort(np.argsort(v))
        ranks[m, j] = order / max(len(v) - 1, 1)
    return ranks, deg, it, scores


def auroc(scores: np.ndarray, y: np.ndarray) -> float:
    return float(evaluate(scores, y, np.array([str(i) for i in range(len(y))]),
                           ks=(10, 50))["auroc"])


def stability(scores: np.ndarray, it: np.ndarray, repeats: np.ndarray):
    """S1 metrics on the shared universe: one score per repeat per gene."""
    n_genes, n_eval = scores.shape
    n_rep = int(repeats.max()) + 1
    per_rep = np.full((n_genes, n_rep), np.nan)
    for r in range(n_rep):
        js = np.where(repeats == r)[0]
        # the genes scored in evaluation j are the test fold of j; a gene
        # is in the test set of exactly one fold per repeat, so the
        # per-repeat score is taken from that evaluation
        for j in js:
            m = it[j]
            per_rep[m, r] = scores[m, j]
    mean = per_rep.mean(axis=1)
    sd = per_rep.std(axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = np.where(np.abs(mean) > 1e-12, sd / np.abs(mean), np.nan)
    order = np.argsort(np.argsort(np.nan_to_num(per_rep, nan=-np.inf), axis=0),
                       axis=0)
    ranks = order / max(n_genes - 1, 1)
    rank_sd = np.nanstd(ranks, axis=1, ddof=1)
    ok = np.isfinite(mean)
    return {
        "score_cv_median": float(np.nanmedian(cv[ok])),
        "rank_sd_median": float(np.nanmedian(rank_sd[ok])),
        "score_sd_median": float(np.nanmedian(sd[ok])),
    }


def repeat_ci(rep: pd.Series) -> tuple[float, float]:
    """95 % t interval of the repeat-level means (the primary unit)."""
    from scipy import stats as sps
    n = len(rep)
    if n < 3:
        return (np.nan, np.nan)
    se = rep.std(ddof=1) / np.sqrt(n)
    h = sps.t.ppf(0.975, n - 1) * se
    return (float(rep.mean() - h), float(rep.mean() + h))


def main() -> int:
    # ---- the shared universe ------------------------------------------
    shared = [l.strip() for l in
              open(NET_DIR / "universe_shared.txt", encoding="utf-8")
              if l.strip()]
    log(f"shared universe: {len(shared):,} genes")

    # ---- per-gene scores of every network on the shared universe ------
    for net in NETWORKS:
        for label in LABELS:
            ensure_scores(net, label)

    # ---- per label: members, consensus variants, evaluation -----------
    member_rows, cons_rows, stab_rows = [], [], []
    for label in LABELS:
        log(f"label {label}")
        ranks, deg, it, raw = {}, {}, {}, {}
        for net in NETWORKS:
            r, dg, i, s = load_shared(net, label, shared)
            ranks[net], deg[net], it[net], raw[net] = r, dg, i, s
        # the split is shared, so all networks share the test mask
        mask = it[NETWORKS[0]]
        repeats = np.asarray(np.load(SCORES_DIR / f"{NETWORKS[0]}__{label}.npz",
                                     allow_pickle=True)["repeat_idx"]).astype(int)

        # label vector per evaluation (from the split file)
        sp = pd.read_csv(SPLITS_DIR / f"{label}_shared__shared.tsv", sep="\t")
        y_by_eval = {}
        for (r, f), g in sp.groupby(["repeat", "fold"])["y"]:
            pass
        y_map = dict(zip(sp["gene"], sp["y"].astype(int)))
        y_shared = np.array([y_map.get(g, 0) for g in shared])

        variants = {"consensus_all6": NETWORKS}
        for net in NETWORKS:
            variants[f"consensus_minus_{net}"] = [n for n in NETWORKS if n != net]

        # ONE degree baseline for every prior on this universe: the mean
        # degree over the six networks. Comparing the consensus against
        # funmap's degree while the members are compared against their
        # own degree would inflate the consensus by the difference
        # between the degree baselines -- so the common baseline is used
        # for the fair comparison, and each network's own degree is kept
        # as the D-12-comparable number.
        deg_common = np.nanmean(np.stack([deg[n] for n in NETWORKS]), axis=0)

        for vname, members in variants.items():
            r_mean = np.nanmean(np.stack([ranks[m] for m in members]), axis=0)
            deltas, deltas_own, per_rep = [], [], []
            for j in range(r_mean.shape[1]):
                m = mask[j]
                if m.sum() == 0 or y_shared[m].sum() == 0:
                    continue
                a_c = auroc(r_mean[m, j], y_shared[m])
                a_d = auroc(deg_common[m], y_shared[m])
                deltas.append(a_c - a_d)
                # against the mean of the members' own degree baselines
                a_do = auroc(np.nanmean(np.stack([deg[x][m] for x in members]),
                                        axis=0), y_shared[m])
                deltas_own.append(a_c - a_do)
                per_rep.append((int(repeats[j]), a_c - a_d))
            st = stability(r_mean, mask, repeats)
            rep = pd.DataFrame(per_rep, columns=["repeat", "d"]).groupby(
                "repeat")["d"].mean()
            ci = repeat_ci(rep)
            cons_rows.append({
                "label": label, "variant": vname,
                "n_members": len(members),
                "members": ";".join(members),
                "delta_auroc": float(np.mean(deltas)),
                "delta_auroc_vs_members_own_degree": float(np.mean(deltas_own)),
                "delta_sd_across_folds": float(np.std(deltas, ddof=1)),
                "n_eval": len(deltas),
                "repeat_mean_delta": float(rep.mean()),
                "repeat_sd_delta": float(rep.std(ddof=1)),
                "repeat_ci95_low": ci[0], "repeat_ci95_high": ci[1],
                "n_repeats": int(len(rep)),
                **st,
            })

        for net in NETWORKS:
            deltas, deltas_own, per_rep = [], [], []
            for j in range(raw[net].shape[1]):
                m = mask[j]
                if m.sum() == 0 or y_shared[m].sum() == 0:
                    continue
                a_c = auroc(raw[net][m, j], y_shared[m])
                a_d = auroc(deg_common[m], y_shared[m])
                deltas.append(a_c - a_d)
                deltas_own.append(a_c - auroc(deg[net][m], y_shared[m]))
                per_rep.append((int(repeats[j]), a_c - a_d))
            st = stability(raw[net], mask, repeats)
            rep = pd.DataFrame(per_rep, columns=["repeat", "d"]).groupby(
                "repeat")["d"].mean()
            ci = repeat_ci(rep)
            member_rows.append({
                "label": label, "network": net,
                "delta_auroc": float(np.mean(deltas)),
                "delta_auroc_vs_own_degree": float(np.mean(deltas_own)),
                "delta_sd_across_folds": float(np.std(deltas, ddof=1)),
                "repeat_mean_delta": float(rep.mean()),
                "repeat_sd_delta": float(rep.std(ddof=1)),
                "repeat_ci95_low": ci[0], "repeat_ci95_high": ci[1],
                **st,
            })
            stab_rows.append({"label": label, "entity": net, "kind": "member",
                              **st})

    cons = pd.DataFrame(cons_rows)
    members = pd.DataFrame(member_rows)
    cons.to_csv(TABLES / "consensus_table.tsv", sep="\t", index=False, lineterminator="\n")
    members.to_csv(TABLES / "consensus_member_table.tsv", sep="\t",
                   index=False, lineterminator="\n")

    # ---- Q1: consensus vs the members --------------------------------
    log("\n=== Q1/Q2/Q3 ===")
    for label in LABELS:
        m = members[members.label == label]
        c = cons[(cons.label == label) & (cons.variant == "consensus_all6")]
        cm = cons[(cons.label == label) &
                  (cons.variant == "consensus_minus_reactome")]
        log(f"{label}: members median {m.delta_auroc.median():+.4f} "
            f"(best {m.delta_auroc.max():+.4f}, worst {m.delta_auroc.min():+.4f})"
            f" | consensus_all6 {c.delta_auroc.iloc[0]:+.4f} "
            f"(rank SD {c.rank_sd_median.iloc[0]:.4f})"
            f" | minus_reactome {cm.delta_auroc.iloc[0]:+.4f}"
            f" | members' own-degree median "
            f"{m.delta_auroc_vs_own_degree.median():+.4f}")
        log(f"    member rank SD: "
            + ", ".join(f"{r.network}={r.rank_sd_median:.4f}"
                        for r in m.itertuples()))

    make_figure(cons, members)
    log(f"wrote consensus_table.tsv ({len(cons)} rows) and "
        f"consensus_member_table.tsv ({len(members)} rows)")
    return 0


def make_figure(cons: pd.DataFrame, members: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    labels = list(dict.fromkeys(members["label"]))
    x = np.arange(len(labels))

    ax = axes[0]
    for net, g in members.groupby("network"):
        gg = g.set_index("label").reindex(labels)
        ax.plot(x, gg["delta_auroc"], "o--", alpha=0.45, lw=0.9,
                color="grey", ms=4)
    c = cons[cons.variant == "consensus_all6"].set_index("label").reindex(labels)
    ax.plot(x, c["delta_auroc"], "s-", color="#C44E52", lw=2.2, ms=7,
            label="consensus (6 networks)")
    ax.axhline(0, color="black", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8.5)
    ax.set_ylabel("delta-AUROC on the shared universe")
    ax.set_title("A  Consensus vs its members (grey = single networks)",
                 fontsize=11)
    ax.legend(fontsize=8.5)

    ax = axes[1]
    for net, g in members.groupby("network"):
        gg = g.set_index("label").reindex(labels)
        ax.plot(x, gg["rank_sd_median"], "o--", alpha=0.45, lw=0.9,
                color="grey", ms=4)
    ax.plot(x, c["rank_sd_median"], "s-", color="#C44E52", lw=2.2, ms=7,
            label="consensus (6 networks)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8.5)
    ax.set_ylabel("median rank SD (lower = more stable)")
    ax.set_title("B  Is the consensus more stable than its members?",
                 fontsize=11)
    ax.legend(fontsize=8.5)

    ax = axes[2]
    sub = cons[cons.label == "intogen_temporal_new"].copy()
    sub = sub[sub.variant.str.startswith("consensus_minus_") |
              (sub.variant == "consensus_all6")]
    sub["label_short"] = sub["variant"].str.replace("consensus_minus_",
                                                    "-")
    sub = sub.sort_values("delta_auroc")
    ax.barh(range(len(sub)), sub["delta_auroc"], color="#4C72B0")
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(sub["label_short"], fontsize=8)
    ax.axvline(0, color="black", lw=1)
    ax.set_xlabel("delta-AUROC on intogen_temporal_new")
    ax.set_title("C  Leave-one-network-out consensus (temporal label)",
                 fontsize=11)

    fig.suptitle("Figure 7  A consensus prior over six networks "
                 "(shared universe, 1,605 genes)", fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG_DIR / "figure7_consensus_prior.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)
    log(f"wrote figure7_consensus_prior.png")


if __name__ == "__main__":
    raise SystemExit(main())
