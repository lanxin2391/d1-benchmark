"""Item 17 (part 2) -- run the protocol on the newly collected networks.

The new networks get their own split directory (the released panels are
never touched) and their own runs, and the result is the per-class
consistency of the margin: the H2 question is whether the networks of one
provenance class behave alike, which needs at least four per class.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(r"D:\Grade3\swxxx\final\D1\d1-benchmark")
sys.path.insert(0, str(ROOT))

from d1.splits import make_splits, write_splits  # noqa: E402

NETS_DIR = ROOT / "data/processed/networks"
LABELS_DIR = ROOT / "data/processed/labels"
SPLIT_DIR = ROOT / "data/processed/splits_extra"
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
NEW_NETS = ["string_full400", "string_full900", "rna_coexp_d150k",
            "rna_coexp_d250k", "biogrid_curated"]
N_REP = 10
K = 5
BASE_SEED = 20260926


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def label_genes(label_id: str) -> list[str]:
    p = LABELS_DIR / f"{label_id}.tsv"
    if not p.exists():
        return []
    import pandas as pd
    return list(pd.read_csv(p, sep="\t", dtype=str)["gene"])


def network_genes(net_id: str) -> list[str]:
    p = NETS_DIR / f"{net_id}.tsv"
    genes = []
    with p.open(encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            a, b = line.rstrip("\n").split("\t")[:2]
            genes.append(a)
            genes.append(b)
    return sorted(set(genes))


def build_splits() -> None:
    SPLIT_DIR.mkdir(parents=True, exist_ok=True)
    for net in NEW_NETS:
        universe = network_genes(net)
        for label in LABELS:
            pos = [g for g in label_genes(label) if g in set(universe)]
            out = SPLIT_DIR / f"{label}__native_{net}.tsv"
            if out.exists():
                continue
            if len(pos) < 20:
                log(f"  {label} x {net}: only {len(pos)} positives, skipped")
                continue
            df = make_splits(universe, pos, base_seed=BASE_SEED, n_rep=N_REP,
                             k=K)
            write_splits(df, label, f"native_{net}", out_dir=SPLIT_DIR)
            log(f"  {label} x {net}: {len(universe):,} genes, {len(pos)} "
                "positives")


def main() -> int:
    log("building the splits for the new networks")
    build_splits()
    log("running the protocol (10 repeats, alpha = 0.5)")
    import subprocess
    nets = [n for n in NEW_NETS
            if (SPLIT_DIR / f"intogen2024__native_{n}.tsv").exists()]
    if not nets:
        log("no network is ready")
        return 1
    cmd = [sys.executable, "scripts/run_a4.py",
           "--alpha", "0.5", "--n-rep", str(N_REP),
           "--splits-dir", str(SPLIT_DIR),
           "--runs-dir", str(ROOT / "results/runs_extra"),
           "--table-suffix", "_extra",
           "--scores-dir", str(ROOT / "results/scores_extra")]
    for n in nets:
        cmd += ["--only-net", n]
    for l in LABELS:
        cmd += ["--only-label", l]
    r = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    out = (r.stdout or "") + (r.stderr or "")
    for line in out.splitlines()[-25:]:
        log(line)
    return 0 if r.returncode == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
