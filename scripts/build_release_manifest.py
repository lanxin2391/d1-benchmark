"""Write the RELEASE manifest: the SHA-256 of the artefacts the release contains.

This is the release's own tamper-evidence record. The pre-registration
manifest (``docs/pre_registration_manifest.json``) is FROZEN -- it records
the state that was filed with Zenodo, and the fingerprint mismatches the
audit reports are the expected consequence of the repairs that recomputed
those tables.

FIX (2026-11-11). The previous manifest (``docs/release_manifest_v1.0-rc2.json``)
recorded the commit ``3192cc5``, which **is not in this repository's history**
(the history was rewritten after the history-rewrite commit 4209064), and
its SHA-256 fingerprints matched only 1 of the 20 tables it listed. The
manifest was therefore not verifiable: it claimed a commit that does not
exist and hashes that did not describe the files in the tree.

This script now:

* takes the release name and the commit as arguments (no hard-coded tag),
* lists every table and figure that exists, recording the SHA-256 of the
  file **as it is in the working tree**,
* records which files were **missing** rather than silently skipping them,
* and writes a manifest that a reader can verify with one command.

Usage
-----
    python scripts/build_release_manifest.py --release v1.0-rc3
    python scripts/build_release_manifest.py --release v1.0-rc3 --commit HEAD
    python scripts/build_release_manifest.py --release v1.0-rc3 --check
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.io import sha256_normalised  # noqa: E402

KEY_TABLES = [
    # the primary endpoint and its statistics
    "results/tables/a4_table1_delta_auroc_v40.tsv",
    "results/tables/a4_table1_rwr_auroc_v40.tsv",
    "results/tables/a4_table1_degree_auroc_v40.tsv",
    "results/tables/results_stats.tsv",
    "results/tables/stats_bootstrap_ci.tsv",
    "results/tables/stats_equivalence.tsv",
    "results/tables/stats_h1_verdict_v40.tsv",
    "results/tables/stats_variance_decomp.tsv",
    # the nulls
    "results/tables/a5_p_values.tsv",
    "results/tables/w5_p_values.tsv",
    # the controls
    "results/tables/d09_table.tsv",
    "results/tables/d11_table.tsv",
    "results/tables/d12_table.tsv",
    "results/tables/d3_split_invariants.tsv",
    "results/tables/d4_v1_vs_v2_comparison.tsv",
    "results/tables/control7_matched_null.tsv",
    "results/tables/pubcount_adj_comparison.tsv",
    "results/tables/control8_seed_degree_summary.tsv",
    "results/tables/control9_size_adjustment.tsv",
    "results/tables/control9_size_model.tsv",
    "results/tables/control10_random_input_summary.tsv",
    # the stability and consensus work
    "results/tables/hub_capped_rwr.tsv",
    "results/tables/stability_by_arm.tsv",
    "results/tables/stability_attribution.tsv",
    "results/tables/network_stability.tsv",
    "results/tables/consensus_table.tsv",
    "results/tables/consensus_robustness.tsv",
    "results/tables/rewiring_mixing_check.tsv",
    # the label panel and the transfer arms
    "results/tables/label_network_overlap.tsv",
    "results/tables/results_stats_transfer.tsv",
    "results/tables/provenance_summary.tsv",
    "results/tables/gnn_courtesy_arm.tsv",
    # the figure data
    "results/tables/figure1_data.tsv",
]

FIGURES = [
    "figures/figure1_pilot_overview.png",
    "figures/figure2_d09_d11_d12.png",
    "figures/figure2_d09_temporal_strip.png",
    "figures/figure5_stability_atlas.png",
    "figures/figure5b_stability_attribution.png",
    "figures/figure6_network_stability.png",
    "figures/figure7_consensus_prior.png",
    "figures/figure8_stability_tradeoff.png",
]


def sha256(path: Path) -> str:
    """SHA-256 of the file with CRLF normalised to LF.

    The repo's `.gitattributes` normalises every `*.tsv` to LF, so this is
    the hash a reader computes from a checkout -- which is the hash a
    tamper-evidence manifest has to record. Hashing the raw working-tree
    bytes instead would fail on Windows, where a regenerated table can
    still carry CRLF endings.
    """
    return sha256_normalised(path)


def git(*args) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=str(ROOT),
                                       text=True, stderr=subprocess.DEVNULL).strip()
    except subprocess.CalledProcessError:
        return ""


def build(release: str, commit: str) -> dict:
    tables: dict[str, str] = {}
    figures: dict[str, str] = {}
    missing: list[str] = []
    for rel in KEY_TABLES:
        p = ROOT / rel
        if p.is_file():
            tables[rel] = sha256(p)
        else:
            missing.append(rel)
    for rel in FIGURES:
        p = ROOT / rel
        if p.is_file():
            figures[rel] = sha256(p)
        else:
            missing.append(rel)
    return {
        "release": release,
        "commit": commit,
        "built_by": "scripts/build_release_manifest.py",
        "note": ("The pre-registration manifest "
                 "(docs/pre_registration_manifest.json) is frozen and "
                 "records the state filed with Zenodo; the fingerprint "
                 "mismatches the audit reports are the expected "
                 "consequence of the repairs that recomputed those "
                 "tables. This file is the release's own record. The "
                 "hashes below describe the files in the working tree "
                 "at the commit named above; a reader can check them "
                 "with scripts/check_release_manifest.py."),
        "n_tables": len(tables),
        "n_figures": len(figures),
        "missing_at_build_time": missing,
        "sha256_tables": tables,
        "sha256_figures": figures,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--release", default="v1.0-rc3")
    ap.add_argument("--commit", default="HEAD")
    ap.add_argument("--out", default=None)
    ap.add_argument("--check", action="store_true",
                    help="verify the existing manifest against the tree")
    args = ap.parse_args()

    commit = git("rev-parse", args.commit) or args.commit
    out = Path(args.out) if args.out else \
        ROOT / "docs" / f"release_manifest_{args.release}.json"

    if args.check:
        if not out.is_file():
            print(f"{out} does not exist")
            return 1
        manifest = json.loads(out.read_text(encoding="utf-8"))
        bad = []
        for rel, want in manifest.get("sha256_tables", {}).items():
            p = ROOT / rel
            if not p.is_file():
                bad.append(f"{rel}: MISSING")
            elif sha256(p) != want:
                bad.append(f"{rel}: hash differs")
        for rel, want in manifest.get("sha256_figures", {}).items():
            p = ROOT / rel
            if not p.is_file():
                bad.append(f"{rel}: MISSING")
            elif sha256(p) != want:
                bad.append(f"{rel}: hash differs")
        print(f"checked {out.name}: "
              f"{len(manifest.get('sha256_tables', {}))} tables, "
              f"{len(manifest.get('sha256_figures', {}))} figures")
        if bad:
            print("FAIL")
            for b in bad:
                print(f"  {b}")
            return 1
        print("PASS -- every fingerprint matches the working tree")
        return 0

    manifest = build(args.release, commit)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, sort_keys=False) + "\n",
                   encoding="utf-8")
    print(f"wrote {out.relative_to(ROOT)}")
    print(f"  tables: {manifest['n_tables']}, "
          f"figures: {manifest['n_figures']}")
    if manifest["missing_at_build_time"]:
        print("  MISSING at build time:")
        for m in manifest["missing_at_build_time"]:
            print(f"    {m}")
    print(f"  commit: {commit}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
