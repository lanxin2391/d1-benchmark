"""Write the RELEASE manifest: the SHA-256 of the artefacts that the
v1.0-rc1 release contains, plus the release zips themselves.

The pre-registration manifest (docs/pre_registration_manifest.json) is
FROZEN -- it records the state that was filed with Zenodo, and the
audit's five fingerprint mismatches are the expected consequence of the
P0 repairs (the tables were recomputed). This file is the release's own
tamper-evidence record.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "release_manifest_v1.0rc1.json"
RELEASE = ROOT / "release"

KEY_TABLES = [
    "results/tables/a4_table1_delta_auroc_v40.tsv",
    "results/tables/a4_table1_rwr_auroc_v40.tsv",
    "results/tables/a4_table1_degree_auroc_v40.tsv",
    "results/tables/stats_bootstrap_ci.tsv",
    "results/tables/stats_equivalence.tsv",
    "results/tables/stats_h1_verdict_v40.tsv",
    "results/tables/stats_variance_decomp.tsv",
    "results/tables/results_stats.tsv",
    "results/tables/a5_p_values.tsv",
    "results/tables/w5_p_values.tsv",
    "results/tables/d09_table.tsv",
    "results/tables/d11_table.tsv",
    "results/tables/d12_table.tsv",
    "results/tables/d1_leakage_summary.tsv",
    "results/tables/d3_split_invariants.tsv",
    "results/tables/d4_v1_vs_v2_comparison.tsv",
    "results/tables/hub_capped_rwr.tsv",
    "results/tables/stability_by_arm.tsv",
    "results/tables/stability_attribution.tsv",
    "results/tables/network_stability.tsv",
    "results/tables/consensus_table.tsv",
    "results/tables/consensus_robustness.tsv",
    "results/tables/rewiring_mixing_check.tsv",
]
FIGURES = [
    "figures/figure1_pilot_overview.png",
    "figures/figure5_stability_atlas.png",
    "figures/figure5b_stability_attribution.png",
    "figures/figure6_network_stability.png",
    "figures/figure7_consensus_prior.png",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git(*args) -> str:
    return subprocess.check_output(["git", *args], cwd=str(ROOT),
                                   text=True).strip()


def main() -> int:
    commit = git("rev-parse", "v1.0-rc1")
    manifest = {
        "release": "v1.0-rc1",
        "commit": commit,
        "built_from_tag": "v1.0-rc1",
        "note": ("The pre-registration manifest "
                 "(docs/pre_registration_manifest.json) is frozen and "
                 "records the state filed with Zenodo; the five "
                 "fingerprint mismatches reported by "
                 "scripts/audit_reproducibility.py are the expected "
                 "consequence of the P0 repairs, which recomputed those "
                 "tables. This file is the release's own record."),
        "sha256_tables": {},
        "sha256_figures": {},
        "sha256_release_zips": {},
    }
    for rel in KEY_TABLES:
        p = ROOT / rel
        if p.is_file():
            manifest["sha256_tables"][rel] = sha256(p)
    for rel in FIGURES:
        p = ROOT / rel
        if p.is_file():
            manifest["sha256_figures"][rel] = sha256(p)
    for z in sorted(RELEASE.glob("d1-benchmark-v1.0-*.zip")):
        manifest["sha256_release_zips"][z.name] = sha256(z)

    OUT.write_text(json.dumps(manifest, indent=2, sort_keys=False) + "\n",
                   encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")
    print(f"  tables: {len(manifest['sha256_tables'])}, "
          f"figures: {len(manifest['sha256_figures'])}, "
          f"zips: {len(manifest['sha256_release_zips'])}")
    print(f"  commit: {commit}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
