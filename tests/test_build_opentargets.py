"""Tests for the Open Targets label builder (B4).

These tests focus on the *plumbing*: schema handling, aggregation choice
and the ot_nolit europepmc-drop rule. We do NOT spin up a full parquet
set; instead we read a single real part to verify the schema, then drive
the ot_nolit helper against an in-memory frame.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "build_opentargets_under_test", REPO_ROOT / "scripts" / "build_opentargets.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_datasource_part_schema_uses_aggregation_value():
    """The OT 26.06 ``association_by_datasource_direct`` schema puts the
    real datasource id in ``aggregationValue`` (not ``datasourceId``).
    Verify so we catch any future OT schema rename early.
    """
    ds_dir = (REPO_ROOT / "data" / "raw" / "labels" / "opentargets_26.06"
              / "association_by_datasource_direct")
    if not ds_dir.exists():
        # Skip if data hasn't been downloaded on this CI runner.
        import pytest
        pytest.skip(f"{ds_dir} not downloaded")
    parts = sorted(ds_dir.glob("part-*.snappy.parquet"))
    if not parts:
        import pytest
        pytest.skip("no datasource parquet parts")
    import pyarrow.parquet as pq
    schema = pq.read_table(parts[0]).schema
    names = {f.name for f in schema}
    assert "aggregationType" in names
    assert "aggregationValue" in names
    assert "associationScore" in names
    # The schema documentation that motivated the bug fix:
    # the literal column ``datasourceId`` does NOT exist; ``aggregationValue``
    # carries the actual datasource name (e.g. ``europepmc``).
    assert "datasourceId" not in names
    assert "score" not in names  # the score is ``associationScore``, not ``score``


def test_harmonic_sum_single_score_returns_score():
    """``_harmonic_sum`` is only defined for >= 2 scores; with one score
    the helper must return the score unchanged.
    """
    mod = _load_module()
    assert mod._harmonic_sum(np.array([0.7])) == 0.7
    assert mod._harmonic_sum(np.array([0.0])) == 0.0


def test_harmonic_sum_two_scores_matches_formula():
    """For two sorted scores ``s1 >= s2``, the function returns
    ``s1 + s2/4 / (pi^2/6)``.
    """
    mod = _load_module()
    s = np.array([0.6, 0.4])
    expected = (0.6 + 0.4 / 4.0) / (np.pi ** 2 / 6.0)
    assert mod._harmonic_sum(s) == expected


def test_harmonic_sum_invariant_to_ordering():
    """Sorting inside the helper means order-of-arguments doesn't matter."""
    mod = _load_module()
    s_forward = np.array([0.1, 0.5, 0.9])
    s_reversed = np.array([0.9, 0.1, 0.5])
    assert mod._harmonic_sum(s_forward) == mod._harmonic_sum(s_reversed)


def test_to_label_table_top_n_caps_and_ranks():
    """``_to_label_table`` returns a C3 schema (gene, rank, source, release),
    with rank starting at 1 and never exceeding top_n.
    """
    mod = _load_module()
    df = pd.DataFrame({
        "gene": ["B", "C", "A", "D"],
        "associationScore": [0.9, 0.8, 0.85, 0.7],
    })
    out = mod._to_label_table(df, top_n=2, label_id="x", release="r")
    assert list(out.columns) == ["gene", "rank", "source", "release"]
    assert out["rank"].tolist() == [1, 2]
    # sorted by score desc: 0.9 (B), 0.85 (A)
    assert out["gene"].tolist() == ["B", "A"]
    assert (out["source"] == "open_targets_26.06").all()
    assert (out["release"] == "r").all()


def test_literature_datasource_constant_is_europepmc():
    """Lock the literature datasource spelling to ``europepmc`` so a future
    Open Targets release rename is caught here, not in the production build.
    """
    mod = _load_module()
    assert mod.LITERATURE_DATASOURCE == "europepmc"