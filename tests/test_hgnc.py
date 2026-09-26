"""Tests for the HGNC mapper (B1). The mapper is frozen at S1."""
from __future__ import annotations

import pandas as pd
import pytest

from d1.hgnc import HGNCMapper, load_hgnc_map_table


@pytest.fixture()
def hgnc_path(tmp_path):
    """Write a minimal HGNC table covering the cases the tests check."""
    p = tmp_path / "hgnc.tsv"
    p.write_text(
        "symbol\tstatus\tlocus_group\tprev_symbol\talias_symbol\t"
        "entrez_id\tensembl_gene_id\tuniprot_ids\n"
        'TP53\tApproved\tprotein-coding gene\t"TP53A|TP53P"\t'
        '"P53|TRP53"\t7157\tENSG00000141510\t"P04637|P04637-2"\n'
        'MDM2\tApproved\tprotein-coding gene\t"HDM2"\t\t'
        '4193\tENSG00000135679\t"Q00987"\n'
        'TESTA\tApproved\tprotein-coding gene\t\t\t\t\t\n'
        'TESTB\tApproved\tprotein-coding gene\t"BADSYM"\t"BADALIAS"\t\t\t\n'
        # Two approved rows with the same previous symbol -> ambiguous
        'TESTC\tApproved\tprotein-coding gene\t"SHARED"\t\t\t\t\n'
        'TESTD\tApproved\tprotein-coding gene\t"SHARED"\t\t\t\t\n'
    )
    return p


def test_maps_approved_symbol_to_itself(hgnc_path):
    m = HGNCMapper(hgnc_path)
    out = m.map(pd.Series(["TP53", "MDM2", "TESTA"]), "symbol")
    assert list(out) == ["TP53", "MDM2", "TESTA"]


def test_is_case_insensitive(hgnc_path):
    m = HGNCMapper(hgnc_path)
    out = m.map(pd.Series(["tp53", "Tp53", "tP53"]), "symbol")
    assert list(out) == ["TP53", "TP53", "TP53"]


def test_accepts_unique_previous_symbol(hgnc_path):
    m = HGNCMapper(hgnc_path)
    # HDM2 was a previous symbol of MDM2 (and not itself approved)
    assert m.map(pd.Series(["HDM2"]), "symbol").iloc[0] == "MDM2"
    # TP53P was a previous symbol of TP53
    assert m.map(pd.Series(["TP53P"]), "symbol").iloc[0] == "TP53"


def test_rejects_ambiguous_previous_symbol(hgnc_path):
    m = HGNCMapper(hgnc_path)
    # SHARED is a previous symbol of both TESTC and TESTD -> ambiguous -> NaN
    assert pd.isna(m.map(pd.Series(["SHARED"]), "symbol").iloc[0])


def test_approved_self_mapping_survives_prev_alias_clash(hgnc_path):
    """Decision D-02: an approved symbol always maps to itself even when its
    upper-cased form also appears as another gene's previous symbol or alias."""
    # Add a row whose prev_symbol points at TESTA (an approved symbol of another
    # gene). The mapper must keep TESTA -> TESTA; the new row's prev_symbol is
    # silently ignored because the key is already approved.
    extra = (
        "TESTX\tApproved\tprotein-coding gene\t\"TESTA\"\t\t\t\t\n"
    )
    import pathlib
    p = pathlib.Path(hgnc_path)
    p.write_text(p.read_text() + extra)
    m = HGNCMapper(p)
    assert m.map(pd.Series(["TESTA"]), "symbol").iloc[0] == "TESTA"


def test_entrez_and_ensembl_and_uniprot(hgnc_path):
    m = HGNCMapper(hgnc_path)
    assert m.map(pd.Series(["7157", "4193"]), "entrez").tolist() == ["TP53", "MDM2"]
    # Ensembl with version suffix is stripped
    assert m.map(pd.Series(["ENSG00000141510.12"]), "ensembl").iloc[0] == "TP53"
    # UniProt isoform suffix is stripped
    assert m.map(pd.Series(["P04637-2", "Q00987"]), "uniprot").tolist() == ["TP53", "MDM2"]


def test_unmapped_returns_nan(hgnc_path):
    m = HGNCMapper(hgnc_path)
    out = m.map(pd.Series(["NOPE", "99999999", ""]), "symbol")
    assert list(out) == [None, None, None] or all(pd.isna(out))


def test_bad_id_type_raises(hgnc_path):
    m = HGNCMapper(hgnc_path)
    with pytest.raises(ValueError, match="unknown id_type"):
        m.map(pd.Series(["TP53"]), "ensembl_protein")  # type: ignore[arg-type]


def test_coverage_helper(hgnc_path):
    m = HGNCMapper(hgnc_path)
    cov = m.coverage(pd.Series(["TP53", "FOO", "MDM2"]), "symbol")
    assert cov["n"] == 3 and cov["n_mapped"] == 2
    assert pytest.approx(cov["frac_mapped"]) == 2 / 3


def test_load_hgnc_map_table_explodes_multi(hgnc_path):
    df = load_hgnc_map_table(hgnc_path)
    # TP53 had two prev_symbols -> two TP53 rows
    tp53 = df[df.symbol == "TP53"]
    assert set(tp53.prev_symbol.dropna()) >= {"TP53A", "TP53P"}
    # TP53 had two UniProt IDs -> two TP53 rows
    assert set(tp53.uniprot_ids.dropna()) >= {"P04637", "P04637-2"}
