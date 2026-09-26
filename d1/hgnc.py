"""HGNC gene-ID mapper (Person B's deliverable B1).

Maps every gene ID the D1 benchmark encounters to a single HGNC approved symbol,
following the rules in docs/decisions.md D-02.

Contract C6: ``from d1.hgnc import HGNCMapper; m = HGNCMapper(path); m.map(s, id_type)``
returns a pandas Series of approved symbols (NaN where the input is unmapped,
ambiguous or rejected).

Design notes
------------
* The mapper is built once and is immutable. ``map()`` is vectorised; millions of
  STRING rows go through it without a Python loop.
* ``id_type`` is one of ``"symbol", "entrez", "ensembl", "uniprot"``. Symbols are
  matched case-insensitively. Ensembl versions and UniProt isoform suffixes are
  stripped before lookup (matches B1's plan in the work protocol).
* An ID that points to more than one approved symbol is *ambiguous* and is mapped
  to NaN. The set of ambiguous keys is kept on the instance as
  ``self.ambiguous_<id_type>`` for inspection.
* Previous symbols are accepted as a fallback for plain ``"symbol"`` lookups when
  (a) the previous symbol maps to exactly one approved symbol and (b) it is not
  itself also an approved symbol of a different gene. Aliases are accepted under
  the same uniqueness rule but only after previous-symbol lookup fails.
* The full approved-symbol table (used for inspection) is exposed as
  ``m.table`` (DataFrame) and ``m.approved`` (set of uppercase symbols).

Locking
-------
The mapper is frozen at S1 (end of W1). Any change after that point needs a
dated entry in docs/decisions.md and a coordinated re-run of every dependent
artefact (network files, label files, splits).
"""
from __future__ import annotations

import os
from typing import Iterable, Literal, Mapping

import numpy as np
import pandas as pd

IDType = Literal["symbol", "entrez", "ensembl", "uniprot"]
_MULTI_SEP = "|"
_QUOTE = '"'


def _split_multi(value):
    """Split an HGNC |-separated, possibly-quoted cell into a list of trimmed tokens.

    HGNC wraps multi-valued cells in double quotes and uses ``|`` as a separator.
    Some Ensembl IDs also appear without quotes; tolerate both.
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    s = str(value).strip()
    if not s or s.upper() in {"NA", "N/A", ""}:
        return []
    if s.startswith(_QUOTE) and s.endswith(_QUOTE):
        s = s[1:-1]
    return [tok.strip() for tok in s.split(_MULTI_SEP) if tok.strip()]


class HGNCMapper:
    """Vectorised HGNC approved-symbol mapper.

    Parameters
    ----------
    path : str or path-like
        Path to ``hgnc_complete_set.txt`` (tab-separated). Only ``status == 'Approved'``
        rows are kept; the rest of the table is discarded.
    """

    REQUIRED_COLUMNS = {
        "symbol", "prev_symbol", "alias_symbol",
        "entrez_id", "ensembl_gene_id", "uniprot_ids",
        "status", "locus_group",
    }
    PROT_ID_COLUMN = "uniprot_ids"

    def __init__(self, path: str | os.PathLike):
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
        missing = self.REQUIRED_COLUMNS - set(df.columns)
        if missing:
            raise ValueError(
                f"HGNC table at {path} is missing expected columns {sorted(missing)}"
            )
        df = df[df["status"] == "Approved"].copy()
        if df.empty:
            raise ValueError(f"No rows with status == 'Approved' in {path}")
        df["symbol"] = df["symbol"].str.strip()
        df = df[df["symbol"].astype(bool)]

        self.table = df.reset_index(drop=True)
        self.approved = set(self.table["symbol"].str.upper())

        self.maps: dict[str, dict[str, str]] = {}
        self.ambiguous: dict[str, set[str]] = {}

        self.maps["symbol"] = self._build_symbol_map()
        self.maps["entrez"] = self._build_simple_map("entrez_id")
        self.maps["ensembl"] = self._build_simple_map("ensembl_gene_id", strip_versions=True)
        self.maps["uniprot"] = self._build_simple_map(self.PROT_ID_COLUMN, strip_isoforms=True)

    # ------------------------------------------------------------------ builders
    def _build_symbol_map(self) -> dict[str, str]:
        """Build the symbol -> symbol map following decision D-02.

        Rules (in this order):
        1. Every approved symbol maps to itself, case-insensitively. If two
           approved rows happen to share an upper-cased symbol we treat the key
           as ambiguous and drop the self-mapping (so neither gene is silently
           merged with the other).
        2. A previous symbol is added as a fallback for ``symbol`` lookups when
           (a) it is not itself an approved symbol and (b) it is not ambiguous
           across multiple previous-symbol rows.
        3. An alias is added only after previous-symbol resolution fails and
           only under the same uniqueness rules.
        """
        m: dict[str, str] = {}
        ambiguous: set[str] = set()

        # 1) approved symbol -> itself
        for sym in self.table["symbol"]:
            key = sym.upper()
            if key in m and m[key] != sym:
                # Two approved rows share the same upper-cased symbol -> ambiguous.
                ambiguous.add(key)
                m.pop(key, None)
            else:
                m[key] = sym

        # 2) previous symbols as a fallback
        for _, row in self.table.iterrows():
            approved = row["symbol"]
            for prev in _split_multi(row.get("prev_symbol", "")):
                key = prev.upper()
                if key in self.approved:
                    # A previous symbol that is *also* a current approved symbol
                    # of a different gene stays a no-op (the approved mapping
                    # above already wins). We never add it to ``ambiguous``
                    # because the approved self-mapping must not be removed.
                    continue
                if key in m and m[key] != approved:
                    ambiguous.add(key)
                else:
                    m[key] = approved

        # 3) aliases as a last-resort fallback
        for _, row in self.table.iterrows():
            approved = row["symbol"]
            for alias in _split_multi(row.get("alias_symbol", "")):
                key = alias.upper()
                if key in self.approved:
                    continue
                if key in m and m[key] != approved:
                    ambiguous.add(key)
                else:
                    m[key] = approved

        # Remove only the truly ambiguous keys (those whose previous/alias
        # entries conflicted); keep all approved self-mappings.
        for key in list(m):
            if key in ambiguous:
                del m[key]
        self.ambiguous["symbol"] = ambiguous
        return m

    def _build_simple_map(self, column: str, *, strip_versions: bool = False,
                          strip_isoforms: bool = False) -> dict[str, str]:
        """Build a 1:many column -> approved-symbol map.

        Keys that resolve to more than one approved symbol are dropped (NaN at
        lookup time). Version suffixes on Ensembl and isoform suffixes on UniProt
        are removed before lookup so the most common input forms work directly.
        """
        m: dict[str, str] = {}
        ambiguous: set[str] = set()
        for _, row in self.table.iterrows():
            approved = row["symbol"]
            for raw in _split_multi(row.get(column, "")):
                key = raw
                if strip_versions and "." in key:
                    key = key.split(".", 1)[0]
                if strip_isoforms and "-" in key:
                    key = key.split("-", 1)[0]
                if not key:
                    continue
                if key in m and m[key] != approved:
                    ambiguous.add(key)
                else:
                    m[key] = approved
        for key in list(m):
            if key in ambiguous:
                del m[key]
        self.ambiguous[column] = ambiguous
        return m

    # ------------------------------------------------------------------ public API
    def map(self, series: pd.Series | Iterable[str], id_type: IDType) -> pd.Series:
        """Map a Series of IDs to approved HGNC symbols.

        Returns a Series the same length as the input, with the original index.
        Unmapped, empty and ambiguous IDs become NaN.
        """
        if id_type not in self.maps:
            raise ValueError(
                f"unknown id_type {id_type!r}; expected one of {sorted(self.maps)}"
            )
        s = pd.Series(series, dtype="object").astype(str).str.strip()
        empty = s == ""
        if id_type == "symbol":
            s = s.str.upper()
        elif id_type == "ensembl":
            s = s.str.split(".").str[0]
        elif id_type == "uniprot":
            s = s.str.split("-").str[0]
        out = s.map(self.maps[id_type])
        # Treat Python's default None as NaN explicitly.
        out = out.where(~empty, np.nan)
        return out

    def is_approved(self, symbol: str) -> bool:
        """True if ``symbol`` is an HGNC approved symbol (case-insensitive)."""
        return str(symbol).strip().upper() in self.approved

    # ------------------------------------------------------------------ diagnostics
    def coverage(self, series: pd.Series, id_type: IDType) -> dict[str, float]:
        """Return a small dict summarising how well ``series`` of given id_type maps.

        Useful at S1 to report how many STRING ENSPs, IntOGent Entrez IDs, etc.
        survive the mapping.
        """
        n = len(series)
        if n == 0:
            return {"n": 0, "n_mapped": 0, "frac_mapped": float("nan")}
        mapped = self.map(series, id_type)
        n_mapped = int(mapped.notna().sum())
        return {"n": int(n), "n_mapped": n_mapped, "frac_mapped": n_mapped / n}


def load_hgnc_map_table(path: str | os.PathLike) -> pd.DataFrame:
    """Read HGNC and return only the columns useful for the human-readable
    ``hgnc_map.tsv`` that the work protocol (B1) asks B to write.

    The full table is large; this strips down to the four ID columns plus symbol,
    status and locus_group, and explodes multi-valued cells (Ensembl, UniProt,
    previous symbol, alias) so every row corresponds to a single mapping.
    """
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    keep = ["symbol", "status", "locus_group", "prev_symbol", "alias_symbol",
            "entrez_id", "ensembl_gene_id", "uniprot_ids"]
    df = df[[c for c in keep if c in df.columns]].copy()
    for col in ("prev_symbol", "alias_symbol", "ensembl_gene_id", "uniprot_ids"):
        if col not in df.columns:
            continue
        df[col] = df[col].apply(_split_multi)
    df = df.explode("prev_symbol", ignore_index=True) if "prev_symbol" in df.columns else df
    df = df.explode("alias_symbol", ignore_index=True) if "alias_symbol" in df.columns else df
    df = df.explode("ensembl_gene_id", ignore_index=True) if "ensembl_gene_id" in df.columns else df
    df = df.explode("uniprot_ids", ignore_index=True) if "uniprot_ids" in df.columns else df
    return df
