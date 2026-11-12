"""Reading and writing the released tables, with platform-stable bytes.

Why this module exists
----------------------
`pandas.DataFrame.to_csv` writes ``os.linesep`` as the row terminator when
it is given a path. On Windows that is ``\\r\\n``; on Linux it is ``\\n``.
The repository's ``.gitattributes`` normalises every ``*.tsv`` to LF, so a
table written on Windows and committed is byte-identical to the same
table written on Linux -- but the file **in the working tree** is not, and
therefore:

* the SHA-256 of the working-tree file differs from the SHA-256 of the
  committed blob, which breaks the release manifest's tamper-evidence
  check (4 of 6 fingerprints in the audit failed for exactly this
  reason), and
* two people regenerating the same table on different platforms get
  different hashes.

The fix is to write the row terminator explicitly as ``\\n`` and to hash
the *normalised* content, so that the hash a reader computes from a
checkout is the hash the manifest records.

The functions here are the single place that knows about this. Every
script that writes a released table should use them.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import pandas as pd

__all__ = [
    "LF",
    "write_table",
    "read_table",
    "sha256_file",
    "sha256_normalised",
]

#: the row terminator every released table uses, on every platform
LF = "\n"


def write_table(df: pd.DataFrame, path: str | Path, sep: str = "\t",
                float_format: str | None = None, **kwargs: Any) -> Path:
    """Write a table with LF line endings, whatever the platform.

    This is the only writer a released table should use.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, sep=sep, index=False, lineterminator=LF,
              float_format=float_format, **kwargs)
    return p


def read_table(path: str | Path, **kwargs: Any) -> pd.DataFrame:
    """Read a released table (tab separated, UTF-8)."""
    return pd.read_csv(path, sep="\t", **kwargs)


def _normalised_bytes(path: str | Path) -> bytes:
    """The file's bytes with every CRLF collapsed to LF."""
    return Path(path).read_bytes().replace(b"\r\n", b"\n")


def sha256_file(path: str | Path) -> str:
    """SHA-256 of the file's bytes exactly as they are on disk."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_normalised(path: str | Path) -> str:
    """SHA-256 of the file's bytes with CRLF normalised to LF.

    This is the hash a reader computes from a git checkout, so it is the
    hash the release manifest must record.
    """
    return hashlib.sha256(_normalised_bytes(path)).hexdigest()
