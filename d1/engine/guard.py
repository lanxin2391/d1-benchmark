"""Shared guard for the control scripts: refuse to silently overwrite a released table.

WHY: the control scripts write their released tables directly
(`results/tables/control8_*.tsv`, `control9_*.tsv`, `control10_*.tsv`).
A run restricted to a subset of the panel -- `--only-net reactome
--only-label ot_all` -- would therefore replace the full table with a
few rows, and a subsequent release-manifest rebuild would bless the
truncated table as the release. This happened twice during the review
pass.

The guard below makes that impossible: a restricted run writes to a
DIFFERENT file (with the restriction appended to the name) and says so,
and a full run overwrites the released table as intended.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def resolve_out_path(released: Path, *, nets: list[str] | None,
                     labels: list[str] | None, extra: str | None = None) -> tuple[Path, bool]:
    """Return (path to write, is_the_released_table).

    A full run (no --only-net / --only-label restriction) writes to
    `released`. A restricted run writes to `released` with a
    `_restricted` suffix, so the released artefact is never truncated
    by accident.
    """
    restricted = bool(nets) or bool(labels)
    if not restricted:
        return released, True
    parts = []
    if nets:
        parts.append("net-" + "-".join(nets))
    if labels:
        parts.append("label-" + "-".join(labels))
    if extra:
        parts.append(extra)
    suffix = "_restricted-" + "_".join(parts)
    out = released.with_name(released.stem + suffix + released.suffix)
    return out, False


def write_table(df: pd.DataFrame, released: Path, *, nets: list[str] | None,
                labels: list[str] | None, extra: str | None = None,
                **kwargs) -> Path:
    """Write a control table, honouring the restricted-run guard.

    Returns the path actually written. Prints a warning when the
    released table was deliberately NOT overwritten.
    """
    out, is_released = resolve_out_path(released, nets=nets, labels=labels,
                                        extra=extra)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False, lineterminator="\n", **kwargs)
    if not is_released:
        print(f"  NOTE: this was a RESTRICTED run, so the released table "
              f"{released.name} was NOT overwritten.")
        print(f"        The partial result was written to {out.name}.")
        print(f"        Re-run without --only-net / --only-label to update "
              f"the released table.")
    return out
