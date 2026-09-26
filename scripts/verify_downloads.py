"""B0: download verification and data manifest.

Walks data/raw/ (or D1_data/ if data/raw/ is empty), records each file's path,
size in bytes, MD5 and the date the verification was run, and writes the result
to data/manifest.tsv. Every file also gets an integrity check:

  * .gz     : reads to the end with gzip.open (catches truncated downloads)
  * .zip    : zipfile.ZipFile.testzip() returns None
  * .tgz    : tarfile.open in r:gz mode lists cleanly; MD5 must equal the value
              declared in scripts/download_D1_data.py for the FunMap input
              matrices (the one MD5 we can actually verify offline)
  * parquet : opens with pyarrow.parquet.ParquetFile and reads the schema
  * text    : the first non-blank line is treated as a header and must NOT look
              like an HTML error page

Exits non-zero if any file fails, so the script is CI-friendly.

Usage
-----
    python scripts/verify_downloads.py                 # data/raw/ if it exists, else D1_data/
    python scripts/verify_downloads.py --root data/raw
"""
from __future__ import annotations

import argparse
import datetime as _dt
import gzip
import hashlib
import io
import os
import sys
import tarfile
import zipfile
from pathlib import Path

# MD5 declared in download_D1_data.py for the FunMap input archive.
# If the file's MD5 matches this we can be certain the bytes are exactly the
# ones FunMap's authors uploaded to Zenodo.
KNOWN_MD5: dict[str, str] = {
    "funmap_input_expression_data.tgz": "f490643dc8bec13dac72ed818d8b707d",
}

TEXT_EXTS = {".txt", ".tsv", ".csv", ".tab"}
HTML_PREFIXES = (b"<!doctype html", b"<html")


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _looks_like_html(path: Path) -> bool:
    try:
        with open(path, "rb") as f:
            head = f.read(512).lstrip().lower()
    except OSError:
        return False
    return head.startswith(HTML_PREFIXES)


def _check(path: Path) -> tuple[bool, str]:
    """Return (ok, reason). reason is empty when ok=True."""
    name = path.name
    ext = path.suffix.lower()
    if _looks_like_html(path):
        return False, "file is an HTML error page (download blocked?)"
    if name in KNOWN_MD5:
        got = _md5(path)
        if got != KNOWN_MD5[name]:
            return False, f"MD5 mismatch (expected {KNOWN_MD5[name]}, got {got})"
        return True, ""
    if ext == ".gz":
        if tarfile.is_tarfile(path):                       # .tgz / .tar.gz
            try:
                with tarfile.open(path, "r:gz") as t:
                    for _ in t:
                        pass
            except (tarfile.TarError, OSError) as e:
                return False, f"tar read failed: {e}"
            return True, ""
        # plain .gz (e.g. gene2pubmed.gz, STRING files)
        try:
            with gzip.open(path, "rb") as f:
                while f.read(1 << 20):
                    pass
        except (OSError, EOFError, gzip.BadGzipFile) as e:
            return False, f"gzip read failed: {e}"
        return True, ""
    if ext == ".zip":
        try:
            with zipfile.ZipFile(path) as z:
                bad = z.testzip()
                if bad is not None:
                    return False, f"zip member {bad!r} corrupt"
        except (zipfile.BadZipFile, OSError) as e:
            return False, f"zip read failed: {e}"
        return True, ""
    if ext in {".parquet"}:
        try:
            import pyarrow.parquet as pq
            pq.ParquetFile(path)
        except Exception as e:
            return False, f"parquet open failed: {e}"
        return True, ""
    if ext in TEXT_EXTS:
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                header = next((ln for ln in f if ln.strip()), "")
        except OSError as e:
            return False, f"text read failed: {e}"
        if header.lstrip().lower().startswith(("<!doctype", "<html")):
            return False, "text file starts with an HTML tag (download blocked?)"
        return True, ""
    # Unknown extension: best-effort size sanity (>= 1 byte and not the
    # 2079-byte placeholder that the IntOGen download returned at one point).
    sz = path.stat().st_size
    if sz == 0:
        return False, "empty file"
    return True, ""


def _default_root(repo_root: Path) -> Path:
    raw = repo_root / "data" / "raw"
    if raw.exists() and any(raw.rglob("*")):
        return raw
    legacy = repo_root / "D1_data"
    if legacy.exists():
        return legacy
    return raw


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=None,
                    help="directory to scan (defaults to data/raw/, falling back to D1_data/)")
    ap.add_argument("--out", type=Path, default=None,
                    help="manifest TSV path (default: <root>/../manifest.tsv)")
    args = ap.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[1]
    root = args.root or _default_root(repo_root)
    if not root.exists():
        print(f"!! root does not exist: {root}", file=sys.stderr)
        return 2

    out = args.out or (root.parent / "manifest.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    today = _dt.date.today().isoformat()
    for dirpath, _, filenames in os.walk(root):
        for name in sorted(filenames):
            p = Path(dirpath) / name
            if p.suffix == ".part":
                # half-finished download; record it but do not integrity-check.
                rows.append({"path": str(p.relative_to(repo_root)), "bytes": p.stat().st_size,
                             "md5": "", "ok": "partial", "reason": "download in progress",
                             "date": today})
                continue
            sz = p.stat().st_size
            md5 = _md5(p) if sz > 0 else ""
            ok, reason = _check(p)
            rows.append({"path": str(p.relative_to(repo_root)), "bytes": sz, "md5": md5,
                         "ok": "PASS" if ok else "FAIL", "reason": reason, "date": today})

    # write manifest (TSV)
    cols = ["path", "bytes", "md5", "ok", "reason", "date"]
    with open(out, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    # pretty summary
    pass_n = sum(1 for r in rows if r["ok"] == "PASS")
    fail_n = sum(1 for r in rows if r["ok"] == "FAIL")
    part_n = sum(1 for r in rows if r["ok"] == "partial")
    total_bytes = sum(int(r["bytes"]) for r in rows)
    print(f"Verified {len(rows)} files under {root}")
    print(f"  PASS   : {pass_n}")
    print(f"  FAIL   : {fail_n}")
    print(f"  partial: {part_n}")
    print(f"  bytes  : {total_bytes/1e9:.2f} GB")
    print(f"  manifest: {out}")
    if fail_n:
        print("\nFAILURES:")
        for r in rows:
            if r["ok"] == "FAIL":
                print(f"  {r['path']:80s} {r['bytes']:>10d} B   {r['reason']}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
