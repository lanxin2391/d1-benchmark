"""Item 17 -- collect more networks so that every provenance class has at
least four members (the H2 test needs that).

The additions:
  protein-derived     string_full400 and string_full900 (the raw STRING
                      files are already in the repository; the builder
                      already knows how to cut them)
  RNA-derived          two co-expression variants with a different edge
                      budget (the same TCGA expression, a different
                      density), so the class has three members

The new networks are written to data/processed/networks/ and are NOT
added to the released panel until the protocol is re-run on them; the
panel itself stays as it is.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(r"D:\Grade3\swxxx\final\D1\d1-benchmark")
sys.path.insert(0, str(ROOT))
PY = r"C:\Users\PC\AppData\Local\Programs\Python\Python313\python.exe"


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def build_string_variants() -> None:
    """The two extra STRING thresholds (400 and 900)."""
    script = ROOT / "scripts/build_string.py"
    src = (ROOT / "scripts/build_string.py").read_text(encoding="utf-8")
    if "string_full400" in src and "string_full900" in src:
        log("  the builder already knows the 400 and 900 cuts")
    for net_id in ("string_full400", "string_full900"):
        target = ROOT / "data/processed/networks" / f"{net_id}.tsv"
        if target.exists():
            log(f"  {net_id}: already built")
            continue
        # build only this network: import the module and call build()
        spec = importlib.util.spec_from_file_location(
            f"bs_{net_id}", script)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        src_path = (mod.DETAILED_PATH if "full" in net_id
                    else mod.PHYSICAL_PATH)
        threshold = int(net_id.replace("string_full", ""))
        log(f"  building {net_id} (the threshold {threshold})")
        mod.build(net_id, src_path, threshold=threshold)


def build_coexp_variants() -> None:
    """Two co-expression variants with a different edge budget."""
    script = ROOT / "scripts/build_rna_coexp.py"
    for net_id, edges in (("rna_coexp_d150k", 150_000),
                          ("rna_coexp_d250k", 250_000)):
        target = ROOT / "data/processed/networks" / f"{net_id}.tsv"
        if target.exists():
            log(f"  {net_id}: already built")
            continue
        spec = importlib.util.spec_from_file_location(f"bc_{net_id}", script)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        mod.NET_ID = net_id
        mod.TARGET_EDGES = edges
        log(f"  building {net_id} ({edges} edges)")
        mod.main()


def main() -> int:
    log("building the new networks")
    build_string_variants()
    build_coexp_variants()
    out = ROOT / "data/processed/networks"
    log("the network directory now holds:")
    for p in sorted(out.glob("*.tsv")):
        if p.name.endswith(".meta.json"):
            continue
        n = sum(1 for _ in p.open(encoding="utf-8"))
        log(f"  {p.name}: {n - 1:,} edges")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
