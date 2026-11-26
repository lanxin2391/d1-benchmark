"""Item 22: reproduce the main tables from the release zips in a clean
environment.

What this does
--------------
1. Creates a fresh virtual environment (a new interpreter, nothing
   inherited from the development environment).
2. Unpacks the three release zips into a scratch directory -- the source
   zip, the splits zip and the results zip, exactly as a reader who
   downloaded them from Zenodo would have them.
3. Installs the package from the unpacked source (`pip install -e .`).
4. Runs the reproduction: the ground-truth recomputation from the
   released tables and the claim verification against the manuscript.
5. Writes the log to `docs/reproduction_log.md` and prints the verdict.

The point is not that the numbers are right (the audit already checks
that) but that **a reader who only has the zips can get the main tables
back**, and that the environment we claim in `environment.yml` is enough.
"""
from __future__ import annotations

import json
import os
import pathlib
import shutil
import subprocess
import sys
import time
import venv

sys.stdout.reconfigure(encoding="utf-8")
REPO = pathlib.Path(r"D:\Grade3\swxxx\final\D1\d1-benchmark")
WORK = pathlib.Path(r"C:\Users\PC\AppData\Local\Temp\opencode\repro")
VENV = WORK / "env"
ZIPDIR = REPO / "release"
PY = VENV / "Scripts" / "python.exe"


def run(cmd, cwd=None, timeout=7200, env=None):
    t0 = time.time()
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=timeout,
                       env=env)
    return r, time.time() - t0


def main() -> int:
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)
    log: list[str] = []
    log.append("# The clean-machine reproduction (plan item 22)\n")
    log.append(f"**Date:** {time.strftime('%Y-%m-%d')}\n")
    log.append("The environment: a fresh virtual environment created from "
               "nothing but the Python standard library, the three release "
               "zips (source, splits, results) and `pip install -e .` from "
               "the unpacked source. The point of the exercise is not that "
               "the numbers are right -- the audit already checks that -- "
               "but that a reader who only has the zips can get the main "
               "tables back.\n")

    # ---- 1. the fresh environment ---------------------------------------
    log.append("## 1. a fresh virtual environment\n")
    # the base interpreter: the conda environment's Python 3.11, which is
    # what environment.yml asks for (a venv inherits the interpreter but
    # none of its packages, so the environment is still clean)
    base = pathlib.Path(r"E:\worksw\miniconda\envs\d1\python.exe")
    if base.exists():
        r0 = subprocess.run([str(base), "-m", "venv", str(VENV)],
                            capture_output=True, text=True)
        if r0.returncode != 0:
            print("the venv creation from the conda interpreter failed:",
                  r0.stderr[-300:])
            venv.EnvBuilder(with_pip=True, clear=True).create(str(VENV))
    else:
        venv.EnvBuilder(with_pip=True, clear=True).create(str(VENV))
    r, dt = run([str(PY), "-c", "import sys; print(sys.version)"])
    log.append(f"```\n$ python -m venv env\n$ python -c 'import sys; "
               f"print(sys.version)'\n{r.stdout.strip()}\n```\n")
    log.append(f"(the environment was created in {dt:.0f}s)\n")

    # ---- 2. unpack the zips ---------------------------------------------
    log.append("## 2. the three release zips\n")
    src_dir = WORK / "src"
    src_dir.mkdir()
    for name in ("d1-benchmark-v1.0-source.zip",
                 "d1-benchmark-v1.0-splits.zip",
                 "d1-benchmark-v1.0-results.zip"):
        z = ZIPDIR / name
        if not z.exists():
            print(f"MISSING: {z}")
            return 1
        import zipfile
        with zipfile.ZipFile(z) as zf:
            zf.extractall(src_dir)
        log.append(f"* `{name}` ({z.stat().st_size / 1e6:.1f} MB) unpacked\n")
    root = src_dir / "d1-benchmark-v1.0"
    log.append(f"The unpacked tree: `{root}`\n")
    n_tables = len(list((root / "results/tables").glob("*.tsv")))
    n_splits = len(list((root / "data/processed/splits_v40").glob("*.tsv")))
    log.append(f"It holds **{n_tables} result tables** and **{n_splits} "
               f"split files** of the 40-repeat design.\n")

    # ---- 3. install ------------------------------------------------------
    log.append("## 3. install the package from the unpacked source\n")
    # the machine's default pip index (a mirror) is not reachable from
    # here, so the reproduction uses the public index explicitly
    INDEX = ["--index-url", "https://pypi.org/simple"]
    r, dt = run([str(PY), "-m", "pip", "install", "--quiet", "--upgrade",
                  "pip"] + INDEX, cwd=str(root), timeout=3600)
    log.append("```\n$ pip install --upgrade pip "
               "--index-url https://pypi.org/simple\n"
               f"{r.stdout.strip()[-400:]}\n```\n")
    log.append("The index is given explicitly because this machine's "
               "default pip index is a mirror that is not reachable from "
               "here; the reproduction must not depend on the reader's "
               "network configuration.\n")
    # the runtime dependencies of environment.yml, installed with pip
    deps = ["numpy>=1.26,<1.27", "scipy", "pandas", "pyarrow", "networkx",
            "scikit-learn", "statsmodels", "matplotlib", "seaborn", "tqdm",
            "requests", "pytest"]
    r, dt = run([str(PY), "-m", "pip", "install", "--quiet"] + INDEX + deps,
                cwd=str(root), timeout=7200)
    log.append("```\n$ pip install " + " ".join(deps) + "\n"
               + (r.stdout.strip()[-600:] if r.returncode == 0 else
                  (r.stdout + r.stderr).strip()[-1200:]) + "\n```\n")
    log.append(f"The dependencies of `environment.yml` were installed with "
               f"pip in {dt / 60:.1f} min.\n")
    r, dt = run([str(PY), "-m", "pip", "install", "--quiet", "-e", "."]
                + INDEX, cwd=str(root), timeout=7200)
    ok = r.returncode == 0
    log.append("```\n$ pip install -e .\n"
               + (r.stdout.strip()[-800:] if ok else
                  (r.stdout + r.stderr).strip()[-1500:]) + "\n```\n")
    log.append(f"The install took {dt / 60:.1f} min "
               f"({'OK' if ok else 'FAILED'}).\n")
    if not ok:
        (REPO / "docs/reproduction_log.md").write_text(
            "\n".join(log), encoding="utf-8")
        print("the install failed; the log is written")
        return 1
    r, _ = run([str(PY), "-m", "pip", "list", "--format=freeze"], cwd=str(root))
    pkgs = [ln for ln in r.stdout.splitlines()
            if ln.split("==")[0].lower() in ("numpy", "scipy", "pandas",
                                             "networkx", "scikit-learn",
                                             "statsmodels", "matplotlib")]
    log.append("The versions that matter:\n\n```\n" + "\n".join(pkgs) +
               "\n```\n")

    # ---- 4. the reproduction ---------------------------------------------
    log.append("## 4. the reproduction: the main tables from the released "
               "data\n")
    r, dt = run([str(PY), "scripts/recompute_ground_truth.py"], cwd=str(root),
                timeout=3600)
    log.append("```\n$ python scripts/recompute_ground_truth.py\n"
               + (r.stdout + r.stderr).strip()[-800:] + "\n```\n")
    n_truth = 0
    gt = root / "review/ground_truth.tsv"
    if gt.exists():
        n_truth = sum(1 for _ in gt.read_text(encoding="utf-8").splitlines()) - 1
    log.append(f"The ground-truth table recomputed in the clean environment "
               f"holds **{n_truth} quantities**.\n")

    r, dt = run([str(PY), "scripts/verify_claims.py"], cwd=str(root),
                timeout=3600)
    tail = (r.stdout + r.stderr)
    verdict = ""
    for line in tail.splitlines():
        if line.startswith("RESULT:"):
            verdict = line.strip()
    log.append("```\n$ python scripts/verify_claims.py\n"
               + "\n".join(tail.strip().splitlines()[-14:]) + "\n```\n")
    log.append(f"**The verdict: {verdict or 'no RESULT line'}** "
               f"(the check ran {dt / 60:.1f} min).\n")

    # a table of the headline numbers as the clean machine computed them
    rows = []
    if gt.exists():
        import csv
        with gt.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["quantity"] in ("delta_auroc", "h1_verdict") and \
                        row["scope"] in ("panel range (max)",
                                          "panel range (min)",
                                          "panel median",
                                          "the arms positive",
                                          "the arms"):
                    rows.append((row["quantity"], row["scope"], row["value"]))
    if rows:
        log.append("The headline numbers as the clean machine recomputed "
                   "them:\n\n| quantity | scope | the value |\n|---|---|---|\n")
        for q, s, v in rows:
            log.append(f"| {q} | {s} | {v} |\n")
        log.append("\n")

    # ---- 5. the verdict --------------------------------------------------
    ok = ("0 FAIL" in verdict) and n_truth > 0
    log.append("### What the clean environment could not recompute, and "
               "why\n\n"
               "The clean environment produces fewer ground-truth "
               "quantities than the repository: the network descriptors "
               "(`network_nodes`, `network_lcc_nodes`, `network_lcc_edges`, "
               "ten networks each) are computed from "
               "`data/processed/networks/`, and the release zips exclude "
               "that directory on purpose -- the networks are rebuilt from "
               "the raw downloads with `scripts/build_*.py`. Everything "
               "else (the primary endpoint, every control table, the "
               "permutation audit, the rewiring grid, the pathway level and "
               "the decomposition) is recomputed from the zips alone.\n\n")
    log.append("## 5. The verdict\n\n")
    log.append(("**A reader who only has the three zips can recompute the "
                "main tables.** The clean environment installed from the "
                "unpacked source, the ground-truth recomputation and the "
                "claim verification all ran, and the verification finished "
                f"with **{verdict}**.\n") if ok else
               ("**The reproduction did not finish cleanly.** The log above "
                "says where it stopped; the environment in "
                "`environment.yml` is the one that has to work.\n"))
    log.append("What this exercise does **not** check: it does not rebuild "
               "the networks, the splits or the result tables from the raw "
                "data (that needs the data downloads and a lot of compute); "
                "it checks that the released artefacts are self-sufficient "
                "for the headline numbers, which is the claim the release "
                "makes.\n")

    (REPO / "docs/reproduction_log.md").write_text("\n".join(log),
                                                    encoding="utf-8")
    print(f"the verdict: {verdict or 'none'}; the ground truth: {n_truth} "
          "quantities")
    print("the log: docs/reproduction_log.md")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
