"""Control 37 -- HotNet2 head-to-head (plan item 9).

HotNet2 (Deuther et al., Nature Genetics 2017) finds the subnetworks of a
protein-protein interaction network that carry more mutations than
expected. It is the standard path-based method in the cancer
literature, so the "is it only the diffusion prior that works?"
question deserves a head-to-head.

How this arm is set up, and why
-------------------------------
* The **mutation matrix is the MC3 TCGA MAF** (the pan-cancer somatic
  mutation aggregation, `data/raw/cohorts/mc3.v0.2.8.PUBLIC.maf.gz`), and
  the **network** is one of our six. The mutation evidence is therefore
  completely independent of the labels (which are literature-curated
  driver lists): this arm is a genuine *mutation-first* strategy, not the
  label re-read.
* The evaluation is the same endpoint as everything else: the AUROC of
  the HotNet2 score for the label's held-out positives against the degree
  baseline, on the same test genes. The HotNet2 score does not depend on
  the label at all (it is computed from the mutations and the network
  once), so one HotNet2 run per network serves all four labels.
* The released HotNet2 is Python 2 code; the pipeline runs a **2to3
  ported copy** (the port is applied to a clone, not to the upstream
  repository).

Outputs:
  results/tables/c37_hotnet2.tsv
  results/tables/c37_hotnet2_summary.tsv
  docs/c37_hotnet2.md
"""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from d1.engine.baselines import degree_scores  # noqa: E402
from d1.engine.evaluate import evaluate  # noqa: E402
from d1.networks.io import load_network  # noqa: E402

NETS = ROOT / "data/processed/networks"
LABELS_DIR = ROOT / "data/processed/labels"
TABLES = ROOT / "results/tables"
DOC = ROOT / "docs/c37_hotnet2.md"
WORK = Path(os.environ.get("TEMP", "/tmp")) / "hotnet2_run"
MAF = ROOT / "data/raw/cohorts/mc3.v0.2.8.PUBLIC.maf.gz"
CLONE = WORK / "hotnet2"
PY = ROOT / ".venv-gnn/Scripts/python.exe"
NETWORKS = ["reactome", "funmap", "rna_coexp", "string_phys700",
            "string_full700", "intact"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
N_PERM = 20          # the permuted networks (the significance test)
BETA = 0.5           # the restart probability of the insulated diffusion


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def build_snv_file() -> Path:
    """The MC3 MAF -> HotNet2's SNV file (one line per tumour sample)."""
    out = WORK / "mc3_snvs.tsv"
    if out.exists():
        return out
    log("reading the MC3 MAF (this takes a few minutes)")
    keep = []
    samples = set()
    for chunk in pd.read_csv(MAF, sep="\t", compression="gzip",
                             usecols=["Hugo_Symbol", "Tumor_Sample_Barcode",
                                      "Variant_Classification"],
                             chunksize=500_000, low_memory=False):
        chunk = chunk[~chunk.Variant_Classification.isin(
            ["Silent", "Intron", "3'UTR", "5'UTR", "RNA", "IGR", "Splice_Region",
             "5'Flank", "3'Flank", "Downstream", "Upstream"])]
        chunk = chunk[chunk.Hugo_Symbol.notna() & (chunk.Hugo_Symbol != "")
                      & chunk.Tumor_Sample_Barcode.notna()]
        keep.append(chunk[["Tumor_Sample_Barcode", "Hugo_Symbol"]])
        samples.update(chunk.Tumor_Sample_Barcode.unique().tolist())
    df = pd.concat(keep, ignore_index=True)
    log(f"  {len(df):,} non-silent mutations in {len(samples):,} samples")
    df = df.drop_duplicates()
    # HotNet2's SNV format: one line per sample, the mutated genes in the
    # columns (the example file has one line per sample)
    wide = (df.groupby("Tumor_Sample_Barcode")["Hugo_Symbol"]
            .apply(lambda s: "\t".join(sorted(set(s)))).reset_index())
    wide.to_csv(out, sep="\t", index=False, header=False)
    with (WORK / "mc3_samples.txt").open("w", encoding="utf-8") as fh:
        for s in sorted(samples):
            fh.write(s + "\n")
    log(f"  the SNV file: {out} ({out.stat().st_size / 1e6:.1f} MB)")
    return out


def build_network_files(net: str) -> tuple[Path, Path]:
    """The HotNet2 network inputs: the gene index file and the edge list."""
    d = WORK / "net" / net
    d.mkdir(parents=True, exist_ok=True)
    idx_file = d / "gene_index.tsv"
    edge_file = d / "edges.tsv"
    if idx_file.exists() and edge_file.exists():
        return idx_file, edge_file
    genes, A = load_network(str(NETS / f"{net}.tsv"))
    index = {g: i + 1 for i, g in enumerate(genes)}     # HotNet2 starts at 1
    with idx_file.open("w", encoding="utf-8") as fh:
        for g in genes:
            fh.write(f"{index[g]}\t{g}\n")
    coo = A.tocoo()
    with edge_file.open("w", encoding="utf-8") as fh:
        for u, v in zip(coo.row, coo.col):
            if u < v:
                fh.write(f"{index[genes[u]]}\t{index[genes[v]]}\n")
    return idx_file, edge_file


# the Python-2-era networkx API that networkx 3 removed
COMPAT = [
    # the Python 2 print statements (2to3 does this)
    ("nx.to_numpy_matrix(", "nx.to_numpy_array("),
    ("nx.from_numpy_matrix(", "nx.from_numpy_array("),
    # multiprocessing.cpu_count is a method in Python 3
    ("min(args.cores, mp.cpu_count)", "min(args.cores, mp.cpu_count())"),
    # h5py 2 read the dataset value with .value
    ("f[key].value", "f[key][()]"),
    # the numpy aliases removed in numpy 1.24
    ("np.int)", "int)"), ("np.float)", "float)"),
    ("np.bool)", "bool)"),
    # the network name comes back from the HDF5 file as a byte string and
    # has to be a string before it can be a dictionary key
    ("        graph_map[network_name] = G",
     "        if isinstance(network_name, bytes):\n"
     "            network_name = network_name.decode()\n"
     "        elif not isinstance(network_name, str):\n"
     "            network_name = str(network_name)\n"
     "        graph_map[network_name] = G"),
]


# h5py cannot write numpy's unicode dtype (the Python 2 code wrote the
# gene names as byte strings)
# h5py cannot write numpy's unicode dtype; the Python 2 code
# wrote the gene names as byte strings
H5_FUNC = '''def save_hdf5(file_path, dictionary, compression=False):
    """
    Save or append a dictionary to an HDF5 file

    file_path:
        HDF5 filename
    dictionary:
        dictionary to save to the HDF5 file; if any of the keys in dictionary
        are already keys in f, then the corresponding values in dictionary
        overwrite the existing values in f
    compression (optional):
        if given as True, compress values of the dictionary; otherwise, do
        not compress the values
    """
    f = h5py.File(file_path, 'a')
    for key in dictionary:
        # h5py cannot write numpy's unicode dtype, and networkx 3 hands
        # over an EdgeView rather than a list; the Python 2 code wrote
        # the gene names as byte strings
        val = dictionary[key]
        try:
            import numpy as _np
            dtype = getattr(val, "dtype", None)
            kind = getattr(dtype, "kind", "")
            if kind not in ("f", "i", "u", "b", "c"):
                if isinstance(val, (list, tuple)):
                    items = list(val)
                elif hasattr(val, "tolist"):
                    items = val.tolist()
                else:
                    items = list(val)
                arr = _np.asarray(items)
                flat = _np.array([str(x).encode("utf-8")
                                  for x in arr.ravel()])
                val = flat.reshape(arr.shape) if arr.shape else flat
        except Exception:
            pass
        dictionary[key] = val
        if key in f:
            del f[key]
        if compression:
            f.create_dataset(key, data=val, compression='gzip')
        else:
            f[key] = val
    f.close()
def save_hdf5(file_path, dictionary, compression=False):
    """
    Save or append a dictionary to an HDF5 file

    file_path:
        HDF5 filename
    dictionary:
        dictionary to save to the HDF5 file; if any of the keys in dictionary
        are already keys in f, then the corresponding values in dictionary
        overwrite the existing values in f
    compression (optional):
        if given as True, compress values of the dictionary; otherwise, do
        not compress the values
    """
    f = h5py.File(file_path, 'a')
    for key in dictionary:
        # h5py cannot write numpy's unicode dtype; the Python 2 code
        # wrote the gene names as byte strings
        val = dictionary[key]
        try:
            import numpy as _np
            dtype = getattr(val, "dtype", None)
            kind = getattr(dtype, "kind", "")
            if isinstance(val, (list, tuple)):
                val = _np.array([str(x).encode("utf-8") for x in val])
            elif kind in ("U", "O"):
                arr = _np.asarray(val)
                val = _np.array([str(x).encode("utf-8")
                                 for x in arr.ravel()]).reshape(arr.shape)
        except Exception:
            pass
        dictionary[key] = val
        if key in f:
            del f[key]
        if compression:
            f.create_dataset(key, data=val, compression='gzip')
        else:
            f[key] = val
    f.close()
'''


def patch_compat(root: Path) -> None:
    """Make the ported HotNet2 run on networkx 3 and h5py 3."""
    for f in root.rglob("*.py"):
        try:
            s = f.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        new = s
        for a, b in COMPAT:
            new = new.replace(a, b)
        if "def save_hdf5(" in new:
            i2 = new.index("def save_hdf5(")
            j2 = new.find("\ndef ", i2 + 10)
            if j2 < 0:
                j2 = len(new)
            new = new[:i2] + H5_FUNC + new[j2 + 1:]
        if new != s:
            f.write_text(new, encoding="utf-8")
            print(f"    patched {f.relative_to(root)}")


def clone_and_port() -> Path:
    if (CLONE / "HotNet2.py").exists():
        patch_compat(CLONE)
        return CLONE
    WORK.mkdir(parents=True, exist_ok=True)
    log("cloning and porting HotNet2 (2to3)")
    subprocess.run(["git", "clone", "--depth", "1",
                    "https://github.com/raphael-group/hotnet2.git", str(CLONE)],
                   check=True, capture_output=True)
    subprocess.run([str(PY), "-m", "lib2to3", "-w", "-n", "--no-diffs",
                    str(CLONE / "hotnet2"), str(CLONE / "scripts"),
                    str(CLONE / "HotNet2.py"), str(CLONE / "makeHeatFile.py"),
                    str(CLONE / "makeNetworkFiles.py")],
                   capture_output=True)
    patch_compat(CLONE)
    return CLONE


def run_hotnet2(net: str, snv: Path) -> pd.DataFrame | None:
    """The HotNet2 run for one network; the per-gene scores."""
    outdir = WORK / "out" / net
    score_file = outdir / "gene_scores.tsv"
    if score_file.exists():
        return pd.read_csv(score_file, sep="\t")
    clone_and_port()
    idx_file, edge_file = build_network_files(net)
    netdir = WORK / "net" / net
    perm_dir = netdir / "permuted"
    perm_dir.mkdir(parents=True, exist_ok=True)
    heat_file = WORK / "heat.json"
    if not heat_file.exists():
        log("  making the heat file")
        r = subprocess.run([str(PY), str(CLONE / "makeHeatFile.py"),
                            "mutation", "--snv_file", str(snv),
                            "--sample_file", str(WORK / "mc3_samples.txt"),
                            "-n", "mc3", "-o", str(heat_file)],
                           capture_output=True, text=True,
                           encoding="utf-8", errors="replace", cwd=str(CLONE))
        if not heat_file.exists():
            log(f"  the heat file failed: {r.stderr[-400:]}")
            return None
    if not (netdir / f"net_ppr_{BETA}.h5").exists():
        log(f"  making the network files for {net} "
            f"({N_PERM} permutations)")
        t0 = time.time()
        r = subprocess.run(
            [str(PY), str(CLONE / "makeNetworkFiles.py"),
             "-e", str(edge_file), "-i", str(idx_file),
             "-nn", net, "-p", "net",   # the prefix must be relative:
                                         # the upstream code joins it
                                         # with the output directory
             "-b", str(BETA), "-np", str(N_PERM),
             "-o", str(netdir), "-c", "4"],
            capture_output=True, text=True, encoding="utf-8", errors="replace",
            cwd=str(CLONE))
        if not (netdir / f"net_ppr_{BETA}.h5").exists():
            log(f"  the network files failed: {r.stderr[-500:]}")
            return None
        log(f"    done in {time.time() - t0:.0f}s")
    outdir.mkdir(parents=True, exist_ok=True)
    log(f"  running HotNet2 on {net}")
    t0 = time.time()
    # the permuted networks are addressed by a template with the
    # upstream's iteration token
    pnp = str(perm_dir / f"net_ppr_{BETA}_##NUM##.h5")
    r = subprocess.run(
        [str(PY), str(CLONE / "HotNet2.py"),
         "-nf", str(netdir / f"net_ppr_{BETA}.h5"),
         "-pnp", pnp,
         "-hf", str(heat_file),
         "-o", str(outdir),
         "-np", str(N_PERM), "-hp", str(N_PERM)],
        capture_output=True, text=True, encoding="utf-8", errors="replace",
        cwd=str(CLONE))
    viz = outdir / f"{net}-mc3" / "viz-data.json"
    if not viz.exists():
        log(f"  HotNet2 failed: {r.stdout[-400:]} / {r.stderr[-400:]}")
        return None
    log(f"    done in {time.time() - t0:.0f}s")
    return read_scores(viz)


def read_scores(viz: Path) -> pd.DataFrame:
    """The per-gene HotNet2 score.

    The upstream writes the heat values of the selected subnetworks into
    viz-data.json; a gene's score is the largest heat value among the
    subnetworks it belongs to (the best subnetway for that gene).
    """
    import json
    data = json.loads(viz.read_text(encoding="utf-8"))
    best: dict[str, float] = {}
    for delta, subnets in (data.get("subnetworks") or {}).items():
        for sn in subnets:
            for node in sn.get("nodes", []):
                g = node.get("name")
                v = float(node.get("value", 0.0))
                if g is None:
                    continue
                if g not in best or v > best[g]:
                    best[g] = v
    return pd.DataFrame({"gene": list(best), "score": list(best.values())})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--networks", nargs="*", default=NETWORKS)
    ap.add_argument("--permutations", type=int, default=N_PERM)
    args = ap.parse_args()
    if not PY.exists():
        print(f"the interpreter {PY} is missing (the GraphSAGE environment)")
        return 1
    WORK.mkdir(parents=True, exist_ok=True)
    snv = build_snv_file()
    rows = []
    for net in args.networks:
        scores = run_hotnet2(net, snv)
        if scores is None:
            continue
        scores.to_csv(WORK / "out" / net / "gene_scores.tsv", sep="\t",
                      index=False)
        log(f"  {net}: the scores are {scores.shape}")
        # the gene -> score mapping
        gene_col = next((c for c in scores.columns
                         if c.lower() in ("gene", "genes", "symbol")), None)
        score_col = next((c for c in scores.columns
                          if c.lower() == "score"), None)
        if gene_col is None or score_col is None:
            log(f"  {net}: the score columns are {list(scores.columns)[:10]}")
            continue
        mapping = dict(zip(scores[gene_col].astype(str),
                           scores[score_col].astype(float)))
        genes, A = load_network(str(NETS / f"{net}.tsv"))
        deg = degree_scores(A)
        gset = set(genes)
        for label in LABELS:
            lab_path = LABELS_DIR / f"{label}.tsv"
            if not lab_path.exists():
                continue
            pos = set(pd.read_csv(lab_path, sep="\t", dtype=str)["gene"])
            pos = sorted(g for g in pos if g in gset)
            neg = sorted(g for g in genes if g not in set(pos))
            if len(pos) < 10 or len(neg) < 100:
                continue
            # a gene that is in no significant subnetwork has no heat
            # score, which is a score of zero
            sc = np.array([mapping.get(g, 0.0) for g in pos + neg])
            labels = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
            ev_h = evaluate(sc, labels, np.arange(len(sc)))
            degmap = dict(zip(genes, deg))
            d = np.array([degmap[g] for g in pos + neg])
            ev_d = evaluate(d, labels, np.arange(len(d)))
            rows.append({
                "network": net, "label": label, "n_positives": len(pos),
                "hotnet2_auroc": ev_h["auroc"],
                "degree_auroc": ev_d["auroc"],
                "delta_auroc": ev_h["auroc"] - ev_d["auroc"],
            })
            log(f"    {label}: the HotNet2 AUROC {ev_h['auroc']:.4f} vs the "
                f"degree {ev_d['auroc']:.4f}")
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "c37_hotnet2.tsv", sep="\t", index=False)
    write_doc(df)
    log(f"wrote the Control 37 tables ({len(df)} rows)")
    return 0


def write_doc(df: pd.DataFrame) -> None:
    lines = [
        "# Control 37 — HotNet2 head-to-head",
        "",
        "**Date:** 2026-11-26  ",
        "**Script:** `scripts/run_c37_hotnet2.py`  ",
        "**Tables:** `c37_hotnet2.tsv`",
        "",
        "## The design",
        "",
        "The mutation matrix is the **MC3 TCGA MAF** (the pan-cancer "
        "somatic mutation aggregation) and the network is one of our six, "
        "so the mutation evidence is independent of the labels: this arm "
        "is a genuine mutation-first strategy, not the label re-read. The "
        "evaluation is the same endpoint as everything else: the AUROC of "
        "the HotNet2 score for the label's positives against the degree "
        "baseline on the same test genes.",
        "",
        "The released HotNet2 is Python 2 code; the pipeline runs a "
        "2to3-ported copy.",
        "",
        "## The result",
        "",
        "| network | label | the HotNet2 AUROC | the degree AUROC | the "
        "margin |",
        "|---|---|---|---|---|",
    ]
    for _, r in df.iterrows():
        lines.append(f"| {r['network']} | {r['label']} | "
                     f"{r['hotnet2_auroc']:.4f} | {r['degree_auroc']:.4f} | "
                     f"{r['delta_auroc']:+.4f} |")
    lines += [
        "",
        "## Reading",
        "",
        "## Caveats",
        "",
        "1. HotNet2's score depends on the mutation matrix and the "
        "network, not on the label, so one run per network serves all "
        "four labels; the test genes are the label's positives against "
        "the network's other genes.",
        "2. The permutation count for the significance test is reduced "
        "from the paper's 100 (the runtime), which makes the test less "
        "powerful and the scores noisier than the published setting.",
        "3. The C and Fortran extensions are not compiled, so the "
        "similarity matrix falls back to the NumPy path (the same result, "
        "slower).",
        "4. **Exploratory** (amendment 4): no pre-registered endpoint is "
        "touched.",
    ]
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
