"""The reproducible HotNet2 porting recipe (and the arm's runner).

Step 1 restores the clone from git and ports it with 2to3.
Step 2 applies the documented compat patches, two of them as whole-function
replacements (the HDF5 reader and writer).
Step 3 runs the driver for the requested networks and evaluates the
per-gene scores against the degree baseline.
"""
from __future__ import annotations

import json
import pathlib
import re
import shutil
import subprocess
import sys
import time

sys.stdout.reconfigure(encoding="utf-8")
ROOT = pathlib.Path(r"D:\Grade3\swxxx\final\D1\d1-benchmark")
sys.path.insert(0, str(ROOT))
CLONE = pathlib.Path(r"C:\Users\PC\AppData\Local\Temp\hotnet2_run\hotnet2")
PY = ROOT / ".venv-gnn/Scripts/python.exe"
WORK = pathlib.Path(r"C:\Users\PC\AppData\Local\Temp\hotnet2_run")
NETS_DIR = ROOT / "data/processed/networks"
LABELS_DIR = ROOT / "data/processed/labels"
TABLES = ROOT / "results/tables"
MAF = ROOT / "data/raw/cohorts/mc3.v0.2.8.PUBLIC.maf.gz"
NETWORKS = ["reactome", "funmap", "string_full700", "string_phys700",
            "rna_coexp"]
LABELS = ["intogen2024", "ot_all", "ot_nolit", "intogen_temporal_new"]
N_PERM = 20
BETA = 0.5

SAVE_FUNC = '''def save_hdf5(file_path, dictionary, compression=False):
    """
    Save or append a dictionary to an HDF5 file
    """
    f = h5py.File(file_path, 'a')
    for key in dictionary:
        # h5py 3 cannot write numpy's unicode dtype, networkx 3 hands
        # over an EdgeView rather than a list, and a plain string must
        # not be exploded into its characters
        val = dictionary[key]
        try:
            import numpy as _np
            dtype = getattr(val, 'dtype', None)
            kind = getattr(dtype, 'kind', '')
            if kind not in ('f', 'i', 'u', 'b', 'c'):
                if isinstance(val, str):
                    val = val.encode('utf-8')
                    dictionary[key] = val
                    if key in f:
                        del f[key]
                    if compression:
                        f.create_dataset(key, data=val, compression='gzip')
                    else:
                        f[key] = val
                    continue
                if isinstance(val, (list, tuple)):
                    items = list(val)
                elif hasattr(val, 'tolist'):
                    items = val.tolist()
                else:
                    items = list(val)
                arr = _np.asarray(items)
                flat = _np.array([str(x).encode('utf-8') for x in arr.ravel()])
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
'''

LOAD_FUNC = '''def load_hdf5(file_path, keys=None):
    """
    Load and return a dictionary from an HDF5 file
    """
    f = h5py.File(file_path, 'r')
    def _dec(v):
        if isinstance(v, bytes):
            return v.decode('utf-8')
        if isinstance(v, np.ndarray) and v.dtype.kind == 'S':
            return np.array([x.decode('utf-8') for x in
                             v.ravel()]).reshape(v.shape)
        return v
    if keys:
        dictionary = {key:_dec(f[key][()]) for key in keys if key in f}
    else:
        dictionary = {key:_dec(f[key][()]) for key in f}
    f.close()
    return dictionary
'''

LINE_PATCHES = [
    ("nx.to_numpy_matrix(", "nx.to_numpy_array("),
    ("nx.from_numpy_matrix(", "nx.from_numpy_array("),
    ("min(args.cores, mp.cpu_count)", "min(args.cores, mp.cpu_count())"),
    ("f[key].value", "f[key][()]"),
    ("np.int)", "int)"), ("np.float)", "float)"), ("np.bool)", "bool)"),
    ("        graph_map[network_name] = G",
     "        if isinstance(network_name, bytes):\n"
     "            network_name = network_name.decode()\n"
     "        elif not isinstance(network_name, str):\n"
     "            network_name = str(network_name)\n"
     "        graph_map[network_name] = G"),
]
SP_ALIASES = ("where", "array", "zeros", "ones", "sum", "mean", "max",
               "min", "abs", "sqrt", "exp", "log", "asarray", "arange",
               "concatenate", "isfinite", "isnan", "argsort", "unique")


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def restore() -> None:
    subprocess.run(["git", "checkout", "--", "."], cwd=str(CLONE),
                   capture_output=True, text=True)
    subprocess.run([str(PY), "-m", "lib2to3", "-w", "-n", "--no-diffs",
                    str(CLONE / "hotnet2"), str(CLONE / "scripts"),
                    str(CLONE / "HotNet2.py"), str(CLONE / "makeHeatFile.py"),
                    str(CLONE / "makeNetworkFiles.py")],
                   capture_output=True, text=True)
    n = 0
    for f in CLONE.rglob("*.py"):
        s = f.read_text(encoding="utf-8", errors="replace")
        new = s
        for a, b in LINE_PATCHES:
            new = new.replace(a, b)
        for a in SP_ALIASES:
            new = re.sub(r"\bsp\." + a + r"\b", "np." + a, new)
        for name, text in (("save_hdf5", SAVE_FUNC), ("load_hdf5", LOAD_FUNC)):
            marker = f"def {name}("
            if marker in new:
                i = new.index(marker)
                j = new.find("\ndef ", i + 10)
                if j < 0:
                    j = len(new)
                new = new[:i] + text + new[j + 1:]
        if new != s:
            f.write_text(new, encoding="utf-8")
            n += 1
    log(f"the clone restored, ported and patched ({n} files)")


def build_snv() -> pathlib.Path:
    out = WORK / "mc3_snvs.tsv"
    if out.exists():
        return out
    import pandas as pd
    log("reading the MC3 MAF")
    keep, samples = [], set()
    for chunk in pd.read_csv(MAF, sep="\t", compression="gzip",
                             usecols=["Hugo_Symbol", "Tumor_Sample_Barcode",
                                      "Variant_Classification"],
                             chunksize=500_000, low_memory=False):
        chunk = chunk[~chunk.Variant_Classification.isin(
            ["Silent", "Intron", "3'UTR", "5'UTR", "RNA", "IGR",
             "Splice_Region", "5'Flank", "3'Flank", "Downstream", "Upstream"])]
        chunk = chunk[chunk.Hugo_Symbol.notna() & (chunk.Hugo_Symbol != "")
                      & chunk.Tumor_Sample_Barcode.notna()]
        keep.append(chunk[["Tumor_Sample_Barcode", "Hugo_Symbol"]])
        samples.update(chunk.Tumor_Sample_Barcode.unique().tolist())
    df = pd.concat(keep, ignore_index=True).drop_duplicates()
    wide = (df.groupby("Tumor_Sample_Barcode")["Hugo_Symbol"]
            .apply(lambda s: "\t".join(sorted(set(s)))).reset_index())
    wide.to_csv(out, sep="\t", index=False, header=False)
    with (WORK / "mc3_samples.txt").open("w", encoding="utf-8") as fh:
        for s in sorted(samples):
            fh.write(s + "\n")
    log(f"the SNV file: {len(wide):,} samples")
    return out


def network_inputs(net: str) -> tuple[pathlib.Path, pathlib.Path]:
    import pandas as pd
    from d1.networks.io import load_network
    d = WORK / "net" / net
    d.mkdir(parents=True, exist_ok=True)
    idx_file, edge_file = d / "gene_index.tsv", d / "edges.tsv"
    if idx_file.exists() and edge_file.exists():
        return idx_file, edge_file
    genes, A = load_network(str(NETS_DIR / f"{net}.tsv"))
    index = {g: i + 1 for i, g in enumerate(genes)}
    idx_file.write_text("".join(f"{index[g]}\t{g}\n" for g in genes),
                        encoding="utf-8")
    coo = A.tocoo()
    with edge_file.open("w", encoding="utf-8") as fh:
        for u, v in zip(coo.row, coo.col):
            if u < v:
                fh.write(f"{index[genes[u]]}\t{index[genes[v]]}\n")
    return idx_file, edge_file


def run_net(net: str, snv: pathlib.Path) -> dict | None:
    import pandas as pd
    outdir = WORK / "out" / net
    viz = outdir / f"{net}-mc3" / "viz-data.json"
    if viz.exists():
        return read_scores(viz)
    idx_file, edge_file = network_inputs(net)
    netdir = WORK / "net" / net
    perm_dir = netdir / "permuted"
    perm_dir.mkdir(parents=True, exist_ok=True)
    heat = WORK / "heat.json"
    if not heat.exists():
        r = subprocess.run([str(PY), str(CLONE / "makeHeatFile.py"),
                            "mutation", "--snv_file", str(snv),
                            "--sample_file", str(WORK / "mc3_samples.txt"),
                            "-n", "mc3", "-o", str(heat)],
                           capture_output=True, text=True, cwd=str(CLONE),
                           encoding="utf-8", errors="replace")
        if not heat.exists():
            log(f"  the heat file failed: {r.stderr[-300:]}")
            return None
    if not (netdir / f"net_ppr_{BETA}.h5").exists():
        log(f"  the network files for {net} ({N_PERM} permutations)")
        t0 = time.time()
        r = subprocess.run(
            [str(PY), str(CLONE / "makeNetworkFiles.py"),
             "-e", str(edge_file), "-i", str(idx_file), "-nn", net,
             "-p", "net", "-b", str(BETA), "-np", str(N_PERM),
             "-o", str(netdir), "-c", "4"],
            capture_output=True, text=True, cwd=str(CLONE),
            encoding="utf-8", errors="replace")
        if not (netdir / f"net_ppr_{BETA}.h5").exists():
            log(f"  the network files failed: {r.stderr[-400:]}")
            return None
        log(f"    done in {time.time() - t0:.0f}s")
    outdir.mkdir(parents=True, exist_ok=True)
    log(f"  running HotNet2 on {net}")
    t0 = time.time()
    r = subprocess.run(
        [str(PY), str(CLONE / "HotNet2.py"),
         "-nf", str(netdir / f"net_ppr_{BETA}.h5"),
         "-pnp", str(perm_dir / f"net_ppr_{BETA}_##NUM##.h5"),
         "-hf", str(heat), "-o", str(outdir),
         "-np", str(N_PERM), "-hp", str(N_PERM)],
        capture_output=True, text=True, cwd=str(CLONE),
        encoding="utf-8", errors="replace")
    if not viz.exists():
        log(f"  HotNet2 failed: {r.stdout[-300:]} / {r.stderr[-300:]}")
        return None
    log(f"    done in {time.time() - t0:.0f}s")
    return read_scores(viz)


def read_scores(viz: pathlib.Path) -> dict:
    """The per-gene HotNet2 score: the largest heat value among the
    significant subnetworks the gene belongs to."""
    data = json.loads(viz.read_text(encoding="utf-8"))
    best: dict[str, float] = {}
    for delta, subnets in (data.get("subnetworks") or {}).items():
        for sn in subnets:
            for node in sn.get("nodes", []):
                g, v = node.get("name"), float(node.get("value", 0.0))
                if g is None:
                    continue
                if g not in best or v > best[g]:
                    best[g] = v
    return best


def main() -> int:
    import numpy as np
    import pandas as pd
    from d1.engine.baselines import degree_scores
    from d1.engine.evaluate import evaluate
    from d1.networks.io import load_network
    restore()
    snv = build_snv()
    rows = []
    for net in NETWORKS:
        scores = run_net(net, snv)
        if scores is None:
            continue
        genes, A = load_network(str(NETS_DIR / f"{net}.tsv"))
        deg = degree_scores(A)
        degmap = dict(zip(genes, deg))
        gset = set(genes)
        for label in LABELS:
            p = LABELS_DIR / f"{label}.tsv"
            if not p.exists():
                continue
            pos = sorted(g for g in set(pd.read_csv(p, sep="\t",
                                                    dtype=str)["gene"])
                        if g in gset)
            posset = set(pos)
            neg = sorted(g for g in genes if g not in posset)
            if len(pos) < 10 or len(neg) < 100:
                continue
            sc = np.array([scores.get(g, 0.0) for g in pos + neg])
            labels = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
            ev_h = evaluate(sc, labels, np.arange(len(sc)))
            d = np.array([degmap[g] for g in pos + neg])
            ev_d = evaluate(d, labels, np.arange(len(d)))
            rows.append({"network": net, "label": label,
                         "n_positives": len(pos),
                         "n_scored": sum(1 for g in pos + neg
                                         if scores.get(g, 0.0) > 0),
                         "hotnet2_auroc": ev_h["auroc"],
                         "degree_auroc": ev_d["auroc"],
                         "delta_auroc": ev_h["auroc"] - ev_d["auroc"]})
            log(f"    {label}: the HotNet2 AUROC {ev_h['auroc']:.4f} vs the "
                f"degree {ev_d['auroc']:.4f}")
    df = pd.DataFrame(rows)
    df.to_csv(TABLES / "c37_hotnet2.tsv", sep="\t", index=False)
    log(f"wrote the Control 37 table ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
