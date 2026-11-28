"""Check whether a full Gurobi licence is active.

The `pip install gurobipy` package ships with a **size-limited** licence
(2,000 variables and 2,000 constraints). The academic (Named User)
licence removes that limit. This script builds a model that is well
beyond the size limit, so it fails with the bundled licence and succeeds
with the academic one -- an unambiguous test.

Usage
-----
    python scripts/check_gurobi_licence.py
"""
from __future__ import annotations

import os
import pathlib
import sys

sys.stdout.reconfigure(encoding="utf-8")


def main() -> int:
    try:
        import gurobipy as gp
    except ImportError:
        print("gurobipy is not installed in this interpreter")
        return 2
    print("the gurobipy version:", ".".join(map(str, gp.gurobi.version())))
    print("the env GRB_LICENSE_FILE:",
          os.environ.get("GRB_LICENSE_FILE", "(unset)"))
    for p in (pathlib.Path.home() / "gurobi.lic",
              pathlib.Path.cwd() / "gurobi.lic"):
        print(f"  the licence at {p}: {p.exists()}")
    n = 6000
    try:
        m = gp.Model("licence_probe")
        m.setParam("OutputFlag", 0)
        x = m.addVars(n, vtype=gp.GRB.BINARY)
        m.setObjective(gp.quicksum((i % 7 + 1) * x[i] for i in range(n)),
                       gp.GRB.MAXIMIZE)
        m.optimize()
        print(f"the model with {n} variables solved: the status "
              f"{m.Status} ({m.ObjVal:.0f})")
        print("=> a full licence is active: the NetMix2 arm can run")
        return 0
    except gp.GurobiError as e:
        print(f"the model with {n} variables was refused: {e}")
        print("=> the size-limited licence is still in force")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
