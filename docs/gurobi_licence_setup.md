# Getting an academic Gurobi licence (for the NetMix2 arm)

**Why:** NetMix2's two subnetwork families (the propagation / edge-dense
one and the cut one) are integer programs solved with Gurobi. The
`pip install gurobipy` package ships with a **size-limited** licence
(2,000 variables and 2,000 constraints), which rejects even our smallest
network (Reactome's LCC, 3,953 nodes). With an academic licence the arm
runs on the panel unchanged.

**Time:** the application takes about five minutes; the approval usually
arrives within a few days (it is an automatic approval for a valid
academic email, but it can take a week if the address is not recognised).

---

## 1. Apply

1. Open **https://www.gurobi.com/academics/** and choose
   **Academic Program** -> **Apply for a licence** (the free
   *Named-User* licence, which is tied to one machine and one user).
2. Fill in: your name, your **academic** email address (the university
   domain matters -- a personal address is rejected or sent to manual
   review), your department and institution, and the country.
3. Leave the "product of interest" as the **Python API** (`gurobipy`).
4. Submit. The confirmation e-mail gives you the licence, and the
   licence file itself is downloaded from the same account page
   ("Download the licence file", a `gurobi.lic`).

## 2. Install the licence file

Put `gurobi.lic` in your home directory:

```
C:\Users\PC\gurobi.lic
```

(gurobipy also accepts it in the current working directory, and the path
can be set explicitly with the environment variable `GRB_LICENSE_FILE`.)

Nothing else is needed: `gurobipy` is already installed in the
project's GraphSAGE environment
(`.venv-gnn\Scripts\python.exe`, version 14.0.0).

## 3. Check that it works

```powershell
D:\Grade3\swxxx\final\D1\d1-benchmark\.venv-gnn\Scripts\python.exe `
    scripts\check_gurobi_licence.py
```

The script builds a model with 6,000 variables. With the size-limited
licence it is refused ("Model too large for size-limited license"); with
the academic licence it solves and the script prints

```
=> a full licence is active: the NetMix2 arm can run
```

## 4. What happens next

Tell us that the check passed, and we run the NetMix2 arm over the
panel (six networks, the MC3 mutation heat, the same AUROC-against-degree
endpoint as every other arm) and write it up as the successor comparison
to HotNet2.

## If you do not want to apply

Two documented fallbacks, in `docs/hotnet2_porting_notes.md`:

* run the arm on a **trimmed network** (about 1,500 nodes fits the
  size-limited licence) and say that the method was measured on a
  smaller universe -- a caveat, not the method as published;
* use NetMix2's **local-FDR component alone** (the per-gene significance
  from the mutation heat) as the mutation-first score, and label it as *a
  component of NetMix2* rather than as NetMix2.

## One note on the HotNet2 arm

HotNet2 has the same scaling limit (its personalised-PageRank matrix is
O(n^2), so it only fits the smallest network) and additionally needs a
Python-2 port, which we have done and documented. The NetMix2 arm
supersedes it: it is the method's own successor, it runs in Python 3, and
it scales. The HotNet2 result stays as the one-arm feasibility check it
is.
