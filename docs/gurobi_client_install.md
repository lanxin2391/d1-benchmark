# Installing the Gurobi client (the full installer) and getting a Named User licence

**Why the client:** the `pip install gurobipy` wheel carries the solver but
**not** the `gurobi_cl` tool, and the Named-User licence is bound to the
machine by an identifier that only the client reports. So one extra
download, and the licence route works without the campus network (which
this machine cannot present: the campus egress is a proxy address, so the
Web Licence Service route is refused).

**What we do NOT need:** a separate Gurobi installation for the solving
(the pip package does that ✓); we need the client only for its ID ✓.

---

## Step 1 — download the full installer

1. Go to **https://www.gurobi.com/downloads**.
2. Under **"Current version"**, download the **Windows 64-bit installer**
   (the file is named `gurobi-<version>-win64.exe`, e.g.
   `gurobi-12.0.2-win64.exe`). The size is roughly 150–200 MB.
3. If the site asks you to sign in first, sign in with the account you
   just registered ✓.

## Step 2 — install it

1. Run the `.exe` (a self-extracting installer ✓).
2. Accept the default destination (`C:\gurobi<version>` ✓).
3. If it offers a licence during the setup, choose the option that leaves
   the licence for later (or simply skip it ✓) — we will place the
   academic licence file ourselves ✓.
4. Finish the wizard. Nothing else is needed; **do not** uninstall or
   replace the pip `gurobipy` package ✓ — the two are independent ✓
   (the client provides the ID; the pip package does the solving ✓).

## Step 3 — read the machine's WLS ID

Open a **new** PowerShell window (so the PATH is re-read ✓) and run:

```powershell
& "C:\gurobi12.0.2\bin\gurobi_cl.exe" --help
```

(adjust the version to what you installed ✓; or add `C:\gurobi<version>\bin`
to the PATH ✓ and just run `gurobi_cl --help` ✓.)

The output contains a block like

```
WLS ID: xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
```

**Copy that WLS ID** ✓ — it is what the licence manager asks for ✓.

## Step 4 — create the Named User licence

1. Go back to **https://license.gurobi.com/manager/licenses** ✓.
2. Click **"Create a license"** ✓ (if the manager shows only the
   "Request a WLS Academic License" link, use the *User Portal* →
   **Licenses → Open the WLS manager** ✓ instead ✓).
3. Fill in:
   - **Name**: anything ✓ (e.g. `d1-netmix2`)
   - **Type**: **Named User** ✓
   - **Product**: **Gurobi Optimizer** ✓ (there is no "Python API" entry
     in the current form ✓ — the Python API is included ✓)
   - **Version**: the latest ✓
   - **WLS ID**: paste the ID from step 3 ✓
4. Create it ✓ → the licence appears in the manager's **Report** page ✓
   (the "Licenses" list of the *User Portal* is a different list ✓).

## Step 5 — download the licence file

1. On the manager's **Report** page ✓, find the licence ✓ → **"Download
   license file"** ✓ → you get `gurobi.lic` ✓.
2. Save it as **`C:\Users\PC\gurobi.lic`** ✓ (the home directory is
   searched by gurobipy ✓; the current working directory works too ✓,
   and `GRB_LICENSE_FILE` overrides both ✓).

## Step 6 — verify

```powershell
D:\Grade3\swxxx\final\D1\d1-benchmark\.venv-gnn\Scripts\python.exe `
    D:\Grade3\swxxx\final\D1\d1-benchmark\scripts\check_gurobi_licence.py
```

Expected: the script prints the licence paths it found ✓ and then

```
the model with 6000 variables solved: the status 2 (...)
=> a full licence is active: the NetMix2 arm can run
```

If it still says *"the size-limited licence is still in force"* ✓, then
the file was not found ✓ — check the path ✓ (step 5 ✓) — or it is not
valid for this machine ✓ (a Named-User licence is machine-bound ✓: if
the WLS ID came from a different machine ✗ → recreate it with this
machine's ID ✓).

## Step 7 — what we do then

We run the NetMix2 arm over the panel (six networks ✓, the MC3 mutation
heat ✓, the same AUROC-against-degree endpoint as every other arm ✓) and
write it up as the successor comparison to HotNet2.

## If you would rather not install it

The two fallbacks are already implemented and documented in
`docs/hotnet2_porting_notes.md`:

* the **trimmed-network** variant (about 1,500 nodes, which the
  size-limited licence can solve) ✓ — the method as published, on a
  smaller universe ✓;
* the **local-FDR component** alone ✓ — the significance step without
  the subnetwork selection ✓, labelled as a component ✓.

Note for the record: on this data the local FDR's null fit is itself
fragile ✗ (the mutation burden's distribution is extremely skewed ✓, so
the "middle of the histogram is nonnormal" error appears ✓) — which is
itself a finding about the mutation-score distribution ✓ and a reason to
prefer the published method ✓.
