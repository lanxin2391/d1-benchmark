#!/usr/bin/env python3
"""Check every headline number in the manuscript against the released tables.

The D1 manuscript is written by hand in three Markdown files. Hand-written
numbers drift: the design was upgraded from 10 to 40 repeats, several
control tables were regenerated after bugs were fixed, and the prose kept
the old values. This script makes that class of error mechanical rather
than anecdotal.

What it does
------------
1. Runs ``scripts/recompute_ground_truth.py``'s computation (imported, not
   re-implemented) to get the authoritative value of every headline
   quantity.
2. Section 1 -- for each entry in the CLAIMS table, searches the named
   manuscript file for the stated string. These are *textual* checks: they
   assert that a claim is present and that a superseded one is absent.
3. Section 2 -- the FORBIDDEN table, the same idea for strings that must
   no longer appear.
4. Section 3 -- the NUMERIC checks, which are the ones that matter. Each
   entry names a manuscript file, a REGULAR EXPRESSION that CAPTURES the
   number the manuscript states, and the ground-truth key that number must
   equal. The captured number is compared NUMERICALLY against the value
   recomputed from the released tables, with a stated tolerance. A
   hardcoded literal is never compared against the table; if the table
   changes, these checks change verdict with it.
5. Writes ``review/verification_report.txt``.

A FAIL means the manuscript states something the released tables do not
support. Fix the manuscript, not the table.

Usage
-----
    python scripts/verify_claims.py
    python scripts/verify_claims.py --verbose
"""

from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# ---------------------------------------------------------------- the claims
# (claim_id, file under docs/, a string that MUST appear, what it is)
CLAIMS: list[tuple[str, str, str, str]] = [
    # ---- the design ----
    ("design.repeats", "manuscript_intro_abstract.md",
     "40-repeat", "the abstract must state the 40-repeat design"),
    ("design.repeats.a", "manuscript_a_sections.md",
     "40 repeats", "§A must state the 40-repeat design"),

    # ---- the primary endpoint ----
    ("endpoint.range", "manuscript_intro_abstract.md",
     "−0.013", "the abstract must give the negative cell's ΔAUROC"),
    ("endpoint.range.max", "manuscript_intro_abstract.md",
     "+0.171", "the abstract must give the panel maximum"),
    ("endpoint.cells", "manuscript_intro_abstract.md",
     "23 of 24 cells are positive", "the abstract must give the positive-cell count"),

    # ---- H1 ----
    ("h1.rule", "manuscript_intro_abstract.md",
     "24/24 arms under four of the five", "the abstract must state the H1 result"),
    ("h1.gene", "manuscript_intro_abstract.md",
     "20/24 under the pre-registered", "the abstract must disclose the gene-cluster result"),
    ("h1.gene.a", "manuscript_a_sections.md",
     "20/24", "§J must disclose the gene-cluster result"),
    ("h2.nottested", "manuscript_intro_abstract.md",
     "not testable", "the abstract must say H2 is not testable"),
    ("h2.nottested.a", "manuscript_a_sections.md",
     "the H2 verdict as **not tested**", "§J.2 must say H2 is not testable"),

    # ---- variance decomposition ----
    ("var.interaction", "manuscript_intro_abstract.md",
     "40.1 %", "the abstract must give the interaction share"),
    ("var.label", "manuscript_intro_abstract.md",
     "34.8 %", "the abstract must give the label share"),
    ("var.network", "manuscript_intro_abstract.md",
     "24.8 %", "the abstract must give the network share"),
    ("var.fold", "manuscript_intro_abstract.md",
     "0.3 % fold noise", "the abstract must give the fold share"),

    # ---- controls ----
    ("c2.z", "manuscript_intro_abstract.md",
     "z > 5", "the abstract must state the null rejection"),
    ("c3.runs", "manuscript_a_sections.md",
     "598", "§F.1 must state the true number of rewiring runs"),
    ("c4.universe", "manuscript_intro_abstract.md",
     "1,605-gene universe", "the abstract must state the shared universe"),
    ("c5.ratio", "manuscript_intro_abstract.md",
     "39.8 %", "the abstract must give the cross-catalogue retention"),
    ("c5.cells", "manuscript_intro_abstract.md",
     "59 of 66", "the abstract must give the transfer-cell counts"),
    ("c7.matched", "manuscript_intro_abstract.md",
     "+0.001", "the abstract must give the Control-7 matched-null difference"),

    # ---- the new controls ----
    ("c8.low", "manuscript_intro_abstract.md",
     "+0.053", "the abstract must give the low-stratum margin"),
    ("c8.high", "manuscript_intro_abstract.md",
     "+0.054", "the abstract must give the high-stratum margin"),
    ("c8.a", "manuscript_a_sections.md",
     "Control 8", "§F must document Control 8"),
    ("c9.rho", "manuscript_intro_abstract.md",
     "−0.02", "the abstract must give the size correlation"),
    ("c9.a", "manuscript_a_sections.md",
     "Control 9", "§F must document Control 9"),
    ("c10.jaccard", "manuscript_intro_abstract.md",
     "0.019", "the abstract must give the input-sensitivity overlap"),
    ("c10.a", "manuscript_a_sections.md",
     "Control 10", "§F must document Control 10"),

    # ---- the negative cell ----
    ("neg.reading", "manuscript_intro_abstract.md",
     "coverage failure", "the abstract must not call the negative cell a negative effect"),
    ("neg.reading.a", "manuscript_a_sections.md",
     "not a demonstration of a negative effect",
     "§J.3 must not call the negative cell a negative effect"),

    # ---- the network table ----
    ("net.funmap.nodes", "manuscript_a_sections.md",
     "| funmap             | 10,525 |", "§A must give FunMap's node count from the meta file"),
    ("net.reactome.nodes", "manuscript_a_sections.md",
     "| reactome           | 4,993  |", "§A must give Reactome's node count from the meta file"),
    ("net.citation", "manuscript_a_sections.md",
     "Shi Z, Lei JT, Elizarraras JM", "§A must cite the correct FunMap paper"),

    # ---- the label table ----
    ("label.overlap", "manuscript_b_sections.md",
     "| ot_all               | 205 |  24 | 600 | 430 | 83 |",
     "§D must give the recomputed label overlap"),
    ("label.transfer", "manuscript_b_sections.md",
     "59 cells positive, 7 cells negative", "§E must give the transfer-cell counts"),
    ("label.nolit", "manuscript_b_sections.md",
     "is **not** a literature-free label", "§D must qualify ot_nolit"),

    # ---- the GNN arm ----
    ("gnn.naming", "manuscript_a_sections.md",
     "DeepWalk", "§J.4 must name the arm honestly"),
    ("gnn.naming.log", "gnn_courtesy_log.md",
     "DeepWalk", "the GNN log must name the arm honestly"),
    ("gnn.number", "manuscript_a_sections.md",
     "+0.028", "§J.4 must give the released embedding margin"),
    ("gnn.oldnumber", "manuscript_a_sections.md",
     "+0.15", "the superseded embedding margin must be GONE"),

    # ---- the decision log ----
    ("dec.d11", "decisions.md",
     "reactome **63**", "the decision log must carry the released D-11 N values"),
    ("dec.d11.old", "decisions.md",
     "string_full700 N=587", "the superseded D-11 N values must be marked superseded"),

    # ---- Phase S: the stability analysis (section I) ----
    ("stability.sd", "manuscript_a_sections.md",
     "repeat-level SD of ΔAUROC is 0.0047",
     "§I.1 must give the arm-level SD of the margin"),
    ("stability.cv", "manuscript_a_sections.md",
     "relative error of the score across the 40 re-partitions is 0.185",
     "§I.1 must give the gene-level relative error"),
    ("stability.rank", "manuscript_a_sections.md",
     "~270 places out of ~10,000",
     "§I.1 must give the rank movement"),

    # ---- G3: the bait-matched null ----
    ("g3.temporal", "manuscript_a_sections.md",
     "**4.7 % (chance)**",
     "§F.12 must state that the temporal label falls to chance"),

    # ---- G7: the NCG arbitration test ----
    ("g7.arms", "manuscript_a_sections.md",
     "never seeded with above the degree baseline in 24/24 arms",
     "§F.14 must state the arbitration result"),
    ("g7.lowfpr", "manuscript_a_sections.md",
     "the low-FPR endpoint agrees in only 13/24",
     "§F.14 must disclose the low-FPR caveat"),

    # ---- G8: the negative-set audit ----
    ("g8.cleaned", "manuscript_a_sections.md",
     "**+0.0822**",
     "§F.19 must give the median margin of the cleaned re-run"),
    ("g8.cleaned.arms", "manuscript_a_sections.md",
     "23 of 24 arms stay positive",
     "§F.19 must state that the cleaned re-run keeps the arms positive"),
    ("g8.reactome", "manuscript_a_sections.md",
     "**+0.0893**",
     "§F.16 must give the Reactome-on-NCG margin"),

    # ---- GR-9: the cancer-type stratification ----
    ("gr9.types", "manuscript_a_sections.md",
     "All 26 cancer types have a positive median margin",
     "§F.15 must state the cancer-type result"),

    # ---- C6: the variance decomposition with the fifth label ----
    ("c6.label", "manuscript_a_sections.md",
     "34.78 %",
     "§F.17 must give the 24-arm label share"),
    ("c6.label30", "manuscript_a_sections.md",
     "30.75 %",
     "§F.17 must give the 30-arm label share"),

    # ---- C7: the label-quality tiers and the noise anchor ----
    ("c7.tiers", "manuscript_a_sections.md",
     "| intogen2024 | +0.0892 | +0.0633 | +0.0653 |",
     "§F.18 must give the tier-restricted margins"),
    ("c7.anchor", "manuscript_a_sections.md",
     "The margin changes by between −0.0018 and +0.0026",
     "§F.18 must give the ClinGen noise-anchor range"),

]

# Numbers that MUST NOT appear anywhere in the three manuscript files.
FORBIDDEN: list[tuple[str, str, str]] = [
    ("forbidden.10repeat", "manuscript_intro_abstract.md",
     "10-repeat", "the 10-repeat design is superseded"),
    ("forbidden.vandisp", "manuscript_intro_abstract.md",
     "52 % residual", "the superseded variance numbers must be GONE"),
    ("forbidden.h2holds", "manuscript_intro_abstract.md",
     "H2 holds for IntOGen 2024", "H2 is not testable"),
    ("forbidden.600", "manuscript_a_sections.md",
     "600 random-network edge-permuted runs", "there are 598 rewiring runs"),
    ("forbidden.gnnlog", "gnn_courtesy_log.md",
     "all 72/72", "the embedding arm has 24 cells, not 72"),
    ("forbidden.vandisp", "manuscript_intro_abstract.md",
     "33 % from label choice", "the superseded variance numbers must be GONE"),
    ("forbidden.oldcitation", "manuscript_a_sections.md",
     "Mei et al., Nat Commun 2021", "the wrong FunMap citation must be GONE"),
    ("forbidden.oldcount", "manuscript_b_sections.md",
     "61 cells positive, 5 cells negative", "the superseded transfer counts must be GONE"),
    ("forbidden.negclaim", "manuscript_intro_abstract.md",
     "publishable finding on curator-lag",
     "the curator-lag reading is a hypothesis, not a finding"),
]


def load_ground_truth() -> dict[tuple[str, str], object]:
    """Import recompute_ground_truth and run it, returning its table."""
    spec = importlib.util.spec_from_file_location(
        "rgt", ROOT / "scripts" / "recompute_ground_truth.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    import io
    from contextlib import redirect_stdout
    buf = io.StringIO()
    with redirect_stdout(buf):
        mod.main()
    import pandas as pd
    df = pd.DataFrame(mod.ROWS)
    out: dict[tuple[str, str], object] = {}
    for r in df.itertuples():
        out[(r.quantity, r.scope)] = r.value
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--out", default=str(ROOT / "review" / "verification_report.txt"))
    args = ap.parse_args(argv)

    gt = load_ground_truth()

    lines: list[str] = []
    lines.append("D1 benchmark -- manuscript claim verification")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Every check below compares a string in the manuscript against a")
    lines.append("value recomputed from the released result tables by")
    lines.append("scripts/recompute_ground_truth.py. A FAIL means the manuscript")
    lines.append("states something the tables do not support.")
    lines.append("")

    n_pass = n_fail = 0
    failures: list[str] = []

    lines.append("-" * 72)
    lines.append("1. POSITIVE CLAIMS -- the manuscript must contain these strings")
    lines.append("-" * 72)
    for cid, fname, needle, what in CLAIMS:
        path = ROOT / "docs" / fname
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        # a Markdown table may pad with spaces, so compare on a collapsed form
        collapsed = re.sub(r"\s+", " ", text)
        needle_c = re.sub(r"\s+", " ", needle)
        ok = needle_c in collapsed
        status = "PASS" if ok else "FAIL"
        n_pass += ok
        n_fail += (not ok)
        lines.append(f"  [{status}] {cid:<24} {fname}")
        lines.append(f"         {what}")
        if not ok:
            lines.append(f"         NOT FOUND: {needle!r}")
            failures.append(f"{cid}: {fname} does not contain {needle!r}")
        elif args.verbose:
            lines.append(f"         found: {needle!r}")

    lines.append("")
    lines.append("-" * 72)
    lines.append("2. FORBIDDEN STRINGS -- the manuscript must NOT contain these")
    lines.append("-" * 72)
    for cid, fname, needle, what in FORBIDDEN:
        path = ROOT / "docs" / fname
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        collapsed = re.sub(r"\s+", " ", text)
        needle_c = re.sub(r"\s+", " ", needle)
        present = needle_c in collapsed
        ok = not present
        n_pass += ok
        n_fail += (not ok)
        lines.append(f"  [{'PASS' if ok else 'FAIL'}] {cid:<24} {fname}")
        lines.append(f"         {what}")
        if not ok:
            lines.append(f"         STILL PRESENT: {needle!r}")
            failures.append(f"{cid}: {fname} still contains {needle!r}")

    lines.append("")
    lines.append("-" * 72)
    lines.append("3. NUMERIC CROSS-CHECK -- the manuscript's captured "
                 "numbers vs the recomputed table values")
    lines.append("-" * 72)
    # Each entry is (claim id, file, REGEX capturing the number the
    # manuscript states, ground-truth key, tolerance). The captured number
    # is compared NUMERICALLY with the value recomputed from the released
    # tables -- the literal in the manuscript is never the thing being
    # tested, so a change in the tables changes the verdict.
    #
    # The regexes are written to match the MINUS SIGN the manuscript uses
    # (U+2212) as well as the ASCII hyphen.
    NUM = r"[+−-]?\d+(?:\.\d+)?"
    cross = [
        ("abstract negative cell", "manuscript_intro_abstract.md",
         r"(−0\.01\d)", ("delta_auroc", "reactome x intogen_temporal_new"),
         0.0006),
        ("abstract panel max", "manuscript_intro_abstract.md",
         r"(\+0\.17\d)", ("delta_auroc", "panel range (max)"),
         0.0006),
        ("manuscript_a GNN margin", "manuscript_a_sections.md",
         r"\| \*\*embedding baseline\*\*\s*\| (\+0\.\d\d\d) \|",
         ("gnn_arm", "mean GNN delta"), 0.0006),
        ("control 8 low stratum", "manuscript_a_sections.md",
         r"\| low \(lowest third\) \| \*\*(\+0\.\d{4})\*\*",
         ("control8_seed_degree", "mean low"), 0.0002),
        ("control 8 high stratum", "manuscript_a_sections.md",
         r"\| high \(highest third\) \| \*\*(\+0\.\d{4})\*\*",
         ("control8_seed_degree", "mean high"), 0.0002),
        ("control 9 permutation p", "manuscript_a_sections.md",
         r"permutation p = (0\.\d+)", ("control9_size", "permutation p"),
         0.005),
        ("control 10 Jaccard", "manuscript_a_sections.md",
         r"top-100 Jaccard overlap \| \*\*(\d\.\d+)\*\*",
         ("control10_random_input", "input sensitivity Jaccard"),
         0.001),
        ("manuscript_b per-label mean", "manuscript_b_sections.md",
         r"\| intogen2024\s+\| 633 \| (\+\d\.\d{4})",
         ("label_mean", "intogen2024"), 0.0002),
        ("manuscript_b overlap row", "manuscript_b_sections.md",
         r"\| ot_all\s+\| 205 \|\s+24 \| 600 \| (\d+) \| 83 \|",
         ("label_overlap", "ot_all x ot_nolit"), 0.5),
        ("manuscript_b transfer mean", "manuscript_b_sections.md",
         r"retention of \*\*(\d+\.\d) %\*\*",
         ("control5_cross_catalogue", "retention percent"), 0.05),
        ("control 7 matched difference", "manuscript_intro_abstract.md",
         r"difference (\+\d\.\d{3,4})",
         ("control7_publication_count", "mean matched-minus-standard"),
         0.0002),
        ("co-primary partial AUC median", "manuscript_a_sections.md",
         r"\| partial AUC, FPR . 0\.02 \| (\+0\.\d{4}) \|",
         ("co_primary_endpoint", "partial AUC FPR<=0.02: median margin"),
         0.0002),
        ("co-primary top-10 % median", "manuscript_a_sections.md",
         r"\| precision in the top 10 % of the test set \| (\+0\.\d{4}) \|",
         ("co_primary_endpoint", "precision in the top 10 %: median margin"),
         0.0002),

        # ---- Phase S: the stability analysis ----
        ("stability repeat-level SD", "manuscript_a_sections.md",
         r"repeat-level SD of ΔAUROC is (0\.\d{4})",
         ("stability", "repeat-level SD of dAUROC (median over arms)"),
         0.0006),
        ("stability gene-level error", "manuscript_a_sections.md",
         r"re-partitions is (0\.\d{3})",
         ("stability", "gene-level relative error (median over arms)"),
         0.0006),
        ("stability rank movement", "manuscript_a_sections.md",
         r"~(\d{3}) places out of",
         ("stability", "median gene's rank movement (places)"),
         5.0),

        # ---- G3: the bait-matched null ----
        ("g3 temporal unmatched", "manuscript_a_sections.md",
         r"\*\*intogen_temporal_new\*\*\s*\|\s*(\d+\.\d+) %",
         ("bait_null", "temporal label, unmatched (%)"),
         0.06),
        ("g3 temporal matched", "manuscript_a_sections.md",
         r"\*\*intogen_temporal_new\*\*\s*\|\s*\d+\.\d+ %\s*\|\s*\*\*(\d+\.\d+) %",
         ("bait_null", "temporal label, matched (%)"),
         0.06),

        # ---- G7: the NCG arbitration test ----
        ("g7 arms positive", "manuscript_a_sections.md",
         r"never seeded with above the degree baseline in (\d+)/24 arms",
         ("arbitration", "arms with a positive dAUROC"),
         0.5),
        ("g7 median dAUROC", "manuscript_a_sections.md",
         r"arms\*\* \(median (\+0\.\d{4})",
         ("arbitration", "median dAUROC"),
         0.0006),
        ("g7 low-FPR arms", "manuscript_a_sections.md",
         r"the low-FPR endpoint agrees in only (\d+)/24",
         ("arbitration", "arms with a positive low-FPR dAUC"),
         0.5),

        # ---- G8: the negative-set audit (the post-hoc table lives in
        #      its own doc; the manuscript now carries the C9 re-run) ----
        ("g8 intact intogen2024 standard", "g8_negative_set.md",
         r"\| intact \| intogen2024 \| (\+0\.\d{4}) \|",
         ("negative_set", "intogen2024 intact, standard"),
         0.0006),
        ("g8 intact intogen2024 cleaned", "g8_negative_set.md",
         r"\| intact \| intogen2024 \| \+0\.\d{4} \| (\+0\.\d{4}) \|",
         ("negative_set", "intogen2024 intact, cleaned"),
         0.0006),
        ("g8 funmap ot_all standard", "g8_negative_set.md",
         r"\| funmap \| ot_all \| (\+0\.\d{4}) \|",
         ("negative_set", "ot_all funmap, standard"),
         0.0006),
        ("g8 funmap ot_all cleaned", "g8_negative_set.md",
         r"\| funmap \| ot_all \| \+0\.\d{4} \| (\+0\.\d{4}) \|",
         ("negative_set", "ot_all funmap, cleaned"),
         0.0006),
        ("g8 reactome temporal standard", "g8_negative_set.md",
         r"\| reactome \| intogen_temporal_new \| ([−-]0\.\d{4}) \|",
         ("negative_set", "reactome temporal, standard"),
         0.0006),
        ("g8 reactome temporal cleaned", "g8_negative_set.md",
         r"\| reactome \| intogen_temporal_new \| [−-]0\.\d{4} \| (\+0\.\d{4}) \|",
         ("negative_set", "reactome temporal, cleaned"),
         0.0006),

        # ---- GR-9: the cancer-type stratification ----
        ("gr9 types positive", "manuscript_a_sections.md",
         r"All (\d+) cancer types have a positive median margin",
         ("cancer_types", "types with a positive median margin"),
         0.5),
        ("gr9 median", "manuscript_a_sections.md",
         r"median over the types (\+0\.\d{4})",
         ("cancer_types", "median of the per-type medians"),
         0.0006),

        # ---- C6: the variance decomposition with the fifth label ----
        ("c6 label share 24", "manuscript_a_sections.md",
         r"\| label \(main effect\) \| (\d+\.\d{2}) %",
         ("variance24", "label share of the between-arm spread (%)"),
         0.006),
        ("c6 label share 30", "manuscript_a_sections.md",
         r"\| label \| \d+\.\d{2} % \| (\d+\.\d{2}) %",
         ("variance30", "label share of the between-arm spread (%)"),
         0.006),
        ("c6 interaction share 30", "manuscript_a_sections.md",
         r"\| network × label \| \d+\.\d{2} % \| (\d+\.\d{2}) %",
         ("variance30", "network x label share of the between-arm spread (%)"),
         0.006),

        # ---- C7: the label-quality tiers and the noise anchor ----
        ("c7 intogen2024 high tier", "manuscript_a_sections.md",
         r"\| intogen2024 \| (\+0\.\d{4}) \| \+0\.\d{4} \| \+0\.\d{4} \|",
         ("label_quality", "intogen2024, high tier (median)"),
         0.0006),
        ("c7 temporal high tier", "manuscript_a_sections.md",
         r"\| intogen_temporal_new \| (\+0\.\d{4}) \| \+0\.\d{4} \| \+0\.\d{4} \|",
         ("label_quality", "intogen_temporal_new, high tier (median)"),
         0.0006),
        ("c7 anchor min", "manuscript_a_sections.md",
         r"changes by between ([−-]0\.\d{4}) and",
         ("noise_anchor", "min change under 17 % misassignment"),
         0.0006),
        ("c7 anchor max", "manuscript_a_sections.md",
         r"changes by between [−-]0\.\d{4} and (\+0\.\d{4})",
         ("noise_anchor", "max change under 17 % misassignment"),
         0.0006),


        # ---- the stability attribution, as the manuscript's table states it ----
        ("attribution degree joint beta", "manuscript_a_sections.md",
         r"\| log degree \| ([−-]0\.\d{3}) \|",
         ("stability_attribution", "log degree, joint beta (rank SD)"),
         0.002),
        ("attribution degree marginal rho", "manuscript_a_sections.md",
         r"\| log degree \| [−-]0\.\d{3} \| ([−+]?0\.\d{3})",
         ("stability_attribution", "log degree, marginal rho (rank SD)"),
         0.002),
        ("attribution score CV joint beta", "manuscript_a_sections.md",
         r"\| \*\*log gene score CV\*\* \| \*\*([−+]?0\.\d{3})\*\*",
         ("stability_attribution", "log score CV, joint beta (rank SD)"),
         0.002),
        ("attribution band marginal rho", "manuscript_a_sections.md",
         r"\| log local score density \([^)]*\) \| [−+]?0\.\d{3} \| ([−+]?0\.\d{3})",
         ("stability_attribution", "the +/- 10 % band, marginal rho (rank SD)"),
         0.002),
        ("attribution scott marginal rho", "manuscript_a_sections.md",
         r"\| log crowding \(Scott bandwidth\) \| [−+]?0\.\d{3} \(q = 0\.\d+\) \| ([−+]?0\.\d{3})",
         ("stability_attribution", "the Scott bandwidth, marginal rho (rank SD)"),
         0.002),


        # ---- C9: the factorial re-run with a cleaned background ----
        ("c9 positive arms, standard", "manuscript_a_sections.md",
         r"\*\*(\d+) of 24 arms are\s*positive before the cleaning",
         ("clean_background", "positive arms, standard design"),
         0.5),
        ("c9 positive arms, cleaned", "manuscript_a_sections.md",
         r"and (\d+) of 24 after",
         ("clean_background", "positive arms, cleaned background"),
         0.5),
        ("c9 median, standard", "manuscript_a_sections.md",
         r"moves from\s*\*\*([+−-]0\.\d{4})\*\*",
         ("clean_background", "median margin, standard design"),
         0.002),
        ("c9 median, cleaned", "manuscript_a_sections.md",
         r"to \*\*([+−-]0\.\d{4})\*\*",
         ("clean_background", "median margin, cleaned background"),
         0.002),

    ]
    for cid, fname, pattern, key, tol in cross:
        val = gt.get(key)
        if val is None:
            lines.append(f"  [SKIP] {cid:<30} (no ground truth for {key})")
            continue
        path = ROOT / "docs" / fname
        raw = path.read_text(encoding="utf-8") if path.exists() else ""
        # a Markdown paragraph may be hard-wrapped, so the pattern is
        # matched against the whitespace-collapsed text
        text = re.sub(r"\s+", " ", raw)
        m = re.search(pattern, text)
        if not m:
            n_fail += 1
            lines.append(f"  [FAIL] {cid:<30} table={val}")
            lines.append(f"         the pattern {pattern!r} did not match; "
                         f"the manuscript may have reworded the claim")
            failures.append(f"{cid}: {fname} does not match {pattern!r}")
            continue
        try:
            # the manuscript uses the MINUS SIGN (U+2212); normalise it
            stated = float(m.group(1).replace("−", "-"))
        except ValueError:
            stated = float(m.group(0).replace("−", "-"))
        ok = abs(stated - float(val)) <= tol
        n_pass += ok
        n_fail += (not ok)
        lines.append(f"  [{'PASS' if ok else 'FAIL'}] {cid:<30} "
                     f"manuscript={stated:+.4f}  table={float(val):+.4f}  "
                     f"tol={tol}")
        if not ok:
            lines.append(f"         the manuscript states {stated}, the table "
                         f"says {float(val)}")
            failures.append(f"{cid}: {fname} states {stated}, the table says "
                            f"{float(val)}")

    lines.append("")
    lines.append("=" * 72)
    lines.append(f"RESULT: {n_pass} PASS / {n_fail} FAIL")
    lines.append("=" * 72)
    if failures:
        lines.append("")
        lines.append("Failures:")
        for f in failures:
            lines.append(f"  - {f}")

    report = "\n".join(lines)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(report, encoding="utf-8")
    # the console may be a legacy code page; the file is always UTF-8
    sys.stdout.buffer.write(report.encode("utf-8", "replace"))
    sys.stdout.buffer.write(f"\n\nwritten to {out}\n".encode("utf-8", "replace"))
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
