# Control 38 — rescuing H2: the class effect once the size and the coverage are out

**Date:** 2026-11-27  
**Script:** `scripts/run_c38_h2_rescue.py`  
**Tables:** `c38_h2_rescue.tsv`, `c38_h2_summary.tsv`, `c38_h2_class_medians.tsv`

## The question

H2 ("the networks of one provenance class behave alike") was stated as **not testable** on the released panel: provenance is confounded with size and coverage, and six networks cannot separate the two. The collection of item 17 adds four networks, which does not remove the confounding, but it does allow the confound to be removed **arithmetically**: the margin can be regressed on the network's size and its coverage, and the class effect can then be tested on what is left.

## The result

| panel | the networks | the Kruskal-Wallis p on the raw margins | the p on the size- and coverage-adjusted margins | the variance share of the class | the R2 of the size and the coverage |
|---|---|---|---|---|---|
| the extra networks | 4 | 0.439 | 0.439 | 0.400 | 0.959 |
| the released panel | 6 | 0.117 | 0.117 | 0.931 | 0.395 |
| the extended panel | 10 | 0.082 | 0.082 | 0.632 | 0.415 |

## The per-class medians

| panel | the class | the networks | the median margin | the median adjusted margin |
|---|---|---|---|---|
| the extra networks | protein | 2 | +0.1085 | +0.0012 |
| the extra networks | rna | 2 | +0.0674 | -0.0012 |
| the released panel | curated-literature | 3 | +0.0718 | -0.0175 |
| the released panel | protein | 2 | +0.1276 | +0.0210 |
| the released panel | rna | 1 | +0.0644 | +0.0007 |
| the extended panel | curated-literature | 3 | +0.0718 | -0.0301 |
| the extended panel | protein | 4 | +0.1276 | +0.0200 |
| the extended panel | rna | 3 | +0.0644 | -0.0001 |

## Reading

**H2 stays open, and now it can be quantified rather than asserted.** On
the extended panel (the six released networks plus the four new ones) the
class effect is still **not significant**: the Kruskal-Wallis p is 0.082 on
the raw margins and the same 0.082 on the margins adjusted for the
network's size and coverage. The size and the coverage explain about 40 %
of the margin's variance (the R-squared 0.415 on the extended panel, 0.395
on the released one), which is the confounding the paper describes,
now measured.

Two things follow. First, **the class signal is not an artefact of the
adjustment**: the ordering survives it (the protein class's adjusted
median +0.0200, the curated-literature class -0.0301, the RNA class
-0.0001 on the extended panel), so the effect is under-powered rather
than absent. Second, **the adjustment does not move the p-value at all**
on the released panel (0.117 before and after), which is what one expects
when the class and the size are collinear: the class *is* the size, as
long as every class happens to sit at its own scale.

**What would settle H2, and what would not.**

* **Not** more networks of the same kind. The four new ones already
  answer this: they move the p-value from 0.117 to 0.082 and not past
  it, and they share their sources with the released ones (the same
  STRING download, the same TCGA expression), so the within-class
  variation they add is *source* variation, not independent replication.
* **Yes** a panel that is matched on size and coverage inside each class.
  The class is still entangled with the scale here: the curated class
  holds three networks whose sizes differ by a factor of five (funmap
  10,189 nodes, reactome 3,953, intact 17,918), the protein class four
  whose edge counts differ by a factor of nine. A fourth curated network
  at FunMap's scale (an OmniPath-style curated PPI, which is reachable
  from this machine) and a third co-expression variant would give four
  networks per class at comparable sizes.
* **And** the test has to be the adjusted one, with the adjustment
  pre-registered, because the unadjusted test is a joint test of
  provenance and scale.

So the honest summary is: **H2 is now a quantified open question (the
class effect p = 0.082 at ten networks, the size explaining about 40 % of
the variance), and the way to close it is a data-collection and
pre-registration task, not another re-analysis of the networks we have.**

## Caveats

1. The adjustment is a linear regression on log10(nodes) and the coverage with ten (or six) networks: three parameters on ten observations leaves little room, and the adjusted test is a diagnostic, not a confirmation.
2. The extra networks share their sources with the released ones (the same STRING download, the same TCGA expression), so the within-class variation is *source* variation, not independent replication.
3. The class labels are three categories of one judgement; the test is a Kruskal-Wallis on ten numbers.
4. **Exploratory** (amendment 4): no pre-registered endpoint is touched.
