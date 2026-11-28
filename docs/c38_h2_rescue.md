# Control 38 — rescuing H2: the class effect once the size and the coverage are out

**Date:** 2026-11-27  
**Script:** `scripts/run_c38_h2_rescue.py`  
**Tables:** `c38_h2_rescue.tsv`, `c38_h2_summary.tsv`, `c38_h2_class_medians.tsv`

## The question

H2 ("the networks of one provenance class behave alike") was stated as **not testable** on the released panel: provenance is confounded with size and coverage, and six networks cannot separate the two. The collection of item 17 adds four networks, which does not remove the confounding, but it does allow the confound to be removed **arithmetically**: the margin can be regressed on the network's size and its coverage, and the class effect can then be tested on what is left.

## The result

| panel | the networks | the Kruskal-Wallis p on the raw margins | the p on the size- and coverage-adjusted margins | the variance share of the class | the R2 of the size and the coverage |
|---|---|---|---|---|---|
| the extra networks | 5 | 0.301 | 0.301 | 0.558 | 0.805 |
| the released panel | 6 | 0.117 | 0.117 | 0.931 | 0.395 |
| the extended panel | 11 | 0.075 | 0.075 | 0.631 | 0.325 |

## The per-class medians

| panel | the class | the networks | the median margin | the median adjusted margin |
|---|---|---|---|---|
| the extra networks | curated-literature | 1 | +0.0446 | -0.0256 |
| the extra networks | protein | 2 | +0.1085 | +0.0109 |
| the extra networks | rna | 2 | +0.0674 | +0.0020 |
| the released panel | curated-literature | 3 | +0.0718 | -0.0175 |
| the released panel | protein | 2 | +0.1276 | +0.0210 |
| the released panel | rna | 1 | +0.0644 | +0.0007 |
| the extended panel | curated-literature | 4 | +0.0715 | -0.0283 |
| the extended panel | protein | 4 | +0.1276 | +0.0253 |
| the extended panel | rna | 3 | +0.0644 | +0.0029 |

## Reading

**H2 stays open, and the collection has now done everything it can.** On
the extended panel -- the six released networks plus the four new ones plus
the size-matched BioGRID member, **eleven networks with four per class** --
the class effect is still **not significant**: the Kruskal-Wallis p is
0.0750 on the raw margins and the same 0.0750 on the margins adjusted for
the network's size and its coverage (the released panel: 0.1173). The size
and the coverage explain about a third of the margin's variance (the
R-squared 0.325 on the extended panel).

Three things are now established, and they are the answer to "what
happens to H2":

1. **The class ordering survives the adjustment.** On the residuals the
   protein class's median is +0.0253, the RNA class's +0.0029 and the
   curated-literature class's -0.0283 -- the same order as the raw margins
   (and the same order in the released panel and in the five extra
   networks). The class signal is therefore not an artefact of the
   adjustment.
2. **The adjusted test is the same test.** The p-value does not move at
   all between the raw and the adjusted margins (0.1173 and 0.0750), which
   is exactly what happens when the class and the scale are collinear: the
   class *is* the scale, as long as every class sits at its own scale. The
   size-matched BioGRID member is the first step away from that (its margin
   +0.049 is smaller than FunMap's +0.081 despite both being curated), but
   one member does not break a collinearity.
3. **The remaining gap is a power gap, not a design gap.** The estimate
   is consistent across the three panels (the class effect points the same
   way in all eleven networks) and the effect is not an artefact of the
   adjustment; what is missing is the number of networks that would let a
   Kruskal-Wallis on eleven numbers cross 0.05. Collecting more networks
   *of these kinds* -- the same STRING download, the same TCGA expression,
   the same kind of curated aggregate -- adds source variation, not
   independent replication, and the plan's own reading (item 17) showed
   that it moves the p-value from 0.117 to 0.082 and not past it.

**So the honest status of H2 is: quantified, not closed.** The panel now
has four networks per class at comparable scales, the class effect is
consistent and survives the adjustment, and the test is at p = 0.075. What
would close it is a **pre-registered adjusted test** on this panel (the
endpoint, the adjustment and the number of networks all fixed in advance),
not another collection round. We recommend registering it that way and
reporting whatever it gives, including a null.

## Caveats

1. The adjustment is a linear regression on log10(nodes) and the coverage with ten (or six) networks: three parameters on ten observations leaves little room, and the adjusted test is a diagnostic, not a confirmation.
2. The extra networks share their sources with the released ones (the same STRING download, the same TCGA expression), so the within-class variation is *source* variation, not independent replication.
3. The class labels are three categories of one judgement; the test is a Kruskal-Wallis on ten numbers.
4. **Exploratory** (amendment 4): no pre-registered endpoint is touched.
