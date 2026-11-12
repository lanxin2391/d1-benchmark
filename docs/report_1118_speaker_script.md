# 1118 终期答辩 · 主讲稿（40 分钟）

**讲者**：本人（主讲 P）· **问答**：另三位成员（含 B），各 1–2 题
**配套材料**：`D1_defence_deck.pptx`（正文 16 页 + 附录 12 页）、`report_1118_qa_bank.md`（29 题）、
`report_1118_facts.md`（全部数字的单一来源）、`datasets_sources.md`（数据来源与登录号）

**怎么用这份稿子**
- 每节开头是**中文提示**（时间点、这一节要落到的点、指哪张图），后面是**逐字英文口播**。
- 括号里的 `【图：…】` 是你要点的图；`【数：…】` 是必须说对的数字。
- **不要背全部段落**：背住"每节三句话"（开场句 / 落点句 / 过渡句），中间用数字和图撑。
- 时间压缩方案在最后一节（若现场只剩 15 分钟）。

---

## 0:00–1:30 · 开场：我们的问题从哪来

【提示：40 秒讲清"领域—问题—我们的切入口—今天要什么"。不要一上来就报数字。】

Good morning. Our project asks a deceptively simple question that the
field has never actually answered: **when a network prior is used to
prioritise cancer-driver genes, how much of the reported gain comes from
the network's topology, how much from plain node degree, and how much
from the fact that well-studied genes have been published about
hundreds of times?** [`【图：Deck slide 2 — research map】`]
>
The motivation is not scepticism for its own sake. Network-based
prioritisation is standard practice — it is the formalisation of
guilt-by-association — and it works in the sense that random walk with
restart usually beats a degree-only ranking. But the *size* of the
reported effect has never been decomposed. And there are three
specific reasons to think the decomposition matters. [`【数：Schaefer
2015: bait usage correlates with degree r = 0.52】`]
>
First, driver genes are hubs and the enrichment is largely an artefact
of research intensity. Second, performance tracks network size and
collapses once coverage is controlled. Third, standard cross-validation
is optimistic when the same protein complexes appear on both sides of
the split.
>
So we built a **pre-registered, leakage-controlled benchmark** on
public data, and we report where the effect survives, where it shrinks
and where it inverts. Today I will show you the design, the control
stack, three result slides, and the limits of what we can claim.

---

## 1:30–6:00 · 我们依据的文献：五篇方法学 + 三篇原始研究

【提示：这是给生信出身的老师听的段落 —— 每篇都必须讲"**它问了什么 / 用了什么数据 / 结论是什么 /
留下什么**"四件事，并且说清"**我们因此做了什么设计决定**"。不要报 DOI 串，报逻辑。】

### 1.1 三篇原始研究（我们三个输入的来源）

**Bailey et al. 2018, Cell 173:371–385 — the recurrence consensus.**
[`【图：slide 3】`]

The first source paper is Bailey and the MC3 Working Group, in Cell.
Their question was: **which genes are mutated recurrently enough,
across tumour types, to be called drivers?** The data were 9,423
TCGA tumour-normal pairs spanning 26 tissue types, and the key
methodological move was the **MC3 consensus** — 26 independent
mutation-calling algorithms vote, and a gene is called a driver when
the methods agree. The result became the field's default gold
standard: the Cancer Gene Census.
>
What it leaves open, and this is the gap we work in: it is a
**purely recurrence-based** catalogue. There is no network, no prior,
and no way to separate "what the network adds" from "the catalogue's
own degree bias". **Our decision**: we use it as our primary label
(`intogen2024`, 633 genes) and we treat the curation choices it
carries as inherited bias that we have to measure rather than assume
away.

**Shi et al. 2025, Nature Cancer 6:205–222 — FunMap.**
[`【图：slide 4】`]

The second source paper is Shi and colleagues' FunMap. Their question
was: **can a pan-cancer proteogenomic functional network prioritise
cancer genes better than the existing interactomes?** The data are the
CPTAC pan-cancer proteogenomics and transcriptomics. The method is
worth one sentence of precision because it changes how we treat the
network: FunMap's edges are **not** raw interactions — each gene pair
has an XGBoost-predicted co-function probability, the pairs are added
in blocks until the **set-level log-likelihood ratio falls below 3.912
(LR = 50)**, and a hard filter requires the pair's **mutual rank to
exceed 0.95 in at least one cohort**; that threshold yields 10,525
nodes and 196,800 edges.
>
The result: FunMap outperforms STRING for cancer-gene prioritisation,
and it ships with a gold standard.
>
The limitation the authors state themselves is the one we care about:
**their gold standard is Reactome-derived** — a pair is positive if
the two genes share a "detailed" Reactome pathway of at most 200
proteins — and Reactome is *also* an input to STRING's `database`
channel. They also note the network is built from 11 cancer types in
bulk tissue, so tumour cells and microenvironment are mixed.
>
**Our decisions**: FunMap is one of our six networks because it is the
most on-topic cancer-specific prior — *and* we treat it as a
**derived, supervised, thresholded object** rather than an
interaction table. We also report its coverage explicitly, because
its 10,525 genes are exactly the protein-coding genes that CPTAC
quantified, so a largest-connected-component analysis runs on a small,
expression-biased, tumour-derived universe.

**Repana et al. 2019, Genome Biology 20:1 — the NCG catalogue.**
[`【图：slide 5】`]

The third source paper is Repana and colleagues' Network of Cancer
Genes. Their question: **which genes have been called drivers across
the published cancer sequencing literature, and how confident is each
call?** They curated 2,372 genes from the published sequencing
literature into three tiers — 708 known tumour suppressors and
oncogenes, 104 multi-screen candidates, 1,560 screening-level genes —
and attached a per-gene confidence level built on curated interaction
channels.
>
What it leaves open: **that confidence level is itself network-derived**,
and the unlabelled genes are implicitly treated as negatives. So NCG
cannot measure anything by itself.
>
**Our decision — and this is the one I would defend hardest**: we
use NCG as an **arbitration catalogue**, not as a label. It shares no
curation with our four label sets, so it can arbitrate them; and
because its confidence level is network-derived, it is also the
right instrument for the **circularity audit** I will show you later.

### 1.2 五篇方法学文献（每一条都变成我们的一个设计决定）

The three source papers give us the inputs. Five further papers give
us the **design**, and I want to be explicit about which decision each
one produced, because every one of them is a place where we could
have gone wrong.

**Picart-Armada et al. 2019, PLoS Comput Biol 15:e1007276** — benchmarking
network propagation for disease-gene identification.
They benchmarked propagation methods and found two things that changed
our design. First, **pure topology baselines against random** must be
included, not just against degree. Second, and critically, **protein
complexes leak across folds**: holding a complex together in the
cross-validation fold drops a strong network classifier from about 12
to fewer than 4.5 hits in the top 20 — a fall as large as the gap
between competing methods and the baseline. They also showed that
**the endpoint changes the ranking of methods**: AUROC and top-k
disagree.
**Our decision**: we implemented **complex-aware cross-validation**
(the whole complex stays in one fold) and we report **two extra
endpoints** — partial AUC at FPR ≤ 0.02 and top-10 % precision — next
to AUROC.

**Guala & Sonnhammer 2017, Scientific Reports 7:46598** — a large-scale
benchmark of gene prioritisation methods.
This is the paper that made us most careful. Its result is a warning:
**MaxLink — a local method with an explicit degree filter, designed so
that highly connected genes cannot rank high just by being connected —
beats the diffusion (RWR-family) methods on the full list.** Diffusion
wins only in the top 1–10 %. They also use partial AUC at FPR ≤ 0.02,
which is where we took our low-FPR endpoint from, and they warn about
the circularity of gold standards that share evidence with the network.
**Our decision**: we never claim RWR is the best method. We freeze it
and measure the **network's contribution relative to degree**, which
is the quantity Guala shows can be small — and we report the top-1–10 %
region separately because that is where diffusion actually wins.

**Hill et al. 2019** — benchmarking 17 network-contextualisation algorithms.
Their finding is the one that makes "did the learner even use its
input?" a mandatory control: **over 10,000 random gene lists, the
algorithms differ enormously in how much of their output is
input-independent.**
**Our decision**: we run a **random-input degeneracy control** — the
learner is fed shuffled inputs and must degrade to the baseline.

**Zhang et al. 2021** — a benchmark of network-based gene prioritisation.
Their contribution is methodological: they show that node2vec's
p/q parameters trade off BFS against DFS, and that node2vec is
skip-gram over random walks treated as sentences. We use this to make
our **GNN courtesy arm** honest: our p = q = 1 makes node2vec exactly
equivalent to DeepWalk-style embeddings, so we report it under that
name rather than borrowing a method's prestige.

**Wright et al. 2025** — the current state of molecular-interaction
omics.
Over 45 human interactomes, **performance tracks network size**, and
the effect largely collapses once gene-set coverage is controlled.
**Our decision**: we built the **shared-gene-universe control**
(D-12) and the **per-network common-size cap** (D-11), so that a
margin cannot be bought with a bigger gene universe.

**Schaefer et al. 2015, Frontiers in Genetics** — correcting for study bias
in protein-interaction networks.
This is the paper our third hypothesis is built on: **bait usage —
how often a gene appears in the published interaction literature —
correlates with its degree (r = 0.52), and correcting for it removes
the cancer-versus-control degree difference.** In other words, what
looks like a biological degree effect can be a publication effect.
**Our decision**: we built a **publication-count-adjusted label** and
a **bait-matched null model** in which each positive is matched to a
negative with the same publication count.

**Boutros et al. 2014, Genome Biology 15:462** — challenge-based
benchmarking.
Their principle is the one we borrowed for the whole design: **the
gold standard dominates any ranking, so the evaluator should be blind
to the algorithm**; a benchmark should be a community challenge with
hidden answers, not a leaderboard people can overfit.
**Our decision**: the NCG arbitration test is our mini-version of
that — an independent catalogue that none of our labels was built
from — and it is also the experiment we name in our validation plan.

---

## 6:00–8:00 · 知识空白与我们的三个可证伪假说

[`【图：slide 6】`]

Putting the literature together, the gap is precise: **no published
study separates the topology contribution from the degree contribution
from the study-visibility contribution, inside one frozen design with
one decision rule.** Each piece exists separately — Zhang for degree,
Schaefer for visibility, Picart-Armada for leakage, Wright for
coverage — but never together, and never on a pre-registered
endpoint.
>
We therefore made three hypotheses, and we wrote down in advance how
each one could fail.
>
**H-vis**: much of the reported margin is study visibility, not
topology. *Falsified* if the margin survives a null model that
matches positives and negatives on their publication count.
>
**H-stable**: the prior's instability is a property of the gene's own
score variance, not of topology, and it cannot be smoothed away.
*Falsified* if any topology descriptor predicts instability, or if a
regularised scorer restores the ranking.
>
**H-lag**: curated networks fail specifically on **recently
discovered** drivers. *Falsified* if the failure appears equally on
older catalogues.
>
[`【数：三个假说分别由 G3、C7/S1–S6、G8/C5 的结果检验】`]

---

## 8:00–13:00 · 我们的数据从哪里来

【提示：这一段必须具体到"哪个文件、哪个版本、我们怎么标注"。老师是生信出身，会问"你的基因符号怎么对齐的"。
**逐条讲清我们的两条数据治理规则**：D-02（HGNC 映射）与 cap。】

### 8.1 网络面板：六个网络，三类来源

[`【图：slide 7】`]

Six networks, deliberately stratified by **molecular provenance** so
that provenance is an experimental factor rather than a nuisance:
three protein-derived, one RNA-derived, two curated/literature.
>
**funmap** — the CPTAC proteogenomic network: 10,189 nodes, 196,605
edges in the largest connected component. **intact** — IntAct
human-human interactions, taxid 9606: 17,918 nodes, 564,706 edges.
**string_phys700** — STRING v12 physical subnetwork, combined score
≥ 700: 9,830 nodes, 85,576 edges. **rna_coexp** — built by us from
FunMap's own published CPTAC mRNA matrices: Pearson correlation per
cohort, Fisher-z averaged across cohorts, the top 196,800 pairs by
mean z to density-match funmap: 9,152 nodes, 196,172 edges.
**string_full700** — STRING v12 all-evidence channels at ≥ 700:
15,882 nodes, 236,712 edges. **reactome** — Reactome v86, human
only: 3,953 nodes, 20,143 edges.
>
Two harmonisation decisions matter. First, **Reactome needs a
pathway-size cap**: Reactome's co-membership graph has no numeric
edge weight, and the pathways signalling and transcription are so
large that the naive graph is nearly complete inside them, so its
connected component carries almost no information as a "prior". We
therefore follow the ≤200-protein rule that FunMap's own gold
standard uses, and we snapshot a specific Reactome version rather
than always taking "current" — because the 3- and 4-star curation
events appeared in June 2023 and historical tools drift.
Second, **the analysis is unweighted** (decision D-01), so the
quantity we measure is topology, not the weight scale.

### 8.2 标签面板：五个目录，四种策展哲学

Five label sets, and they are not interchangeable:
**intogen2024** (633 genes) — mutation recurrence, from IntOGen, and
it is the most defensible pure-positive source we have, because a
positive is defined by positive-selection evidence in mutation data,
not by literature co-occurrence or pathway co-membership — that is
exactly the independence a network-prior benchmark needs.
**ot_all** and **ot_nolit** (600 each) — Open Targets 26.06, the
second with the literature channel removed; their difference is
control 7 in label form.
**intogen_temporal_new** (152 genes) — the newest IntOGen drivers, a
deliberately harder label.
**clingen_label** (84 genes) — clinical gene-disease validity, the
hardest evidence of all: 79 Definitive, 1 Strong, 4 Moderate.
**ncg_known** (708 genes) — the NCG known tier, added later as a
pre-planned exploratory extension.
>
Two curation caveats we had to handle. IntOGen is a **gene × cancer
type** structure — 99 breast drivers mixed with 20 neuroblastoma
drivers — so we stratify by cancer type rather than treating it as a
flat list. And **ranking a label by Cancer Gene Census membership
would measure literature curation intensity, not driver prediction**,
because 75 % of the IntOGen compendium is CGC-annotated.

### 8.3 符号对齐与数据治理

Everything is mapped to HGNC approved symbols under decision D-02:
the approved symbol maps to itself; a previous symbol or alias maps
only when the mapping is one-to-one; ambiguous symbols are **dropped
and counted** — for NCG that was 83 symbols, and we keep them
visible in the release rather than silently discarding them.
>
[`【数：NCG 9 个改名基因：CARS→CARS1、MKL1→MRTFA、FGFR1OP→CEP43 等】`]
>
Two further provenance items. The **protein-complex map** used for
the complex-aware splits was built from the UniProt REST API by
cross-referencing ComplexPortal identifiers — 2,569 complexes over
3,708 genes. And every external dataset carries an explicit
identifier and release: the pre-registration is frozen at Zenodo
DOI 10.5281/zenodo.23121762, the FunMap matrices are Zenodo record
7948944, IntOGen releases 2020/2023/2024, Open Targets 26.06,
STRING v12.0, Reactome v86, IntAct 2024-09, the ClinGen file of
2026-09-26, and DepMap's CRISPR gene-effect matrix.

---

## 13:00–17:00 · 我们的设计与它的纪律

[`【图：slide 8、slide 9】`]

Now the design, and I want to spend a minute on the **discipline**,
because that is what makes the result interpretable.
>
**The learner is frozen.** RWR at α = 0.5, no tuning, the same in all
24 arms and all controls. The reason is not that RWR is best — Guala
shows a degree-filtered local method can beat it on a full list — it
is that once the learner is tunable, the quantity of interest becomes
"we tuned a method" instead of "the network contributes".
>
**The resampling unit is the repeat.** 5 folds × 40 repeats = 200
evaluations per arm; 6 networks × 4 labels × 3 alphas × 200 = the full
design. We moved from 10 to 40 repeats because the pre-registered
interval must resample repeats, and a coverage simulation showed the
fold-level interval is over-conservative (coverage 1.000) while the
repeat-level percentile is calibrated (0.936–0.944).
>
**Splits are fixed and released** with SHA-256 manifests, and the
first 10 repeats of the 40-repeat split files are byte-identical to
the original 10-repeat files, so the new design is a strict
superset.
>
**The common-size cap** is a control, not the primary design: the
primary arms use each label's own positive set, and the cap is what
equalises prevalence inside an arm. Note the order of operations: the
cap must be applied **after** restricting to the universe, otherwise
funmap holds 78 positives instead of 104 — a defect we found and
recorded as amendment A2.3i.
>
**Everything is pre-registered and every departure is written down.**
The pre-registration was frozen on 2026-09-27; there are now
amendments A2.3a to A2.3k plus a register of the fourteen analysis
families, which states for each one whether it is in the verification
family or exploratory. The statistical layer is a mixed-effects
variance decomposition, a **gene-cluster bootstrap** (the
pre-registered §9.3 unit: resample genes, not folds), a TOST
equivalence test for the borderline arm, and BH-FDR **inside the
pre-declared family only** — we do not run BH across 72 correlated
cells and call that multiplicity control.
>
And we should mention the self-audit: an independent review of the
repository found **ten P0-class defects** — among them a bootstrap
that accepted an alpha parameter and never used it, a TOST that was
not TOST, a cap file whose name collided so each network silently
overwrote the previous cap, a shared-universe control that did not use
the shared universe, a temporal control with zero replication, a
variance decomposition that treated folds as a noise bucket, and a
manuscript that quietly relaxed the H1 decision rule. All of them are
fixed and recorded in the amendments.

---

## 17:00–19:30 · 十一层控制栈

[`【图：slide 10】`]

Eleven controls, one artefact each. In order: **degree-matched random
seeds** — does the prior just re-rank hubs? **Degree-preserving
rewiring** — is the signal carried by the wiring? **Shared gene
universe** — gene-set size differences? **Temporal split** — is the
effect only in the past? **Common-size cap** — unequal positive-set
sizes? **Complex-aware cross-validation** — complex leakage?
**Hub-capped RWR** — is the gain a hub effect? **Publication-count
adjusted label** — study visibility? **Bait-matched null** — does the
gain track publication count? **Label-quality tiers** — does the gain
live in the weak positives? **ClinGen 17 % noise anchor** — is the
gain an artefact of label noise?
>
Every one of these has a number, and I will show the numbers that
changed our conclusions.

---

## 19:30–22:30 · 结果一：主端点与 H1 的准确表述

[`【图：slide 11】`]

The primary endpoint is the margin of RWR over a degree-only ranking
on identical folds: ΔAUROC = AUROC(RWR) − AUROC(degree), at α = 0.5.
>
[`【数：24 臂中 23 个为正；中位 +0.0809；最好 string_phys700 × ot_nolit +0.1711；
唯一负 reactome × intogen_temporal_new −0.0125】`]
>
The single negative cell is outlined in the figure. We do **not** claim
it is a negative effect: the equivalence test says it cannot be
distinguished from zero, so the honest statement is "on this arm there
is no detectable contribution". And it is not a network failure — it
is a **label-time** failure, which I will show in the third result.
>
Now the H1 verdict, stated exactly as the pre-registration defines it.
The rule is: the 95 % confidence interval excludes zero in **every**
arm. It is met under **four of the five** interval constructions —
fold-percentile, repeat-percentile, repeat-t and repeat-BCa all give
**24/24** — and it is **not** met under the fifth, the conservative
gene-cluster bootstrap, which gives **20/24**. The four exceptions
are **all** the newest-driver label, the one with the smallest
effects. [`【数：gene-cluster 区间宽度约 12×，因为它问的是"换一批基因来测"】`]
>
I want to be careful here, because this is exactly the kind of place
where a paper quietly relaxes a rule. We did not relax it. The
pre-registered endpoint is what it is, and we report the verdict under
all five constructions, and if you force a single hard sentence, it
is: **in 20 of 24 arms the increment over degree excludes zero at 95 %;
in the remaining four, the effect is too small for us to distinguish
it from zero.** That is what "label-dependent" means literally.
>
One more result on this slide, from the G4 work: if we change the
endpoint, the margin shrinks. Median ΔAUROC is +0.081 on AUROC, but
**+0.031 on partial AUC at FPR ≤ 0.02 and +0.041 on top-10 %
precision**, and the arm order is broadly stable (ρ = 0.80–0.90).
So the advantage is **smallest exactly where practitioners use the
method** — in the top of the list.

---

## 22:30–26:30 · 结果二：稳定性从哪里来，以及它为什么抹不掉

[`【图：slide 12】`]

The second result is the one the teacher asked us to make
first-class: **where does the stability of a network prior come from?**
We took the approach of Klingenhuber and colleagues — treat stability
as a first-class quantity, measure it, attribute it, and publish it.
>
The picture has two layers. At the **arm level** the mean margin is
stable — SD of ΔAUROC across the 200 evaluations is 0.0047. At the
**gene level** the picture is the opposite: the median per-gene score
CV is **18 %**, and a median gene moves about **270 ranks** between
repeats. Both are true at once, and that is the point: the endpoint
is stable, the individual gene calls are not.
>
[`【图：slide 12 面板 B】`]
>
Now the attribution, and here we had to **correct our own first
explanation**. Our initial hypothesis was crowding — that genes in
dense regions have unstable scores. The data do not support it: the
two crowding definitions have **opposite signs** (a ±10 % band gives
+0.468 on 24 of 24 arms, a Scott bandwidth gives −0.246).
>
What does explain it is the gene's **own score variance**: β = +0.616
with Spearman ρ = **1.00** and q = 2e-137, while **no classical
topology descriptor predicts it at all**. And once the gene's own
score CV is in the model, the degree coefficient collapses from −0.995
to −0.083 — which tells us the degree effect we first measured was a
**suppression artefact**, not a mechanism.
>
[`【数：网络层：稳定性由覆盖度与 2-core 比例预测（ρ = −0.905, q = 0.012），与规模/密度/模块度无关】`]
>
And the mechanism is a theorem, not a correlation: **RWR weights
seeds linearly**, so the score of a gene is a linear function of the
number of seeds near it; averaging the prior over bootstrapped seeds
— bagging — is therefore a provable no-op. **Regularising the prior
cannot restore the ranking.** We tested that (a stability-regularised
scorer that shrinks each gene's score toward a confidence-weighted
neighbourhood median) and it does not work, which is why we report
the negative result instead of hiding it. [`【数：I.6 负结果】`]
>
One more piece of independent evidence: the whole picture replicates
on **ClinGen** as an independent system, and the transfer experiment
reproduces the ranking only at the two extremes (ρ = +0.37).

---

## 26:30–32:00 · 结果三：外部证据，以及我们自己的数字有多脆

[`【图：slide 13 + 附录 A6–A9】`]

Third result, and it is two-sided: an independent line of support, and
an uncomfortable finding about our own numbers.
>
**The support: DepMap.** We took 1,361 cancer cell lines with their
CRISPR gene-effect scores and their expression, defined a candidate as
"essential in that cell line" (CERES below −0.5 and expressed above
log2(TPM+1) > 1), and asked, per cell line, whether the prior's
top-K candidates are enriched for essential genes relative to the
cell line's own gene set (hypergeometric test).
>
[`【数：K = 30 时 6 % 的细胞系 p < 0.05（随机水平 5 %，**没有富集**）；K = 100 时 48 %；
中位 p 从 0.52 变成 0.058】`]
>
So the honest reading is: **the prior carries set-level signal, not
gene-level calls.** At the scale where practitioners work — a handful
of candidates — it is indistinguishable from chance; at the scale of a
hundred, it is clearly informative. And the negative result at K = 30
is the more important half, because K ≈ 30 is roughly the number of
genuinely essential genes in one cell line.
>
**The second support: the NCG arbitration test.** This is the
mini-challenge in the Boutros sense. We seed the walk with our own
label's positives, and then ask: among the drivers **our labels do not
contain** — 27 % of the NCG known tier — does the prior rank them
above the rest of the universe? The answer is **24/24 arms positive**
(median +0.0413), and 26 of 26 cancer types positive. But on the
low-FPR endpoint it is only 13/24, and the split is systematic: the
four protein networks pass in every label while funmap and RNA
co-expression fail in all four. So the prior transfers to drivers we
never labelled — but the advantage in the very top of the list is a
property of the **interaction network**, not of the prior alone.
>
**Now the uncomfortable part: our own negative set is not clean.**
Using NCG as an arbiter, we measured that **10–22 % of our held-out
"negatives" are curated drivers we simply do not label**, and 1.3–9.5 %
are NCG-known high-confidence ones. Reactome's arms are the worst,
the dense protein networks the least.
>
[`【数：清洗后 IntOGen margin 腰斩 0.065→0.030；Open Targets margin 上升 0.085→0.103；
reactome × temporal 翻号 −0.0125→+0.0073】`]
>
That single table is why our abstract says the numbers depend on how
the negative set is built, and why we think **every benchmark in this
area should audit its own negative set**.
>
**Two more controls to close this section.** The **bait-matched null**
— matching each positive to a negative with the same publication
count — shows how much of the margin is visibility: the top-1000
enrichment falls from 17.6 % to 11.6 %, and on the temporal label it
falls to **4.7 %, which is chance (5 %)**. And the **label-quality
tiering** — ClinGen's review status, IntOGen's number of detecting
methods, NCG's literature count, and the literature/non-literature
split for Open Targets — shows the **high-confidence tier carries the
larger margin in four of five labels**, while for NCG the order
inverts, because NCG's tier *is* a publication count, so its high tier
is exactly the set the degree baseline already ranks well. Same
statement, two sides. And the **ClinGen 17 % noise anchor** — demoting
a random 17 % of positives in every fold — moves the margin by
−0.0018 to +0.0026, while a perfect scorer would lose 8.5 AUROC
points under the same misassignment. So the margin is **not** an
artefact of the label noise that ClinGen's own curation process
documents.

---

## 32:00–35:00 · 怎么解释 Reactome 失效，以及第五个标签告诉我们什么

[`【图：slide 16 / 附录 A10】`]

Let me close the loop on the single negative cell, because it is where
our interpretation is most interesting and where we were most careful
not to overclaim.
>
The **curator-lag** hypothesis says: a curated pathway database lags
the literature, so on the **newest** drivers it fails. Our temporal
label was designed as the test — 152 newly discovered drivers — and
Reactome is negative there (−0.0125) and on IntOGen 2024 (−0.0110).
>
But that evidence is entangled: the temporal label is also the label
most affected by publication-count effects, and its margin collapses
to +0.002 under the bait-matched null. So we needed a test that
separates **catalogue age** from **visibility**.
>
The fifth label gives exactly that: **NCG is a 2019 catalogue**. If
the failure were about the network, Reactome would fail on every
catalogue. If it is about the age of the catalogue, Reactome should
be fine on the older one.
>
[`【数：reactome × NCG 2019 = +0.0893（正）；reactome × IntOGen 2024 = −0.0110；
reactome × temporal = −0.0125】`]
>
**The Reactome failure is label-specific, not network-specific.** That
is the first test of the curator-lag hypothesis that does not rely on
the temporal split. We report it as *supported* rather than proven,
because a 2019 catalogue and an IntOGen 2024 catalogue differ in more
than age.
>
The fifth label also let us re-run the variance decomposition on a
30-cell grid. Note what we could and could not conclude: because NCG
is 65 % identical to IntOGen (Jaccard 0.522 — above the 50 %
threshold we had set in advance for calling two labels non-independent),
the six NCG arms are a near-replication, not new information. The
label main effect's share of the between-arm spread therefore falls
from 34.8 % to 30.8 %, the interaction rises from 40.1 % to 42.0 %,
and the within-cell noise rises from 40.4 % to 41.6 % of the total.
The honest reading: the label main effect is real but not dominant —
across five catalogues of different provenance, a quarter of the
between-arm spread is the network, a third is the label, and the
largest single piece is the interaction.

---

## 35:00–38:00 · 结论、局限与验证计划

[`【图：slide 14、slide 15】`]

Five conclusions.
>
**One**: a network prior contributes signal beyond node degree — 23 of
24 arms positive, the pre-registered rule met under four of five
interval constructions. **Two**: the contribution is real but modest
and label-dependent — the median margin is +0.081, and cleaning our own
negative set halves it on the strictest label. **Three**: much of the
published margin is study visibility — a bait-matched null removes
most of the enrichment and the temporal label's advantage collapses to
+0.002. **Four**: the prior's instability is intrinsic — the gene's
own score variance explains it, and bagging the prior is a provable
no-op, so no learner can remove it. **Five**: curated networks fail on
recently discovered drivers — Reactome is +0.089 on the 2019 catalogue
and negative on the two newest ones.
>
The limitations we know about, each with the experiment that would
settle it: the gene-cluster interval misses four arms — re-run on an
independent cell-line-derived label; our negative set is contaminated —
re-derive it with explicit NCG/ClinGen exclusion; the evaluation is
self-referential — run a blind expert-panel challenge; we have no
perturbation experiment, so we claim no causality — edit the seed set
and measure the effect on the predicted neighbourhood; and we test one
learner — re-run the locked protocol with a second fixed learner.
>
[`【图：slide 15（局限 → 验证实验）】`]

---

## 38:00–40:00 · 交付物、AI 使用声明与结束

[`【图：slide 16】`]

Everything is public and re-runnable: the code, the tables, the figures
and the split manifests are on GitHub at the v1.0-rc3 tag; the
pre-registration is frozen at Zenodo DOI 10.5281/zenodo.23121762 with
its amendments; two independent audit scripts pass 50/50 and 63/63
checks; and every figure in this talk is generated from a released
table by a script in the repository.
>
**AI use statement.** AI tools assisted with code drafting, debugging,
literature triage and slide layout. Every number in this talk comes
from our own scripts and tables, verified by the two audit scripts;
every figure is generated from a released table; **no simulated result
is presented as a measured one**. The scientific decisions — the
panel, the frozen learner, the endpoint, the decision rule, and the
interpretation of every negative cell — were made by the group and are
recorded in `docs/decisions.md`.
>
So the one sentence I would like you to take away is this: **network
priors do help, but the honest number is much smaller than the
headline one, and the difference is mostly about how often a gene has
been published.** Our benchmark, the pre-registration and the
amendments are public — please re-run it and check us.
>
Thank you. We would like to take your questions.

---

# 附录 A · 必须记住的数字（考前 10 分钟过一遍）

| 数字 | 含义 | 出处 |
|---|---|---|
| **+0.0809** | 24 臂 ΔAUROC 中位 | `a4_table1_delta_auroc_v40.tsv` |
| 23 / 24 | ΔAUROC 为正的臂 | 同上 |
| +0.1711 / −0.0125 | 最好 / 最差臂 | 同上 |
| **24/24 (×4) 与 20/24** | H1 在 5 种区间下的判定 | `stats_h1_verdict_v40.tsv` |
| 12× | gene-cluster 区间相对 repeat 区间的宽度 | amendment A2.3b |
| +0.031 / +0.041 | 低 FPR pAUC / top-10 % 的中位 margin | G4 |
| 0.0047 / 18 % / ~270 | 臂级 SD / 基因级 CV 中位 / 中位基因位移 | Phase S |
| **ρ = 1.00, q = 2e-137** | 自身 score CV 对不稳定性的解释 | `stability_attribution.tsv` |
| −0.995 → −0.083 | 加入 score CV 后 degree 系数的变化 | 同上 |
| ρ = −0.905, q = 0.012 | 网络级稳定性的预测量（覆盖度、2-core） | `network_stability.tsv` |
| 6 % / 48 % | DepMap K=30 / K=100 的细胞系富集率 | `depmap_validation.tsv` |
| **24/24, +0.0413** | NCG 仲裁测试的 ΔAUROC | `g7_ncg_arbitration.tsv` |
| 13/24 | 仲裁测试的低 FPR 端点 | 同上 |
| 26/26, +0.0271 | 癌种分层的中位 margin | `gr9_cancer_type_summary.tsv` |
| **10–22 %** | 阴性集中策展 driver 的比例 | `g8_negative_contamination.tsv` |
| 0.065→0.030 / 0.085→0.103 | 清洗背景后的 margin 变化 | `g8_cleaned_background.tsv` |
| **+0.0893 vs −0.0110** | reactome × NCG 2019 vs × IntOGen 2024 | `d11_table_v2.tsv` |
| 17.6 %→11.6 % / 6.8 %→4.7 % | bait 匹配前后的 top-1000 富集（随机 5 %） | `g3_bait_matched_null.tsv` |
| −0.0018..+0.0026 | ClinGen 17 % 噪声锚点下的 margin 变化 | `clingen_noise_anchor.tsv` |
| 34.8 %→30.8 % | label 主效应占比（24→30 臂） | `stats_variance_decomp_comparison.tsv` |
| Jaccard 0.522 | NCG 与 IntOGen 的重叠 | `c2_label_jaccard_matrix.tsv` |
| 70.2 vs 70.2 | NCG 与 IntOGen 在 intact 的平均度百分位 | `c2_ncg_circularity.tsv` |
| ρ = 0.47–0.72 | NCG 自身 ppin_degree 与我们网络度的相关 | 同上 |
| 134 / 50 / 63 | 测试通过数 / 复现审计 / 声明核对 | pytest / 两个审计脚本 |

---

# 附录 B · 最可能被问到的十个问题（一句话版）

1. **这是新发现还是边界刻画？** —— 边界刻画；但每条结论都有可证伪的实验与具体会变的表。
2. **为什么 RWR 不用 GNN？** —— learner 必须冻结，否则测的是"我们调出了好方法"；我们做了 GNN 礼节性对照（RWR +0.079 vs 嵌入 +0.028）。
3. **阴性集干净吗？** —— 不干净，10–22 % 是策展 driver；这是我们自曝的，且已给清洗后的两套数。
4. **H1 到底成立吗？** —— 5 种区间里 4 种 24/24，gene-cluster 20/24，四个例外全是 temporal 标签。
5. **会不会是"网络越大越好"？** —— D-11 cap 与 D-12 共享宇宙把这件事压下来了（−0.0360），但覆盖度确实解释一部分。
6. **你的 ΔAUROC 能读成"多找出 8 % 的 driver"吗？** —— 不能，低 FPR 端点只有 +0.031。
7. **为什么用 M1/M2/M3/M4 而不是同一个标签？** —— 因为策展口味不同；跨癌种 26/26 为正说明不是某一份目录的怪癖。
8. **能不能说因果？** —— 不能，没有扰动实验；验证计划里有种子扰动实验。
9. **AI 会不会编数字？** —— 全部数字来自已发布表格 + 两个审计脚本（50/50、63/63）；PPT 的数字由脚本从表格导出。
10. **如果核心差异不成立会怎样？** —— 见题库 C2：三个假说各自对应哪张表会变。

---

# 附录 C · 时间压缩方案（若只剩 15 分钟）

| 保留 | 时长 | 理由 |
|---|---|---|
| 开场 + 研究地图 | 1.5 min | 建立问题 |
| 三篇来源论文（压缩成"一问/一数/一局限"各一句） | 3 min | 评分要求 |
| 知识空白 + 三个假说 | 1.5 min | 评分要求 |
| 数据与设计（合并，只讲 provenance 与 D-02） | 3 min | 评分要求 |
| 主端点 + H1 的准确表述 | 2 min | 核心结果 |
| 稳定性机制 | 1.5 min | 老师点名 |
| 外部证据 + 阴性集污染 | 1.5 min | 诚实性得分点 |
| 结论 + 局限 + 交付物 | 1 min | 收尾 |

**砍掉**：控制栈逐条（只说一句"11 层，每层一个数"）、第五个标签的方差分解细节、
参数表的逐行（只在被问时展开）。**绝不能砍**：H1 的准确表述、阴性集污染、稳定性机制。

---

# 附录 D · 彩排清单

- [ ] 第一遍：按本稿计时，目标 **14 分钟**（40 分钟版本讲不满会显得内容不足，15 分钟版本讲超时会被打断）。
- [ ] 第二遍：合上稿子讲，只看 `report_1118_facts.md` 的数字表。
- [ ] 第三遍：三人一起过 `report_1118_qa_bank.md`，每人抽 3 题互问。
- [ ] 检查：PPT 在**另一台机器**上能打开（字体 Times New Roman 已嵌入要求；图片是 PNG，不依赖外部链接）。
- [ ] 检查：讲稿里的每个数字都能在 `report_1118_facts.md` 里找到出处（口头数字错一个就全盘不可信）。
- [ ] 备份：`D1_defence_deck.pptx` 与本稿各存一份到 U 盘。
