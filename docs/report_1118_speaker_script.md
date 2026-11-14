# 1118 终期答辩 · 主讲稿（单人 40 分钟版）

**讲者**：本人（唯一主讲，讲满 40 分钟）· 其余三人在讲完后回答问题
**配套**：`D1_defence_deck.pptx`（正文 16 页 + 附录 10 页，附录页在讲稿里被用到）、
`report_1118_qa_bank.md`（29 题）、`report_1118_facts.md`（全部数字的单一来源）、
`datasets_sources.md`（数据来源与登录号）

**这份稿子的写法**
- 每节开头：**中文提示**（时间点、这一节要落到的点、指哪页）。
- 正文：**逐字英文口播**。凡出现数字，口播里都必须交代四件事——
  **它是什么**（量纲/分母）、**怎么算出来的**（哪张表、哪个比较）、
  **它说明什么**、**它不能说明什么**（这一条最容易丢分，也最容易显得专业）。
- `【数：…】` 是**只给你自己看的提醒**（数字 + 出处），不要念。
- 每节末尾有一句**过渡句**，用来接下一节。
- 附录 A 是**必背数字卡**（23 张），附录 C 是**15 分钟压缩方案**，附录 D 是**彩排清单**。
- **语速**：正文英文口播共约 **6,840 词**，按 140–170 词/分是 **40–49 分钟**；配上下面的
  逐页配时表走，正常语速下正好落在 40 分钟上下。讲快了就用附录 C 砍。


---

## 逐页配时表（正文 16 页 = 30.5 分钟 + 10 分钟的"深讲"，合计 40 分钟）

括号里是**该页要用的分钟数**；附录页（17–26）**不在讲稿里主动讲**，它们是回答追问时的
证据页（可随时跳转）。40 分钟的主讲按下面的配时走；**如果只剩 15 分钟，见附录 C**。

| 时间 | 页 | 讲什么 | 分钟 |
|---|---|---|---|
| 0:00 | 1 | 标题、四个成员与分工 | 0.3 |
| 0:20 | 2 | 研究地图：问题 → 三个假说 → 设计 → 结论 | 1.2 |
| 1:30 | 3 | Bailey 2018（复现共识） | 1.2 |
| 2:45 | 4 | Shi 2025（FunMap，含 LLR/MR 方法细节与 Reactome 金标准） | 1.5 |
| 4:15 | 5 | Repana 2019（NCG，三层与网络置信度） | 1.3 |
| 5:40 | — | 五篇方法学文献 → 十一个设计决定（可指附录 A1） | 2.3 |
| 8:00 | 6 | 知识空白 + 三个可证伪假说 | 1.5 |
| 9:30 | 7 + A2 | 数据从哪来：6 个网络、5 个标签、D-02、cap 顺序、复合物映射、每个输入的版本与标识符 | 6.5 |
| 16:00 | 8 | 流程与复现：冻结 learner、40×5、切分发布、预注册与 13 个 P0 自审修复 | 2.0 |
| 18:00 | 9 | 参数与稳健性：为什么锁、怎么检验、诚实声明 | 2.0 |
| 20:00 | 10 | 十一层控制栈（每层一句 + 一个数） | 2.0 |
| 22:00 | 11 + A3–A5 | 主端点：23/24、H1 在五种区间下的精确表述、低 FPR 端点 | 4.0 |
| 26:00 | 12 | 稳定性：两层结构、拥挤假说的自我更正、种子数定律 | 5.0 |
| 31:00 | 13 + A6–A10 | DepMap、NCG 仲裁、阴性集污染、bait 匹配、标签质量、第五标签与 curator lag | 7.0 |
| 38:00 | 14 | 五条结论 | 1.0 |
| 39:00 | 15 | 局限 → 验证实验 | 0.8 |
| 39:50 | 16 | 交付物、AI 声明、结束语 | 0.7 |

---

## 0:00–1:30 · 开场（正文 slide 1–2）

【提示：40 秒。不要报任何数字，先把问题立住。】

Good morning. Our project asks a question that sounds simple and that
the field has never actually answered: **when a network prior is used
to prioritise cancer-driver genes, how much of the reported gain comes
from the network's topology, how much from plain node degree, and how
much from the fact that well-studied genes have been published about
hundreds of times?** [`【图：slide 2 研究地图】`]
>
The motivation is not scepticism for its own sake. Network-based
prioritisation is standard practice — it is the formalisation of
guilt-by-association — and it works in the sense that random walk
with restart usually beats a degree-only ranking. What has never been
done is to **decompose** that gain.
>
There are three concrete reasons to think the decomposition matters,
and each one is a published result. **First**, driver genes are hubs,
but the enrichment is largely an artefact of research intensity: bait
usage — how often a gene appears in the published interaction
literature — correlates with its degree at **r = 0.52**, and correcting
for it removes the cancer-versus-control degree difference entirely.
That is Schaefer and colleagues, 2015, and it is the single most
important paper for us. **Second**, performance tracks network size:
across more than **45 human interactomes**, accuracy follows the size
of the network, and the effect largely collapses once gene-set
coverage is controlled. That is Wright and colleagues, 2025. **Third**,
standard cross-validation is optimistic when the same protein
complexes appear on both sides of the split: holding a complex
together in the fold drops a strong network classifier from about **12**
hits in the top 20 to **fewer than 4.5** — a fall as large as the gap
between competing methods and the baseline. That is Picart-Armada and
colleagues, 2019.
>
So we built a **pre-registered, leakage-controlled benchmark on public
data**, and in the next forty minutes I will show you the design, the
control stack, three result slides, and — just as important — the
limits of what we are allowed to claim.

---

## 1:30–8:00 · 我们依据的文献（正文 slide 3–5 + 附录 A1）

【提示：这是给生信出身的老师听的段落。每篇论文讲四件事：**问了什么 / 用了什么数据 / 结论是什么 /
留下什么**，并说清**我们因此做了什么设计决定**。数字必须带量纲和分母。不要报 DOI 串，报逻辑。
附录 A1 把这八篇与十一个设计决定一一对应，讲到这里可以指一下。】

### 1.1 Bailey et al. 2018, Cell 173:371–385 — 复现共识

[`【图：slide 3】` `【数：9,423 对肿瘤/正常；26 个组织类型；26 种突变检测方法；我们的 intogen2024 = 633 个基因】`]

The first source paper is Bailey and the MC3 Working Group, in Cell.
Their question was: **which genes are mutated recurrently enough,
across tumour types, to be called drivers?** The data were **9,423
TCGA tumour-normal pairs** across **26 tissue types**. The
methodological core is the **MC3 consensus**: **26 independent
mutation-calling algorithms vote**, and a gene becomes a driver when
the methods agree — not when one method is confident. The result
became the field's default gold standard, the Cancer Gene Census.
>
The numbers I want you to hold: **9,423 samples** is a large
cohort, which is why the consensus is trustworthy; **26 voters** is
what makes it a consensus rather than another predictor; and the
mutation-recurrence panel we inherit from this lineage has **633
positive genes** in our build.
>
What it leaves open, and this is the gap we work in: it is **purely
recurrence-based**. There is no network, no prior, and no way to
separate "what a network adds" from "the catalogue's own degree
bias". And note the consequence for us: because we use this catalogue
as our primary label, **we inherit its curation choices** — which is
why we measure label quality rather than assume it (I will show you
that measurement in result three).

### 1.2 Shi et al. 2025, Nature Cancer 6:205–222 — FunMap

[`【图：slide 4】` `【数：XGBoost 逐对共功能概率；集合级 LLR 3.912（= LR 50）；mutual rank > 0.95；10,525 节点 / 196,800 边；11 个癌种；≤200 蛋白通路规则；FunMap 的金标准来自 Reactome】`]

The second source paper is Shi and colleagues' FunMap. Their question:
**can a pan-cancer proteogenomic functional network prioritise cancer
genes better than the existing interactomes?** The data are the
**CPTAC pan-cancer proteogenomics and transcriptomics**.
>
The method deserves thirty seconds of precision, because it changes how
we treat the network. **FunMap's edges are not raw interactions.**
Every gene pair gets an XGBoost-predicted co-function probability;
the pairs are then added **in blocks** until the **set-level
log-likelihood ratio falls below 3.912 — that is, a likelihood ratio
of 50**, meaning "these two genes are 50 times more likely to be
functionally associated than not"; and before scoring, a hard filter
requires the pair's **mutual rank to exceed 0.95 in at least one
cohort**. That threshold is what yields **10,525 nodes and 196,800
edges**. Two consequences for us: the published edge list is a
**derived, supervised, thresholded object**, not an interaction table,
so we treat it as a prior network and we report its coverage; and the
network is built from **11 cancer types** in **bulk tissue**, so tumour
cells and microenvironment are mixed inside one graph.
>
The result: FunMap outperforms STRING for cancer-gene prioritisation,
and it ships with a gold standard. And here is the limitation the
authors state themselves, which is the one we care about most:
**their gold standard is Reactome-derived** — a gene pair counts as
positive if the two genes share a "detailed" Reactome pathway of at
most **200 proteins** — and Reactome is *also* an input to STRING's
`database` channel. So the network inherits Reactome's notion of
function, and the comparison networks "incorporate the knowledge used
for evaluation".
>
**Our design decisions, three of them.** First, FunMap is one of our
six networks because it is the most on-topic cancer-specific prior,
and we report its coverage explicitly: its **10,525** genes are the
protein-coding genes CPTAC actually quantified, so a
largest-connected-component analysis runs on a small,
expression-biased, tumour-derived universe. Second, we **audit
circularity per network** rather than in general. Third, because the
gold standard is Reactome-derived, we do **not** treat "FunMap beats
STRING" as independent evidence — we treat it as one more reason the
whole field needs the decomposition we are doing.

### 1.3 Repana et al. 2019, Genome Biology 20:1 — NCG

[`【图：slide 5】` `【数：2,372 个基因；708 known / 104 candidate / 1,560 screening；置信度来自 HuRI、GTRD、Reactome、STRING 实验通道；我们的 ncg_known = 708】`]

The third source paper is Repana and colleagues' Network of Cancer
Genes. Their question: **which genes have been called drivers across
the published cancer sequencing literature, and how confident is each
call?** They curated **2,372 genes** from that literature into three
tiers: **708 known** tumour suppressors and oncogenes, **104
multi-screen candidates**, and **1,560 screening-level genes**, and
attached to each gene a confidence level computed from curated
interaction channels — HuRI for binary interactions, GTRD for
transcription-factor binding, Reactome for pathways, STRING's
experimental channel.
>
What it leaves open, and it is a sharp point: **that confidence level
is itself network-derived**, and every gene the catalogue does not
mention is implicitly a negative. So NCG cannot measure anything by
itself — it can only arbitrate.
>
**Our decision**, which is the one I would defend hardest: we use
NCG as an **arbitration catalogue and as an audit instrument**, never
as a confirmatory label. Two things follow. First, because it shares
no curation with our four label sets, it can arbitrate them — I will
show you that test, and it is our strongest positive result. Second,
because its confidence level is network-derived, it is exactly the
right instrument for the **circularity audit**: it tells us whether
our network panel and our label panel overlap in their evidence.
And note the arithmetic that makes this work: **708 known genes** in
the catalogue against **633 positives** in our primary label — so
they are not the same set, and that is the whole point.

### 1.4 五篇方法学文献：每一条都变成一个设计决定

[`【图：附录 A1（文献 → 决定对照表）】`]

The three source papers give us the inputs. **Five** further papers
give us the **design**, and I want to be explicit about which decision
each one produced, because each of them is a place where we could
have gone wrong.

**Picart-Armada et al. 2019, PLoS Comput Biol 15:e1007276.**
They benchmarked network propagation for disease-gene identification
and found three things we adopted. First, **pure topology baselines
against random** have to be in the benchmark, not only degree.
Second, the **complex leakage** number I quoted: holding a complex
together in the fold drops a strong network classifier from about
**12** hits in the top 20 to **fewer than 4.5**. Third, **the endpoint
changes the ranking of methods** — AUROC and top-k disagree.
**Our decisions**: complex-aware cross-validation, where the whole
protein complex stays in one fold; and **two extra endpoints** —
partial AUC at FPR ≤ 0.02 and top-10 % precision — reported next to
AUROC, so that "the prior wins" cannot be an artefact of the metric.

**Guala & Sonnhammer 2017, Scientific Reports 7:46598.**
This is the paper that made us most careful, and I want to state its
result plainly because it is the sharpest warning in this field:
**MaxLink — a local method with an explicit degree filter, built so
that highly connected genes cannot rank high just by being connected —
beats the diffusion (RWR-family) methods on the full list.** Diffusion
wins only in the **top 1–10 %**. They also use partial AUC at
**FPR ≤ 0.02**, which is where we took our low-FPR endpoint from, and
they warn about gold standards that share evidence with the network.
**Our decisions**: we never claim RWR is the best method — we freeze
it and measure the **network's contribution relative to degree**,
which is precisely the quantity Guala shows can be small; and we
report the top-1–10 % region separately, because that is where
diffusion actually wins.

**Hill et al. 2019.**
They benchmarked **17 network-contextualisation algorithms** and found
that over **10,000 random gene lists**, the algorithms differ
enormously in how much of their output is **input-independent**. In
other words: a learner that ignores its input still returns the hubs.
**Our decision**: a **random-input degeneracy control** — feed the
learner shuffled inputs, and it must degrade to the baseline.

**Wright et al. 2025.**
Over **45 human interactomes**, performance tracks network size, and
the effect largely collapses once gene-set coverage is controlled.
**Our decisions**: the **shared-gene-universe control** (all arms
evaluated on the same 1,605 genes) and the **per-network common-size
cap**, so that a margin cannot be bought with a bigger gene universe.

**Schaefer et al. 2015, Frontiers in Genetics.**
This is the paper our third hypothesis is built on, and I will come
back to it in the results because **the number r = 0.52** is the hinge
of the whole talk: **bait usage correlates with degree, and correcting
for it removes the cancer-versus-control degree difference.** What
looks like a biological degree effect can be a publication effect.
**Our decisions**: a **publication-count-adjusted label** and a
**bait-matched null model** in which each positive is paired with a
negative carrying the same publication count.

**Boutros et al. 2014, Genome Biology 15:462.**
Their principle is the one we borrowed for the whole design: **the
gold standard dominates any ranking, so the evaluator must be blind
to the algorithm**; a benchmark should be a community challenge with
hidden answers, not a leaderboard people can overfit.
**Our decisions**: the **NCG arbitration test** is our mini-version
of that — an independent catalogue none of our labels was built from
— and it is also the experiment we name first in our validation plan.


---

## 8:00–9:30 · 知识空白与三个可证伪假说（正文 slide 6）

[`【图：slide 6】` `【数：三个假说各有一条"如何被证伪"的判据】`]

Putting the literature together, the gap is precise and it is not a
generic "more research is needed". **No published study separates the
topology contribution from the degree contribution from the
study-visibility contribution inside one frozen design with one
decision rule.** Each piece exists separately — Zhang's degree
stratification, Schaefer's bait correction, Picart-Armada's leakage,
Wright's coverage — but never together, and never on a pre-registered
endpoint with the resampling unit fixed in advance.
>
So we made three hypotheses, and for each one we wrote down in advance
how it could fail.
>
**H-vis**: much of the reported margin is **study visibility**, not
topology. *Falsified* if the margin survives a null model that matches
positives and negatives on their publication count.
>
**H-stable**: the prior's instability is a property of the **gene's
own score variance**, not of topology, and it cannot be smoothed away.
*Falsified* if any topology descriptor predicts instability, or if a
regularised scorer restores the ranking.
>
**H-lag**: curated networks fail specifically on **recently
discovered** drivers. *Falsified* if the failure appears equally on
older catalogues.
>
Each of the three has a number attached, and I will now show you the
numbers.

---

## 9:30–16:00 · 数据从哪来（正文 slide 7 + 附录 A2）

[`【图：slide 7（面板与标签）→ 附录 A2（每个输入的版本与标识符）】`]

**The network panel: six networks, three provenance classes.**
Provenance is an experimental factor here, not a nuisance to be
removed. Three are protein-derived, one is RNA-derived, two are
curated or literature-derived.
>
**funmap**, the CPTAC proteogenomic network: **10,189 nodes and
196,605 edges** in the largest connected component. **intact**,
IntAct human-human interactions at taxid 9606: **17,918 nodes and
564,706 edges** — by far the densest, and remember that density is
what the coverage control is about. **string_phys700**, the STRING v12
physical subnetwork at combined score ≥ 700: **9,830 nodes and 85,576
edges**. **rna_coexp**, which we built ourselves from FunMap's
published CPTAC mRNA matrices: Pearson correlation per cohort, Fisher-z
averaged across cohorts, and then the top **196,172** pairs by mean z
to density-match funmap — note that we matched the **density**, not
the size, so the RNA network has comparable edge density but its own
gene universe: **9,152 nodes**. **string_full700**, STRING v12 all
evidence channels at ≥ 700: **15,882 nodes and 236,712 edges**.
**reactome**, version 86, human only: **3,953 nodes and 20,143 edges**
— the smallest network by an order of magnitude, and that turns out to
matter.
>
**Two harmonisation decisions.** First, **Reactome needs a pathway
size cap**. Reactome provides **no numeric edge weight at all** — its
only quantitative layer is the curation star rating — and the
signalling and transcription pathways are so large that a naive
co-membership graph is nearly complete inside them, so its connected
component carries almost no information as a "prior". We therefore
follow the **≤ 200-protein** rule that FunMap's own gold standard
uses, and we snapshot a **specific** Reactome version rather than
always taking "current" — the 3- and 4-star curation events appeared
in June 2023, and historical tools drift.
Second, the primary analysis is **unweighted** — decision D-01 — so
what we measure is topology, not the weight scale.
>
**The label panel: five catalogues, four curation philosophies.**
**intogen2024**, **633** positive genes, mutation recurrence from
IntOGen — and this is the most defensible pure-positive source we
have, because a positive is defined by positive-selection evidence in
mutation data, not by literature co-occurrence or pathway
co-membership, which is exactly the independence a network-prior
benchmark needs. **ot_all** and **ot_nolit**, **600** genes each,
Open Targets 26.06, the second with the literature channel removed;
the difference between them is control 7 in label form.
**intogen_temporal_new**, **152** genes — the newest IntOGen drivers,
a deliberately harder label. **clingen_label**, **84** genes, clinical
gene-disease validity: **79 Definitive, 1 Strong and 4 Moderate** —
the hardest evidence of all, and we use it twice, as an extra label
and as the external noise anchor. **ncg_known**, **708** genes, the
NCG known tier, added later as a pre-planned exploratory extension.
>
**Two curation caveats we had to handle.** IntOGen is a **gene ×
cancer type** structure — 99 breast drivers mixed with 20
neuroblastoma drivers — so we stratify by cancer type rather than
treating it as a flat list; that stratification is what gives us
**26 cancer types with a positive median margin** later on. And
**ranking a label by Cancer Gene Census membership would measure
literature curation intensity, not driver prediction**, because 75 %
of the IntOGen compendium is CGC-annotated — so we never use CGC
membership as a ranking tier.
>
**Symbol alignment and provenance.** Everything maps to HGNC approved
symbols under decision D-02: the approved symbol maps to itself, a
previous symbol or alias maps only when the mapping is one-to-one, and
ambiguous symbols are **dropped and counted** — for NCG that was
**83** symbols, and we keep them visible in the release rather than
discarding them silently. Nine NCG genes were renamed outright —
CARS → CARS1, MKL1 → MRTFA, FGFR1OP → CEP43 and six others — and we
carry the old names in the label file for traceability.
>
Two more provenance items. The **protein-complex map** for the
complex-aware splits was built from the UniProt REST API by
cross-referencing ComplexPortal identifiers: **2,569 complexes** over
**3,708 genes**. And every external input carries an explicit release
and identifier, which is in appendix A2: the pre-registration is
frozen at Zenodo DOI **10.5281/zenodo.23121762**, the FunMap matrices
are Zenodo record **7948944**, IntOGen releases **2020, 2023 and
2024**, Open Targets **26.06**, STRING **v12.0**, Reactome **v86**,
IntAct **2024-09**, the ClinGen file of **2026-09-26**, and DepMap's
CRISPR gene-effect matrix in two public releases.

---

## 16:00–20:00 · 设计与纪律（正文 slide 8–9）

[`【图：slide 8（流程）→ slide 9（参数表）】` `【数：α=0.5；5 折 × 40 repeat = 200 次/臂；6×4×3×200；覆盖模拟 1.000 vs 0.936–0.944；13 个 P0；10 次 repeat 与原文件逐字节相同】`]

Now the design, and I want to spend a few minutes on the discipline,
because that is what makes the result interpretable rather than
merely interesting.
>
**The learner is frozen.** RWR at **α = 0.5**, no tuning, the same
value in all 24 arms and all controls. The reason is not that RWR is
best — Guala just showed a degree-filtered local method can beat it
on a full list. The reason is that once the learner is tunable, the
quantity of interest silently becomes "we tuned a method" instead of
"the network contributes". We ran a **GNN courtesy arm** anyway: with
p = q = 1, node2vec is exactly DeepWalk-style embeddings, so we report
it under that name — RWR **+0.079** against **+0.028** for the
embedding, both as margins over degree.
>
**The resampling unit is the repeat, not the fold.** 5 folds × 40
repeats = **200 evaluations per arm**; 6 networks × 4 labels × 3
alphas × 200 evaluations is the whole primary design. We moved from
10 to 40 repeats because the pre-registered interval must resample
repeats, and a coverage simulation showed why: at 40 repeats the
fold-level percentile interval has coverage **1.000** — it is
over-conservative — while the repeat-level percentile is calibrated at
**0.936–0.944**. And the design is a strict superset: the first 10
repeats of the 40-repeat split files are **byte-identical** to the
original 10-repeat files.
>
**Splits are fixed and released** with SHA-256 manifests, so another
group re-runs on the same folds.
>
**The common-size cap is a control, not the primary design.** The
primary arms use each label's own positive set; the cap is what
equalises prevalence inside an arm. One detail that matters and that
we got wrong first: the cap must be applied **after** restricting to
the network's universe, otherwise funmap holds **78** positives
instead of **104** — we found that as a defect and recorded it as
amendment A2.3i.
>
**Everything is pre-registered, and every departure is written
down.** The pre-registration was frozen on **2026-09-27**; we now
have amendments **A2.3a to A2.3k** plus a register of the **14**
analysis families, each marked as belonging to the verification family
or the exploratory family.
>
**The statistical layer.** A mixed-effects variance decomposition; a
**gene-cluster bootstrap** — the pre-registered §9.3 unit, which
resamples **genes**, not folds; a **TOST** equivalence test for the
borderline arm; and **BH-FDR inside the pre-declared family only** —
we do not run BH across 72 correlated cells and call that multiplicity
control.
>
**And we should mention the self-audit**, because a benchmark that
reports its own defects is more credible than one that does not. An
independent review of this repository found **ten P0-class defects**:
a bootstrap that accepted an alpha parameter and never used it, a
TOST that was not TOST, a cap file whose name collided so each
network silently overwrote the previous cap, a shared-universe control
that did not use the shared universe, a temporal control with zero
replication, a variance decomposition that treated folds as a noise
bucket, a released results table with 4 rows instead of 24, and a
manuscript that quietly relaxed the H1 decision rule. All are fixed,
and all are recorded in the amendments.

---

## 20:00–22:00 · 十一层控制栈（正文 slide 10）

[`【图：slide 10】` `【数：11 层；每一层一个数】`]

Eleven controls, one artefact each. In order, and I will give the
number for each as we go: **degree-matched random seeds** — does the
prior only re-rank hubs? **Degree-preserving rewiring** — is the
signal carried by the wiring? **Shared gene universe** — gene-set
size differences? **Temporal split** — is the effect only in the
past? **Common-size cap** — unequal positive-set sizes?
**Complex-aware cross-validation** — complex leakage? **Hub-capped
RWR** — is the gain a hub effect? **Publication-count adjusted
label** — study visibility? **Bait-matched null** — does the gain
track publication count? **Label-quality tiers** — does the gain live
in the weak positives? **ClinGen 17 % noise anchor** — is the gain an
artefact of label noise?
>
Every one of these has a number, and four of them changed what we are
allowed to conclude.


---

## 22:00–26:00 · 结果一：主端点（正文 slide 11 + 附录 A3–A5）

[`【图：slide 11（30 格热图）→ 附录 A3（24 臂全表）→ 附录 A4（H1 五种区间）→ 附录 A5（低 FPR 端点）】`]

The primary endpoint is the margin of RWR over a degree-only ranking
on identical folds: **ΔAUROC = AUROC(RWR) − AUROC(degree)**, at
α = 0.5, with the degree baseline computed **on the same folds, in
the same universe**, so the only thing that differs between the two
numbers is the scorer.
>
**The headline numbers.** Across the **24 pre-registered cells** — six
networks times four labels — the margin is positive in **23**, and the
**median is +0.0809**. The best cell is STRING physical against
Open Targets without literature, at **+0.1711**; the worst is Reactome
against the newest IntOGen drivers, at **−0.0125**. Each cell is the
mean over **200 evaluations** — 40 repeats times 5 folds — and each
has its own standard deviation, which is in appendix A3.
>
**What +0.0809 does and does not mean.** It is a difference of two
AUROC values on the same test folds. It does **not** mean we find 8 %
more drivers; an AUROC difference of 0.08 is a modest, real
discrimination difference. And it does **not** mean the prior is
*better* than degree — it means the *information the prior adds
beyond degree* is worth about 0.08 AUROC on average. Which is why
the next two slides matter more than this one.
>
**The one negative cell, stated honestly.** Reactome against the
newest drivers is **−0.0125**. We do **not** claim this is a negative
effect. A confidence interval excluding zero says an effect is
*distinguishable* from zero; it does not say the effect is *real*.
The equivalence test answers the second question, and with a
corrected TOST the answer is that this cell is **not distinguishable
from "no effect"**. So the correct sentence is: **on this arm there is
no detectable contribution.** We had this inference wrong in an early
draft of the amendment and corrected it in A2.3c.
>
**Now the H1 verdict, exactly as the pre-registration defines it.**
The rule: the 95 % confidence interval excludes zero **in every arm**.
It is met under **four of the five** interval constructions —
fold-percentile, repeat-percentile, repeat-t and repeat-BCa all give
**24 of 24** — and it is **not** met under the fifth, the
gene-cluster bootstrap, which gives **20 of 24**. The four exceptions
are **all** the newest-driver label, the one with the smallest
effects. Why is that interval wider — **about 12 times** wider than
the repeat-level one? Because it asks a different question: resampling
genes asks "what if a different set of genes had been measured",
while resampling repeats asks "what if the same universe were
re-partitioned". Appendix A4 has the five rows.
>
I want to be careful here, because this is precisely where a paper
quietly relaxes a rule, and our own review flagged exactly that in our
first draft. We did not relax it. If you force one hard sentence:
**in 20 of the 24 arms the increment over degree excludes zero at
95 % confidence; in the remaining four, the effect is too small for us
to distinguish it from zero.** That is what "label-dependent" means
literally, not rhetorically.
>
**And if we change the endpoint, the margin shrinks.** Median
**+0.0809** on AUROC, but **+0.0308** on partial AUC at FPR ≤ 0.02
and **+0.0412** on top-10 % precision, and the number of arms with a
positive margin falls from 23 to **20** for both. The arm *order* is
broadly stable, with Spearman ρ between **0.80 and 0.90** depending
on the endpoint pair — so the ranking of networks is robust, but the
size of the advantage is not, and it is **smallest exactly where
practitioners use the method**, in the top of the list. That is the
honest reading, and it is in appendix A5.

---

## 26:00–31:00 · 结果二：稳定性从哪来，为什么抹不掉（正文 slide 12）

[`【图：slide 12（面板 A 网络级 + 面板 B 基因级）】`]

The second result is the one our teacher asked us to make
first-class: **where does the stability of a network prior come from?**
We took the approach of Klingelhuber and colleagues — treat stability
as a first-class quantity, measure it, attribute it, and publish it.
>
**The two layers, which are both true at once.** At the **arm level**
the mean margin is stable: the standard deviation of ΔAUROC across
the 200 evaluations of an arm is **0.0047**. At the **gene level** the
picture is the opposite: the median per-gene score coefficient of
variation is **18 %**, and a median gene moves about **270 ranks**
between repeats. These are not in conflict — they say the *endpoint*
is stable while the *individual gene calls* are not, and that is
exactly what a benchmark user should be told.
>
**The attribution — and the correction of our own first explanation.**
Our initial hypothesis was **crowding**: genes in dense regions have
unstable scores. **The data do not support it.** The two crowding
definitions have **opposite signs** — a ±10 % band gives **+0.478** in 24 of 24 arms, while a Scott bandwidth gives **-0.175** (negative in 23 of 24 arms). A
hypothesis with opposite signs under two reasonable definitions is
not a mechanism.
>
**What does explain it is the gene's own score variance**: in the
gene-level model, the coefficient on the gene's own score CV is
**β = +0.616**, with Spearman **ρ = 1.00** and **q = 2e-137** — and
**no classical topology descriptor predicts the instability at all**:
degree, k-core, clustering, betweenness, PageRank and community size
are all non-significant.
>
**And here is the diagnostic that convinced us the degree effect we
first measured was an artefact.** Put the gene's own score CV into the
model and degree's *marginal* association is **+0.038** — that is, nothing — while its coefficient in the model is **−0.083**. That is the signature of **suppression**: the topology
appeared to matter only because it was standing in for the score
variance.
>
**At the network level, something does predict stability**: coverage
and the 2-core fraction, at Spearman **ρ = −0.905**, **q = 0.012**;
network size, density and modularity do **not** predict it. So a
network is stable when it covers the label well and has a
well-connected core — which is a statement about curation, not about
graph theory.
>
**The mechanism is a theorem, not a correlation.** RWR weights seeds
**linearly**, so a gene's score is a linear function of how many
seeds sit near it. Averaging the prior over bootstrapped seeds —
**bagging** — is therefore a **provable no-op**. We tested the
consequence: a stability-regularised scorer that shrinks each gene's
score toward a confidence-weighted neighbourhood median **cannot
restore the ranking**. That is why we report it as a **negative
result** in section I.6 of the manuscript instead of hiding it.
>
**Two independent checks.** The whole picture replicates on **ClinGen**
as an independent system, and in a transfer experiment the ranking
reproduces only at the two extremes, with **ρ = +0.37** — which is
itself the honest statement of how much of this transfers.

---

## 31:00–38:00 · 结果三：外部证据与我们自己数字的脆弱性（正文 slide 13 + 附录 A6–A10）

[`【图：slide 13（DepMap 图 + 右栏要点）→ 附录 A6（bait 匹配）→ A7（阴性集审计）→ A8（等阳性率面板）→ A9（标签质量）→ A10（方差分解）】`]

The third result is deliberately two-sided: an independent line of
support, and an uncomfortable finding about our own numbers.
>
**The support, first: DepMap.** We took **1,361 cancer cell lines**
with their CRISPR gene-effect scores and their expression. A gene
counts as essential in a cell line if its CERES score is below **−0.5**
and it is expressed above **log2(TPM+1) > 1**. Then, per cell line,
we asked whether the prior's top-K candidates are enriched for
essential genes relative to that cell line's own gene set, with a
hypergeometric test.
>
**The numbers.** At **K = 30** candidates — and 30 is roughly the
number of genuinely essential genes in one cell line — **6 % of cell
lines** reach p < 0.05. Chance is **5 %**. So at that scale the
prior is **indistinguishable from chance**: the median p-value is
**0.52**. At **K = 100**, **48 %** of cell lines reach p < 0.05 and
the median p is **0.058**. The result is robust to the CERES
threshold — **−0.3, −0.5 and −0.8** give the same picture.
>
**What that means, and what it does not.** It means the prior carries
**set-level signal, not gene-level calls**: it can tell you something
about the composition of a hundred candidates, and nothing reliable
about which single candidate is a driver. The negative half — nothing
at K = 30 — is the more important half, because that is the scale
people actually work at.
>
**The second support: the NCG arbitration test.** This is the
mini-challenge in Boutros's sense. We seed the walk with our own
label's positives, and then ask: among the drivers **our labels do not
contain** — **27 %** of the NCG known tier — does the prior rank them
above the rest of the universe? The answer is **24 of 24 arms
positive**, median **+0.0413**. And stratified by cancer type,
**26 of 26 cancer types** have a positive median margin, **+0.0271**
over the types. So the prior transfers to drivers we never labelled,
and the advantage is not a quirk of one tumour class.
>
**But** on the low-FPR endpoint it is only **13 of 24**, and the
split is systematic: the four protein networks pass in every label,
while funmap and RNA co-expression fail in all four. So the advantage
*in the very top of the list* is a property of the **interaction
network**, not of the prior alone. One technical note so nobody
questions it: the background was sub-sampled to **20,000** genes for
tractability, and re-computing one arm with **200,000** and with the
**full** background gives **exactly** the same value, +0.1269 — so
the sub-sampling is immaterial.
>
**Now the uncomfortable part: our own negative set is not clean.**
Using NCG as an arbiter, we measured that **10 % to 22 %** of our
held-out "negatives" are curated drivers we simply do not label, and
**1.3 % to 9.5 %** are NCG-*known* high-confidence ones. Reactome's
arms are the worst, the dense protein networks the least.
>
**The consequence, and the correction we had to make.** Our first
answer here was wrong, and the way it was wrong is instructive. The
post-hoc version removes every NCG gene from the *evaluation set* — and
because NCG and IntOGen overlap heavily, that also removes **87 % of
the IntOGen positives** (633 down to 80), which is why the margin
appeared to halve, from +0.0651 to +0.0304. That is a statement about a
much smaller evaluation set, not about the background.

The clean experiment removes the curated drivers from the **background
only** — the positives are untouched — and re-runs the entire benchmark
on those splits. The result: **23 of 24 arms stay positive, no sign
flips, the median margin moves from +0.0809 to +0.0822** (slightly up
in almost every arm), and the conservative gene-cluster interval now
passes in **21 of 24** arms instead of 20. So the contamination of our
negative set is real, and it is **not** what carries the result — which
is a stronger statement than the one we started with, and it is the one
we will defend.
>
> **And four controls on the design itself**, because the headline is
> one number on one learner on one panel. **First, the margin survives
> losing seeds**: dropping 10 % of the training-fold seeds costs a median
> of **1.9 %** of the margin and only one of the 24 arms loses more than
> 10 % of it; halving them costs **11.7 %**. What predicts the loss is the
> *number* of seeds, not the gene-level noise — the arm with about 50
> seeds loses 29 % of its margin when half are dropped, the arm with
> about 470 loses 8 %.
> **Second, the four arms that fail the conservative gene-cluster
> interval fail for a reason of power, not of size**: on the equal-
> positive-size panel the failures are no longer only the newest-driver
> label — every arm with a small effect fails, so what the four share is
> a small effect.
> **Third, five labels are not five independent opinions**: each label's
> per-network means are **84 % to 99.8 %** predictable from the other
> four, and the overlap matrix has a participation ratio of **1.54** — so
> the "label main effect" in the variance decomposition is mostly one
> shared difficulty factor, and it is not five pieces of evidence about
> driver biology.
> **Fourth, and this is the one we have to say out loud: the increment is
> learner-specific.** Under a second, equally locked learner — the
> distance to the nearest seed — the margin is positive in **16 of 24**
> arms with a median of **+0.011**, against +0.081 for RWR, and on the
> densest network it is negative in every label. That is exactly the
> warning Guala and Sonnhammer raise, now quantified. Our claim is
> therefore only the pre-registered one: *under a fixed, pre-registered
> learner, the network contributes beyond degree* — and we do not
> generalise it to arbitrary priors. Removing STRING's `database`
> channel, by contrast, changes little, so the concern that our panel
> shares an upstream knowledge base with our labels is real in
> principle and small in this measurement.

**Three controls close this section.**
**The bait-matched null**, where each positive is paired with a
negative carrying the same publication count: the fraction of
positives above the null falls from a median **17.6 % to 11.6 %**
across all arms, and on the temporal label from **6.8 % to 4.7 %** —
and chance is **5 %**. On the other three labels the fall is from
**24.8 % to 13.7 %**. So a large part of the apparent margin is
visibility.
**The label-quality tiering**, where each gene gets a tier from its
own catalogue — ClinGen's review status, IntOGen's number of detecting
methods, NCG's literature count, and the literature/non-literature
split for Open Targets: the **high-confidence tier carries the larger
margin in four of five labels** (IntOGen **+0.0892** against +0.0633
and +0.0653; the temporal label **+0.0743** against +0.0123 and
+0.0232). NCG **inverts** the order, and that is not a
contradiction: NCG's tier *is* a publication count, so its high tier
is exactly the set the degree baseline already ranks well. The two
orderings are the same statement from two sides: **the prior's
increment is largest where the catalogue's own evidence is thinnest.**
**And the ClinGen 17 % noise anchor**: demoting a random 17 % of the
positives in every fold — 17 % being ClinGen's own inter-laboratory
disagreement rate — moves the margin by between **−0.0018** and
**+0.0026**, while a *perfect* scorer would lose **8.5** AUROC points
under the same misassignment. So the margin is **not** an artefact of
the label noise that ClinGen's curation process documents.
>
**One more diagnostic worth a sentence**: the gene-cluster bootstrap
interval width correlates with the arm's positive-set size at Spearman
**ρ = −0.966** (p = 2e-14, 24 arms), and with the high-tier share
only at **ρ = −0.475**. So the wide intervals are a
**small-positive-set phenomenon**, not a label-quality phenomenon.
>
**Finally, the fifth label and the curator lag.** On the
equal-prevalence panel, where all five labels have the same number of
positives inside each arm, **Reactome is positive on the 2019 NCG
catalogue at +0.0893** and negative on the two newest catalogues, at
**−0.0110** and **−0.0125**. The Reactome failure is therefore
**label-specific, not network-specific** — the first test of the
curator-lag hypothesis that does not rely on the temporal split. We
report it as *supported*, not proven, because a 2019 catalogue
differs from IntOGen 2024 in more than its age. And re-running the
variance decomposition on the 30-cell grid: the label main effect's
share of the between-arm spread falls from **34.78 % to 30.75 %**, the
interaction rises from **40.07 % to 41.95 %**, and the within-cell
noise rises from **40.39 % to 41.56 %** of the total variance. But
remember NCG is **65 % identical** to IntOGen — Jaccard **0.522** —
so this measures the effect of adding a **redundant** label, not an
independent one.

---

## 38:00–40:00 · 结论、局限与交付物（正文 slide 14–16）

[`【图：slide 14（结论）→ slide 15（局限 → 验证实验）→ slide 16（交付物 + AI 声明）】`]

Five conclusions.
**One**: a network prior contributes signal beyond node degree — 23 of
24 arms positive, the pre-registered rule met under four of five
interval constructions. **Two**: the contribution is real but modest
and label-dependent — median +0.081, and re-running the benchmark with the curated drivers removed from the negative set leaves it at +0.082 (23/24 arms, no sign flip) — the contamination of the background is not what carries the result. **Three**: much of the
published margin is study visibility — a bait-matched null removes
most of the enrichment, and the temporal label's advantage collapses
to +0.002. **Four**: the prior's instability is intrinsic — the
gene's own score variance explains it, and bagging the prior is a
provable no-op, so no learner can remove it. **Five**: curated
networks fail on recently discovered drivers — Reactome is +0.089 on
the 2019 catalogue and negative on the two newest.
>
**The limits, and what would change our mind** — every limitation we
know of has a named experiment: the gene-cluster interval misses four
arms → re-run on an independent cell-line-derived label; our negative
set is contaminated → re-derive it with explicit NCG and ClinGen
exclusion; the evaluation is self-referential → run a blind
expert-panel challenge; we have no perturbation experiment, so we
claim no causality → edit the seed set and measure the effect on the
predicted neighbourhood; we test one learner → re-run the locked
protocol with a second fixed learner.
>
**Everything is public and re-runnable**: the code, tables, figures
and split manifests are on GitHub at the v1.0-rc3 tag; the
pre-registration is frozen at Zenodo DOI 10.5281/zenodo.23121762 with
its amendments; two independent audit scripts pass **50 of 50** and
**63 of 63** checks; and every figure in this talk is generated from a
released table by a script in the repository.
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
So the one sentence I would like you to take away: **network priors do
help, but the honest number is much smaller than the headline one, and
the difference is mostly about how often a gene has been published.**
The benchmark, the pre-registration and the amendments are public —
please re-run it and check us.
>
Thank you.

---

# 附录 A · 必背数字卡（23 张）

每张卡：**数字 → 它是什么 → 怎么算的 → 说明什么 → 不能说明什么 → 一句话**。

| # | 数字 | 它是什么 | 怎么算的 | 说明什么 | 不能说明什么 | 一句话 |
|---|---|---|---|---|---|---|
| 1 | **+0.0809** | 24 臂 ΔAUROC 的中位 | 每臂 200 次 (repeat, fold) 评估的均值，再取 24 个臂的中位 | 先验相对度基线有约 0.08 AUROC 的平均增量 | 不等于"多找出 8 % 的 driver" | "这是增量的量级，不是召回率" |
| 2 | **23/24** | ΔAUROC 为正的臂数 | 同上 | 效应方向稳健 | 不等于"统计显著"（显著性看 CI） | "方向稳健，显著性看区间" |
| 3 | **+0.1711** | 最好臂（string_phys700 × ot_nolit） | 同上 | 蛋白网络 + 去文献标签这一组合最优 | 不能推广到别的 learner | "最好的一格，不是一般水平" |
| 4 | **−0.0125** | 唯一负臂（reactome × temporal） | 同上 | Reactome 在最新驱动上没有贡献 | **不能说它是"负效应"**（等效性检验：与 0 不可区分） | "没有可检出的贡献" |
| 5 | **24/24（×4）** | 四种区间构造下 H1 达成 | fold/repeat-percentile、repeat-t、repeat-BCa 的 95 % CI 是否排除 0 | 结论在四种重采样单位下都成立 | 不代表在第五种（gene-cluster）下成立 | "四种区间下成立，第五种下 20/24" |
| 6 | **20/24** | gene-cluster 区间下 H1 达成数 | 以基因为聚类单位重采样 2,000 次 | 保守界：4 个臂不够 | 失败的 4 个全是 temporal 标签（效应最小） | "保守界下有 4 个臂不够" |
| 7 | **12×** | gene-cluster 区间宽度 ÷ repeat 区间宽度 | 两者平均宽度之比 | 两个区间问的是不同问题 | 不是"gene-cluster 更有信息" | "换了问题，不是换了答案" |
| 8 | **0.936–0.944** | repeat 区间在 40 repeat 下的覆盖率 | 覆盖模拟（400 次迭代/格/模拟器） | repeat 区间是标定的 | 不是"95 % 置信等于 95 % 覆盖"的保证 | "标定过的区间" |
| 9 | **+0.0308 / +0.0412** | 低 FPR pAUC、top-10 % precision 的中位 margin | 同臂同折换端点 | 顶端的优势比 AUROC 上小得多 | 不能说"换端点结论就变了"（臂序 ρ=0.80–0.90 稳定） | "顶端的优势更小" |
| 10 | **0.0047** | 臂级 ΔAUROC 的标准差 | 200 次评估的 SD | 端点是稳的 | 与基因级不稳定不矛盾 | "均值稳，基因不稳" |
| 11 | **18 %** | 每基因分数 CV 的中位 | 同一基因在 200 次评估中的分数 SD/均值 | 单基因打分很抖 | 不等于"排名随机"（中位基因只挪约 270 位） | "单基因打分抖，端点不抖" |
| 12 | **≈270 位** | 中位基因的排名位移 | 两次评估之间该基因名次差 | 不稳定是真实的、不是幻觉 | 不能说"每次都乱" | "中位基因挪 270 位" |
| 13 | **+0.478 / −0.175** | 两种拥挤定义的效应 | 两种拥挤度量的回归系数 | 拥挤假说被**否证**（符号相反） | 不能说"拥挤完全无关" | "拥挤假说不成立" |
| 14 | **+0.616 / ρ=1.00 / q=2e-137** | 自身 score CV 对不稳定性的解释 | 基因级模型系数、Spearman、BH-FDR | 不稳定主要是基因自身分数方差 | 不能说"拓扑无关"（网络层仍有覆盖度效应） | "不稳定来自自身，不是拓扑" |
| 15 | **+0.038 → −0.083** | 加入 score CV 后 degree 系数的变化 | 同一模型的嵌套比较 | 原来的度效应是**抑制假象** | 不能说"度完全无影响" | "度效应原来是假象" |
| 16 | **ρ=−0.905, q=0.012** | 网络级稳定性的预测量 | 8 个网络的 Spearman + BH-FDR | 覆盖度与 2-core 预测稳定性 | 规模/密度/模块度不预测 | "稳定取决于覆盖与核心" |
| 17 | **6 % vs 48 %** | DepMap 细胞系富集率（K=30 / K=100） | 每细胞系超几何检验，1,361 个细胞系 | 集合层面有信号，基因层面没有 | 不能说"先验在单基因水平有用" | "集合有效，基因无效" |
| 18 | **24/24, +0.0413** | NCG 仲裁测试的 ΔAUROC | 阳性=我们标签没有的 NCG driver | 先验能迁移到没见过的 driver | 低 FPR 端点只有 13/24 | "能迁移，但顶端取决于网络" |
| 19 | **26/26, +0.0271** | 癌种分层的中位 margin | 按 NCG 癌种注释拆分仲裁阳性 | 优势与癌种无关 | 不能说"每个癌种都很强" | "与癌种无关" |
| 20 | **10–22 %** | 阴性集中策展 driver 的比例 | 用 NCG 仲裁我们的 holdout 阴性 | 基准的阴性集不干净 | **但把它清掉并重跑，结论不变**（23/24 臂为正、中位 +0.082） | "阴性集脏，但结论不靠它" |
| 21 | **17.6→11.6 %；6.8→4.7 %** | bait 匹配前后阳性超过零模型的比例 | 每个阳性配一个发表量相同的阴性 | 一部分 margin 是可见度 | 不能说"全是发表量"（其他标签仍留有信号） | "temporal 标签归零到随机水平" |
| 22 | **−0.0018..+0.0026** | 17 % 误标下的 margin 变化 | 每折随机降级 17 % 阳性 | margin 不是标签噪声的产物 | 不能说"我们的标签无噪声" | "17 % 噪声不伤 margin" |
| 23 | **+0.0893 vs −0.0110** | reactome × NCG 2019 vs × IntOGen 2024 | 等阳性率面板 | Reactome 失效是标签特异的 | 不能说"curator lag 已被证明" | "是标签年龄，不是网络" |
| 24 | **1.9 % / 11.7 %** | 掉 10 % / 50 % 种子后 margin 的损失 | 每个臂随机抽掉一部分训练种子重跑 | margin 对丢种子稳健 | 不能说丢一半种子也没事（temporal 臂损失 29 %） | 掉 10 % 种子只掉 2 % |
| 25 | **20/24 → 18/24** | gene-cluster 区间通过数（等阳性数面板） | D-11 capped 面板上重算（阳性数相同） | 失败的是小效应臂，不是小标签 | 不能说把阳性数等量化就能救回 | 是功效问题 |
| 26 | **84 %–99.8 %** | 每个标签的每网络均值可由其他标签预测的 R² | 5×6 arm 均值矩阵的回归 | 五个标签≈一个共同维度 | 不能说五种策展哲学 | 标签几乎是同一个维度 |
| 27 | **16/24, +0.011** | 第二个 learner（邻近度）下的 margin | 同切分同基线，score = −到最近种子的距离 | 增量是 learner-specific 的 | 不能说换个方法也成立 | 只在我们锁定的方法下成立 |

---

# 附录 B · 最可能被问到的十个问题（一句话版）

1. **新发现还是边界刻画？** —— 边界刻画；每条结论都配了可证伪实验与"哪张表会变"。
2. **为什么不用 GNN？** —— learner 必须冻结，否则测的是"调出了好方法"；GNN 礼节性对照 RWR +0.079 vs 嵌入 +0.028。
3. **阴性集干净吗？** —— 不干净，10–22 % 是策展 driver；我们自曝并给了清洗后两套数。
4. **H1 到底成立吗？** —— 5 种区间里 4 种 24/24，gene-cluster 20/24，4 个例外全是 temporal。
5. **是不是"网络越大越好"？** —— D-11 cap 与 D-12 共享宇宙压住了这个解释（−0.0360），但覆盖度确实解释一部分（ρ=−0.905）。
6. **ΔAUROC 能读成"多找出 8 % driver"吗？** —— 不能，低 FPR 端点只有 +0.031。
7. **为什么用五个目录？** —— 策展哲学不同；26/26 癌种为正说明不是某一份目录的怪癖。
8. **能说因果吗？** —— 不能，没有扰动实验；验证计划里有种子扰动实验。
9. **AI 会不会编数字？** —— 全部数字来自已发布表格 + 两个审计脚本（50/50、63/63）；PPT 的数字由脚本导出。
10. **核心差异不成立会怎样？** —— 见题库 C2：三个假说各自对应哪张表会变。

---

# 附录 C · 时间压缩方案（若只剩 15 分钟）

| 保留 | 时长 | 理由 |
|---|---|---|
| 开场 + 研究地图 | 1.5 min | 立住问题 |
| 三篇来源论文（每篇"一问/一数/一局限"一句） | 3 min | 评分硬要求 |
| 知识空白 + 三个假说 | 1.5 min | 评分硬要求 |
| 数据与设计（只讲 provenance、D-02、cap 顺序） | 3 min | 评分硬要求 |
| 主端点 + H1 的准确表述 | 2 min | 核心结果 |
| 稳定性机制（含"拥挤假说被否证"） | 1.5 min | 老师点名 |
| 外部证据 + 阴性集污染 | 1.5 min | 诚实性得分点 |
| 结论 + 局限 + 交付物 | 1 min | 收尾 |

**砍掉**：控制栈逐条（只说"11 层，每层一个数"）、第五标签的方差分解细节、参数表逐行（被问再展开）。
**绝不能砍**：H1 的精确表述、阴性集污染、稳定性机制与拥挤假说的自我更正。

---

# 附录 D · 彩排清单

- [ ] 第一遍：按本稿计时，目标 **14 分钟**（40 分钟讲不满会显得内容不足；15 分钟讲超时会被打断）。
- [ ] 第二遍：合上稿子讲，只看 `report_1118_facts.md` 的数字表 + 附录 A 的 23 张卡。
- [ ] 第三遍：三人互问，每人在附录 B 之外再抽 3 题。
- [ ] 数字核对：讲稿里的每个数字都能在 `report_1118_facts.md` 里找到出处（**口头数字错一个，全盘不可信**）。
- [ ] 检查：PPT 在另一台机器上能打开（字体 Times New Roman，图为 PNG，不依赖外链）。
- [ ] 备份：`D1_defence_deck.pptx`、本稿、题库各一份到 U 盘。
- [ ] 提醒自己：被问到不会的，**说"我们没有做这个实验，这是我们下一步要做的"**，不要临场编数。
