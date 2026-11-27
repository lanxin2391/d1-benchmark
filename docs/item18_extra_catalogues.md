# Item 18 — more independent label catalogues

**Date:** 2026-11-26

## The question

The paper says the four catalogues are **one curation dimension, not four
independent votes** (the label factor explains 97.8 % of the label main
effect, F.23). The direct check needs more catalogues, in particular
catalogues that are curated by *different* mechanisms.

## What is reachable from this machine (probed 2026-11-26)

| source | the status | what it would give |
|---|---|---|
| the **ICGC / PCAWG** driver list | **reachable** (the API answers) | a somatic-mutation-driven catalogue, independent of the literature curation |
| the **COSMIC** gene census | the site answers, **a login is needed** for the download | the classic somatic census |
| the **NCG** | the site is **down for maintenance** (the DNS does not resolve) | the literature census (we use the 2018 supplementary material) |
| the **CIViC** | the API returns 410 (gone) | the curated clinical variants |
| **OncoKB** | the API answers 401 (a **token** is needed) | the therapeutic-curated list |
| the cBioPortal studies | the study id I tried returns 404 | the Karchin 2019 census |

So: **the ICGC catalogue can be built here**; the COSMIC, OncoKB and NCG
catalogues need the user's account or download.

## What the analysis needs

With two or three more catalogues the label-dimension analysis can be
re-run: the variance decomposition with six or seven labels instead of
four, the per-label agreement matrix, and the question the paper states
as H4 — whether "the four catalogues are not four independent pieces of
evidence" still holds with more catalogues of a different kind.

## The honest dependency

Two of the four additional sources need something this machine does not
have (a COSMIC login, an OncoKB token) and one needs the user to download
it while the NCG site is down. The ICGC one needs no credential and can
be started immediately.
