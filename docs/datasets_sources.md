# Public data sources, releases and identifiers

**What this file is for.** The course requires a "reference and public
dataset source list" with dataset identifiers (GEO / SRA / PRIDE /
MetaboLights). This project re-analyses **already published, already
public** data — it produces no new sequencing data — so the honest
identifier for most of our inputs is the **release identifier and the
exact download URL**, which is what `scripts/download_D1_data.py`
actually fetched. Where a DOI-style accession exists we give it.
Anything we could not verify offline is marked **TO CONFIRM** and must be
checked in a browser before submission.

**Rule followed.** Every URL below was read out of our own download
script or our own data directory; none is written from memory.

---

## 1. The network panel

| # | Dataset | Release / version | Source URL (as downloaded) | Local file |
|---|---|---|---|---|
| 1.1 | **FunMap** functional network (protein-derived; CPTAC pan-cancer proteogenomics) | FunMap v1.0 | `https://funmap.linkedomics.org/data/download/funmap.tsv` | `data/raw/networks/funmap.tsv` |
| 1.2 | **FunMap** published CPTAC mRNA matrices (the RNA-derived network) | Zenodo record **7948944** | `https://zenodo.org/records/7948944/files/funmap_input_expression_data.tgz?download=1` | `data/raw/networks/funmap_expression.tgz` |
| 1.3 | **IntAct** human-human interactions (taxid 9606) | 2024-09 (the `intact.zip` in the current release) | `https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip` | `data/raw/networks/intact.zip` |
| 1.4 | **STRING v12** all-evidence links (combined_score ≥ 700) | 12.0 | `https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz` | `data/raw/networks/9606.protein.links.detailed.v12.0.txt.gz` |
| 1.5 | **STRING v12** physical subnetwork (≥ 700) | 12.0 | `https://stringdb-downloads.org/download/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz` | `data/raw/networks/9606.protein.physical.links.v12.0.txt.gz` |
| 1.6 | **STRING v12** protein info / aliases (the symbol crosswalk) | 12.0 | `https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz` and `.../protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz` | `data/raw/networks/` |
| 1.7 | **Reactome** human interactions | v86 (content version 86, September 2023) | `https://reactome.org/download/current/interactors/reactome.homo_sapiens.interactions.tab-delimited.txt` | `data/raw/networks/reactome.tsv` |

## 2. The label panel

| # | Dataset | Release | Source URL (as downloaded) | Local file |
|---|---|---|---|---|
| 2.1 | **IntOGen** driver catalogues (MC3 / recurrence) | 2024-09-20, 2023-05-31, 2020-02-01 | `https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip` (and the 20230531 / 20200201 zips) | `data/raw/labels/intogen/*.zip` |
| 2.2 | **Open Targets Platform** 26.06 | 26.06 | `https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.06/output` (association_overall_direct, association_by_datasource_direct, association_by_datatype_direct, target, disease) | `data/raw/labels/opentargets_26.06/` |
| 2.3 | **ClinGen** gene-disease validity | the file created **2026-09-26** | `https://search.clinicalgenome.org/kb/gene-validity/download` | `data/raw/labels/clingen/clingen_gene_validity.csv` |
| 2.4 | **NCG** driver catalogue (the arbitration set) | the 2019 supplementary files of the **Genome Biology** paper (see §5) | the article's Additional file 2 / 3 | `data/raw/labels/ncg/13059_2018_1612_MOESM{2,3}_ESM.xlsx` |
| 2.5 | **HGNC** complete set (the symbol authority, decision D-02) | the file fetched with the release | `https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt` | `data/raw/reference/hgnc_complete_set.txt` |
| 2.6 | **PubMed** gene→publication counts (the publication-count adjustment and the bait-matched null) | the `gene2pubmed` dump fetched with the release | `https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz` | `data/raw/reference/gene2pubmed.gz` |

## 3. The external validations

| # | Dataset | Release | Source | Local file |
|---|---|---|---|---|
| 3.1 | **DepMap** CRISPR gene effect (the primary, 1,361 cell lines) | the release whose file we hold | DepMap portal (download needs a free account) | `data/raw/depmap/CRISPRGeneEffect (1).csv` |
| 3.2 | **DepMap** CRISPR gene effect (the earlier release, 1,150 cell lines) | earlier | DepMap portal | `data/raw/depmap/CRISPRGeneEffect.csv` |
| 3.3 | **DepMap** expression (log2(TPM+1), protein-coding) | as above | DepMap portal | `data/raw/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv` |
| 3.4 | **DepMap** expression (log2(TPM+1), human protein-coding) | as above | DepMap portal | `data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` |
| 3.5 | **UniProtKB / ComplexPortal** (the protein-complex map used for the complex-aware splits) | the REST snapshot fetched on 2026-11-05 | `https://rest.uniprot.org/uniprotkb/search` (cross-references to ComplexPortal) | `data/raw/complexportal/` → `data/processed/complexes/complex_map.tsv` |

## 4. Our own artefacts (the identifiers a reader needs)

| Artefact | Identifier |
|---|---|
| The pre-registration (frozen 2026-09-27) | **Zenodo DOI 10.5281/zenodo.23121762** (`https://doi.org/10.5281/zenodo.23121762`) |
| The amendments to it | `docs/pre_registration_amendments.md` in the repository (A2.3a–A2.3k + the analysis-family register) |
| The code, tables, figures and split manifests | **GitHub `github.com/lanxin2391/d1-benchmark`**, tag `v1.0-rc3` |
| The release bundles | three zips (source / splits / results) + `docs/release_manifest_v1.0-rc2.json` |
| The v1.0 Zenodo DOI | **TO CONFIRM** — the v1.0 record is minted in a browser (the placeholder `10.5281/zenodo.d1-benchmark-v1.0` is not yet live) |
| The bio.tools entry | **TO CONFIRM** — registered in a browser |
| The leaderboard placeholder | `docs/leaderboard.md` (**TO CONFIRM** — not created yet) |

## 5. A correction we must carry into the reference list

The NCG catalogue has **two** associated papers, and the supplementary
files we actually downloaded come from the **Genome Biology** one:

* Repana D, Nulsen J, Dressler L, Bortolomeazzi M, Kuppili Venkata SK,
  Tourna A, Yakovleva A, Palmieri T, Ciccarelli FD. **The Network of
  Cancer Genes (NCG): a comprehensive catalogue of known and candidate
  cancer genes from cancer sequencing screens.** *Genome Biology*
  2019;20:1. doi:**10.1186/s13059-018-1612-0** — **this is the source of
  Additional file 2 and 3, i.e. of our 2,372-gene catalogue.**
* Repana D, D'Antonio M, Pendino V, Ratti E, Gentilini M, Ciccarelli
  FD. The network of cancer genes. *Nucleic Acids Res.*
  2019;47(D1):D671-D675 — the NCG 6.0 methods paper.

The raw-data note in `data/raw/labels/ncg/README.md` cited only the NAR
paper; it has been corrected to name the Genome Biology paper as the
source of the files we used.

## 6. TO CONFIRM before the submission (needs a browser)

1. **CPTAC proteogenomics accessions (PRIDE).** The FunMap network is
   built from the CPTAC pan-cancer proteogenomics; the upstream raw
   mass-spectrometry files have PRIDE / ProteomeXchange accessions.
   We did **not** re-download them (we used FunMap's published matrices
   and edge list), so we do not record a PXD id we cannot verify.
   **Action:** copy the PXD accession(s) from the FunMap paper's data
   availability section (Shi Z et al., *Nat Cancer* 2025;6:205-222) into
   this list.
2. **DepMap release version.** Our two `CRISPRGeneEffect` files come
   from two DepMap public releases (1,150 and 1,361 cell lines); the
   release tag (e.g. 24Q2 / 25Q2) is in the download portal's file
   listing. **Action:** write the exact release tag next to each file.
3. **IntOGen release identifiers.** The three driver zips are named by
   release date (2020-02-01, 2023-05-31, 2024-09-20); the compendium
   inside the 2024-09-20 archive still carries a "RELEASE 20230531"
   header. **Action:** confirm with the IntOGen portal which release each
   archive corresponds to.
4. **Reactome v86 vs "current".** We downloaded the file from the
   `current` interactor export but cite content version 86
   (September 2023). **Action:** re-download the archived v86 export
   (or state the version we actually used) so the file and the citation
   agree.
5. **The clean-machine reproduction.** Not an identifier, but the release
   protocol's hard gate: clone → install → `pytest` → re-run the two
   audit scripts on a clean machine, and record the result in
   `docs/reproducibility_audit.md`.

## 7. Why there are no GEO / SRA / MetaboLights identifiers here

The course asks for "public dataset identifiers (GEO, SRA, PRIDE,
MetaboLights)" **for the data you actually used**. Our project uses
**aggregated releases and curated portals** — a network edge list, a
driver catalogue, a gene-disease validity table, a CRISPR screen matrix
and a publication-count table. None of these is a sequencing run
produced by us or deposited by us, and none of them has a GEO/SRA
accession; the ones that do have a repository accession (the CPTAC raw
proteomics behind FunMap, Zenodo 7948944) are listed above with the
identifier we can verify, and the one we cannot verify offline is
explicitly flagged for the browser check. Inventing an accession that
was never in our download record would be the exact failure the course's
academic-integrity clause forbids.
