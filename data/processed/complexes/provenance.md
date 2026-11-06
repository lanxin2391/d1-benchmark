# Provenance of the complex map

**Built:** 2026-11-06 21:09 by `scripts/build_complex_map.py`

## Source

ComplexPortal complex identifiers (CPX-...) as curated in **UniProtKB**,
retrieved through the UniProt REST search API
(`https://rest.uniprot.org/uniprotkb/search`). UniProt release header:
`unknown`.

**Why not the canonical ComplexPortal dump:** this machine has no
external FTP or website access. `ftp.ebi.ac.uk/pub/databases/complexportal/`
returns 404 and `www.ebi.ac.uk/complexportal` is a JavaScript
single-page application with no static data files; CORUM
(`mips.helmholtz-muenchen.de`) is unreachable. UniProt's REST API is
reachable and carries the ComplexPortal cross-references, so the
identifiers are ComplexPortal identifiers and the underlying curation is
ComplexPortal's. The difference from the canonical dump is that only the
complexes annotated to a UniProt entry are seen, and the complex
*member list* is not used at all (see below).

## Why the member list is not needed

Two genes are complex-linked iff their CPX-ID lists intersect. The
leakage audit (D1) and the complex-aware splitter (D2) both work on
gene-level complex links, so only `gene -> CPX-IDs` is required.

## Gene identifier mapping

`data/raw/reference/hgnc_complete_set.txt` (the HGNC complete set,
already in the repo) supplies the exact `symbol -> UniProt accession`
mapping from its `uniprot_ids` column, preferring the reviewed
(Swiss-Prot) accession when several are listed. No symbol was guessed.

## Counts

| quantity | value |
|---|---|
| panel genes (union of the six networks' LCCs) | 19,442 |
| mapped to a UniProt accession | 18,908 |
| unmapped | 534 |
| accessions queried | 18,778 |
| complexes found | 2,569 |
| distinct genes in at least one complex | 3,708 |
| (gene, complex) pairs | 9,619 |

## Raw artefact

`data/raw/complexes/uniprot_complexportal.tsv` -- SHA-256
`46a84588cc46ae3c5dddca248a1ea15635146b069549b3ce4aebf539d0ced551`

## Caveat to state in the paper

The complex set is UniProt's view of ComplexPortal, not ComplexPortal's
complete dump, and the complexes are the *curated* ones. It is used
only to decide which genes may share a fold; a gene in no complex is
treated as unconstrained, exactly as the protocol specifies.
