# AlzKG

An Alzheimer's disease (AD) knowledge graph, built with the same pipeline shape
as [PrimeKG](../PrimeKG) (Chandak, Huang & Zitnik, *Scientific Data* 2023 --
[`../Paper#7-s41597-023-01960-3.pdf`](../Paper%237-s41597-023-01960-3.pdf)):
primary data resources -> processing scripts -> a harmonized edge list -> a
case-study notebook. PrimeKG's own `case_study/autism.ipynb` does exactly this
for autism spectrum disorder; `case_study/alzheimers_disease.ipynb` here does
the AD equivalent.

This folder does **not** duplicate PrimeKG. It reuses PrimeKG's generic,
disease-agnostic infrastructure (MONDO, HPO, GO, Reactome, UBERON, NCBI Gene,
CTD, SIDER, UMLS, DisGeNET, DrugBank, Drug Central...) by relative path, and
adds only what's specific to Alzheimer's disease: AD GWAS/variant evidence,
AD clinical trials, AD pharmacogenomics, and an AD literature-mining feed.

## How this maps to the project's other files

| File | Role |
|---|---|
| [`../Data_source_for_research_AD-KG.xlsx`](../Data_source_for_research_AD-KG.xlsx) | The full catalog of 217 data sources found across 30+ papers, each mapped to what it contributes to the KG. Exported here as [`vocab/AD_KG_Data_Sources.csv`](vocab/AD_KG_Data_Sources.csv). |
| `Node_and_Edge_Types` sheet (same workbook) | The node/edge schema derived from that catalog. Exported here as [`vocab/Node_and_Edge_Types.csv`](vocab/Node_and_Edge_Types.csv) -- this is the schema `knowledge_graph/build_ad_kg.ipynb` targets. |
| [`../AD_Knowledge_Graphs_Bibliography.md`](../AD_Knowledge_Graphs_Bibliography.md) | 131-paper bibliography. Part I motivates *which* sources matter for AD; Part II ("Dynamic Knowledge Integration") motivates *why* `pubtator_literature_ad.py` and `clinicaltrials_ad.py` are written as re-runnable/incremental jobs rather than one-shot downloads. |
| [`../PrimeKG`](../PrimeKG) | The disease-agnostic KG this pipeline is structurally modeled on, and whose processed ontology outputs this pipeline reads directly. |

## Why not script all 217 sources?

Most of the catalog's 217 rows are either (a) generic ontologies/interactomes
PrimeKG already processes (MONDO, GO, HPO, Reactome, UBERON, STRING, ...), (b)
sources requiring restricted/manual access that PrimeKG itself documents as
manual steps rather than scripts (DrugBank, UMLS, ADNI, NIAGADS Qualified
Access, DrugCentral's Postgres dump), or (c) one-off literature citations, not
databases with a stable bulk-download endpoint. This repo scripts the sources
that are (1) AD-specific and (2) programmatically retrievable, and documents
the rest inline (see each script's docstring and `datasets/primary_data_resources.sh`).

## Structure

```
AlzKG/
├── datasets/
│   ├── primary_data_resources.sh      # curl/API commands for every scriptable AD-specific source
│   ├── processing_scripts/
│   │   ├── gwas_catalog.py            # GWAS Catalog -> variant/gene-disease edges
│   │   ├── clinvar_ad.py              # ClinVar -> variant-disease, variant-gene edges
│   │   ├── disgenet_ad.py             # DisGeNET (curated) filtered to AD -> gene-disease edges
│   │   ├── clinicaltrials_ad.py       # ClinicalTrials.gov API v2 -> drug-disease-via-trial edges
│   │   ├── pharmgkb_ad.py             # PharmGKB relationships filtered to AD -> pharmacogenomic edges
│   │   ├── niagads_ad.py              # NIAGADS GWAS summary stats (manual download) -> variant/gene-disease edges
│   │   └── pubtator_literature_ad.py  # PubTator3 API -> literature co-mention triples (re-runnable/incremental)
│   └── data/                          # gitignored; created by primary_data_resources.sh
├── knowledge_graph/
│   ├── build_ad_kg.ipynb              # harmonizes all of the above (+ PrimeKG's generic outputs) into ad_kg.csv
│   ├── visualize_ad_kg.py             # writes an interactive HTML explorer of the core AD subgraph
│   └── ad_kg_explorer_template.html   # page template used by visualize_ad_kg.py
├── case_study/
│   └── alzheimers_disease.ipynb       # queries ad_kg.csv, mirrors PrimeKG/case_study/autism.ipynb
├── scripts/
│   └── utils.py                       # shared helpers (copied from PrimeKG/scripts/utils.py)
├── vocab/
│   ├── AD_KG_Data_Sources.csv         # exported from ../Data_source_for_research_AD-KG.xlsx
│   └── Node_and_Edge_Types.csv        # exported from the same workbook -- the target schema
└── requirements.txt
```

## Node & edge schema

Every processing script emits rows with the same 10 columns PrimeKG's
`build_graph.ipynb` uses, so AlzKG can be diffed against or merged with PrimeKG:

```
relation, display_relation, x_id, x_type, x_name, x_source, y_id, y_type, y_name, y_source
```

`knowledge_graph/build_ad_kg.ipynb` harmonizes entities across sources (genes ->
`NCBIGene:<id>` via NCBI `gene_info`, variants -> dbSNP `rs<id>`, every spelling of
Alzheimer's disease -> `MONDO:0004975`, drugs/other diseases merged by normalized
name), links AD subtypes to the AD node, restricts the graph to the 2-hop AD
neighborhood, then adds `x_index`/`y_index` and one extra column, `edge_source`
(the dataset that asserted the edge, since merged nodes no longer say). See `vocab/Node_and_Edge_Types.csv` for the full catalog of node types
(Disease, Gene/Protein, Drug/Chemical/Compound, Variant/Mutation, Clinical
Trial/Study, Publication/Article, ...) and edge types (Drug-Target,
Disease-Gene, Variant-Disease, Drug-Disease via clinical trial, literature
co-mention, ...) this schema is meant to cover.

## Building the graph

```bash
pip install -r requirements.txt

# 1. (Optional) Generic ontology/pathway/interaction sources from PrimeKG.
#    build_ad_kg.ipynb skips them if ../PrimeKG is not present.
cd ../PrimeKG/datasets && bash primary_data_resources.sh && cd ../../AlzKG

# 2. AD-specific sources (~700 MB of downloads on first run; cached afterwards).
cd datasets && bash primary_data_resources.sh && cd ..

# 3. Harmonize into ad_kg.csv / ad_kg_nodes.csv.
python -m nbconvert --to notebook --execute --inplace knowledge_graph/build_ad_kg.ipynb

# 4. Sanity-check the result.
python -m nbconvert --to notebook --execute --inplace case_study/alzheimers_disease.ipynb

# 5. Interactive explorer -> datasets/data/kg/ad_kg_explorer.html
python knowledge_graph/visualize_ad_kg.py
```

(`python -m nbconvert` works even when the `jupyter` launcher isn't on PATH, as
with the Windows Store Python.)

Steps 1 and 2 hit live public APIs/FTP servers. DisGeNET and NIAGADS require a
manual download first (see the relevant script's docstring); until then their
scripts print a `[skip]` message and the rest of the pipeline runs without them,
exactly as PrimeKG's own pipeline requires manual steps for DrugBank, UMLS,
and Drug Central.

Source endpoints changed since this pipeline was first written: the GWAS
Catalog's `/api/search/downloads/full` endpoint is gone (releases are on
`ftp.ebi.ac.uk`), and PharmGKB is now ClinPGx (`api.clinpgx.org`).

## Keeping it current (dynamic integration)

`clinicaltrials_ad.py` and `pubtator_literature_ad.py` are written to be
re-run periodically rather than once: the former re-queries the live
ClinicalTrials.gov API (new/updated AD trials appear automatically), and the
latter tracks its own checkpoint file so repeated runs only fetch and append
new literature since the last run. This is the "dynamic, incrementally
updated" pattern that Part 11 of `AD_Knowledge_Graphs_Bibliography.md`
argues for, in contrast to AlzKB/PrimeKG/ADKG's static, periodically-rebuilt
snapshots.
