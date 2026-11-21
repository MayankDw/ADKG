# AlzKG: Abstract and Summary

> The figures below come from the harmonized build committed in `2757855` (15,813 nodes, 27,537 edges).
> The uncommitted `knowledge_graph/build_ad_kg.ipynb` in the working copy is an older version (no gene-ID
> harmonization, no 2-hop filter, no `edge_source` column), and the `ad_kg.csv` currently on disk was written
> by it (17,488 nodes, 29,615 edges, with Alzheimer's disease split across several IDs). Before citing these
> numbers, run `git checkout knowledge_graph/build_ad_kg.ipynb` and rebuild.

---

## Abstract

Alzheimer's disease (AD) is the leading cause of dementia. Evidence about its genetic architecture, pathogenic
variants, drug–gene interactions and therapeutic pipeline keeps growing, but it is spread across GWAS
repositories, clinical-variant archives, pharmacogenomic databases, trial registries and the literature. Each
uses its own identifiers, vocabulary and update cycle. Disease-agnostic biomedical knowledge graphs treat AD as
one node among thousands of diseases and are rebuilt as static snapshots, so they leave out much of the
AD-specific genetic evidence, the active trial pipeline and recent findings.

Here, we present **AlzKG**, a knowledge graph centered on Alzheimer's disease. It is built with PrimeKG's
pipeline architecture and edge schema. AlzKG integrates five AD-specific, programmatically retrievable
resources: the NHGRI-EBI GWAS Catalog, ClinVar, ClinicalTrials.gov, PharmGKB/ClinPGx and PubTator3 literature
mining. They were chosen from a catalog of 217 data sources surveyed across 30 AD knowledge-graph publications.
Genes are harmonized to NCBI Gene, variants to dbSNP, and every spelling and subtype of Alzheimer's disease to
MONDO. The result has **15,813 nodes** of six types (genetic variants, genes/proteins, clinical interventions,
drugs/chemicals, diseases and phenotypes) and **27,537 edges** of 13 relation types, limited to the two-hop
neighborhood of AD.

The AD node alone links to 12,817 distinct entities. These include 5,543 risk variants, 3,571 genes, and 3,221
drug and intervention links taken from 4,299 registered AD trials. Every edge records the dataset that asserted
it, so you can query how many sources support a claim. For example, APOE is supported by GWAS, ClinVar,
PharmGKB and literature evidence together. AlzKG also includes drug–disease edges labeled by trial
participation and drug–disease correlations mined from the literature, both positive and negative. It connects
329 drugs to AD through 52 AD-associated genes, which gives candidate paths for drug-repurposing and mechanism
studies.

Unlike static snapshots, the trial and literature layers are built to be re-run incrementally, so the graph can
follow new trials and publications as they appear. AlzKG keeps PrimeKG's ten-column edge format, so it can be
merged with PrimeKG's ontology, pathway and interactome layers or compared with it directly.

---

## Summary

### 1. What this is

AlzKG is a disease-focused counterpart to PrimeKG. PrimeKG (Chandak, Huang & Zitnik, *Sci. Data* 2023) covers
17,080 diseases broadly. AlzKG covers one disease in depth. It copies PrimeKG's structure: primary data
resources → processing scripts → a harmonized edge list → a case-study notebook.
`case_study/alzheimers_disease.ipynb` is the AD version of PrimeKG's `autism.ipynb`. The project reuses PrimeKG
for generic ontologies such as MONDO, GO, HPO and Reactome instead of rebuilding them, and adds only the
evidence that is specific to AD.

### 2. How it is built

1. **Source selection.** 217 data sources found in 30+ AD-KG papers were cataloged
   (`vocab/AD_KG_Data_Sources.csv`). A target schema of 21 node types and 20 edge types was derived from them
   (`vocab/Node_and_Edge_Types.csv`). Only sources that are both AD-specific and downloadable by script were
   automated.
2. **Acquisition.** `datasets/primary_data_resources.sh` downloads everything from public APIs and FTP; files
   already downloaded are skipped.
3. **Per-source processing.** Each script in `datasets/processing_scripts/` filters its source to AD and writes
   edges in PrimeKG's 10-column format:
   `relation, display_relation, x_id, x_type, x_name, x_source, y_id, y_type, y_name, y_source`.
4. **Harmonization** (`knowledge_graph/build_ad_kg.ipynb`):
   - Gene symbols → `NCBIGene:<id>`, using NCBI `gene_info` (193,752 official symbols and 65,003 unambiguous
     synonyms).
   - Variants → dbSNP rsIDs.
   - Drugs, diseases and phenotypes merged by normalized name (426 merges).
   - 152 AD subtype and trait nodes (early- and late-onset AD, "AD biomarker measurement", and so on) linked to
     `MONDO:0004975`.
5. **Scoping.** Only edges within two hops of AD are kept: 27,582 of 29,087. This drops literature co-mentions
   unrelated to AD but keeps second-hop paths such as gene → drug.
6. **Provenance.** An `edge_source` column is added, because merged nodes no longer show where an edge came
   from.
7. **Validation and exploration.** The case-study notebook walks outward from the AD node.
   `knowledge_graph/visualize_ad_kg.py` writes an interactive HTML explorer.

### 3. Data sources and what each contributes

| Source | What it is | What it adds to AlzKG | Edges |
|---|---|---|---|
| **GWAS Catalog** (NHGRI-EBI) | Curated genome-wide association studies | Variants and genes statistically linked to AD, late- or early-onset AD, family history, age of onset and AD biomarkers: the population-genetics layer | 14,815 |
| **ClinicalTrials.gov** (API v2) | Registry of clinical studies | Drug and intervention → AD edges ("studied in clinical trial") from 4,299 trials, about 2,300 of them completed and 600 recruiting, across phases 1–4: the therapeutic-pipeline layer | 9,404 |
| **ClinVar** (NCBI) | Archive of clinically interpreted variants | Pathogenic or likely pathogenic variants (e.g. PSEN1, PSEN2, APP) with clinical significance kept on the edge label, plus variant → gene edges: the Mendelian, familial-AD layer | 1,575 |
| **PubTator3** (NCBI) | Entity and relation extraction from PubMed | Literature-derived drug/gene/variant–disease relations (association, positive or negative correlation, binding, cotreatment) from about 300 papers per run: the emerging-evidence layer | 1,199 |
| **PharmGKB / ClinPGx** | Curated pharmacogenomics | Gene–drug and variant–drug relations (e.g. CYP2D6, APOE, BCHE) that tie AD genes to drugs: the precision-medicine layer | 392 |
| AlzKG harmonization | Generated by the pipeline | Subtype → AD links | 152 |
| *DisGeNET, NIAGADS* | Curated gene–disease data; AD GWAS summary statistics | Scripts are written but need a manual login download, so they are currently skipped | — |
| *PrimeKG layer* | MONDO, HPO, GO, Reactome, UBERON, CTD, SIDER, DrugBank… | Optional; adds pathways, PPIs, anatomy and side effects. It is not in the current build because `../PrimeKG` is absent | — |

### 4. What AlzKG offers that general KGs don't

- **AD depth instead of breadth.** In PrimeKG, AD is one node with generic edges. In AlzKG, AD has about 12,900
  edges covering population genetics, clinical genetics, pharmacogenomics, trials and literature.
- **Graded genetic evidence.** Variant edges keep their ClinVar significance (pathogenic, likely pathogenic,
  risk factor, protective, conflicting), so a model can tell causal Mendelian variants apart from GWAS
  associations.
- **The trial pipeline as edges.** Failed, ongoing and approved interventions (donepezil, memantine, lecanemab,
  donanemab, aducanumab…) are first-class nodes. This gives the negative and in-progress therapeutic evidence
  that indication-only KGs leave out.
- **Signed literature relations.** Positive and negative drug–disease correlations give direction to
  drug–disease edges, which is useful for repurposing.
- **Provenance-aware edges.** `edge_source` lets you rank entities by how many independent sources support
  them. APOE is supported by four sources, APP/ABCA7/SORL1 by three.
- **Dynamic updates.** The trial and literature scripts can be re-run, and PubTator keeps a checkpoint so each
  run fetches only new papers. Existing AD KGs (AlzKB, PrimeKG, ADKG) are static snapshots.
- **PrimeKG compatibility.** Because the schema is identical, AlzKG can be merged into PrimeKG or diffed
  against it without remapping.

### 5. Why it is useful

- **Drug repurposing.** 329 drugs are two hops from AD through 52 AD genes (e.g. CYP2D6 → 228 drugs, CHRNA7 →
  donepezil and galantamine). These are ready-made candidate paths for link-prediction models such as
  TxGNN-style GNNs.
- **Target prioritization.** Genes supported by GWAS, ClinVar, PharmGKB and literature together (APOE, APP,
  ABCA7, SORL1, TREM2, CLU, BIN1) can be ranked by convergent evidence.
- **Trial landscape analysis.** You can ask which mechanisms or targets are crowded and which are
  under-explored, or which genes have variants but no drugs in trials.
- **Pharmacogenomic stratification.** CYP2D6 and APOE drug edges support analyses of how treatment response
  might differ by genotype.
- **Grounding LLMs and RAG.** Typed edges with provenance can be used as retrievable, citable facts for AD
  question answering.

### 6. Limitations to state in a paper

- The current release has six node types, not the ten biological scales in PrimeKG. Pathway, GO, anatomy, PPI
  and side-effect layers only arrive once the PrimeKG layer is attached.
- DisGeNET and NIAGADS aren't integrated yet.
- Some noise remains:
  - Free-text trial conditions become disease nodes.
  - Some PubTator genes are mouse orthologs (`App` vs `APP`).
  - Drugs are matched by name, not to DrugBank or ChEMBL IDs.
- The literature layer is small (about 300 papers per run) unless it is run with more pages.
- ClinVar benign rows and PharmGKB "not associated" rows were left out by design.
