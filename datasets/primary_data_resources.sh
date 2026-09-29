# PRIMARY DATA RESOURCE RETRIEVAL AND PROCESSING -- AlzKG
#
# Mirrors PrimeKG's datasets/primary_data_resources.sh, but split into two parts:
#   (A) Generic ontology/pathway/interaction sources that AD shares with every
#       other disease. These are NOT duplicated here -- run PrimeKG's own
#       ../../PrimeKG/datasets/primary_data_resources.sh first (or at least the
#       MONDO/HPO/GO/Reactome/UBERON/CTD/SIDER/NCBI-Gene/UMLS/DisGeNET sections
#       of it) so their processed outputs exist under ../../PrimeKG/datasets/data/.
#       ../../knowledge_graph/build_ad_kg.ipynb reads those outputs by relative path.
#   (B) AD-specific primary sources (this file), which feed the new
#       processing_scripts/*.py in this folder.
#
# Run this script from the datasets/ directory.

echo "Making required directories..."
mkdir -p data/gwas_catalog data/clinvar data/disgenet data/aact data/pharmgkb \
         data/literature data/niagads data/kg data/kg/auxillary

# ---------------------------------------------------------------------------
# GWAS CATALOG
# Database: NHGRI-EBI GWAS Catalog, Script: gwas_catalog.py
# Output: gwas_variant_disease.csv, gwas_gene_disease.csv
echo "Downloading GWAS Catalog associations..."
curl -L "https://www.ebi.ac.uk/gwas/api/search/downloads/full" -o data/gwas_catalog/gwas_catalog_associations.tsv
cd processing_scripts
python gwas_catalog.py
cd ..

# ---------------------------------------------------------------------------
# CLINVAR
# Database: NCBI ClinVar, Script: clinvar_ad.py
# Output: clinvar_ad_variant_disease.csv
echo "Downloading ClinVar variant_summary..."
curl -L "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz" -o data/clinvar/variant_summary.txt.gz
cd processing_scripts
python clinvar_ad.py
cd ..

# ---------------------------------------------------------------------------
# DISGENET
# Database: DisGeNET, Script: disgenet_ad.py
# Output: disgenet_ad_gene_disease.csv
# NOTE: DisGeNET now requires a free account/API key for bulk downloads
# (https://www.disgenet.org/downloads). Either reuse the file PrimeKG already
# fetched (../../PrimeKG/datasets/data/disgenet/curated_gene_disease_associations.tsv)
# by symlinking/copying it into data/disgenet/, or download your own copy there.
echo "Skipping DisGeNET auto-download (requires login) -- place"
echo "curated_gene_disease_associations.tsv in data/disgenet/ manually, then:"
echo "  cd processing_scripts && python disgenet_ad.py && cd .."

# ---------------------------------------------------------------------------
# CLINICALTRIALS.GOV (AACT)
# Database: ClinicalTrials.gov API v2, Script: clinicaltrials_ad.py
# Output: clinicaltrials_drug_disease.csv, clinicaltrials_trials.csv
# No bulk file to curl -- the script queries the live API directly (AD-specific
# filtering happens server-side via query.cond). Re-run any time to pick up
# newly registered/updated trials (this is the "dynamic ingestion" pattern).
echo "Fetching Alzheimer's disease trials from ClinicalTrials.gov API v2..."
cd processing_scripts
python clinicaltrials_ad.py
cd ..
# Bulk alternative: download the full AACT Postgres snapshot instead
# (https://aact.ctti-clinicaltrials.org/download), the same way PrimeKG dumps
# DrugCentral's Postgres database in its own primary_data_resources.sh.

# ---------------------------------------------------------------------------
# PHARMGKB
# Database: PharmGKB, Script: pharmgkb_ad.py
# Output: pharmgkb_ad_edges.csv
echo "Downloading PharmGKB relationships..."
curl -L "https://api.pharmgkb.org/v1/download/file/data/relationships.zip" -o data/pharmgkb/relationships.zip
unzip -o data/pharmgkb/relationships.zip -d data/pharmgkb/
cd processing_scripts
python pharmgkb_ad.py
cd ..

# ---------------------------------------------------------------------------
# NIAGADS
# Database: NIAGADS Open Access / Qualified Access, Script: niagads_ad.py
# Output: niagads_variant_disease.csv, niagads_gene_disease.csv
# Manual step required -- see processing_scripts/niagads_ad.py docstring for
# how to request/download a specific accession's GWAS summary statistics.
echo "Skipping NIAGADS auto-download (manual catalog request required) -- see"
echo "processing_scripts/niagads_ad.py docstring, then:"
echo "  cd processing_scripts && python niagads_ad.py --accession NG00075 --in-path ../data/niagads/NG00075_summary_stats.tsv && cd .."

# ---------------------------------------------------------------------------
# PUBTATOR3 LITERATURE MINING (dynamic ingestion layer)
# Database: PubTator3 API, Script: pubtator_literature_ad.py
# Output: literature_triples.csv (append-only; tracks its own checkpoint file)
# No bulk file to curl -- run directly, and re-run periodically (e.g. via cron)
# to keep the literature-derived slice of the graph current.
echo "Mining AD literature via PubTator3 API..."
cd processing_scripts
python pubtator_literature_ad.py
cd ..

echo "Done. See ../knowledge_graph/build_ad_kg.ipynb to harmonize these outputs"
echo "together with PrimeKG's generic ontology outputs into ad_kg.csv."
