# PRIMARY DATA RESOURCE RETRIEVAL AND PROCESSING -- AlzKG
#
# Mirrors PrimeKG's datasets/primary_data_resources.sh, but split into two parts:
#   (A) Generic ontology/pathway/interaction sources that AD shares with every
#       other disease. These are NOT duplicated here -- run PrimeKG's own
#       ../../PrimeKG/datasets/primary_data_resources.sh first (or at least the
#       MONDO/HPO/GO/Reactome/UBERON/CTD/SIDER/NCBI-Gene/UMLS/DisGeNET sections
#       of it) so their processed outputs exist under ../../PrimeKG/datasets/data/.
#       ../knowledge_graph/build_ad_kg.ipynb reads those outputs by relative path
#       and simply skips them if PrimeKG is not present.
#   (B) AD-specific primary sources (this file), which feed the
#       processing_scripts/*.py in this folder.
#
# Run this script from the datasets/ directory:  bash primary_data_resources.sh
# Large downloads are skipped if the file already exists; set REFRESH=1 to re-download.
# Python is invoked as $PYTHON (default: python).

set -e
PYTHON=${PYTHON:-python}

download() {  # download <url> <dest>
    if [ -s "$2" ] && [ -z "$REFRESH" ]; then
        echo "  already have $2 (set REFRESH=1 to re-download)"
    else
        curl -fL --retry 3 "$1" -o "$2"
    fi
}

echo "Making required directories..."
mkdir -p data/gwas_catalog data/clinvar data/disgenet data/aact data/pharmgkb \
         data/literature data/niagads data/ncbi_gene data/kg

# ---------------------------------------------------------------------------
# NCBI GENE (gene identifier harmonization)
# Used by ../knowledge_graph/build_ad_kg.ipynb to map gene symbols / PharmGKB ids
# / PubTator ids onto NCBI Gene ids so all sources share one gene node.
echo "Downloading NCBI Gene human gene_info..."
download "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz" \
         data/ncbi_gene/Homo_sapiens.gene_info.gz

# ---------------------------------------------------------------------------
# GWAS CATALOG
# Database: NHGRI-EBI GWAS Catalog, Script: gwas_catalog.py
# Output: gwas_variant_disease.csv, gwas_gene_disease.csv
# (The old https://www.ebi.ac.uk/gwas/api/search/downloads/full endpoint now 404s;
# releases are published on the EBI FTP instead.)
echo "Downloading GWAS Catalog associations (~200 MB zip)..."
download "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip" \
         data/gwas_catalog/gwas_catalog_associations.zip
$PYTHON processing_scripts/gwas_catalog.py

# ---------------------------------------------------------------------------
# CLINVAR
# Database: NCBI ClinVar, Script: clinvar_ad.py
# Output: clinvar_ad_variant_disease.csv
echo "Downloading ClinVar variant_summary (~450 MB)..."
download "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz" \
         data/clinvar/variant_summary.txt.gz
$PYTHON processing_scripts/clinvar_ad.py

# ---------------------------------------------------------------------------
# DISGENET
# Database: DisGeNET, Script: disgenet_ad.py
# Output: disgenet_ad_gene_disease.csv
# NOTE: DisGeNET now requires an account/API key for bulk downloads
# (https://www.disgenet.org/downloads). Either reuse the file PrimeKG already
# fetched (../../PrimeKG/datasets/data/disgenet/curated_gene_disease_associations.tsv)
# by copying it into data/disgenet/, or download your own copy there.
# The script prints a [skip] message and exits cleanly if the file is absent.
$PYTHON processing_scripts/disgenet_ad.py

# ---------------------------------------------------------------------------
# CLINICALTRIALS.GOV
# Database: ClinicalTrials.gov API v2, Script: clinicaltrials_ad.py
# Output: clinicaltrials_drug_disease.csv, clinicaltrials_trials.csv
# No bulk file to curl -- the script queries the live API directly (AD-specific
# filtering happens server-side via query.cond). Re-run any time to pick up
# newly registered/updated trials (this is the "dynamic ingestion" pattern).
echo "Fetching Alzheimer's disease trials from ClinicalTrials.gov API v2..."
$PYTHON processing_scripts/clinicaltrials_ad.py
# Bulk alternative: download the full AACT Postgres snapshot instead
# (https://aact.ctti-clinicaltrials.org/download), the same way PrimeKG dumps
# DrugCentral's Postgres database in its own primary_data_resources.sh.

# ---------------------------------------------------------------------------
# PHARMGKB (ClinPGx)
# Database: PharmGKB, now ClinPGx (api.pharmgkb.org no longer resolves), Script: pharmgkb_ad.py
# Output: pharmgkb_ad_edges.csv   (the script reads the zip directly)
echo "Downloading PharmGKB/ClinPGx relationships..."
download "https://api.clinpgx.org/v1/download/file/data/relationships.zip" \
         data/pharmgkb/relationships.zip
$PYTHON processing_scripts/pharmgkb_ad.py

# ---------------------------------------------------------------------------
# NIAGADS
# Database: NIAGADS Open Access / Qualified Access, Script: niagads_ad.py
# Output: niagads_variant_disease.csv, niagads_gene_disease.csv
# Manual step required -- see processing_scripts/niagads_ad.py docstring for
# how to request/download a specific accession's GWAS summary statistics, then:
#   python processing_scripts/niagads_ad.py --accession NG00075 --in-path data/niagads/NG00075_summary_stats.tsv
# The script prints a [skip] message and exits cleanly if the file is absent.
$PYTHON processing_scripts/niagads_ad.py

# ---------------------------------------------------------------------------
# PUBTATOR3 LITERATURE MINING (dynamic ingestion layer)
# Database: PubTator3 API, Script: pubtator_literature_ad.py
# Output: literature_triples.csv (append-only; tracks its own checkpoint file)
# No bulk file to curl -- run directly, and re-run periodically (e.g. via cron)
# to keep the literature-derived slice of the graph current. Each run fetches the
# next 30 result pages (300 articles); use --pages N for more.
echo "Mining AD literature via PubTator3 API..."
$PYTHON processing_scripts/pubtator_literature_ad.py

echo "Done. Next: python -m nbconvert --to notebook --execute --inplace ../knowledge_graph/build_ad_kg.ipynb"
