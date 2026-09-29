"""
DisGeNET -> Alzheimer's disease Gene-Disease edges.

Source: DisGeNET curated_gene_disease_associations.tsv (same file PrimeKG uses;
reuse ../../../PrimeKG/datasets/data/disgenet/curated_gene_disease_associations.tsv
if you already ran PrimeKG's pipeline, or fetch a fresh copy per
../primary_data_resources.sh -> data/disgenet/curated_gene_disease_associations.tsv).
Note: DisGeNET now requires a free API key/login for bulk downloads
(https://www.disgenet.org/downloads) -- this script only parses a file you
already have locally.

Filters to UMLS CUIs / disease names corresponding to Alzheimer's disease and its
subtypes (early-onset/familial AD, late-onset AD, AD with cerebral amyloid
angiopathy, etc.) and emits a gene-disease edge list in the harmonized KG schema.
"""
import re
import pandas as pd

IN_PATH = "../data/disgenet/curated_gene_disease_associations.tsv"
OUT_PATH = "../data/disgenet/disgenet_ad_gene_disease.csv"

# Known UMLS CUIs for Alzheimer's disease and closely related concepts.
AD_CUIS = {
    "C0002395",  # Alzheimer's Disease
    "C0494463",  # Alzheimer Disease, Late Onset
    "C1863051",  # Alzheimer Disease, Early Onset
    "C0276496",  # Alzheimer Disease, Familial
    "C0338451",  # Dementia of the Alzheimer type
}
AD_NAME_PATTERN = re.compile(r"alzheimer", flags=re.IGNORECASE)


def load_associations(path=IN_PATH):
    return pd.read_csv(path, sep="\t", low_memory=False)


def filter_ad(df):
    mask = df["diseaseId"].isin(AD_CUIS)
    if "diseaseName" in df.columns:
        mask |= df["diseaseName"].astype(str).str.contains(AD_NAME_PATTERN, na=False)
    return df.loc[mask].copy()


def build_edges(ad_df):
    df = ad_df.dropna(subset=["geneId", "diseaseId"]).copy()

    edges = pd.DataFrame({
        "relation": "gene_disease",
        "display_relation": "associated with",
        "x_id": df["geneId"].astype(int).astype(str),
        "x_type": "gene/protein",
        "x_name": df["geneSymbol"],
        "x_source": "DisGeNET",
        "y_id": df["diseaseId"],
        "y_type": "disease",
        "y_name": df["diseaseName"],
        "y_source": "DisGeNET",
    })
    if "score" in ad_df.columns:
        edges["evidence_score"] = ad_df["score"].values
    return edges.drop_duplicates()


if __name__ == "__main__":
    associations = load_associations()
    ad_associations = filter_ad(associations)
    print(f"{len(ad_associations)} Alzheimer's-related DisGeNET associations found "
          f"out of {len(associations)} total curated associations.")

    edges = build_edges(ad_associations)
    edges.to_csv(OUT_PATH, index=False)
    print(f"Wrote {len(edges)} gene-disease edges -> {OUT_PATH}")
