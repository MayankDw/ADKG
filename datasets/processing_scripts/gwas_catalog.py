"""
GWAS Catalog -> Alzheimer's disease Variant/Gene-Disease edges.

Source: NHGRI-EBI GWAS Catalog "all associations" TSV.
Download (see ../primary_data_resources.sh):
    https://www.ebi.ac.uk/gwas/api/search/downloads/full -> data/gwas_catalog/gwas_catalog_associations.tsv

This script filters the full associations file down to Alzheimer's-disease-relevant
traits (matched on DISEASE/TRAIT and MAPPED_TRAIT text) and emits two edge lists:
  - gwas_variant_disease.csv : SNPS (rsID) -> DISEASE/TRAIT
  - gwas_gene_disease.csv    : MAPPED_GENE -> DISEASE/TRAIT

Output columns match the harmonized KG edge schema used in ../../knowledge_graph/build_ad_kg.ipynb:
    relation, display_relation, x_id, x_type, x_name, x_source, y_id, y_type, y_name, y_source
"""
import re
import numpy as np
import pandas as pd

IN_PATH = "../data/gwas_catalog/gwas_catalog_associations.tsv"
OUT_VARIANT_DISEASE = "../data/gwas_catalog/gwas_variant_disease.csv"
OUT_GENE_DISEASE = "../data/gwas_catalog/gwas_gene_disease.csv"

# Trait strings that identify Alzheimer's-disease-relevant GWAS Catalog entries.
# Matches published EFO trait labels such as "Alzheimer's disease",
# "Alzheimer's disease biomarkers", "late-onset Alzheimer's disease",
# "Alzheimer's disease and age at onset", etc.
AD_TRAIT_PATTERN = re.compile(r"alzheimer", flags=re.IGNORECASE)


def load_associations(path=IN_PATH):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    return df


def filter_ad(df):
    trait_cols = [c for c in ["DISEASE/TRAIT", "MAPPED_TRAIT"] if c in df.columns]
    mask = np.zeros(len(df), dtype=bool)
    for col in trait_cols:
        mask |= df[col].astype(str).str.contains(AD_TRAIT_PATTERN, na=False)
    return df.loc[mask].copy()


def build_variant_disease_edges(ad_df):
    df = ad_df.dropna(subset=["SNPS", "DISEASE/TRAIT"]).copy()
    df = df.get(["SNPS", "DISEASE/TRAIT", "MAPPED_TRAIT_URI", "PVALUE_MLOG", "PUBMEDID"]).drop_duplicates()

    edges = pd.DataFrame({
        "relation": "variant_disease",
        "display_relation": "associated with",
        "x_id": df["SNPS"],
        "x_type": "variant",
        "x_name": df["SNPS"],
        "x_source": "GWAS_Catalog",
        "y_id": df["MAPPED_TRAIT_URI"].fillna(df["DISEASE/TRAIT"]),
        "y_type": "disease",
        "y_name": df["DISEASE/TRAIT"],
        "y_source": "GWAS_Catalog",
    })
    return edges.drop_duplicates()


def build_gene_disease_edges(ad_df):
    df = ad_df.dropna(subset=["MAPPED_GENE", "DISEASE/TRAIT"]).copy()
    df = df.get(["MAPPED_GENE", "DISEASE/TRAIT", "MAPPED_TRAIT_URI"]).drop_duplicates()

    # MAPPED_GENE can contain multiple genes separated by ", " or " - " (intergenic).
    rows = []
    for _, r in df.iterrows():
        genes = re.split(r",\s*|\s*-\s*", str(r["MAPPED_GENE"]))
        for gene in genes:
            gene = gene.strip()
            if gene:
                rows.append((gene, r["DISEASE/TRAIT"], r["MAPPED_TRAIT_URI"]))
    gene_df = pd.DataFrame(rows, columns=["gene", "disease", "trait_uri"]).drop_duplicates()

    edges = pd.DataFrame({
        "relation": "gene_disease",
        "display_relation": "associated with",
        "x_id": gene_df["gene"],
        "x_type": "gene/protein",
        "x_name": gene_df["gene"],
        "x_source": "GWAS_Catalog",
        "y_id": gene_df["trait_uri"].fillna(gene_df["disease"]),
        "y_type": "disease",
        "y_name": gene_df["disease"],
        "y_source": "GWAS_Catalog",
    })
    return edges.drop_duplicates()


if __name__ == "__main__":
    associations = load_associations()
    ad_associations = filter_ad(associations)
    print(f"{len(ad_associations)} Alzheimer's-related GWAS Catalog associations found "
          f"out of {len(associations)} total.")

    variant_disease = build_variant_disease_edges(ad_associations)
    gene_disease = build_gene_disease_edges(ad_associations)

    variant_disease.to_csv(OUT_VARIANT_DISEASE, index=False)
    gene_disease.to_csv(OUT_GENE_DISEASE, index=False)
    print(f"Wrote {len(variant_disease)} variant-disease edges -> {OUT_VARIANT_DISEASE}")
    print(f"Wrote {len(gene_disease)} gene-disease edges -> {OUT_GENE_DISEASE}")
