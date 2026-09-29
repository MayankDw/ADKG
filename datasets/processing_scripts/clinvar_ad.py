"""
ClinVar -> Alzheimer's disease Variant-Disease edges.

Source: NCBI ClinVar variant_summary.
Download (see ../primary_data_resources.sh):
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
    -> data/clinvar/variant_summary.txt.gz

Filters variant_summary.txt for rows whose PhenotypeList mentions Alzheimer's
disease (covers early-onset familial AD genes such as APP, PSEN1, PSEN2, as well
as APOE and other risk/protective variants curated in ClinVar), restricted to the
GRCh38 assembly to avoid duplicate GRCh37 rows, and emits a variant-disease edge
list in the harmonized KG edge schema.
"""
import re
import pandas as pd

IN_PATH = "../data/clinvar/variant_summary.txt.gz"
OUT_PATH = "../data/clinvar/clinvar_ad_variant_disease.csv"

AD_PHENOTYPE_PATTERN = re.compile(r"alzheimer", flags=re.IGNORECASE)

USE_COLS = [
    "VariationID", "Name", "GeneSymbol", "ClinicalSignificance",
    "PhenotypeIDS", "PhenotypeList", "Assembly", "RS# (dbSNP)",
]


def load_variant_summary(path=IN_PATH):
    return pd.read_csv(path, sep="\t", usecols=USE_COLS, low_memory=False, compression="infer")


def filter_ad(df):
    df = df.query('Assembly == "GRCh38"').copy()
    mask = df["PhenotypeList"].astype(str).str.contains(AD_PHENOTYPE_PATTERN, na=False)
    return df.loc[mask].copy()


def build_edges(ad_df):
    df = ad_df.dropna(subset=["VariationID", "PhenotypeList"]).copy()

    edges = pd.DataFrame({
        "relation": "variant_disease",
        "display_relation": ad_df["ClinicalSignificance"].str.lower(),
        "x_id": "ClinVar:" + df["VariationID"].astype(str),
        "x_type": "variant",
        "x_name": df["Name"],
        "x_source": "ClinVar",
        "y_id": df["PhenotypeIDS"],
        "y_type": "disease",
        "y_name": df["PhenotypeList"],
        "y_source": "ClinVar",
    })

    gene_edges = df.dropna(subset=["GeneSymbol"]).copy()
    gene_edges = gene_edges.query('GeneSymbol != "-"')
    variant_gene = pd.DataFrame({
        "relation": "variant_gene",
        "display_relation": "located in",
        "x_id": "ClinVar:" + gene_edges["VariationID"].astype(str),
        "x_type": "variant",
        "x_name": gene_edges["Name"],
        "x_source": "ClinVar",
        "y_id": gene_edges["GeneSymbol"],
        "y_type": "gene/protein",
        "y_name": gene_edges["GeneSymbol"],
        "y_source": "ClinVar",
    })

    return pd.concat([edges, variant_gene], ignore_index=True).drop_duplicates()


if __name__ == "__main__":
    variant_summary = load_variant_summary()
    ad_variants = filter_ad(variant_summary)
    print(f"{len(ad_variants)} Alzheimer's-related ClinVar records found "
          f"out of {len(variant_summary)} GRCh38 records.")

    edges = build_edges(ad_variants)
    edges.to_csv(OUT_PATH, index=False)
    print(f"Wrote {len(edges)} edges -> {OUT_PATH}")
