"""
PharmGKB -> Alzheimer's disease pharmacogenomic edges (Gene-Drug-Disease).

Source: PharmGKB bulk downloads (no login required for the core TSVs).
Download (see ../primary_data_resources.sh):
    https://api.pharmgkb.org/v1/download/file/data/relationships.zip
    -> data/pharmgkb/relationships.zip -> data/pharmgkb/relationships.tsv

    https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip
    -> data/pharmgkb/clinicalAnnotations.zip -> data/pharmgkb/clinical_annotations.tsv

Filters PharmGKB's relationships file to rows where either entity mentions
Alzheimer's disease (e.g. APOE-donepezil response annotations, CYP2D6-related
cholinesterase-inhibitor metabolism) and emits Gene-Drug and Drug-Disease edges
in the harmonized KG schema.
"""
import re
import pandas as pd

IN_RELATIONSHIPS = "../data/pharmgkb/relationships.tsv"
OUT_PATH = "../data/pharmgkb/pharmgkb_ad_edges.csv"

AD_PATTERN = re.compile(r"alzheimer", flags=re.IGNORECASE)

ENTITY_TYPE_MAP = {
    "Gene": "gene/protein",
    "Chemical": "drug/chemical/compound",
    "Disease": "disease",
    "Variant": "variant",
}


def load_relationships(path=IN_RELATIONSHIPS):
    return pd.read_csv(path, sep="\t", low_memory=False)


def filter_ad(df):
    mask = (
        df["Entity1_name"].astype(str).str.contains(AD_PATTERN, na=False)
        | df["Entity2_name"].astype(str).str.contains(AD_PATTERN, na=False)
    )
    return df.loc[mask].copy()


def build_edges(ad_df):
    df = ad_df.copy()

    def node_type(pharmgkb_type):
        return ENTITY_TYPE_MAP.get(pharmgkb_type, "generic_entity")

    edges = pd.DataFrame({
        "relation": df["Association"].str.lower().fillna("associated_with"),
        "display_relation": df["Association"],
        "x_id": df["Entity1_id"],
        "x_type": df["Entity1_type"].map(node_type),
        "x_name": df["Entity1_name"],
        "x_source": "PharmGKB",
        "y_id": df["Entity2_id"],
        "y_type": df["Entity2_type"].map(node_type),
        "y_name": df["Entity2_name"],
        "y_source": "PharmGKB",
    })
    return edges.drop_duplicates()


if __name__ == "__main__":
    relationships = load_relationships()
    ad_relationships = filter_ad(relationships)
    print(f"{len(ad_relationships)} Alzheimer's-related PharmGKB relationships found "
          f"out of {len(relationships)} total.")

    edges = build_edges(ad_relationships)
    edges.to_csv(OUT_PATH, index=False)
    print(f"Wrote {len(edges)} edges -> {OUT_PATH}")
