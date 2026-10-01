"""
PharmGKB (now ClinPGx) -> Alzheimer's disease pharmacogenomic edges (Gene/Variant/Drug-Disease).

Source: ClinPGx bulk downloads (no login required for the core TSVs). PharmGKB
was rebranded ClinPGx in 2025 and the old api.pharmgkb.org host no longer resolves.
Download (see ../primary_data_resources.sh):
    https://api.clinpgx.org/v1/download/file/data/relationships.zip
    -> data/pharmgkb/relationships.zip   (read directly; no unzip step needed)

Filters the relationships file to rows where either entity is Alzheimer's
disease (e.g. APOE / CYP2D6 haplotypes associated with AD or with
cholinesterase-inhibitor response) and emits edges in the harmonized KG schema.
It then adds the Gene-Chemical (drug response / metabolism) relationships of
every gene that is itself linked to AD, giving the drug layer of the graph.
Rows whose Association is "not associated" are dropped -- they record a tested,
negative result and should not become graph edges.
"""
import re
import zipfile
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "pharmgkb"
IN_ZIP = DATA_DIR / "relationships.zip"
IN_RELATIONSHIPS = DATA_DIR / "relationships.tsv"
OUT_PATH = DATA_DIR / "pharmgkb_ad_edges.csv"

AD_PATTERN = re.compile(r"alzheimer", flags=re.IGNORECASE)

ENTITY_TYPE_MAP = {
    "Gene": "gene/protein",
    "Chemical": "drug/chemical/compound",
    "Disease": "disease",
    "Variant": "variant",
    "Haplotype": "variant",
}
# Short tags used to build PrimeKG-style relation names, e.g. gene_disease.
RELATION_TAG = {
    "gene/protein": "gene", "drug/chemical/compound": "drug",
    "disease": "disease", "variant": "variant",
}


def load_relationships():
    if IN_RELATIONSHIPS.exists():
        return pd.read_csv(IN_RELATIONSHIPS, sep="\t", low_memory=False)
    with zipfile.ZipFile(IN_ZIP) as z:
        with z.open("relationships.tsv") as f:
            return pd.read_csv(f, sep="\t", low_memory=False)


def filter_ad(df):
    mask = (
        df["Entity1_name"].astype(str).str.contains(AD_PATTERN, na=False)
        | df["Entity2_name"].astype(str).str.contains(AD_PATTERN, na=False)
    )
    mask &= df["Association"].astype(str).str.lower() != "not associated"
    return df.loc[mask].copy()


def ad_gene_drug_rows(relationships, ad_df):
    """Gene-Chemical relationships for the genes that have an AD edge."""
    ad_gene_ids = (set(ad_df.loc[ad_df["Entity1_type"] == "Gene", "Entity1_id"])
                   | set(ad_df.loc[ad_df["Entity2_type"] == "Gene", "Entity2_id"]))
    is_gene_drug = [{a, b} == {"Gene", "Chemical"} for a, b in
                    zip(relationships["Entity1_type"], relationships["Entity2_type"])]
    df = relationships[is_gene_drug]
    df = df[df["Entity1_id"].isin(ad_gene_ids) | df["Entity2_id"].isin(ad_gene_ids)]
    return df[df["Association"].astype(str).str.lower() != "not associated"]


def build_edges(ad_df):
    df = ad_df.copy()
    x_type = df["Entity1_type"].map(lambda t: ENTITY_TYPE_MAP.get(t, "generic_entity"))
    y_type = df["Entity2_type"].map(lambda t: ENTITY_TYPE_MAP.get(t, "generic_entity"))

    # Orient edges as gene/variant -> disease and gene -> drug, matching the other scripts.
    flip = (((x_type == "disease") & (y_type != "disease"))
            | ((x_type == "drug/chemical/compound") & (y_type == "gene/protein")))
    ent1 = ["Entity1_id", "Entity1_name"]
    ent2 = ["Entity2_id", "Entity2_name"]
    df.loc[flip, ent1 + ent2] = df.loc[flip, ent2 + ent1].values
    x_type, y_type = x_type.where(~flip, y_type), y_type.where(~flip, x_type)

    edges = pd.DataFrame({
        "relation": [f"{RELATION_TAG.get(a, 'entity')}_{RELATION_TAG.get(b, 'entity')}"
                     for a, b in zip(x_type, y_type)],
        "display_relation": df["Association"].fillna("associated").str.lower(),
        "x_id": df["Entity1_id"],
        "x_type": x_type,
        "x_name": df["Entity1_name"],
        "x_source": "PharmGKB",
        "y_id": df["Entity2_id"],
        "y_type": y_type,
        "y_name": df["Entity2_name"],
        "y_source": "PharmGKB",
    })
    return edges.drop_duplicates()


if __name__ == "__main__":
    relationships = load_relationships()
    ad_relationships = filter_ad(relationships)
    print(f"{len(ad_relationships)} Alzheimer's-related PharmGKB relationships "
          f"(excluding 'not associated') found out of {len(relationships)} total.")

    gene_drug = ad_gene_drug_rows(relationships, ad_relationships)
    print(f"{len(gene_drug)} Gene-Chemical relationships for those AD-linked genes.")

    edges = build_edges(pd.concat([ad_relationships, gene_drug], ignore_index=True))
    edges.to_csv(OUT_PATH, index=False)
    print(f"Wrote {len(edges)} edges -> {OUT_PATH}")
    print(edges["relation"].value_counts().to_string())
