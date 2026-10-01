"""
GWAS Catalog -> Alzheimer's disease Variant/Gene-Disease edges.

Source: NHGRI-EBI GWAS Catalog "all associations (ontology-annotated)" release.
Download (see ../primary_data_resources.sh):
    https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip
    -> data/gwas_catalog/gwas_catalog_associations.zip
(The old https://www.ebi.ac.uk/gwas/api/search/downloads/full endpoint now returns 404.)

This script filters the full associations file down to Alzheimer's-disease-relevant
traits and emits two edge lists:
  - gwas_variant_disease.csv : SNPS (rsID) -> mapped AD trait (MONDO/EFO term)
  - gwas_gene_disease.csv    : MAPPED_GENE -> mapped AD trait (MONDO/EFO term)

A row's MAPPED_TRAIT / MAPPED_TRAIT_URI can list several ontology terms (e.g.
"Alzheimer disease, educational attainment" for a pleiotropy study). Only the
terms that are themselves Alzheimer's terms are kept as edge targets, so the
graph stays AD-anchored.

Output columns match the harmonized KG edge schema used in ../../knowledge_graph/build_ad_kg.ipynb:
    relation, display_relation, x_id, x_type, x_name, x_source, y_id, y_type, y_name, y_source
"""
import re
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "gwas_catalog"
IN_PATH = DATA_DIR / "gwas_catalog_associations.zip"
OUT_VARIANT_DISEASE = DATA_DIR / "gwas_variant_disease.csv"
OUT_GENE_DISEASE = DATA_DIR / "gwas_gene_disease.csv"

USE_COLS = ["DISEASE/TRAIT", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "SNPS", "MAPPED_GENE",
            "PVALUE_MLOG", "PUBMEDID"]

# Matches published trait labels such as "Alzheimer disease", "late-onset Alzheimers
# disease", "family history of Alzheimer's disease", "Alzheimer's disease biomarker measurement".
AD_TRAIT_PATTERN = re.compile(r"alzheimer", flags=re.IGNORECASE)


def load_associations(path=IN_PATH):
    # The release file contains a few non-UTF-8 bytes (curly apostrophes).
    return pd.read_csv(path, sep="\t", usecols=USE_COLS, low_memory=False,
                       encoding_errors="replace")


def filter_ad(df):
    mask = (df["DISEASE/TRAIT"].astype(str).str.contains(AD_TRAIT_PATTERN, na=False)
            | df["MAPPED_TRAIT"].astype(str).str.contains(AD_TRAIT_PATTERN, na=False))
    return df.loc[mask].copy()


def uri_to_curie(uri):
    """http://purl.obolibrary.org/obo/MONDO_0004975 -> MONDO:0004975"""
    last = str(uri).strip().rstrip("/").split("/")[-1]
    return last.replace("_", ":", 1) if "_" in last else last


def clean_label(label):
    return re.sub(r"[�’`]", "'", str(label)).strip()


def ad_traits(row):
    """Return [(curie, label, node_type)] for the AD-specific mapped traits of one row."""
    labels = [clean_label(t) for t in str(row["MAPPED_TRAIT"]).split(", ")]
    uris = [u.strip() for u in str(row["MAPPED_TRAIT_URI"]).split(",")]
    if len(labels) != len(uris):  # a label itself contained ", " -- can't pair reliably
        labels = [clean_label(row["MAPPED_TRAIT"])] * len(uris)
    out = []
    for label, uri in zip(labels, uris):
        if AD_TRAIT_PATTERN.search(label) or uri.endswith(("MONDO_0004975", "EFO_0000249")):
            node_type = "disease" if "MONDO_" in uri or label.lower().endswith("disease") else "effect/phenotype"
            out.append((uri_to_curie(uri), label, node_type))
    return out


def split_snps(snps):
    """'rs1 x rs2' (interaction) / 'rs1; rs2' (haplotype) -> ['rs1', 'rs2']"""
    return [s.strip() for s in re.split(r";|,|\s+x\s+", str(snps)) if s.strip() and s.strip() != "nan"]


def split_genes(mapped_gene):
    """'A, B' (multiple) / 'A - B' (intergenic, flanking genes) / 'A; A' (haplotype) /
    'A x B' (interaction) -> ['A', 'B'].
    Hyphenated symbols such as HLA-DRB1 or MMADHC-DT are kept intact."""
    genes = [g.strip() for g in re.split(r",\s*|;\s*|\s+-\s+|\s+x\s+", str(mapped_gene))]
    return list(dict.fromkeys(g for g in genes if g and g != "nan"))


def explode(ad_df):
    rows = []
    for _, r in ad_df.iterrows():
        for curie, label, node_type in ad_traits(r):
            rows.append({"snps": split_snps(r["SNPS"]), "genes": split_genes(r["MAPPED_GENE"]),
                         "trait_id": curie, "trait": label, "trait_type": node_type})
    return pd.DataFrame(rows)


def _edges(df, entity_col, relation, x_type):
    df = df.explode(entity_col).dropna(subset=[entity_col])
    df = df.drop_duplicates(subset=[entity_col, "trait_id"])
    return pd.DataFrame({
        "relation": relation,
        "display_relation": "associated with",
        "x_id": df[entity_col],
        "x_type": x_type,
        "x_name": df[entity_col],
        "x_source": "GWAS_Catalog",
        "y_id": df["trait_id"],
        "y_type": df["trait_type"],
        "y_name": df["trait"],
        "y_source": "GWAS_Catalog",
    })


def build_variant_disease_edges(exploded):
    return _edges(exploded, "snps", "variant_disease", "variant")


def build_gene_disease_edges(exploded):
    return _edges(exploded, "genes", "gene_disease", "gene/protein")


if __name__ == "__main__":
    associations = load_associations()
    ad_associations = filter_ad(associations)
    print(f"{len(ad_associations)} Alzheimer's-related GWAS Catalog associations found "
          f"out of {len(associations)} total.")

    exploded = explode(ad_associations)
    variant_disease = build_variant_disease_edges(exploded)
    gene_disease = build_gene_disease_edges(exploded)

    variant_disease.to_csv(OUT_VARIANT_DISEASE, index=False)
    gene_disease.to_csv(OUT_GENE_DISEASE, index=False)
    print(f"Wrote {len(variant_disease)} variant-disease edges -> {OUT_VARIANT_DISEASE}")
    print(f"Wrote {len(gene_disease)} gene-disease edges -> {OUT_GENE_DISEASE}")
