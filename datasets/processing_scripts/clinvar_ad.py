"""
ClinVar -> Alzheimer's disease Variant-Disease and Variant-Gene edges.

Source: NCBI ClinVar variant_summary.
Download (see ../primary_data_resources.sh):
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
    -> data/clinvar/variant_summary.txt.gz

Filters variant_summary.txt for rows whose PhenotypeList mentions Alzheimer's
disease (covers early-onset familial AD genes such as APP, PSEN1, PSEN2, as well
as APOE and other risk/protective variants curated in ClinVar), restricted to the
GRCh38 assembly to avoid duplicate GRCh37 rows, and emits edge lists in the
harmonized KG edge schema.

Notes on the ClinVar format handled here:
  - PhenotypeList / PhenotypeIDS hold several phenotypes separated by '|' and ';'
    (same layout in both columns). They are split and paired, and only the AD
    phenotypes become edge targets (identified by MONDO id when ClinVar has one).
  - Variants classified Benign / Likely benign are dropped: they are evidence of
    *no* disease effect and should not become "associated with AD" edges.
  - Variants are identified by dbSNP rsID when ClinVar has one (so they join with
    GWAS Catalog / NIAGADS / PubTator variant nodes), else by ClinVar VariationID.
"""
import re
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "clinvar"
IN_PATH = DATA_DIR / "variant_summary.txt.gz"
OUT_PATH = DATA_DIR / "clinvar_ad_variant_disease.csv"

AD_PHENOTYPE_PATTERN = re.compile(r"alzheimer", flags=re.IGNORECASE)
BENIGN_PATTERN = re.compile(r"^(benign|likely benign|benign/likely benign)$", flags=re.IGNORECASE)

USE_COLS = [
    "VariationID", "Name", "GeneSymbol", "ClinicalSignificance",
    "PhenotypeIDS", "PhenotypeList", "Assembly", "RS# (dbSNP)",
]


def load_variant_summary(path=IN_PATH, chunksize=500_000):
    """Stream the ~450 MB gzip in chunks, keeping only GRCh38 AD rows (keeps memory low)."""
    total, kept = 0, []
    for chunk in pd.read_csv(path, sep="\t", usecols=USE_COLS, low_memory=False,
                             compression="infer", chunksize=chunksize):
        chunk = chunk[chunk["Assembly"] == "GRCh38"]
        total += len(chunk)
        kept.append(filter_ad(chunk))
    return pd.concat(kept, ignore_index=True), total


def filter_ad(df):
    df = df[df["Assembly"] == "GRCh38"]
    mask = df["PhenotypeList"].astype(str).str.contains(AD_PHENOTYPE_PATTERN, na=False)
    mask &= ~df["ClinicalSignificance"].astype(str).str.match(BENIGN_PATTERN)
    return df.loc[mask].copy()


def preferred_id(id_group):
    """'MONDO:MONDO:0004975,MeSH:D000544,MedGen:C0002395' -> 'MONDO:0004975' (else first id)."""
    ids = [i.strip() for i in id_group.split(",") if i.strip()]
    for i in ids:
        if i.startswith("MONDO:MONDO:"):
            return i.replace("MONDO:MONDO:", "MONDO:")
    return ids[0] if ids else None


def ad_phenotypes(phenotype_list, phenotype_ids):
    names = re.split(r"[|;]", str(phenotype_list))
    ids = re.split(r"[|;]", str(phenotype_ids))
    pairs = zip(names, ids) if len(names) == len(ids) else ((n, None) for n in names)
    out = []
    for name, id_group in pairs:
        if AD_PHENOTYPE_PATTERN.search(name):
            out.append((preferred_id(id_group) if id_group else name.strip(), name.strip()))
    return out


def variant_id(row):
    rs = row["RS# (dbSNP)"]
    return f"rs{int(rs)}" if pd.notna(rs) and int(rs) > 0 else f"ClinVar:{row['VariationID']}"


def build_edges(ad_df):
    df = ad_df.dropna(subset=["VariationID", "PhenotypeList"]).copy()
    df["variant_id"] = df.apply(variant_id, axis=1)

    vd_rows = []
    for _, r in df.iterrows():
        for pid, pname in ad_phenotypes(r["PhenotypeList"], r["PhenotypeIDS"]):
            vd_rows.append({
                "relation": "variant_disease",
                "display_relation": str(r["ClinicalSignificance"]).lower(),
                "x_id": r["variant_id"],
                "x_type": "variant",
                "x_name": r["Name"],
                "x_source": "ClinVar",
                "y_id": pid,
                "y_type": "disease",
                "y_name": pname,
                "y_source": "ClinVar",
            })
    variant_disease = pd.DataFrame(vd_rows)

    gene_df = df.dropna(subset=["GeneSymbol"]).copy()
    gene_df["GeneSymbol"] = gene_df["GeneSymbol"].str.split(";")
    gene_df = gene_df.explode("GeneSymbol").query('GeneSymbol != "-"')
    variant_gene = pd.DataFrame({
        "relation": "variant_gene",
        "display_relation": "located in",
        "x_id": gene_df["variant_id"],
        "x_type": "variant",
        "x_name": gene_df["Name"],
        "x_source": "ClinVar",
        "y_id": gene_df["GeneSymbol"],
        "y_type": "gene/protein",
        "y_name": gene_df["GeneSymbol"],
        "y_source": "ClinVar",
    })

    return pd.concat([variant_disease, variant_gene], ignore_index=True).drop_duplicates()


if __name__ == "__main__":
    ad_variants, n_grch38 = load_variant_summary()
    print(f"{len(ad_variants)} Alzheimer's-related (non-benign) ClinVar records found "
          f"out of {n_grch38} GRCh38 records.")

    edges = build_edges(ad_variants)
    edges.to_csv(OUT_PATH, index=False)
    print(f"Wrote {len(edges)} edges -> {OUT_PATH}")
    print(edges["relation"].value_counts().to_string())
