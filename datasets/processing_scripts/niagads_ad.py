"""
NIAGADS -> Alzheimer's disease GWAS summary-statistics Variant/Gene-Disease edges.

Source: NIAGADS Open Access GWAS summary statistics (e.g. Kunkle et al. 2019
stage 1 meta-analysis, NIAGADS accession NG00075; Bellenguez et al. 2022,
NG00115). These are catalog entries, not a single stable URL, so -- like
DrugBank and UMLS in PrimeKG's own pipeline -- the file must be requested/
downloaded manually first:
    1. Browse https://www.niagads.org/opendata (or apply for Qualified Access at
       https://www.niagads.org/adgc/data for individual-level ADGC/ADSP data).
    2. Download the summary-statistics file for your chosen accession to
       data/niagads/<accession>_summary_stats.tsv (gzip is also handled).
    3. Set ACCESSION and IN_PATH below (or pass --accession/--in-path) and run
       this script.

Expects the common GWAS summary-stats column layout used by IGAP/NIAGADS
releases: Chromosome, Position, MarkerName (rsID), Effect_allele,
Non_Effect_allele, Beta, SE, Pvalue, and (when present) Nearest_Gene.
Column names are matched case-insensitively and with light fuzzing since exact
headers vary slightly by NIAGADS release.
"""
import argparse
import re
import sys
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "niagads"
ACCESSION = "NG00075"  # Kunkle et al. 2019 IGAP stage 1, by default
IN_PATH = DATA_DIR / "NG00075_summary_stats.tsv"
OUT_VARIANT_DISEASE = DATA_DIR / "niagads_variant_disease.csv"
OUT_GENE_DISEASE = DATA_DIR / "niagads_gene_disease.csv"

# Anchor disease node for this whole dataset (genome-wide AD case-control GWAS).
AD_DISEASE_ID = "MONDO:0004975"
AD_DISEASE_NAME = "Alzheimer's disease"

PVALUE_GENOME_WIDE_SIGNIFICANT = 5e-8

COLUMN_ALIASES = {
    "marker": "MarkerName", "snp": "MarkerName", "rsid": "MarkerName",
    "pvalue": "Pvalue", "p_value": "Pvalue", "p": "Pvalue",
    "nearest_gene": "Nearest_Gene", "gene": "Nearest_Gene", "mapped_gene": "Nearest_Gene",
}


def _normalize_columns(df):
    rename = {}
    for col in df.columns:
        key = re.sub(r"[^a-z]", "", col.lower())
        for alias, canonical in COLUMN_ALIASES.items():
            if key == re.sub(r"[^a-z]", "", alias.lower()):
                rename[col] = canonical
    return df.rename(columns=rename)


def load_summary_stats(path):
    df = pd.read_csv(path, sep=None, engine="python", compression="infer")
    return _normalize_columns(df)


def filter_genome_wide_significant(df, threshold=PVALUE_GENOME_WIDE_SIGNIFICANT):
    pvalues = pd.to_numeric(df["Pvalue"], errors="coerce")
    return df.loc[pvalues <= threshold].assign(Pvalue=pvalues).copy()


def build_variant_disease_edges(sig_df, accession):
    df = sig_df.dropna(subset=["MarkerName"]).copy()
    edges = pd.DataFrame({
        "relation": "variant_disease",
        "display_relation": "associated with (GWAS)",
        "x_id": df["MarkerName"],
        "x_type": "variant",
        "x_name": df["MarkerName"],
        "x_source": f"NIAGADS:{accession}",
        "y_id": AD_DISEASE_ID,
        "y_type": "disease",
        "y_name": AD_DISEASE_NAME,
        "y_source": f"NIAGADS:{accession}",
        "pvalue": df["Pvalue"].values,
    })
    return edges.drop_duplicates()


def build_gene_disease_edges(sig_df, accession):
    if "Nearest_Gene" not in sig_df.columns:
        return pd.DataFrame()
    df = sig_df.dropna(subset=["Nearest_Gene"]).copy()
    edges = pd.DataFrame({
        "relation": "gene_disease",
        "display_relation": "associated with (GWAS, nearest gene)",
        "x_id": df["Nearest_Gene"],
        "x_type": "gene/protein",
        "x_name": df["Nearest_Gene"],
        "x_source": f"NIAGADS:{accession}",
        "y_id": AD_DISEASE_ID,
        "y_type": "disease",
        "y_name": AD_DISEASE_NAME,
        "y_source": f"NIAGADS:{accession}",
    })
    return edges.drop_duplicates()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--accession", default=ACCESSION)
    parser.add_argument("--in-path", default=IN_PATH)
    args = parser.parse_args()

    if not Path(args.in_path).exists():
        print(f"[skip] {args.in_path} not found. NIAGADS summary statistics must be requested and\n"
              "       downloaded manually (see this script's docstring), then re-run with --in-path.")
        sys.exit(0)

    summary_stats = load_summary_stats(args.in_path)
    sig = filter_genome_wide_significant(summary_stats)
    print(f"{len(sig)} genome-wide-significant (p<={PVALUE_GENOME_WIDE_SIGNIFICANT}) variants "
          f"out of {len(summary_stats)} total rows in {args.accession}.")

    variant_disease = build_variant_disease_edges(sig, args.accession)
    gene_disease = build_gene_disease_edges(sig, args.accession)

    variant_disease.to_csv(OUT_VARIANT_DISEASE, index=False)
    print(f"Wrote {len(variant_disease)} variant-disease edges -> {OUT_VARIANT_DISEASE}")
    if len(gene_disease):
        gene_disease.to_csv(OUT_GENE_DISEASE, index=False)
        print(f"Wrote {len(gene_disease)} gene-disease edges -> {OUT_GENE_DISEASE}")
