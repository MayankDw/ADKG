"""
Technical validation (c): recovery of established AD genes and drugs in AlzKG.

Reference sets
  * GWAS loci: Bellenguez et al., Nat Genet 2022 (doi:10.1038/s41588-022-01024-z),
    Table 1 (previously reported loci) and Table 2 (42 new loci), gene names as
    reported. [VERIFY against the published tables before submission.]
  * Mendelian / major-effect genes: APP, PSEN1, PSEN2 (autosomal-dominant AD), APOE.
  * Approved AD drugs (FDA): cholinesterase inhibitors, memantine, anti-amyloid
    antibodies, and drugs approved for AD-related symptoms.

For each reference gene the script reports whether it is a node, whether it links to
the AD set (MONDO:0004975 plus its linked subtype/trait nodes) directly or through an
AD-associated ClinVar variant, and which datasets support that link. Because the GWAS
Catalog is itself an input, recovery is also reported with GWAS evidence removed.
It then tests whether AlzKG's evidence-convergence score (number of independent
datasets supporting a gene-AD link) ranks the reference genes above other AD-linked
genes (ROC AUC, Mann-Whitney U).

Usage:
    python known_ad_recovery.py [path/to/ad_kg.csv] [output_dir]
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu  # U statistic only; p from permutation

HERE = Path(__file__).resolve().parent
KG_PATH = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent / "datasets/data/kg/ad_kg.csv"
OUT_DIR = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "results"
AD_ANCHOR = "MONDO:0004975"

BELLENGUEZ_KNOWN = """CR1 BIN1 INPP5D CLNK HLA-DQA1 UNC5CL TREM2 TREML2 CD2AP EPDR1 SPDYE3 EPHA1
PTK2B CLU USP6NL SPI1 MS4A4A EED SORL1 FERMT2 SLC24A4 SPPL2A MINDY2 APH1B BCKDK IL34
PLCG2 SCIMP WNT3 ABI3 TSPOAP1 ACE ABCA7 CASS4 ADAMTS1""".split()
# "IGH gene cluster" (locus 27) is not a single gene and is excluded.
BELLENGUEZ_NEW = """SORT1 ADAM17 PRKD3 NCK2 WDR12 MME IDUA RHOH ANKH COX7C TNIP1 RASGEF1C HS3ST5
UMAD1 ICA1 TMEM106B JAZF1 SEC61G CTSB SHARPIN ABCA1 ANK3 TSPAN14 BLNK PLEKHA1 TPCN1
SNX1 CTSH DOC2A MAF FOXF1 PRDM7 WDR81 MYO15A GRN KLF16 SIGLEC11 LILRB2 RBCK1 SLC2A4RG
APP""".split()
MENDELIAN = ["APP", "PSEN1", "PSEN2", "APOE"]

# Lead variants as reported in Bellenguez et al. 2022 Tables 1-2 [VERIFY].
LEAD_VARIANTS = {
    "CR1": ["rs679515"], "BIN1": ["rs6733839"], "INPP5D": ["rs10933431"], "CLNK": ["rs6846529"],
    "HLA-DQA1": ["rs6605556"], "UNC5CL": ["rs10947943"], "TREM2": ["rs143332484", "rs75932628"],
    "TREML2": ["rs60755019"], "CD2AP": ["rs7767350"], "EPDR1": ["rs6966331"], "SPDYE3": ["rs7384878"],
    "EPHA1": ["rs11771145"], "PTK2B": ["rs73223431"], "CLU": ["rs11787077"], "USP6NL": ["rs7912495"],
    "SPI1": ["rs10437655"], "MS4A4A": ["rs1582763"], "EED": ["rs3851179"],
    "SORL1": ["rs74685827", "rs11218343"], "FERMT2": ["rs17125924"],
    "SLC24A4": ["rs7401792", "rs12590654"], "SPPL2A": ["rs8025980"], "MINDY2": ["rs602602"],
    "APH1B": ["rs117618017"], "BCKDK": ["rs889555"], "IL34": ["rs4985556"],
    "PLCG2": ["rs12446759", "rs72824905"], "SCIMP": ["rs7225151"], "WNT3": ["rs199515"],
    "ABI3": ["rs616338"], "TSPOAP1": ["rs2526377"], "ACE": ["rs4277405"], "ABCA7": ["rs12151021"],
    "CASS4": ["rs6014724"], "ADAMTS1": ["rs2830489"],
    "SORT1": ["rs141749679"], "ADAM17": ["rs72777026"], "PRKD3": ["rs17020490"],
    "NCK2": ["rs143080277"], "WDR12": ["rs139643391"], "MME": ["rs16824536", "rs61762319"],
    "IDUA": ["rs3822030"], "RHOH": ["rs2245466"], "ANKH": ["rs112403360"], "COX7C": ["rs62374257"],
    "TNIP1": ["rs871269"], "RASGEF1C": ["rs113706587"], "HS3ST5": ["rs785129"], "UMAD1": ["rs6943429"],
    "ICA1": ["rs10952097"], "TMEM106B": ["rs13237518"], "JAZF1": ["rs1160871"], "SEC61G": ["rs76928645"],
    "CTSB": ["rs1065712"], "SHARPIN": ["rs34173062"], "ABCA1": ["rs1800978"], "ANK3": ["rs7068231"],
    "TSPAN14": ["rs6586028"], "BLNK": ["rs6584063"], "PLEKHA1": ["rs7908662"], "TPCN1": ["rs6489896"],
    "IGH gene cluster": ["rs7157106", "rs10131280"], "SNX1": ["rs3848143"], "CTSH": ["rs12592898"],
    "DOC2A": ["rs1140239"], "MAF": ["rs450674"], "FOXF1": ["rs16941239"], "PRDM7": ["rs56407236"],
    "WDR81": ["rs35048651"], "MYO15A": ["rs2242595"], "GRN": ["rs5848"], "KLF16": ["rs149080927"],
    "SIGLEC11": ["rs9304690"], "LILRB2": ["rs587709"], "RBCK1": ["rs1358782"],
    "SLC2A4RG": ["rs6742"], "APP": ["rs2154481"],
}

APPROVED_DRUGS = {
    "donepezil": "cholinesterase inhibitor",
    "galantamine": "cholinesterase inhibitor",
    "rivastigmine": "cholinesterase inhibitor",
    "tacrine": "cholinesterase inhibitor (withdrawn)",
    "memantine": "NMDA receptor antagonist",
    "aducanumab": "anti-amyloid antibody (discontinued)",
    "lecanemab": "anti-amyloid antibody",
    "donanemab": "anti-amyloid antibody",
    "brexpiprazole": "agitation in AD dementia",
    "suvorexant": "insomnia in AD",
}


def norm(s):
    return re.sub(r"[^a-z0-9]+", " ", str(s).lower()).strip()


def load(path):
    kg = pd.read_csv(path, dtype=str)
    if "edge_source" not in kg.columns:
        sys.exit(f"{path} has no edge_source column: it is not the harmonized release build.")
    sub = kg[(kg.edge_source == "AlzKG_harmonization") & (kg.y_id == AD_ANCHOR)]
    ad_set = {AD_ANCHOR} | set(sub.x_id)
    return kg, ad_set


def gene_support(kg, ad_set):
    """gene symbol -> set of datasets supporting its link to the AD set."""
    rows = []
    direct = kg[(kg.x_type == "gene/protein") & kg.y_id.isin(ad_set)]
    rows += list(zip(direct.x_name, direct.edge_source))
    direct = kg[(kg.y_type == "gene/protein") & kg.x_id.isin(ad_set)]
    rows += list(zip(direct.y_name, direct.edge_source))
    # Indirect: gene hosts a variant that is itself linked to the AD set.
    ad_variants = set(kg.loc[(kg.x_type == "variant") & kg.y_id.isin(ad_set), "x_id"])
    vg = kg[(kg.relation == "variant_gene") & kg.x_id.isin(ad_variants)]
    rows += [(g, "ClinVar") for g in vg.y_name]
    df = pd.DataFrame(rows, columns=["gene", "source"])
    return df.groupby("gene")["source"].agg(lambda s: set(s))


def gene_table(kg, support):
    gene_nodes = set(kg.loc[kg.x_type == "gene/protein", "x_name"]) | set(
        kg.loc[kg.y_type == "gene/protein", "y_name"])
    out = []
    for group, genes in [("Bellenguez 2022 known", BELLENGUEZ_KNOWN),
                         ("Bellenguez 2022 new", BELLENGUEZ_NEW),
                         ("Mendelian / APOE", MENDELIAN)]:
        for g in genes:
            src = support.get(g, set())
            out.append({"group": group, "gene": g, "in_graph": g in gene_nodes,
                        "linked_to_AD": bool(src), "n_sources": len(src),
                        "sources": ", ".join(sorted(src)),
                        "linked_without_GWAS": bool(src - {"GWAS_Catalog"})})
    return pd.DataFrame(out)


def drug_table(kg, ad_set):
    drug_types = {"drug/chemical/compound", "clinical_intervention"}
    out = []
    for drug, cls in APPROVED_DRUGS.items():
        is_x = kg.x_type.isin(drug_types) & kg.x_name.map(norm).str.contains(rf"\b{drug}\b")
        is_y = kg.y_type.isin(drug_types) & kg.y_name.map(norm).str.contains(rf"\b{drug}\b")
        exact_nodes = set(kg.loc[is_x & (kg.x_name.map(norm) == drug), "x_id"]) | set(
            kg.loc[is_y & (kg.y_name.map(norm) == drug), "y_id"])
        to_ad = kg[(is_x & kg.y_id.isin(ad_set)) | (is_y & kg.x_id.isin(ad_set))]
        # Two-hop: drug - gene - AD (through PharmGKB gene-drug edges).
        genes = set(kg.loc[is_y & (kg.x_type == "gene/protein"), "x_name"])
        out.append({"drug": drug, "class": cls, "canonical_node": len(exact_nodes) > 0,
                    "matching_nodes": int(len(set(kg.loc[is_x, "x_id"]) | set(kg.loc[is_y, "y_id"])))
                    , "direct_AD_edges": len(to_ad),
                    "direct_AD_sources": ", ".join(sorted(set(to_ad.edge_source))),
                    "linked_genes": ", ".join(sorted(genes))})
    return pd.DataFrame(out)


def variant_table(kg, ad_set):
    """Is each locus's lead variant a node linked to the AD set, and which gene does AlzKG map it to?"""
    ad_var = kg[(kg.x_type == "variant") & kg.y_id.isin(ad_set)]
    gwas_gene = kg[(kg.relation == "gene_disease") & (kg.edge_source == "GWAS_Catalog")]
    out = []
    for locus, rsids in LEAD_VARIANTS.items():
        hit = ad_var[ad_var.x_id.isin(rsids)]
        out.append({"locus": locus, "group": "new" if locus in BELLENGUEZ_NEW or locus == "IGH gene cluster" else "known",
                    "lead_variants": ";".join(rsids), "any_lead_linked_to_AD": len(hit) > 0,
                    "AD_nodes": "; ".join(sorted(set(hit.y_name))),
                    "sources": ", ".join(sorted(set(hit.edge_source))),
                    "gene_named_in_graph": locus in set(gwas_gene.x_name)})
    return pd.DataFrame(out)


def permutation_p(n, is_ref, reps=100_000, seed=0):
    rng = np.random.default_rng(seed)
    v, m = n.values, int(is_ref.sum())
    obs = v[is_ref].mean() - v[~is_ref].mean()
    total = v.sum()
    hits = 0
    for _ in range(reps):
        s = v[rng.choice(len(v), m, replace=False)].sum()
        hits += (s / m - (total - s) / (len(v) - m)) >= obs
    return (hits + 1) / (reps + 1)


def convergence_auc(support, reference):
    """ROC AUC of n_sources for reference genes vs all other AD-linked genes."""
    n = support.map(len)
    is_ref = n.index.isin(reference)
    pos, neg = n[is_ref].values, n[~is_ref].values
    u, p = mannwhitneyu(pos, neg, alternative="greater")
    return u / (len(pos) * len(neg)), permutation_p(n, is_ref), len(pos), len(neg), n


if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    kg, ad_set = load(KG_PATH)
    print(f"{KG_PATH}: {len(kg):,} edges; AD set = anchor + {len(ad_set) - 1} subtype/trait nodes")

    support = gene_support(kg, ad_set)
    genes = gene_table(kg, support)
    genes.to_csv(OUT_DIR / "known_gene_recovery.csv", index=False)
    summary = genes.groupby("group").agg(n=("gene", "size"), in_graph=("in_graph", "sum"),
                                         linked_to_AD=("linked_to_AD", "sum"),
                                         linked_without_GWAS=("linked_without_GWAS", "sum"),
                                         median_sources=("n_sources", "median"))
    print("\nGene recovery\n" + summary.to_string())
    print("\nReference genes not linked to AD:",
          ", ".join(genes.loc[~genes.linked_to_AD, "gene"]) or "none")
    print("Reference genes with >=2 supporting datasets:")
    print(genes[genes.n_sources >= 2].drop_duplicates("gene")[["gene", "n_sources", "sources"]].to_string(index=False))

    reference = set(BELLENGUEZ_KNOWN) | set(BELLENGUEZ_NEW) | set(MENDELIAN)
    auc, p, npos, nneg, n = convergence_auc(support, reference)
    print(f"\nEvidence convergence: {npos} reference vs {nneg} other AD-linked genes; "
          f"AUC = {auc:.3f}, one-sided permutation p = {p:.1e} (100,000 label permutations)")
    dist = pd.DataFrame({"reference": n[n.index.isin(reference)].value_counts(),
                         "other": n[~n.index.isin(reference)].value_counts()}).fillna(0).astype(int)
    dist.index.name = "n_sources"
    print(dist.sort_index().to_string())
    multi = n[n >= 2]
    print(f"Genes with >=2 sources: {len(multi)}, of which {multi.index.isin(reference).sum()} "
          f"are reference genes (precision {multi.index.isin(reference).mean():.2f})")

    variants = variant_table(kg, ad_set)
    variants.to_csv(OUT_DIR / "lead_variant_recovery.csv", index=False)
    vs = variants.groupby("group").agg(loci=("locus", "size"), lead_linked=("any_lead_linked_to_AD", "sum"),
                                       gene_named=("gene_named_in_graph", "sum"))
    print("\nLead-variant recovery (Bellenguez 2022 loci)\n" + vs.to_string())
    print("Loci whose lead variant is absent:",
          ", ".join(variants.loc[~variants.any_lead_linked_to_AD, "locus"]) or "none")
    print("Loci recovered by variant but whose reported gene is not a GWAS gene node:",
          ", ".join(variants.loc[variants.any_lead_linked_to_AD & ~variants.gene_named_in_graph, "locus"]) or "none")

    drugs = drug_table(kg, ad_set)
    drugs.to_csv(OUT_DIR / "approved_drug_recovery.csv", index=False)
    print("\nApproved-drug recovery\n" + drugs.to_string(index=False))
    pd.DataFrame([{"auc": auc, "p": p, "n_reference": npos, "n_other": nneg}]).to_csv(
        OUT_DIR / "convergence_auc.csv", index=False)
