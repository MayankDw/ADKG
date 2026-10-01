"""
Render an interactive, self-contained HTML explorer for the AlzKG core subgraph.

Reads ../datasets/data/kg/ad_kg.csv + ad_kg_nodes.csv (produced by build_ad_kg.ipynb)
and the per-source edge files, picks the most informative ~150 nodes around
Alzheimer's disease (multi-source genes, key variants, AD subtypes, AD drugs and their
pharmacogenomic genes), and writes ../datasets/data/kg/ad_kg_explorer.html.

    python knowledge_graph/visualize_ad_kg.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "datasets" / "data"
KG_DIR = DATA / "kg"
OUT = KG_DIR / "ad_kg_explorer.html"
AD_ID = "MONDO:0004975"

SOURCE_FILES = {
    "GWAS Catalog": ["gwas_catalog/gwas_variant_disease.csv", "gwas_catalog/gwas_gene_disease.csv"],
    "ClinVar": ["clinvar/clinvar_ad_variant_disease.csv"],
    "ClinicalTrials.gov": ["aact/clinicaltrials_drug_disease.csv"],
    "PharmGKB": ["pharmgkb/pharmgkb_ad_edges.csv"],
    "PubTator3": ["literature/literature_triples.csv"],
    "DisGeNET": ["disgenet/disgenet_ad_gene_disease.csv"],
    "NIAGADS": ["niagads/niagads_variant_disease.csv", "niagads/niagads_gene_disease.csv"],
}
KEY_GENES = ["APOE", "APP", "PSEN1", "PSEN2", "TREM2", "CLU", "BIN1", "ABCA7", "SORL1", "CR1",
             "PICALM", "CD33", "MS4A6A", "EPHA1", "CD2AP", "MAPT", "BCHE", "ACHE", "CHRNA7", "CYP2D6"]
KEY_DRUGS = ["donepezil", "rivastigmine", "galantamine", "memantine", "lecanemab",
             "aducanumab", "donanemab", "tacrine"]
KEY_VARIANTS = ["rs429358", "rs7412"]


def neighbors_of(kg, idx):
    return pd.concat([kg.loc[kg.x_index == idx, "y_index"], kg.loc[kg.y_index == idx, "x_index"]])


def select_nodes(kg, nodes):
    by_id = nodes.set_index("node_id")
    ad_idx = int(by_id.loc[AD_ID, "node_index"])
    keep = {ad_idx}
    ad_edges = kg[(kg.x_index == ad_idx) | (kg.y_index == ad_idx)].copy()
    ad_edges["nbr"] = ad_edges.x_index.where(ad_edges.y_index == ad_idx, ad_edges.y_index)
    ad_edges = ad_edges.merge(nodes[["node_index", "node_type", "node_name", "node_source"]],
                              left_on="nbr", right_on="node_index")

    # Genes: named AD genes + genes backed by the most independent datasets.
    genes = ad_edges[ad_edges.node_type == "gene/protein"]
    support = genes.groupby("nbr")["edge_source"].nunique().sort_values(ascending=False)
    keep |= set(support.head(18).index)
    named = nodes[(nodes.node_type == "gene/protein") & nodes.node_name.isin(KEY_GENES)
                  & (nodes.node_source == "NCBIGene")]
    keep |= set(named.node_index)

    # AD subtypes with ontology ids (MONDO), not free-text trial conditions.
    subtypes = ad_edges[(ad_edges.relation == "disease_disease")
                        & ad_edges.node_source.isin(["MONDO", "ClinVar", "GWAS_Catalog"])]
    keep |= set(subtypes[subtypes.node_name.str.len() < 60].nbr.head(16))

    # Variants: APOE e4/e2 SNPs + pathogenic ClinVar variants in the familial AD genes.
    keep |= set(nodes[nodes.node_id.isin(KEY_VARIANTS)].node_index)
    path = kg[(kg.edge_source == "ClinVar") & (kg.display_relation == "pathogenic")]
    for gene in ["APP", "PSEN1", "PSEN2"]:
        gv = path[path.x_name.str.contains(f"({gene})", regex=False)].drop_duplicates("x_index")
        keep |= set(gv.x_index.head(3)) | set(gv.y_index.head(3))
    # Variants supported by several datasets.
    vsup = ad_edges[ad_edges.node_type == "variant"].groupby("nbr")["edge_source"].nunique()
    keep |= set(vsup[vsup >= 3].index[:6])

    # Drugs: approved/late-stage AD drugs + drugs most often trialled for AD.
    drugs = ad_edges[ad_edges.node_type == "drug/chemical/compound"]
    keep |= set(drugs[drugs.node_name.str.lower().isin(KEY_DRUGS)].nbr)
    ct = pd.read_csv(DATA / "aact/clinicaltrials_drug_disease.csv", usecols=["x_name", "trial_id"])
    top_trial = ct.groupby(ct.x_name.str.lower())["trial_id"].nunique().sort_values(ascending=False)
    lower = nodes.assign(l=nodes.node_name.str.lower())
    trial_drugs = lower[(lower.node_type == "drug/chemical/compound") & lower.l.isin(top_trial.head(25).index)]
    keep |= set(trial_drugs.node_index.head(14))

    # Literature: strongest PubTator chemical/gene co-relations with AD.
    lit = ad_edges[ad_edges.edge_source == "PubTator3"].groupby("nbr").size().sort_values(ascending=False)
    keep |= set(lit.head(8).index)
    return keep, ad_idx


def build_payload(kg, nodes):
    keep, ad_idx = select_nodes(kg, nodes)
    sub = kg[kg.x_index.isin(keep) & kg.y_index.isin(keep)]
    # Collapse parallel edges into one link listing every supporting dataset.
    links = (sub.groupby(["x_index", "y_index"])
             .agg(relation=("relation", "first"),
                  labels=("display_relation", lambda s: sorted(set(map(str, s)))[:4]),
                  sources=("edge_source", lambda s: sorted(set(s))))
             .reset_index())
    used = set(links.x_index) | set(links.y_index)
    deg_full = pd.concat([kg.x_index, kg.y_index]).value_counts()
    n = nodes[nodes.node_index.isin(used)]
    node_list = [{
        "id": int(r.node_index), "label": r.node_name if len(str(r.node_name)) < 34 else str(r.node_name)[:32] + "…",
        "name": r.node_name, "cid": r.node_id, "type": r.node_type, "vocab": r.node_source,
        "degree": int(deg_full.get(r.node_index, 0)),
    } for r in n.itertuples()]
    link_list = [{"from": int(r.x_index), "to": int(r.y_index), "relation": r.relation,
                  "labels": r.labels, "sources": r.sources} for r in links.itertuples()]

    per_source = {}
    for src, files in SOURCE_FILES.items():
        total = 0
        for f in files:
            p = DATA / f
            if p.exists():
                total += sum(1 for _ in open(p, encoding="utf-8")) - 1
        per_source[src] = total
    summary = {
        "per_source": per_source,
        "raw_edges": sum(per_source.values()),
        "final_edges": int(len(kg)),
        "final_nodes": int(len(nodes)),
        "node_types": nodes.node_type.value_counts().to_dict(),
        "edge_sources": kg.edge_source.value_counts().to_dict(),
        "ad_degree": int(((kg.x_index == ad_idx) | (kg.y_index == ad_idx)).sum()),
        "subgraph_nodes": len(node_list), "subgraph_links": len(link_list),
    }
    return {"nodes": node_list, "links": link_list, "summary": summary, "ad": ad_idx}


def main():
    kg = pd.read_csv(KG_DIR / "ad_kg.csv", low_memory=False)
    nodes = pd.read_csv(KG_DIR / "ad_kg_nodes.csv", low_memory=False)
    payload = build_payload(kg, nodes)
    template = (Path(__file__).parent / "ad_kg_explorer_template.html").read_text(encoding="utf-8")
    html = template.replace("/*__DATA__*/null", json.dumps(payload, ensure_ascii=False))
    head = ('<!doctype html>\n<html lang="en">\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1">\n')
    OUT.write_text(head + html, encoding="utf-8")
    s = payload["summary"]
    print(f"Explorer subgraph: {s['subgraph_nodes']} nodes, {s['subgraph_links']} links "
          f"(full graph {s['final_nodes']:,} nodes / {s['final_edges']:,} edges) -> {OUT}")


if __name__ == "__main__":
    main()
