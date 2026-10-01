"""
PubTator3 -> Alzheimer's disease literature-derived triples (the "dynamic
ingestion" layer).

Source: NCBI PubTator3 API (no auth required).
Docs: https://www.ncbi.nlm.nih.gov/research/pubtator3-api/

Unlike the other processing_scripts, this one is designed to be re-run
periodically (cron/scheduled job) rather than once against a static bulk
download -- see AD_Knowledge_Graphs_Bibliography.md Part 11 ("Dynamic, Temporal
& Incremental Knowledge Graph Methods") and Part 14 (#126 PubTator 3.0, #128
LLM-Empowered KG Construction). Each run:
  1. Searches PubTator3 for PMIDs that have at least one machine-extracted
     relation involving Alzheimer's disease (relation query, see QUERY).
  2. Pulls BioC-JSON for those PMIDs (export returns {"PubTator3": [doc, ...]}).
  3. Reads each document's `relations` list (document-level, with the two
     normalized entities in infons.role1 / infons.role2) into triples, and
     appends them to literature_triples.csv.
  4. Advances a checkpoint file so the next run continues from the next page.

This keeps the literature-mined slice of the graph continuously current instead
of rebuilding a one-shot static snapshot, which is the central limitation the
bibliography's Part II literature calls out in AlzKB/PrimeKG/ADKG.

Entity/relation types follow PubTator3's BioREx vocabulary (Gene, Disease,
Chemical, Variant, Species, CellLine; Association, Positive_Correlation,
Negative_Correlation, Bind, Cotreatment, Comparison, Drug_Interaction,
Conversion) and are mapped onto the harmonized KG node/edge schema in
../../vocab/Node_and_Edge_Types.csv.

Usage:
    python pubtator_literature_ad.py               # next MAX_PAGES pages (10 PMIDs/page)
    python pubtator_literature_ad.py --pages 50    # fetch more this run
    python pubtator_literature_ad.py --reset       # forget checkpoint, start from page 1
"""
import argparse
import json
import time
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

SEARCH_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/"
EXPORT_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson"

# PubTator3 relation search: any relation whose one side is the Alzheimer Disease
# concept. A plain text query ("alzheimer disease") mostly returns papers with
# no extracted relations at all.
QUERY = "relations:ANY|@DISEASE_Alzheimer_Disease|ANY"
MAX_PAGES = 30   # search pages per run (PubTator3 returns 10 PMIDs per page)
BATCH_SIZE = 100  # PMIDs per export call

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "literature"
OUT_TRIPLES = DATA_DIR / "literature_triples.csv"
CHECKPOINT = DATA_DIR / ".checkpoint_last_pmid_page.json"

PUBTATOR_TYPE_MAP = {
    "Gene": "gene/protein",
    "Disease": "disease",
    "Chemical": "drug/chemical/compound",
    "Variant": "variant",
    "Species": "species/taxon",
    "CellLine": "cell/cell_line",
}

SESSION = requests.Session()
SESSION.mount("https://", HTTPAdapter(max_retries=Retry(
    total=5, backoff_factor=2, status_forcelist=(429, 500, 502, 503, 504))))


def load_checkpoint():
    if CHECKPOINT.exists():
        state = json.loads(CHECKPOINT.read_text())
        if state.get("query") == QUERY:
            return state
    return {"query": QUERY, "last_page": 0}


def save_checkpoint(state):
    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    CHECKPOINT.write_text(json.dumps(state))


def search_pmids(query=QUERY, start_page=1, max_pages=MAX_PAGES):
    """Page through PubTator3 search results starting after the last checkpointed page."""
    pmids = []
    page = start_page
    for _ in range(max_pages):
        resp = SESSION.get(SEARCH_URL, params={"text": query, "page": page}, timeout=60)
        resp.raise_for_status()
        payload = resp.json()
        if page == start_page:
            print(f"  PubTator3 reports {payload.get('count')} matching articles")
        results = payload.get("results", [])
        if not results:
            break
        pmids.extend(str(r["pmid"]) for r in results if r.get("pmid"))
        page += 1
        time.sleep(0.34)  # NCBI rate limit: ~3 req/sec
    return pmids, page


def fetch_biocjson(pmids, batch_size=BATCH_SIZE):
    docs = []
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]
        resp = SESSION.get(EXPORT_URL, params={"pmids": ",".join(batch)}, timeout=120)
        resp.raise_for_status()
        payload = resp.json()
        docs.extend(payload.get("PubTator3", []) if isinstance(payload, dict) else payload)
        print(f"  fetched annotations for {len(docs)}/{len(pmids)} articles")
        time.sleep(0.34)
    return docs


def extract_triples(docs):
    rows = []
    for doc in docs:
        pmid = doc.get("pmid") or doc.get("id")
        for relation in doc.get("relations", []):
            infons = relation.get("infons", {})
            node1, node2 = infons.get("role1") or {}, infons.get("role2") or {}
            if not node1.get("identifier") or not node2.get("identifier"):
                continue
            rel_type = infons.get("type", "Association")
            rows.append({
                "relation": "literature_" + rel_type.lower(),
                "display_relation": rel_type.replace("_", " ").lower(),
                "x_id": node1["identifier"],
                "x_type": PUBTATOR_TYPE_MAP.get(node1.get("type"), "generic_entity"),
                "x_name": node1.get("name") or node1["identifier"],
                "x_source": "PubTator3",
                "y_id": node2["identifier"],
                "y_type": PUBTATOR_TYPE_MAP.get(node2.get("type"), "generic_entity"),
                "y_name": node2.get("name") or node2["identifier"],
                "y_source": "PubTator3",
                "pmid": pmid,
                "score": infons.get("score"),
            })
    return pd.DataFrame(rows).drop_duplicates()


def append_triples(new_triples, path=OUT_TRIPLES):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        existing = pd.read_csv(path, dtype=str)
        combined = pd.concat([existing, new_triples.astype(str)], ignore_index=True).drop_duplicates()
    else:
        combined = new_triples
    combined.to_csv(path, index=False)
    return combined


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pages", type=int, default=MAX_PAGES)
    parser.add_argument("--reset", action="store_true")
    args = parser.parse_args()

    state = {"query": QUERY, "last_page": 0} if args.reset else load_checkpoint()
    pmids, next_page = search_pmids(start_page=state["last_page"] + 1, max_pages=args.pages)
    print(f"Found {len(pmids)} candidate PMIDs (pages {state['last_page'] + 1}..{next_page - 1}).")

    if pmids:
        docs = fetch_biocjson(pmids)
        triples = extract_triples(docs)
        combined = append_triples(triples)
        print(f"Extracted {len(triples)} new triples this run; {len(combined)} total in {OUT_TRIPLES}.")
    else:
        print("No new PMIDs since last checkpoint; nothing to do.")

    save_checkpoint({"query": QUERY, "last_page": next_page - 1})
