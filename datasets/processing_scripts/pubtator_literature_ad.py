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
  1. Searches PubTator3 for new/updated PMIDs matching the AD query.
  2. Pulls BioC-JSON annotations + relations for those PMIDs.
  3. Appends newly found triples to literature_triples.csv.
  4. Advances a checkpoint file so the next run only fetches what changed.

This keeps the literature-mined slice of the graph continuously current instead
of rebuilding a one-shot static snapshot, which is the central limitation the
bibliography's Part II literature calls out in AlzKB/PrimeKG/ADKG.

Entity/relation types follow PubTator3's BioREx vocabulary (Gene, Disease,
Chemical, Variant, Species, CellLine) and are mapped onto the harmonized KG
node/edge schema in ../../vocab/Node_and_Edge_Types.csv.
"""
import json
import time
from pathlib import Path

import requests
import pandas as pd

SEARCH_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/"
EXPORT_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson"

QUERY = "alzheimer disease"
BATCH_SIZE = 100  # PMIDs per export call
OUT_TRIPLES = Path("../data/literature/literature_triples.csv")
CHECKPOINT = Path("../data/literature/.checkpoint_last_pmid_page.json")

PUBTATOR_TYPE_MAP = {
    "Gene": "gene/protein",
    "Disease": "disease",
    "Chemical": "drug/chemical/compound",
    "Variant": "variant",
    "Species": "species/taxon",
    "CellLine": "cell/cell_line",
}


def load_checkpoint():
    if CHECKPOINT.exists():
        return json.loads(CHECKPOINT.read_text())
    return {"last_page": 0}


def save_checkpoint(state):
    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    CHECKPOINT.write_text(json.dumps(state))


def search_pmids(query=QUERY, start_page=1, max_pages=5):
    """Page through PubTator3 search results starting after the last checkpointed page."""
    pmids = []
    page = start_page
    for _ in range(max_pages):
        resp = requests.get(SEARCH_URL, params={"text": query, "page": page}, timeout=30)
        resp.raise_for_status()
        payload = resp.json()
        results = payload.get("results", [])
        if not results:
            break
        pmids.extend(str(r["pmid"]) for r in results if r.get("pmid"))
        page += 1
        time.sleep(0.34)  # NCBI E-utilities-style rate limit: ~3 req/sec
    return pmids, page


def fetch_biocjson(pmids, batch_size=BATCH_SIZE):
    docs = []
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]
        resp = requests.get(EXPORT_URL, params={"pmids": ",".join(batch)}, timeout=60)
        resp.raise_for_status()
        text = resp.text.strip()
        for line in text.splitlines():
            if line.strip():
                docs.append(json.loads(line))
        time.sleep(0.34)
    return docs


def extract_triples(docs):
    rows = []
    for doc in docs:
        pmid = doc.get("pmid") or doc.get("id")
        passages = doc.get("passages", [])
        annotations_by_id = {}
        for passage in passages:
            for ann in passage.get("annotations", []):
                infons = ann.get("infons", {})
                ann_id = infons.get("identifier") or ann.get("id")
                annotations_by_id[ann.get("id")] = {
                    "id": ann_id,
                    "type": infons.get("type"),
                    "text": ann.get("text"),
                }
            for relation in passage.get("relations", []):
                infons = relation.get("infons", {})
                node1 = annotations_by_id.get(infons.get("entity1"))
                node2 = annotations_by_id.get(infons.get("entity2"))
                if not node1 or not node2:
                    continue
                rows.append({
                    "relation": "literature_co_mention",
                    "display_relation": infons.get("type", "associated_with"),
                    "x_id": node1["id"],
                    "x_type": PUBTATOR_TYPE_MAP.get(node1["type"], "generic_entity"),
                    "x_name": node1["text"],
                    "x_source": "PubTator3",
                    "y_id": node2["id"],
                    "y_type": PUBTATOR_TYPE_MAP.get(node2["type"], "generic_entity"),
                    "y_name": node2["text"],
                    "y_source": "PubTator3",
                    "pmid": pmid,
                })
    return pd.DataFrame(rows).drop_duplicates()


def append_triples(new_triples, path=OUT_TRIPLES):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        existing = pd.read_csv(path)
        combined = pd.concat([existing, new_triples], ignore_index=True).drop_duplicates()
    else:
        combined = new_triples
    combined.to_csv(path, index=False)
    return combined


if __name__ == "__main__":
    state = load_checkpoint()
    pmids, next_page = search_pmids(start_page=state["last_page"] + 1)
    print(f"Found {len(pmids)} candidate PMIDs (pages {state['last_page'] + 1}..{next_page - 1}).")

    if pmids:
        docs = fetch_biocjson(pmids)
        triples = extract_triples(docs)
        combined = append_triples(triples)
        print(f"Extracted {len(triples)} new triples this run; {len(combined)} total in {OUT_TRIPLES}.")
    else:
        print("No new PMIDs since last checkpoint; nothing to do.")

    save_checkpoint({"last_page": next_page - 1})
