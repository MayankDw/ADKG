"""
ClinicalTrials.gov -> Alzheimer's disease Drug/Intervention-Disease-Trial edges.

Source: ClinicalTrials.gov API v2 (no auth required, JSON, paginated).
Docs: https://clinicaltrials.gov/data-api/api
Endpoint used: https://clinicaltrials.gov/api/v2/studies

This is a live network call (unlike the other processing_scripts, which parse a
file you downloaded via ../primary_data_resources.sh) because the API already
does the AD-specific filtering server-side via query.cond -- there is no bulk
file to curl first. This also doubles as the "dynamic ingestion" entry point
described in AD_Knowledge_Graphs_Bibliography.md Part II: re-running this script
picks up newly registered/updated AD trials without rebuilding anything else.

Equivalent bulk alternative: download the AACT Postgres snapshot
(https://aact.ctti-clinicaltrials.org/download) the same way PrimeKG dumps
DrugCentral's Postgres database (see primary_data_resources.sh), then query the
`conditions`, `interventions`, and `studies` tables directly.

Emits:
  - clinicaltrials_drug_disease.csv : Intervention (drug/biological/device/...) -> Disease, via a Clinical Trial
  - clinicaltrials_trials.csv       : one row per trial (NCT ID, phase, status, title)

Placebo/sham arms and "healthy volunteer" pseudo-conditions are dropped so they
don't become hub nodes in the graph.
"""
import html
import re
import time
from pathlib import Path

import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

API_URL = "https://clinicaltrials.gov/api/v2/studies"
DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "aact"
OUT_DRUG_DISEASE = DATA_DIR / "clinicaltrials_drug_disease.csv"
OUT_TRIALS = DATA_DIR / "clinicaltrials_trials.csv"

PLACEBO_PATTERN = re.compile(r"placebo|sham|vehicle|no intervention|usual care|standard of care",
                             flags=re.IGNORECASE)
HEALTHY_PATTERN = re.compile(r"^healthy|healthy (volunteer|subject|control|adult|older)", flags=re.IGNORECASE)

QUERY_COND = "Alzheimer Disease"
PAGE_SIZE = 200
FIELDS = [
    "NCTId", "BriefTitle", "OverallStatus", "Phase",
    "Condition", "InterventionName", "InterventionType",
]


def _session():
    session = requests.Session()
    retry = Retry(total=5, backoff_factor=1, status_forcelist=(429, 500, 502, 503, 504))
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def fetch_studies(query_cond=QUERY_COND, page_size=PAGE_SIZE, max_pages=None):
    studies = []
    params = {
        "query.cond": query_cond,
        "fields": ",".join(FIELDS),
        "pageSize": page_size,
        "countTotal": "true",
    }
    session = _session()
    page = 0
    while True:
        resp = session.get(API_URL, params=params, timeout=60)
        resp.raise_for_status()
        payload = resp.json()
        studies.extend(payload.get("studies", []))
        if page == 0 and "totalCount" in payload:
            print(f"  ClinicalTrials.gov reports {payload['totalCount']} matching studies")
        page += 1
        print(f"  page {page}: {len(studies)} studies fetched", end="\r")
        next_token = payload.get("nextPageToken")
        if not next_token or (max_pages and page >= max_pages):
            break
        params["pageToken"] = next_token
        time.sleep(0.2)  # be polite to the public API
    return studies


def flatten(studies):
    rows = []
    for s in studies:
        protocol = s.get("protocolSection", {})
        ident = protocol.get("identificationModule", {})
        status = protocol.get("statusModule", {})
        design = protocol.get("designModule", {})
        conditions_mod = protocol.get("conditionsModule", {})
        arms = protocol.get("armsInterventionsModule", {})

        # Registry free text sometimes carries HTML entities, e.g. "Alzheimer&#39;s Disease".
        nct_id = ident.get("nctId")
        title = html.unescape(ident.get("briefTitle") or "")
        overall_status = status.get("overallStatus")
        phases = design.get("phases", [])
        conditions = [html.unescape(c) for c in conditions_mod.get("conditions", [])]
        interventions = [{**i, "name": html.unescape(i["name"])} if i.get("name") else i
                         for i in arms.get("interventions", [])]

        rows.append({
            "nct_id": nct_id,
            "title": title,
            "status": overall_status,
            "phase": ";".join(phases) if phases else None,
            "conditions": conditions,
            "interventions": interventions,
        })
    return rows


def build_trial_table(rows):
    return pd.DataFrame([{
        "nct_id": r["nct_id"],
        "title": r["title"],
        "status": r["status"],
        "phase": r["phase"],
    } for r in rows]).drop_duplicates()


def build_drug_disease_edges(rows):
    edges = []
    for r in rows:
        nct_id = r["nct_id"]
        conditions = [c for c in r["conditions"] or [] if not HEALTHY_PATTERN.search(c)]
        interventions = [i for i in r["interventions"] or []
                         if i.get("name") and not PLACEBO_PATTERN.search(i["name"])]
        for cond in conditions:
            for interv in interventions:
                edges.append({
                    "relation": "drug_disease",
                    "display_relation": "studied in clinical trial",
                    "x_id": interv.get("name"),
                    "x_type": "drug/chemical/compound" if interv.get("type") in
                              ("DRUG", "BIOLOGICAL") else "clinical_intervention",
                    "x_name": interv.get("name"),
                    "x_source": "ClinicalTrials.gov",
                    "y_id": cond,
                    "y_type": "disease",
                    "y_name": cond,
                    "y_source": "ClinicalTrials.gov",
                    "trial_id": nct_id,
                })
    return pd.DataFrame(edges).drop_duplicates()


if __name__ == "__main__":
    studies = fetch_studies()
    print()
    print(f"Fetched {len(studies)} ClinicalTrials.gov studies for condition='{QUERY_COND}'.")

    rows = flatten(studies)
    trials = build_trial_table(rows)
    edges = build_drug_disease_edges(rows)

    trials.to_csv(OUT_TRIALS, index=False)
    edges.to_csv(OUT_DRUG_DISEASE, index=False)
    print(f"Wrote {len(trials)} trials -> {OUT_TRIALS}")
    print(f"Wrote {len(edges)} drug/intervention-disease edges -> {OUT_DRUG_DISEASE}")
