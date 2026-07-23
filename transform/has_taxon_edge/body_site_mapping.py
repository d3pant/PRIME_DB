from pathlib import Path

import pandas as pd
import requests
from node_synonymizer import NodeSynonymizer

SAMPLES_PARQUET = Path(__file__).parent.parent.parent / "data/extract/samples.parquet"
OUT_PATH = Path(__file__).parent / "data/body_site_mapping.parquet"

NAMERES_BASE = "https://name-resolution-sri.renci.org"

df = pd.read_parquet(SAMPLES_PARQUET)

# Body_Site is 0% null. Values like "Stool/feces" put the specific term on the
# right of the slash -- take that side, not the left (confirmed earlier: "Stool"
# is the generic/redundant half, "feces" is the specific one).
unique_body_sites = set()
for site in df["Body_Site"].dropna():
    unique_body_sites.add(site.split("/")[-1].strip().lower())

body_site_terms = sorted(unique_body_sites)
print(f"unique body site terms: {len(body_site_terms)}")

# Stage 1: NodeSynonymizer, accept only UBERON hits.
syn = NodeSynonymizer()
synonymizer_results = syn.get_canonical_curies(names=body_site_terms)

accepted = {}
needs_fallback = []
for term in body_site_terms:
    result = synonymizer_results.get(term)
    if result and result["preferred_curie"].startswith("UBERON:"):
        accepted[term] = {"id": result["preferred_curie"], "name": result["preferred_name"]}
    else:
        needs_fallback.append(term)

print(f"accepted from NodeSynonymizer directly (UBERON): {len(accepted)}")
print(f"flagged for NameRes fallback: {len(needs_fallback)}")

# Stage 2: NameRes fallback, only accepting a candidate whose CURIE prefix is
# UBERON (checked against ranked candidates, not just top-1) -- strict, no other
# prefix accepted, matching the body-site pipeline built earlier in samples.ipynb.
resp = requests.post(
    f"{NAMERES_BASE}/bulk-lookup",
    json={"strings": needs_fallback, "autocomplete": False, "limit": 10},
)
resp.raise_for_status()
nameres_raw = resp.json()

fallback_accepted = {}
unresolved = []
for term in needs_fallback:
    candidates = nameres_raw.get(term) or []
    chosen = next((c for c in candidates if c["curie"].startswith("UBERON:")), None)
    if chosen:
        fallback_accepted[term] = {"id": chosen["curie"], "name": chosen["label"]}
    else:
        unresolved.append(term)

print(f"NameRes fallback resolved (UBERON): {len(fallback_accepted)}")
print(f"still unresolved: {len(unresolved)}")

rows = []
for term in body_site_terms:
    if term in accepted:
        r = accepted[term]
        rows.append({"original_name": term, "id": r["id"], "name": r["name"]})
    elif term in fallback_accepted:
        r = fallback_accepted[term]
        rows.append({"original_name": term, "id": r["id"], "name": r["name"]})
    else:
        rows.append({"original_name": term, "id": None, "name": None})

mapping_df = pd.DataFrame(rows)
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
mapping_df.to_parquet(OUT_PATH, index=False)

print(f"\nwrote {len(mapping_df)} body site mappings to {OUT_PATH}")
print(f"resolved: {mapping_df['id'].notna().sum()} / {len(mapping_df)}")
