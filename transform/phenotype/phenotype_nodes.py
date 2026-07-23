import re
from pathlib import Path

import pandas as pd
import requests
from node_synonymizer import NodeSynonymizer

SAMPLES_PARQUET = Path(__file__).parent.parent.parent / "data/extract/samples.parquet"
MAPPING_PARQUET = Path(__file__).parent.parent.parent / "data/extract/phenotype_mapping.parquet"
OUT_PATH = Path(__file__).parent.parent.parent / "data/transformed/nodes/phenotype.parquet"

ACCEPTED_PREFIXES = ("MONDO:", "HP:")
ACCEPTED_CATEGORIES = {"biolink:Disease", "biolink:PhenotypicFeature", "biolink:DiseaseOrPhenotypicFeature"}
NAMERES_BASE = "https://name-resolution-sri.renci.org"


def normalize_phenotype_cell(raw: str) -> list[str]:
    cleaned = re.sub(r"\([^)]*\)", "", raw)
    cleaned = cleaned.replace("'", "").replace("’", "").replace("‘", "")
    return [t.strip().lower() for t in cleaned.split(",") if t.strip()]


samples_df = pd.read_parquet(SAMPLES_PARQUET)
mapping_df = pd.read_parquet(MAPPING_PARQUET)

unique_phenotype = set()
for raw in samples_df["Phenotype"].dropna():
    unique_phenotype.update(normalize_phenotype_cell(raw))

# phenotype_mapping.parquet already reflects the retain_flag=False drops decided
# during manual curation -- a term missing from it was rejected, not merely unseen.
retained_terms = set(mapping_df["original_name"])
phenotype_terms = sorted(unique_phenotype & retained_terms)

print(f"unique phenotype terms (post-normalization): {len(unique_phenotype)}")
print(f"terms retained after curation drop:          {len(phenotype_terms)}")

# --- 1st API call: NodeSynonymizer, accept only MONDO/HP hits ---
syn = NodeSynonymizer()
synonymizer_results = syn.get_canonical_curies(names=phenotype_terms)

accepted = {}
needs_fallback = {}

for term in phenotype_terms:
    result = synonymizer_results.get(term)
    if not result:
        needs_fallback[term] = None
        continue
    curie = result["preferred_curie"]
    category = result.get("preferred_category")
    if category not in ACCEPTED_CATEGORIES or not curie.startswith(ACCEPTED_PREFIXES):
        needs_fallback[term] = result
        continue
    accepted[term] = {"preferred_curie": curie, "preferred_name": result["preferred_name"]}

print(f"accepted from NodeSynonymizer directly (MONDO/HP): {len(accepted)}")
print(f"flagged for NameRes fallback:                       {len(needs_fallback)}")

# --- 2nd API call: NameRes fallback, accept only MONDO/HP hits ---
fallback_terms = list(needs_fallback.keys())
resp = requests.post(
    f"{NAMERES_BASE}/bulk-lookup",
    json={"strings": fallback_terms, "autocomplete": False, "limit": 10},
)
resp.raise_for_status()
nameres_raw = resp.json()

fallback_accepted = {}
unresolved_terms = []

for term in fallback_terms:
    candidates = nameres_raw.get(term) or []
    chosen = next((c for c in candidates if c["curie"].startswith(ACCEPTED_PREFIXES)), None)
    if chosen:
        fallback_accepted[term] = {"preferred_curie": chosen["curie"], "preferred_name": chosen["label"]}
    else:
        unresolved_terms.append(term)

print(f"NameRes fallback resolved (MONDO/HP): {len(fallback_accepted)}")
print(f"still unresolved after both API calls: {len(unresolved_terms)}")

# --- Assemble final curie/name per term, then apply phenotype_mapping.parquet as
# an override layer -- both calls still run for every term; the mapping only
# replaces the result where curation recorded a different final_curie. ---
mapping_lookup = mapping_df.set_index("original_name")[["final_curie", "final_curie_name"]].to_dict("index")

rows = []
for term in phenotype_terms:
    if term in accepted:
        curie, name = accepted[term]["preferred_curie"], accepted[term]["preferred_name"]
    elif term in fallback_accepted:
        curie, name = fallback_accepted[term]["preferred_curie"], fallback_accepted[term]["preferred_name"]
    else:
        curie, name = None, None

    override = mapping_lookup.get(term)
    if override:
        curie, name = override["final_curie"], override["final_curie_name"]

    if curie is None:
        continue

    rows.append({"id": curie, "original_name": term, "category": "biolink:DiseaseOrPhenotypicFeature", "_final_curie_name": name})

phenotype_df = pd.DataFrame(rows)

# A node table must have exactly one row per distinct id -- multiple source terms
# can resolve (or be overridden) onto the same CURIE.
before_dedup = len(phenotype_df)
collapsed = {
    curie: group["original_name"].tolist()
    for curie, group in phenotype_df.groupby("id")
    if len(group) > 1
}
phenotype_df = phenotype_df.drop_duplicates(subset="id", keep="first").drop(columns="_final_curie_name")
after_dedup = len(phenotype_df)

phenotype_df.to_parquet(OUT_PATH, index=False)

print(f"\nresolved terms (pre-dedup): {before_dedup}")
print(f"wrote {after_dedup} Phenotype nodes to {OUT_PATH}")
if collapsed:
    print(f"\n{len(collapsed)} CURIE(s) absorbed multiple source terms:")
    for curie, terms in collapsed.items():
        kept, dropped = terms[0], terms[1:]
        print(f"  {curie}: kept {kept!r}, dropped {dropped}")

if not phenotype_df["id"].is_unique:
    raise AssertionError("phenotype node table has duplicate ids after dedup")
