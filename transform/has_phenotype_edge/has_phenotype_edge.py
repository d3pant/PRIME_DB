import re
from pathlib import Path

import pandas as pd

SAMPLES_PARQUET = Path(__file__).parent.parent.parent / "data/extract/samples.parquet"
PHENOTYPE_MAPPING_PARQUET = Path(__file__).parent.parent.parent / "data/extract/phenotype_mapping.parquet"
OUT_PATH = Path(__file__).parent.parent.parent / "data/transformed/edges/has_phenotype.parquet"

samples_df = pd.read_parquet(SAMPLES_PARQUET)
mapping_df = pd.read_parquet(PHENOTYPE_MAPPING_PARQUET)
term_to_curie = mapping_df.set_index("original_name")[["final_curie", "final_curie_name"]].to_dict("index")


def normalize_phenotype_cell(raw: str) -> list[str]:
    cleaned = re.sub(r"\([^)]*\)", "", raw)
    cleaned = cleaned.replace("'", "").replace("’", "").replace("‘", "")
    return [t.strip().lower() for t in cleaned.split(",") if t.strip()]


# BioSample is the edge subject (deliberately BioSample-only, not per-Run) --
# BioProject and Doi are confirmed identical across every run of a given BioSample
# (0 conflicts verified), safe to take the first row's value. Phenotype itself is
# NOT always consistent across runs of the same BioSample (3 real conflicts found,
# e.g. one BioSample's runs disagreeing between "kidney transplant donor" vs
# "recipient") -- those BioSamples are dropped entirely rather than guessing.
first_per_biosample = samples_df.groupby("BioSample", sort=False).first()

phenotype_sets_by_biosample = (
    samples_df.dropna(subset=["Phenotype"])
    .groupby("BioSample")["Phenotype"]
    .apply(lambda s: {t for raw in s for t in normalize_phenotype_cell(raw)})
)

conflicting = set()
for biosample, terms in phenotype_sets_by_biosample.items():
    per_run_sets = samples_df.loc[samples_df["BioSample"] == biosample, "Phenotype"].dropna()
    if per_run_sets.nunique() > 1:
        conflicting.add(biosample)

rows = []
for biosample, terms in phenotype_sets_by_biosample.items():
    if biosample in conflicting:
        continue

    info = first_per_biosample.loc[biosample]
    publication = None
    if pd.notna(info["Doi"]):
        publication = {"doi": f"doi:{info['Doi']}", "category": "biolink:Publication"}

    for term in terms:
        match = term_to_curie.get(term)
        if match is None:
            continue

        curie = match["final_curie"]
        numeric_id = curie.split(":", 1)[1]

        rows.append(
            {
                "id": f"{biosample}_has_phenotype_{numeric_id}",
                "association": {
                    "category": "biolink:MaterialSampleToDiseaseOrPhenotypicFeatureAssociation",
                    "predicate": "biolink:has_phenotype",
                    "BioProject": info["BioProject"],
                    "primary_knowledge_source": ["infores:PRIMEDB"],
                    "aggregator_knowledge_source": None,
                    "agent_type": "biolink:manual_agent",
                    "publication": publication,
                },
                "subject": biosample,
                "object": curie,
            }
        )

edge_df = pd.DataFrame(rows)
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
edge_df.to_parquet(OUT_PATH, index=False)

print(f"wrote {len(edge_df)} HAS_PHENOTYPE edges to {OUT_PATH}")
print(f"dropped {len(conflicting)} BioSamples with conflicting Phenotype values across runs")
print(f"unique BioSamples with >=1 has_phenotype edge: {edge_df['subject'].nunique()}")

if not edge_df["id"].is_unique:
    raise AssertionError("has_phenotype edge table has duplicate ids")
