from pathlib import Path

import pandas as pd

IN_PATH = Path(__file__).parent.parent.parent / "data/extract/samples.parquet"
OUT_PATH = Path(__file__).parent.parent.parent / "data/transformed/nodes/sample.parquet"

# Sample node identity: BioSample ID, 0% null in source data. One BioSample can
# span multiple SRA runs; demographics are identical across those rows for
# every BioSample (verified: 0 conflicts across 1,574 multi-run BioSamples) --
# safe to take the first row's value with no exception handling needed.
df = pd.read_parquet(IN_PATH)

TAKE_FIRST_COLS = [
    "Host_Age",
    "Host_Sex",
    "Host_Height",
    "Host_Weight",
    "Host_BMI",
    "Race_or_Ethnicity",
    "Country",
    "Continent",
    "Smoke_Status",
    "Drinking_Status",
    "Diet_Type",
]

# Rename to snake_case at the source so PRIMEDB's own Sample schema is
# internally consistent (previously only the 5 fields overlapping with
# GMREPO_V3's naming got renamed, and only at merge time -- the other 6
# PRIMEDB-only fields stayed Title_Case forever, a real inconsistency).
# host_age/sex/BMI/country/diet additionally match GMREPO_V3's own naming for
# the same concepts; BMI stays uppercase since it's an acronym, not a casing
# violation (matches GMREPO_V3's own "BMI" column too).
RENAME_TO_SNAKE_CASE = {
    "Host_Age": "host_age",
    "Host_Sex": "sex",
    "Host_Height": "host_height",
    "Host_Weight": "host_weight",
    "Host_BMI": "BMI",
    "Race_or_Ethnicity": "race_or_ethnicity",
    "Country": "country",
    "Continent": "continent",
    "Smoke_Status": "smoke_status",
    "Drinking_Status": "drinking_status",
    "Diet_Type": "diet",
}

grouped = df.groupby("BioSample", sort=False)

sample_df = grouped[TAKE_FIRST_COLS].first()
sample_df = sample_df.reset_index()
sample_df = sample_df.rename(columns={"BioSample": "id", **RENAME_TO_SNAKE_CASE})
sample_df["category"] = "biolink:MaterialSample"

sample_df.to_parquet(OUT_PATH, index=False)

print(f"wrote {len(sample_df)} Sample nodes to {OUT_PATH}")
