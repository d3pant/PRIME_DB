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

grouped = df.groupby("BioSample", sort=False)

sample_df = grouped[TAKE_FIRST_COLS].first()
sample_df = sample_df.reset_index()
sample_df = sample_df.rename(columns={"BioSample": "id"})
sample_df["category"] = "biolink:MaterialSample"

sample_df.to_parquet(OUT_PATH, index=False)

print(f"wrote {len(sample_df)} Sample nodes to {OUT_PATH}")
