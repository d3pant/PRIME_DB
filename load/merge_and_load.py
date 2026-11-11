import json
from pathlib import Path

import pandas as pd

# Hardcoded for now -- reproducible cross-machine repo discovery (env var /
# config file / git submodule) was left as an open question during design and
# still needs a real answer before this script works anywhere but this machine.
GMREPO_ROOT = Path("/Users/devangpant/GMREPO_V3")
PRIMEDB_ROOT = Path(__file__).parent.parent

OUT_DIR = PRIMEDB_ROOT / "data/load"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GRAPH_NAME = "primedb_gmrepo_merged"


def dedup_by_id_keep_more_complete(df: pd.DataFrame) -> pd.DataFrame:
    """For duplicate ids (post-union), keep the row with the most non-null fields."""
    df = df.copy()
    df["_non_null_count"] = df.notna().sum(axis=1)
    df = df.sort_values("_non_null_count", ascending=False)
    df = df.drop_duplicates(subset="id", keep="first")
    return df.drop(columns="_non_null_count")


# ===========================================================================
# NODES
# ===========================================================================

# --- Sample nodes -----------------------------------------------------------
# transform/sample/sample_nodes.py now emits snake_case column names directly
# (host_age, sex, BMI, country, diet, host_height, ...), already matching
# GMREPO_V3's naming for the concepts both track -- no renaming needed here
# anymore. Everything else (from both sides) is unioned as-is, nothing dropped.
prime_sample = pd.read_parquet(PRIMEDB_ROOT / "data/transformed/nodes/sample.parquet")
gmrepo_sample = pd.read_parquet(GMREPO_ROOT / "data/transformed/nodes/sample.parquet")

sample_df = pd.concat([prime_sample, gmrepo_sample], ignore_index=True)
sample_df = dedup_by_id_keep_more_complete(sample_df)
print(f"Sample nodes: {len(prime_sample)} (PRIMEDB) + {len(gmrepo_sample)} (GMREPO_V3) -> {len(sample_df)} merged")

# --- Phenotype nodes --------------------------------------------------------
# Identical schema on both sides (id, original_name, category) -- straight union + dedup.
prime_phenotype = pd.read_parquet(PRIMEDB_ROOT / "data/transformed/nodes/phenotype.parquet")
gmrepo_phenotype = pd.read_parquet(GMREPO_ROOT / "data/transformed/nodes/phenotype.parquet")

phenotype_df = pd.concat([prime_phenotype, gmrepo_phenotype], ignore_index=True)
phenotype_df = dedup_by_id_keep_more_complete(phenotype_df)
print(f"Phenotype nodes: {len(prime_phenotype)} (PRIMEDB) + {len(gmrepo_phenotype)} (GMREPO_V3) -> {len(phenotype_df)} merged")

# --- Taxon nodes -------------------------------------------------------------
# Schemas match except PRIMEDB's extra original_name (list) column -- unions in,
# null for GMREPO_V3 rows.
prime_taxon = pd.read_parquet(PRIMEDB_ROOT / "data/transformed/nodes/taxon.parquet")
gmrepo_taxon = pd.read_parquet(GMREPO_ROOT / "data/transformed/nodes/taxon.parquet")

taxon_df = pd.concat([prime_taxon, gmrepo_taxon], ignore_index=True)
taxon_df = dedup_by_id_keep_more_complete(taxon_df)
print(f"Taxon nodes: {len(prime_taxon)} (PRIMEDB) + {len(gmrepo_taxon)} (GMREPO_V3) -> {len(taxon_df)} merged")

# ===========================================================================
# EDGES
# ===========================================================================

# --- HAS_TAXON ---------------------------------------------------------------
# GMREPO_V3 has zero columns PRIMEDB lacks (verified) -- straight union.
# 214 edge ids DO overlap between the two repos: neither side has any internal
# duplicate ids (verified 0 and 0), so these are 214 real (BioSample, Run, taxid)
# triples that both pipelines captured independently, not an id-collision bug.
# Both sides already spell this edge's id with uppercase HAS_TAXON, so the ids
# collide as intended and the dedup below catches them.
prime_has_taxon = pd.read_parquet(PRIMEDB_ROOT / "data/transformed/edges/has_taxon.parquet")
gmrepo_has_taxon = pd.read_parquet(GMREPO_ROOT / "data/transformed/edges/has_taxon.parquet")

has_taxon_df = pd.concat([prime_has_taxon, gmrepo_has_taxon], ignore_index=True)
n_taxon_overlap = has_taxon_df.duplicated(subset="id").sum()
has_taxon_df = dedup_by_id_keep_more_complete(has_taxon_df)
print(
    f"HAS_TAXON edges: {len(prime_has_taxon)} (PRIMEDB) + {len(gmrepo_has_taxon)} (GMREPO_V3) "
    f"-> {len(has_taxon_df)} merged ({n_taxon_overlap} overlapping ids deduped)"
)

# --- HAS_PHENOTYPE ------------------------------------------------------------
# 125 edge ids DO overlap (same BioSample + same resolved phenotype, captured
# independently by both repos) -- confirmed real duplicates, not a bug. Same
# dedup rule applies: keep whichever row has more non-null fields.
#
# This dedup depends on both sides spelling the id's infix identically. They now
# do (uppercase _HAS_PHENOTYPE_), fixed in each pipeline's own transform rather
# than patched up here -- the merge stays a pure union + dedup and does no value
# rewriting. GMREPO_V3 previously emitted a lowercase _has_phenotype_ infix; with
# that mismatch the 125 stop colliding and dedup silently reports 0, which is the
# signal that one side's transform has drifted.
prime_has_phenotype = pd.read_parquet(PRIMEDB_ROOT / "data/transformed/edges/has_phenotype.parquet")
gmrepo_has_phenotype = pd.read_parquet(GMREPO_ROOT / "data/transformed/edges/has_phenotype.parquet")

has_phenotype_df = pd.concat([prime_has_phenotype, gmrepo_has_phenotype], ignore_index=True)
n_overlap = has_phenotype_df.duplicated(subset="id").sum()
has_phenotype_df = dedup_by_id_keep_more_complete(has_phenotype_df)
print(
    f"HAS_PHENOTYPE edges: {len(prime_has_phenotype)} (PRIMEDB) + {len(gmrepo_has_phenotype)} (GMREPO_V3) "
    f"-> {len(has_phenotype_df)} merged ({n_overlap} overlapping ids deduped)"
)

# ===========================================================================
# TYPED CSV GENERATION -- FalkorDB bulk-loader conventions (--enforce-schema),
# matching GMREPO_V3's load/gen_csv/gen_csv.py exactly: one :ID column per node
# CSV, :START_ID/:END_ID per edge CSV, every other column name:TYPE. ARRAY
# cells must be a JSON-array string (loader parses via ast.literal_eval; JSON
# double-quoted strings are valid Python literal syntax).
# ===========================================================================


def array_col_to_literal_str(series: pd.Series) -> pd.Series:
    def to_plain(x):
        if hasattr(x, "item"):  # numpy scalar (int64, float64, ...) -> plain Python type
            return x.item()
        return x

    def encode(val):
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return "[]"
        items = list(val) if not isinstance(val, (list, tuple)) else val
        # json.dumps (not str(list(...))) -- source strings containing apostrophes
        # (e.g. dbGaP study descriptions) would break single-quote Python-repr-style
        # serialization. Filter None out BEFORE dumping: json.dumps([None]) produces
        # "[null]", which is valid JSON but not valid Python literal syntax, and the
        # bulk-loader parses ARRAY cells with ast.literal_eval (which needs `None`,
        # not `null`).
        return json.dumps([to_plain(x) for x in items if x is not None])

    return series.apply(encode)


def write_csv(df: pd.DataFrame, header: list[str], out_path: Path) -> None:
    out_df = df.copy()
    out_df.columns = header

    # Proactive check, before writing: a literal "nan"/"NaN" *string* cell (as
    # opposed to an actual missing value, which na_rep below converts to "")
    # means some column was pre-stringified (e.g. via str()/astype(str) or an
    # f-string) upstream of this function, bypassing na_rep entirely -- this is
    # exactly the failure mode that crashes falkordb-bulk-insert on a
    # DOUBLE/LONG column ("Could not parse 'nan' as a double"). Vectorized
    # DataFrame check (fast, in-memory) rather than re-reading the written CSV,
    # which would be slow and fragile for the multi-GB edge files here.
    for col in out_df.columns:
        if out_df[col].dtype == object:
            # A real missing value (real np.nan/None/pd.NA) is correctly handled
            # by na_rep below -- pd.isna() still recognizes it even inside a
            # mixed-type merged column (verified: concatenating a string-dtype
            # column with pd.NA against a float64 column with np.nan keeps both
            # forms individually na-detectable). Only flag a cell that renders as
            # "nan"/"NaN" text AND is NOT recognized as missing -- i.e. a literal
            # 3-character string baked in upstream (e.g. via str()/f-string),
            # which na_rep can't catch because it's already text, not an actual NaN.
            looks_like_nan_text = out_df[col].astype(str).str.strip().isin(["nan", "NaN"])
            bad = looks_like_nan_text & ~out_df[col].isna()
            if bad.any():
                raise AssertionError(f"{out_path.name} column {col!r} has {bad.sum()} literal 'nan' string cell(s)")

    # na_rep="" converts real missing values (NaN/None/NaT) to an empty CSV cell
    # (NULL to the bulk loader) for every column, numeric or not.
    out_df.to_csv(out_path, index=False, na_rep="")
    print(f"wrote {len(out_df)} rows to {out_path}")


# --- Sample nodes CSV ---------------------------------------------------------
sample_df["collection_date"] = pd.to_datetime(sample_df["collection_date"], errors="coerce").dt.strftime("%Y-%m-%d")
sample_df["project_id"] = array_col_to_literal_str(sample_df["project_id"])
sample_df["original_sample_description"] = array_col_to_literal_str(sample_df["original_sample_description"])

sample_header = [
    "id:ID",
    "host_age:STRING",
    "sex:STRING",
    "BMI:DOUBLE",
    "country:STRING",
    "diet:STRING",
    "host_height:DOUBLE",
    "host_weight:DOUBLE",
    "race_or_ethnicity:STRING",
    "continent:STRING",
    "smoke_status:STRING",
    "drinking_status:STRING",
    "category:STRING",
    "sample_name:STRING",
    "is_disease_stage_available:STRING",
    "disease_stage:STRING",
    "more:STRING",
    "more_info:STRING",
    "collection_date:STRING",
    "longitude:DOUBLE",
    "latitude:DOUBLE",
    "Recent_Antibiotics_Use:STRING",
    "antibiotics_used:STRING",
    "Antibiotics_Dose:STRING",
    "Days_Without_Antibiotics_Use:STRING",
    "project_id:ARRAY",
    "original_sample_description:ARRAY",
]
sample_df = sample_df[
    [
        "id", "host_age", "sex", "BMI", "country", "diet", "host_height", "host_weight",
        "race_or_ethnicity", "continent", "smoke_status", "drinking_status", "category",
        "sample_name", "is_disease_stage_available", "disease_stage", "more", "more_info",
        "collection_date", "longitude", "latitude", "Recent_Antibiotics_Use", "antibiotics_used",
        "Antibiotics_Dose", "Days_Without_Antibiotics_Use", "project_id", "original_sample_description",
    ]
]
write_csv(sample_df, sample_header, OUT_DIR / "sample.csv")

# --- Phenotype nodes CSV -------------------------------------------------------
phenotype_header = ["id:ID", "original_name:STRING", "category:STRING"]
write_csv(phenotype_df, phenotype_header, OUT_DIR / "phenotype.csv")

# --- Taxon nodes CSV ------------------------------------------------------------
taxon_df["parent_taxid"] = taxon_df["parent_taxid"].apply(lambda v: "" if pd.isnull(v) else str(int(v)))
taxon_df["lineage"] = array_col_to_literal_str(taxon_df["lineage"])
taxon_df["original_name"] = array_col_to_literal_str(taxon_df["original_name"])

taxon_header = [
    "id:ID",
    "taxid:LONG",
    "name:STRING",
    "parent_taxid:LONG",
    "lineage:ARRAY",
    "original_name:ARRAY",
    "rank:STRING",
    "category:STRING",
]
taxon_df = taxon_df[["id", "taxid", "name", "parent_taxid", "lineage", "original_name", "rank", "category"]]
write_csv(taxon_df, taxon_header, OUT_DIR / "taxon.csv")

# --- HAS_PHENOTYPE edges CSV -----------------------------------------------------
has_phenotype_df["primary_knowledge_source"] = array_col_to_literal_str(has_phenotype_df["primary_knowledge_source"])

has_phenotype_header = [
    "id:STRING",
    ":START_ID",
    ":END_ID",
    "category:STRING",
    "predicate:STRING",
    "BioProject:STRING",
    "primary_knowledge_source:ARRAY",
    "aggregator_knowledge_source:STRING",
    "agent_type:STRING",
    "publication_doi:STRING",
    "publication_category:STRING",
]
has_phenotype_df = has_phenotype_df[
    [
        "id", "subject", "object", "category", "predicate", "BioProject",
        "primary_knowledge_source", "aggregator_knowledge_source", "agent_type",
        "publication_doi", "publication_category",
    ]
]
write_csv(has_phenotype_df, has_phenotype_header, OUT_DIR / "has_phenotype.csv")

# --- HAS_TAXON edges CSV -----------------------------------------------------------
has_taxon_df["primary_knowledge_source"] = array_col_to_literal_str(has_taxon_df["primary_knowledge_source"])

has_taxon_header = [
    "id:STRING",
    ":START_ID",
    ":END_ID",
    "category:STRING",
    "predicate:STRING",
    "BioProject:STRING",
    "sra_run_id:STRING",
    "project_name:STRING",
    "relative_abundance:DOUBLE",
    "absolute_abundance:DOUBLE",
    "instrument:STRING",
    "primary_knowledge_source:ARRAY",
    "aggregator_knowledge_source:STRING",
    "agent_type:STRING",
    "publication_doi:STRING",
    "publication_category:STRING",
    "anatomical_entity_id:STRING",
    "anatomical_entity_original_name:STRING",
    "anatomical_entity_category:STRING",
]
has_taxon_df = has_taxon_df[
    [
        "id", "subject", "object", "category", "predicate", "BioProject", "sra_run_id",
        "project_name", "relative_abundance", "absolute_abundance", "instrument",
        "primary_knowledge_source", "aggregator_knowledge_source", "agent_type",
        "publication_doi", "publication_category", "anatomical_entity_id",
        "anatomical_entity_original_name", "anatomical_entity_category",
    ]
]
write_csv(has_taxon_df, has_taxon_header, OUT_DIR / "has_taxon.csv")

print(f"\nAll CSVs written to {OUT_DIR}.")
print(f"Next: verify against the existing 'gmrepo' graph, then run load/load_db.py to bulk-insert into new graph '{GRAPH_NAME}'.")
