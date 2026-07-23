from pathlib import Path

import pandas as pd

SAMPLES_PARQUET = Path(__file__).parent.parent.parent / "data/extract/samples.parquet"
ABSOLUTE_PARQUET = Path(__file__).parent.parent.parent / "data/extract/silva_genus_absolute.parquet"
RELATIVE_PARQUET = Path(__file__).parent.parent.parent / "data/extract/silva_genus_relative.parquet"
LINEAGE_LOOKUP_PARQUET = Path(__file__).parent.parent.parent / "data/transformed/nodes/taxon_lineage_lookup.parquet"
BODY_SITE_MAPPING_PARQUET = Path(__file__).parent / "data/body_site_mapping.parquet"
OUT_PATH = Path(__file__).parent.parent.parent / "data/transformed/edges/has_taxon.parquet"

samples_df = pd.read_parquet(SAMPLES_PARQUET)
absolute_df = pd.read_parquet(ABSOLUTE_PARQUET)
relative_df = pd.read_parquet(RELATIVE_PARQUET)
lineage_lookup = pd.read_parquet(LINEAGE_LOOKUP_PARQUET)
body_site_mapping = pd.read_parquet(BODY_SITE_MAPPING_PARQUET)

lineage_to_id = dict(zip(lineage_lookup["lineage_string"], lineage_lookup["id"]))
body_site_lookup = body_site_mapping.set_index("original_name")[["id", "name"]].to_dict("index")

sample_info = samples_df.set_index("Run")[
    ["BioSample", "BioProject", "Project_name", "Instrument", "Doi", "Body_Site"]
].to_dict("index")

lineage_cols = [c for c in absolute_df.columns if c != "Sample"]

# absolute_df and relative_df have identical columns and identical row order
# (verified: same Sample values in the same order) -- safe to index positionally
# by Run without a merge.
absolute_by_run = absolute_df.set_index("Sample")
relative_by_run = relative_df.set_index("Sample")

rows = []
for run in absolute_by_run.index:
    info = sample_info.get(run)
    if info is None:
        continue

    biosample = info["BioSample"]
    body_site_raw = info["Body_Site"]
    body_site_term = body_site_raw.split("/")[-1].strip().lower() if body_site_raw else None
    body_site_match = body_site_lookup.get(body_site_term) if body_site_term else None

    anatomical_entity = None
    if body_site_term:
        anatomical_entity = {
            "id": body_site_match["id"] if body_site_match else None,
            "original_name": body_site_term,
            "category": "biolink:AnatomicalEntity",
        }

    publication = None
    if info["Doi"]:
        publication = {"doi": f"doi:{info['Doi']}", "category": "biolink:Publication"}

    absolute_row = absolute_by_run.loc[run]
    relative_row = relative_by_run.loc[run]

    # Multiple SILVA lineage columns can resolve to the same taxon id (e.g. many
    # unresolved genera falling back to the same shared family/class ancestor --
    # confirmed 48 different genus-level lineages collapsing onto one Clostridia
    # node). Sum abundance across all such lineages per (Run, taxon) before
    # emitting an edge, so id stays a genuine unique key rather than colliding.
    per_taxon_abundance = {}
    for lineage_col in lineage_cols:
        absolute_abundance = absolute_row[lineage_col]
        if pd.isna(absolute_abundance) or absolute_abundance == 0:
            continue

        taxon_id = lineage_to_id.get(lineage_col)
        if taxon_id is None:
            continue

        relative_abundance = relative_row[lineage_col]
        if taxon_id not in per_taxon_abundance:
            per_taxon_abundance[taxon_id] = {"absolute": 0, "relative": 0.0}
        per_taxon_abundance[taxon_id]["absolute"] += absolute_abundance
        per_taxon_abundance[taxon_id]["relative"] += relative_abundance

    for taxon_id, abundance in per_taxon_abundance.items():
        taxid_suffix = taxon_id.split(":", 1)[1]

        rows.append(
            {
                "id": f"{biosample}_{run}_has_taxon_{taxid_suffix}",
                "association": {
                    "category": "MaterialSampleToOrganismTaxonAssociation",
                    "predicate": "has_taxon",
                    "BioProject": info["BioProject"],
                    "SRA Run ID": run,
                    "project_name": info["Project_name"],
                    "relative_abundance": abundance["relative"],
                    "absolute_abundance": abundance["absolute"],
                    "instrument": info["Instrument"],
                    "primary_knowledge_source": ["infores:PRIMEDB"],
                    "aggregator_knowledge_source": None,
                    "agent_type": "biolink:manual_agent",
                    "publication": publication,
                    "anatomical_entity": anatomical_entity,
                },
                "subject": biosample,
                "object": taxon_id,
            }
        )

edge_df = pd.DataFrame(rows)
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
edge_df.to_parquet(OUT_PATH, index=False)

print(f"wrote {len(edge_df)} HAS_TAXON edges to {OUT_PATH}")
