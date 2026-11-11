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

    anatomical_entity_id = None
    anatomical_entity_original_name = None
    anatomical_entity_category = None
    if body_site_term:
        anatomical_entity_id = body_site_match["id"] if body_site_match else None
        anatomical_entity_original_name = body_site_term
        anatomical_entity_category = "biolink:AnatomicalEntity"

    publication_doi = None
    publication_category = None
    if info["Doi"]:
        publication_doi = f"doi:{info['Doi']}"
        publication_category = "biolink:Publication"

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
                # HAS_TAXON (SCREAMING_SNAKE_CASE) matches the relationship type used
                # at bulk-load time and GMREPO_V3's own id convention for this edge.
                "id": f"{biosample}_{run}_HAS_TAXON_{taxid_suffix}",
                "subject": biosample,
                "object": taxon_id,
                # Deliberately NOT biolink:-prefixed, unlike has_phenotype's category
                # ("biolink:MaterialSampleToDiseaseOrPhenotypicFeatureAssociation").
                # The prefix standardization below applies to `predicate` only; both
                # this repo and GMREPO_V3 agree on this exact unprefixed category value
                # for HAS_TAXON, so changing it would break that agreement for no gain.
                "category": "MaterialSampleToOrganismTaxonAssociation",
                # Prefixed style, consistent with has_phenotype's predicate below --
                # PRIMEDB previously mixed "has_taxon" (bare) and "biolink:has_phenotype"
                # (prefixed) across its two edge types; standardized on prefixed.
                "predicate": "biolink:has_taxon",
                # BioProject: the study ACCESSION (e.g. "PRJDB13875") -- matches
                # GMREPO_V3's own edge-level "BioProject" property name/casing exactly.
                # project_name: a separate, human-readable study name (e.g.
                # "MUSC-JP-2024"), not a duplicate/alias of BioProject -- intentionally
                # kept as its own property, GMREPO_V3's edges have no equivalent.
                "BioProject": info["BioProject"],
                "sra_run_id": run,
                "project_name": info["Project_name"],
                "relative_abundance": abundance["relative"],
                "absolute_abundance": abundance["absolute"],
                "instrument": info["Instrument"],
                "primary_knowledge_source": ["infores:PRIMEDB"],
                "aggregator_knowledge_source": None,
                "agent_type": "biolink:manual_agent",
                "publication_doi": publication_doi,
                "publication_category": publication_category,
                "anatomical_entity_id": anatomical_entity_id,
                "anatomical_entity_original_name": anatomical_entity_original_name,
                "anatomical_entity_category": anatomical_entity_category,
            }
        )

edge_df = pd.DataFrame(rows)
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
edge_df.to_parquet(OUT_PATH, index=False)

print(f"wrote {len(edge_df)} HAS_TAXON edges to {OUT_PATH}")
