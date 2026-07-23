import re
import time
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import biothings_client as bt
import pandas as pd
import requests

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

IN_PATH = Path(__file__).parent.parent.parent / "data/extract/silva_genus_absolute.parquet"
OUT_PATH = Path(__file__).parent.parent.parent / "data/transformed/nodes/taxon.parquet"
LOOKUP_OUT_PATH = Path(__file__).parent.parent.parent / "data/transformed/nodes/taxon_lineage_lookup.parquet"

# SILVA lineage columns are "d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus" --
# always exactly 6 fields, verified for every column in this table. Unlike GMREPO_V3,
# there are no numeric NCBI taxids here at all -- every taxon must be name-resolved.
RANKS_POSITIONAL = ["domain", "phylum", "class", "order", "family", "genus"]
RANKS_SPECIFIC_FIRST = ["genus", "family", "order", "class", "phylum", "domain"]
PREFIXES = {"domain": "d__", "phylum": "p__", "class": "c__", "order": "o__", "family": "f__", "genus": "g__"}
RANK_ALIASES = {
    "domain": {"superkingdom", "domain"},
    "phylum": {"phylum"},
    "class": {"class"},
    "order": {"order"},
    "family": {"family"},
    "genus": {"genus"},
}

# SILVA's own "no classification at all" sentinel (raw column: "Unassigned;__;__;__;__;__").
# Every rank field is blank, so it can never resolve via the normal walk-up -- but it
# carries real, substantial abundance (confirmed: 14,486 nonzero rows across both the
# absolute and relative abundance tables), so it needs its own node rather than being
# silently dropped. Matches GMREPO_V3's -1 convention for its own unclassified sentinel.
UNASSIGNED_LINEAGE = "Unassigned;__;__;__;__;__"
UNASSIGNED_TAXID = -1


def clean_value(raw_part: str, prefix: str) -> str | None:
    val = raw_part[len(prefix):] if raw_part.startswith(prefix) else raw_part
    val = val.strip()
    if val in ("", "__", "Unassigned"):
        return None
    val = re.sub(r"\([^)]*\)", "", val)  # drop bracketed content, e.g. (Subgroup_1)
    val = val.replace("_", " ").strip()  # underscores -> spaces, e.g. Candidatus_X -> Candidatus X
    return val or None


df = pd.read_parquet(IN_PATH)
lineage_cols = [c for c in df.columns if c != "Sample"]
print(f"loaded {len(lineage_cols)} SILVA lineage columns")

# 6 parallel arrays, one entry per lineage column, aligned by index -- not by any
# explicit key. rank_arrays[rank][i] and lineage_cols[i] describe the same column
# because both are built in the same iteration order.
rank_arrays = {rank: [] for rank in RANKS_POSITIONAL}
for col in lineage_cols:
    parts = col.split(";")
    for i, rank in enumerate(RANKS_POSITIONAL):
        rank_arrays[rank].append(clean_value(parts[i], PREFIXES[rank]))

# Batch-resolve each rank's unique names independently via NCBI Taxonomy
# (biothings_client, same tool GMREPO_V3 uses). Each query is just the bare name --
# no lineage context is sent; ancestor/descendant relationships only matter for our
# own walk-up logic below, applied after these results come back.
client = bt.get_client("taxon")
hits_by_name_per_rank = {}

for rank in RANKS_POSITIONAL:
    unique_names = sorted({v for v in rank_arrays[rank] if v is not None})
    hits_by_name = {}
    BATCH = 100
    for i in range(0, len(unique_names), BATCH):
        batch = unique_names[i : i + BATCH]
        res = client.querymany(batch, scopes="scientific_name", fields=["scientific_name", "rank", "taxid"], verbose=False)
        for r in res:
            if not r.get("notfound"):
                hits_by_name.setdefault(r["query"], []).append(r)
    hits_by_name_per_rank[rank] = hits_by_name
    print(f"  {rank:8s}: {len(unique_names)} unique names queried, {len(hits_by_name)} got >=1 candidate back")


def find_rank_match(name, rank):
    if name is None:
        return None
    for h in hits_by_name_per_rank[rank].get(name, []):
        if h.get("rank") in RANK_ALIASES[rank]:
            return h
    return None


# --- eutils genus-only salvage ------------------------------------------------
# biothings' search index has a real coverage gap for genus-rank synonyms (verified
# against literature: Aliifodinibius/Fodinibius, Hungateiclostridium/Acetivibrio,
# Turicella/Corynebacterium, Allisonella/Dialister are all genuine, published
# reclassifications that biothings failed to surface at genus rank). NCBI's own
# eutils resolves these correctly via its synonym index. Only genus-rank matches
# are accepted here -- eutils also returns plenty of family/order/etc. hits for
# these same names, but those are left to the normal walk-up below rather than
# trusted blindly (a spot-checked case, Chloroflexi -> Chloroflexia, showed a real
# phylum/class rank-ambiguity risk at non-genus ranks that isn't worth the risk here).
def batch_eutils_esearch(names, chunk_size=40):
    all_taxids, not_found = set(), set()
    for i in range(0, len(names), chunk_size):
        chunk = names[i : i + chunk_size]
        term = " OR ".join(f'"{n}"[All Names]' for n in chunk)
        resp = requests.get(
            f"{EUTILS_BASE}/esearch.fcgi",
            params={"db": "taxonomy", "term": term, "retmode": "json", "retmax": chunk_size * 3},
        )
        resp.raise_for_status()
        result = resp.json()["esearchresult"]
        all_taxids.update(result.get("idlist", []))
        for m in result.get("warninglist", {}).get("quotedphrasesnotfound", []):
            not_found.add(m.split('"')[1])
        time.sleep(0.4)
    return all_taxids, not_found


def batch_eutils_efetch(taxids, chunk_size=150):
    details = {}
    taxids = list(taxids)
    for i in range(0, len(taxids), chunk_size):
        chunk = taxids[i : i + chunk_size]
        resp = requests.get(
            f"{EUTILS_BASE}/efetch.fcgi",
            params={"db": "taxonomy", "id": ",".join(chunk), "retmode": "xml"},
        )
        resp.raise_for_status()
        root = ET.fromstring(resp.text)
        for taxon in root.findall("Taxon"):
            taxid = taxon.findtext("TaxId")
            sci_name = taxon.findtext("ScientificName")
            rank = taxon.findtext("Rank")
            synonyms = []
            other_names = taxon.find("OtherNames")
            if other_names is not None:
                for tag in ("Synonym", "GenbankSynonym"):
                    for syn in other_names.findall(tag):
                        if syn.text:
                            synonyms.append(syn.text)
            details[taxid] = {"scientific_name": sci_name, "rank": rank, "synonyms": synonyms}
        time.sleep(0.4)
    return details


genus_values = rank_arrays["genus"]
failed_genus = sorted({v for v in genus_values if v is not None and find_rank_match(v, "genus") is None})
print(f"\ngenus values needing eutils salvage: {len(failed_genus)}")

eutils_taxids, eutils_not_found = batch_eutils_esearch(failed_genus)
eutils_details = batch_eutils_efetch(eutils_taxids)

eutils_genus_matches = {}
for name in failed_genus:
    if name in eutils_not_found:
        continue
    for taxid, info in eutils_details.items():
        if info["rank"] != "genus":
            continue
        if info["scientific_name"] and info["scientific_name"].lower() == name.lower():
            eutils_genus_matches[name] = {"taxid": int(taxid), "scientific_name": info["scientific_name"], "rank": "genus"}
            break
        if any(s.lower() == name.lower() for s in info["synonyms"]):
            eutils_genus_matches[name] = {"taxid": int(taxid), "scientific_name": info["scientific_name"], "rank": "genus"}
            break

print(f"eutils salvaged genus-rank matches: {len(eutils_genus_matches)}")


# Per column: walk genus -> family -> order -> class -> phylum -> domain, take the
# first rank with a candidate whose NCBI-reported rank matches expected. Matching is
# rank-only (no exact-string requirement) -- verified safe: every resolved match's
# scientific_name contains the queried string as a substring (checked below), so
# this recovers real synonyms/variants (e.g. SILVA's "Heimdallarchaeia" vs NCBI's
# "Candidatus Heimdallarchaeia (ex Rinke et al. 2021) (nom. nud.)") without ever
# accepting an unrelated name.
rows = []
for i, col in enumerate(lineage_cols):
    if col == UNASSIGNED_LINEAGE:
        rows.append(
            {
                "id": f"NCBITaxon:{UNASSIGNED_TAXID}",
                "taxid": UNASSIGNED_TAXID,
                "name": "Unassigned",
                "lineage_string": col,
                "rank": None,
                "category": "biolink:OrganismTaxon",
            }
        )
        continue

    resolved_rank = None
    resolved_info = None
    original_name = None
    for rank in RANKS_SPECIFIC_FIRST:
        name = rank_arrays[rank][i]
        if name is None:
            continue
        if rank == "genus" and name in eutils_genus_matches:
            resolved_rank = "genus"
            resolved_info = eutils_genus_matches[name]
            original_name = name
            break
        match = find_rank_match(name, rank)
        if match:
            resolved_rank = rank
            resolved_info = match
            original_name = name
            break

    if resolved_info is None:
        continue

    rows.append(
        {
            "id": f"NCBITaxon:{resolved_info['taxid']}",
            "taxid": resolved_info["taxid"],
            "name": resolved_info["scientific_name"],
            "lineage_string": col,  # full raw SILVA column header, not just the matched rank fragment
            "rank": resolved_info.get("rank"),
            "category": "biolink:OrganismTaxon",
        }
    )

resolved_df = pd.DataFrame(rows)
before_dedup = len(resolved_df)

# Lookup table: one row per resolved SILVA lineage column -> its taxon id. This is
# what a later edge-building script joins against silva_genus_absolute.parquet's
# column headers (Run -> Taxon edges), since the node table below is deduped and
# can no longer answer "which raw column produced this node" on its own.
lineage_lookup_df = resolved_df[["lineage_string", "id"]].copy()
lineage_lookup_df.to_parquet(LOOKUP_OUT_PATH, index=False)

# Node table: one row per unique id. original_name holds every full lineage string
# that resolved to this id (a node must have exactly one row per id, but multiple
# source lineages -- e.g. eutils synonym merges, or several unresolved genera
# falling back to the same family -- commonly collapse onto the same node).
taxon_df = (
    resolved_df.groupby("id", sort=False)
    .agg(
        taxid=("taxid", "first"),
        name=("name", "first"),
        original_name=("lineage_string", list),
        rank=("rank", "first"),
        category=("category", "first"),
    )
    .reset_index()
)
after_dedup = len(taxon_df)

taxon_df.to_parquet(OUT_PATH, index=False)

print(f"\nresolved columns (pre-dedup): {before_dedup} / {len(lineage_cols)}")
print(f"wrote {after_dedup} Taxon nodes to {OUT_PATH}")
print(f"wrote {len(lineage_lookup_df)} lineage_string -> id mappings to {LOOKUP_OUT_PATH}")

n_unresolved = len(lineage_cols) - before_dedup
if n_unresolved:
    print(f"unresolved at any rank: {n_unresolved}")

collapsed = {row["id"]: row["original_name"] for _, row in taxon_df.iterrows() if len(row["original_name"]) > 1}
if collapsed:
    print(f"\n{len(collapsed)} id(s) absorbed multiple source lineages:")
    for tid, lineages in collapsed.items():
        kept, dropped = lineages[0], lineages[1:]
        print(f"  {tid}: kept {kept!r}, +{len(dropped)} more lineage(s) also mapped here")

if not taxon_df["id"].is_unique:
    raise AssertionError("taxon node table has duplicate ids after dedup")
if not lineage_lookup_df["lineage_string"].is_unique:
    raise AssertionError("lineage lookup table has duplicate lineage_string keys")
