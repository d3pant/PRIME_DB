import subprocess
import sys
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data/load"
BULK_INSERT_BIN = Path(sys.executable).parent / "falkordb-bulk-insert"

GRAPH_NAME = "primedb_gmrepo_merged"
SERVER_URL = "redis://127.0.0.1:6379"

cmd = [
    str(BULK_INSERT_BIN),
    GRAPH_NAME,
    "--server-url", SERVER_URL,
    "--enforce-schema",
    # Default escapechar (\) makes the CSV reader treat backslash specially,
    # corrupting the \" sequences inside JSON-encoded ARRAY cells that contain
    # an embedded double-quote (e.g. a dbGaP study title with an internal
    # quotation). This is the exact failure mode documented in GMREPO_V3's own
    # load_db.py -- disabling it relies on plain CSV "" quote-doubling instead,
    # which is what pandas.to_csv already writes.
    "--escapechar", "none",
    # Relationship types are SCREAMING_SNAKE_CASE (standard Cypher/Neo4j-style
    # convention, matching GMREPO_V3's own --relations-with-type invocation) --
    # NOT the same thing as the "predicate" property value stored on each edge
    # (that's a biolink-model data property, e.g. "biolink:has_taxon"; see
    # transform/has_taxon_edge/has_taxon_edge.py and has_phenotype_edge.py).
    "--nodes-with-label", "Sample", str(DATA_DIR / "sample.csv"),
    "--nodes-with-label", "Taxon", str(DATA_DIR / "taxon.csv"),
    "--nodes-with-label", "Phenotype", str(DATA_DIR / "phenotype.csv"),
    "--relations-with-type", "HAS_PHENOTYPE", str(DATA_DIR / "has_phenotype.csv"),
    "--relations-with-type", "HAS_TAXON", str(DATA_DIR / "has_taxon.csv"),
    "--index", "Sample:id",
    "--index", "Taxon:id",
    "--index", "Phenotype:id",
]

if __name__ == "__main__":
    # Deliberately requires an explicit --yes flag: this creates a NEW graph
    # (GRAPH_NAME) on the shared FalkorDB instance and should only run after the
    # regenerated CSVs (data/load/*.csv) have been verified against the
    # existing 'gmrepo' graph -- see load/merge_and_load.py's printed next step.
    if "--yes" not in sys.argv:
        print("Refusing to run without --yes -- verify data/load/*.csv against 'gmrepo' first.")
        print("Command that would run:")
        print(" ".join(cmd))
        sys.exit(1)
    subprocess.run(cmd, check=True)
