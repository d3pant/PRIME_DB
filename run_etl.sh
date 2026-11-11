#!/usr/bin/env bash
# Run the PRIMEDB ETL, unchanged, one script after another, IN PLACE in this repo.
#
#   stage_inputs.sh DATA_FOLDER     # first: links the downloaded files into data/raw/ and
#   run_etl.sh                      # then this: runs the scripts; results land in
#                                   #   data/transformed/
#
# The scripts read and write relative to where they live (data/extract, data/transformed),
# so this overwrites this repo's data/ -- which is why it refuses to run unless
# stage_inputs.sh has prepared data/ (it leaves a marker in data/transformed/).
#
# Set PRIMEDB_PYTHON to pick the interpreter (default: python3).
set -euo pipefail

REPO="$(cd "$(dirname "$0")" && pwd -P)"
PY="${PRIMEDB_PYTHON:-python3}"

if [ ! -f "$REPO/data/transformed/.primedb_etl_scratch" ]; then
    echo "refusing: data/ was not staged -- run stage_inputs.sh DATA_FOLDER first" >&2
    exit 2
fi

run() {
    echo
    echo "=== $1 ==="
    "$PY" "$REPO/$1"
}

run extract/csv_to_parquet.py
run transform/sample/sample_nodes.py
run transform/taxon/taxon_nodes.py
run transform/phenotype/phenotype_nodes.py
run transform/has_taxon_edge/body_site_mapping.py
run transform/has_taxon_edge/has_taxon_edge.py
run transform/has_phenotype_edge/has_phenotype_edge.py

echo
echo "=== ETL finished: outputs in $REPO/data/transformed ==="
