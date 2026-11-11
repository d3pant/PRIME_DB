#!/usr/bin/env bash
# Stage the inputs for the PRIMEDB ETL, which runs IN PLACE from this repo (see run_etl.sh).
#
#   stage_inputs.sh DATA_FOLDER
#
# DATA_FOLDER  folder holding the files the dumper downloaded:
#                samples.csv, silva_genus_absolute.csv, silva_genus_relative.csv
#
# The ETL scripts read and write relative to where they live, so this script makes
# <repo>/data/ ready for them, from scratch, on every run:
#   data/raw/        symlinks to the three downloaded files (nothing is copied; the links
#                    are rebuilt each run, so the dated hub folder can be anything)
#   data/extract/    empty, plus the curated reference/phenotype_mapping.parquet
#   data/transformed/{nodes,edges}/   empty
# and clears transform/has_taxon_edge/data/ (body_site_mapping.parquet is written there,
# inside the scripts folder, and has_taxon_edge.py would read a stale one).
# data/load/ and anything else under data/ is left alone.
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "usage: $0 DATA_FOLDER" >&2
    exit 2
fi

REPO="$(cd "$(dirname "$0")" && pwd -P)"
DATA_FOLDER="$(cd "$1" && pwd -P)"
DATA="$REPO/data"
BODY_SITE_DIR="$REPO/transform/has_taxon_edge/data"
MARKER="$DATA/transformed/.primedb_etl_scratch"
CURATED="$REPO/reference/phenotype_mapping.parquet"
INPUTS="samples.csv silva_genus_absolute.csv silva_genus_relative.csv"

# --- checks first: nothing is touched until all of them pass ------------------------------

# Curated, hand-made input (not produced by any script).
if [ ! -f "$CURATED" ]; then
    echo "missing curated input: $CURATED" >&2
    exit 1
fi
for f in $INPUTS; do
    if [ ! -f "$DATA_FOLDER/$f" ]; then
        echo "missing input: $DATA_FOLDER/$f" >&2
        exit 1
    fi
done

# The folders below are wiped, so the downloads must not live inside them.
case "$DATA_FOLDER/" in
    "$DATA"/raw/*|"$DATA"/extract/*|"$DATA"/transformed/*|"$BODY_SITE_DIR"/*)
        echo "refusing: DATA_FOLDER $DATA_FOLDER is inside a folder this script wipes" >&2
        exit 2 ;;
esac

# Only wipe data/ that an earlier staging run made. Anything else in data/raw, data/extract
# or data/transformed (hand-placed files, earlier outputs worth keeping) stops the run.
has_content() {  # the directory holds something besides .gitkeep / .DS_Store
    [ -d "$1" ] && [ -n "$(find "$1" -mindepth 1 ! -name .gitkeep ! -name .DS_Store -print | head -n 1)" ]
}
if [ ! -f "$MARKER" ]; then
    for d in raw extract transformed; do
        if has_content "$DATA/$d"; then
            echo "refusing: $DATA/$d already holds files that were not staged by this script." >&2
            echo "Back them up, then delete data/raw, data/extract and data/transformed and rerun." >&2
            exit 2
        fi
    done
fi

# --- stage -----------------------------------------------------------------------------------

for d in raw extract; do
    mkdir -p "$DATA/$d"
    find "$DATA/$d" -mindepth 1 ! -name .gitkeep -delete   # keeps the tracked .gitkeep
done
rm -rf "$DATA/transformed" "$BODY_SITE_DIR"
mkdir -p "$DATA/transformed/nodes" "$DATA/transformed/edges"
touch "$MARKER"

for f in $INPUTS; do
    ln -s "$DATA_FOLDER/$f" "$DATA/raw/$f"
done
cp "$CURATED" "$DATA/extract/phenotype_mapping.parquet"

echo "staged $DATA from $DATA_FOLDER"
