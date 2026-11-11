"""BioThings parser for PRIMEDB.

load_data(data_folder) does not re-implement the ETL: it runs stage_inputs.sh (links the
files the dumper downloaded into the repo's data/raw/ and resets the rest of data/), then
run_etl.sh (the repo's unchanged ETL scripts, run in place), then streams the ETL's parquet
outputs from <repo>/data/transformed as documents -- HAS_TAXON and HAS_PHENOTYPE edges in
one stream, each edge carrying its subject (sample) and object (taxon / phenotype) node
properties.

Inputs in data_folder (written by the dumper, named from Content-Disposition):
    samples.csv, silva_genus_absolute.csv, silva_genus_relative.csv
data_folder must not be the repo's own data/raw (stage_inputs.sh wipes that folder).

This module expects to live inside the PRIMEDB repo checkout (plugin/primedb/) so it can
find the two scripts; set PRIMEDB_REPO to point at the checkout otherwise. The ETL's own
dependencies (pandas, pyarrow, node-synonymizer, biothings_client, requests) must be
installed in the interpreter running the hub.

Dev switch: PRIMEDB_SKIP_ETL=1 reuses the existing <repo>/data/transformed output instead
of re-running the ETL (which calls live APIs and takes minutes).
"""

from __future__ import annotations

import collections
import logging
import os
import re
import subprocess
import sys
from pathlib import Path

import pyarrow.parquet as pq

log = logging.getLogger("primedb_parser")

REPO = Path(os.environ.get("PRIMEDB_REPO") or Path(__file__).resolve().parents[2])
STAGE_SCRIPT = REPO / "stage_inputs.sh"
ETL_SCRIPT = REPO / "run_etl.sh"
BATCH_ROWS = 50_000

# Node properties copied into each document. Taxon nodes also hold `lineage` and
# `original_name` (every SILVA lineage merged into the node, up to dozens of long strings),
# which would be repeated on millions of edges, so they are left out.

# NodeSynonymizer (used by phenotype_nodes.py / body_site_mapping.py) logs a dead batch and
# returns "no match" instead of raising, so the ETL exits 0 with silently missing values.
# The scripts can't be changed, so their output is watched for these messages instead.
SILENT_API_FAILURE = re.compile(r"giving up after|finished with \d+ failed batches")


class EtlFailure(RuntimeError):
    """A staging/ETL script failed, or an API silently dropped values during the ETL."""
    ##if youre wrapping around errors in this class, can we add specific error types here?


def _run_script(script: Path, args: list[str]) -> list[str]:
    """Run a bash script, logging each output line. Returns the lines that match
    SILENT_API_FAILURE; raises EtlFailure if the script is missing or exits nonzero."""
    if not script.exists():
        raise EtlFailure(f"{script} not found -- run from a PRIMEDB checkout or set PRIMEDB_REPO")

    env = {**os.environ, "PRIMEDB_PYTHON": sys.executable}
    proc = subprocess.Popen(
        ["bash", str(script), *args],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
    )
    tail = collections.deque(maxlen=40)
    silent = []
    for line in proc.stdout:
        line = line.rstrip("\n")
        tail.append(line)
        log.info("etl: %s", line)
        if SILENT_API_FAILURE.search(line):
            silent.append(line)
    if proc.wait() != 0:
        raise EtlFailure(f"{script.name} exited {proc.returncode}; last output:\n" + "\n".join(tail))
    return silent


def run_etl(data_folder: Path) -> None:
    _run_script(STAGE_SCRIPT, [str(data_folder)])
    silent = _run_script(ETL_SCRIPT, [])
    if silent:
        raise EtlFailure("an API dropped values silently during the ETL:\n" + "\n".join(silent[:5]))


def _read_nodes(path: Path, fields: list[str] | None = None) -> dict:
    """{id: node dict without null properties}"""
    nodes = {}
    for row in pq.read_table(path, columns=fields).to_pylist():
        nodes[row["id"]] = {k: v for k, v in row.items() if v is not None}
    return nodes


def _present(value) -> bool:
    return value is not None and value == value  # drops None and NaN


def _iter_edges(path: Path, subjects: dict, objects: dict):
    for batch in pq.ParquetFile(path).iter_batches(batch_size=BATCH_ROWS):
        for row in batch.to_pylist():
            subject_id, object_id = row.pop("subject"), row.pop("object")
            yield {
                "_id": row.pop("id"),
                "association": {k: v for k, v in row.items() if _present(v)},
                "subject": dict(subjects.get(subject_id)),
                "object": dict(objects.get(object_id)),
            }


def load_data(data_folder):
    data_folder = Path(data_folder).resolve()
    out = REPO / "data" / "transformed"

    if os.environ.get("PRIMEDB_SKIP_ETL") and (out / "edges" / "has_taxon.parquet").exists():
        log.warning("PRIMEDB_SKIP_ETL set: reusing ETL output in %s", out)
    else:
        run_etl(data_folder)

    samples = _read_nodes(out / "nodes" / "sample.parquet")
    taxa = _read_nodes(out / "nodes" / "taxon.parquet")
    phenotypes = _read_nodes(out / "nodes" / "phenotype.parquet")

    counts = collections.Counter()
    # phenotype edges are the small table; taxon edges are millions of rows
    for doc in _iter_edges(out / "edges" / "has_phenotype.parquet", samples, phenotypes):
        counts[doc["association"]["predicate"]] += 1
        yield doc
    for doc in _iter_edges(out / "edges" / "has_taxon.parquet", samples, taxa):
        counts[doc["association"]["predicate"]] += 1
        yield doc
    log.info("PRIMEDB parse finished: %s", dict(counts))
  