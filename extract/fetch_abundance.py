"""Fetch SILVA genus-level abundance data from PrimeDB into data/raw.

Downloads the absolute and relative genus-level abundance tables from the
PrimeDB cache endpoint and saves them as CSV in data/raw.
"""

from pathlib import Path

import requests

BASE_URL = "https://primedb.sjtu.edu.cn/api/v1/abundances/cache/silva/genus"
SAMPLES_URL = "https://primedb.sjtu.edu.cn/api/download/sample?format=csv"
ABUNDANCE_TYPES = ["absolute", "relative"]

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = REPO_ROOT / "data" / "raw"


def fetch_csv(abundance_type: str, dest: Path) -> None:
    url = f"{BASE_URL}/{abundance_type}"
    print(f"Fetching {url}")
    response = requests.get(url, timeout=300)
    response.raise_for_status()
    dest.write_bytes(response.content)
    print(f"Saved {dest} ({dest.stat().st_size:,} bytes)")

def fetch_samples(url: str, dest: Path)-> None:
    print(f"Fetching {url}")
    response = requests.get(url, timeout=300)
    response.raise_for_status()
    dest.write_bytes(response.content)
    print(f"Saved {dest} ({dest.stat().st_size:,} bytes)")


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    for abundance_type in ABUNDANCE_TYPES:
        filename = f"silva_genus_{abundance_type}.csv"
        fetch_csv(abundance_type, RAW_DIR / filename)

    filename = "samples.csv"
    fetch_samples(SAMPLES_URL, RAW_DIR / filename)

if __name__ == "__main__":
    main()
