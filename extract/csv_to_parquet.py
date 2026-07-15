"""Convert every CSV in data/raw to Parquet in data/extract."""

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = REPO_ROOT / "data" / "raw"
EXTRACT_DIR = REPO_ROOT / "data" / "extract"


def convert(csv_path: Path) -> None:
    parquet_path = EXTRACT_DIR / (csv_path.stem + ".parquet")
    df = pd.read_csv(csv_path, low_memory=False)

    object_columns = df.select_dtypes(include="object").columns
    df[object_columns] = df[object_columns].astype("string")

    df.to_parquet(parquet_path, index=False)
    print(f"Wrote {parquet_path} ({parquet_path.stat().st_size:,} bytes)")


def main() -> None:
    EXTRACT_DIR.mkdir(parents=True, exist_ok=True)

    csv_files = sorted(RAW_DIR.glob("*.csv"))
    if not csv_files:
        print(f"No CSV files found in {RAW_DIR}")
        return

    for csv_path in csv_files:
        convert(csv_path)


if __name__ == "__main__":
    main()
