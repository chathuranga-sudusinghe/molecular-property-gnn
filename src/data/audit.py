"""Audit the Delaney ESOL CSV without cleaning, filtering, or writing data."""

import argparse
import csv
import hashlib
import io
import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase

SMILES_COLUMN = "smiles"
TARGET_COLUMN = "measured log solubility in mols per litre"
EXPECTED_COLUMNS = (
    "Compound ID",
    "ESOL predicted log solubility in mols per litre",
    "Minimum Degree",
    "Molecular Weight",
    "Number of H-Bond Donors",
    "Number of Rings",
    "Number of Rotatable Bonds",
    "Polar Surface Area",
    TARGET_COLUMN,
    SMILES_COLUMN,
)


def _missing(series: pd.Series) -> pd.Series:
    """Only nulls and empty/whitespace-only text are missing (not literal 'NA')."""
    return series.isna() | series.map(lambda value: isinstance(value, str) and not value.strip())


def validate_smiles(value: Any) -> str | None:
    """Return an invalidity reason, or None for a nonempty sanitized molecule.

    Strip surrounding whitespace in memory only. Reject RDKit name suffixes and
    CXSMILES extensions so a valid prefix cannot conceal trailing invalid text.
    """
    if not isinstance(value, str):
        return "missing or non-string SMILES"
    if not value.strip():
        return "empty SMILES"
    parameters = Chem.SmilesParserParams()
    parameters.parseName = False
    parameters.allowCXSMILES = False
    molecule = Chem.MolFromSmiles(value.strip(), parameters)
    if molecule is None:
        return "RDKit parsing or sanitization failed"
    if molecule.GetNumAtoms() == 0:
        return "molecule has no atoms"
    return None


def _duplicate_summary(series: pd.Series) -> dict[str, Any]:
    """Count extra occurrences and all affected records, excluding missing keys."""
    populated = series.loc[~_missing(series)]
    groups = [
        {
            "smiles": value,
            "count": int(count),
            "record_numbers": [
                int(position + 1) for position in np.flatnonzero(series.eq(value).fillna(False))
            ],
        }
        for value, count in populated.value_counts(sort=False).items()
        if count > 1
    ]
    return {
        "extra_rows": int(populated.duplicated().sum()),
        "affected_rows": int(populated.duplicated(keep=False).sum()),
        "groups": groups,
    }


def audit_dataframe(frame: pd.DataFrame) -> dict[str, Any]:
    """Profile an ESOL table without mutating it; record numbers start at one.

    A schema mismatch raises ValueError. Content problems remain in the report;
    no rows are dropped. Target statistics use only finite numeric values.
    """
    if list(frame.columns) != list(EXPECTED_COLUMNS):
        raise ValueError(f"Expected columns in this order: {list(EXPECTED_COLUMNS)!r}")

    smiles = frame[SMILES_COLUMN]
    invalid_smiles = []
    invalid_targets = []
    numeric_target = pd.to_numeric(frame[TARGET_COLUMN], errors="coerce")
    finite = numeric_target.notna() & np.isfinite(numeric_target)
    for position, (_, row) in enumerate(frame.iterrows(), start=1):
        reason = validate_smiles(row[SMILES_COLUMN])
        if reason:
            invalid_smiles.append(
                {
                    "record_number": position,
                    "compound_id": str(row["Compound ID"]),
                    "value": None if pd.isna(row[SMILES_COLUMN]) else str(row[SMILES_COLUMN]),
                    "reason": reason,
                }
            )
        if not finite.iloc[position - 1]:
            value = row[TARGET_COLUMN]
            invalid_targets.append(
                {
                    "record_number": position,
                    "compound_id": str(row["Compound ID"]),
                    "value": None if pd.isna(value) else str(value),
                    "reason": "target must be a finite numeric value",
                }
            )

    target = numeric_target.loc[finite].astype(float)
    statistics = {
        "count": int(target.size),
        "mean": float(target.mean()) if len(target) else None,
        "std_sample": float(target.std(ddof=1)) if len(target) > 1 else None,
        "min": float(target.min()) if len(target) else None,
        "q25": float(target.quantile(0.25)) if len(target) else None,
        "median": float(target.median()) if len(target) else None,
        "q75": float(target.quantile(0.75)) if len(target) else None,
        "max": float(target.max()) if len(target) else None,
        "skew": float(target.skew()) if len(target) > 2 else None,
    }
    # Automatic edges span all finite observations; the final bin includes its right edge.
    counts, edges = np.histogram(target, bins=10) if len(target) else ([], [])
    histogram = {
        "counts": [int(count) for count in counts],
        "bin_edges": [float(edge) for edge in edges],
        "convention": "[left, right), except the final bin is [left, right]",
    }
    missing_values = {column: int(_missing(frame[column]).sum()) for column in frame}
    duplicate_rows = {
        "extra_rows": int(frame.duplicated().sum()),
        "affected_rows": int(frame.duplicated(keep=False).sum()),
    }
    trimmed = smiles.map(lambda value: value.strip() if isinstance(value, str) else value)
    return {
        "row_count": len(frame),
        "column_count": len(frame.columns),
        "columns": list(frame.columns),
        "dtypes": {column: str(frame[column].dtype) for column in frame},
        "missing_values": missing_values,
        "duplicate_rows": duplicate_rows,
        "duplicate_smiles_exact": _duplicate_summary(smiles),
        "duplicate_smiles_trimmed": _duplicate_summary(trimmed),
        "smiles_column": SMILES_COLUMN,
        "smiles_surrounding_whitespace_count": int(
            sum(isinstance(value, str) and value != value.strip() for value in smiles)
        ),
        "valid_smiles_count": len(frame) - len(invalid_smiles),
        "invalid_smiles": invalid_smiles,
        "target_column": TARGET_COLUMN,
        "target_definition": "log10(S / (mol/L)), where S is measured aqueous solubility",
        "invalid_targets": invalid_targets,
        "target_statistics": statistics,
        "target_histogram": histogram,
        "passes_validation": bool(
            len(frame)
            and not invalid_smiles
            and not invalid_targets
            and not any(missing_values.values())
        ),
    }


def audit_csv(path: str | Path) -> dict[str, Any]:
    """Read one byte snapshot, reject malformed CSV structure, and attach provenance."""
    path = Path(path)
    content = path.read_bytes()
    if b"\x00" in content:
        raise ValueError("CSV contains a NUL byte; refusing parser-dependent truncation")
    rows = csv.reader(io.StringIO(content.decode("utf-8-sig"), newline=""), strict=True)
    try:
        header = next(rows, None)
        if header != list(EXPECTED_COLUMNS):
            raise ValueError(f"Expected CSV header: {list(EXPECTED_COLUMNS)!r}; got {header!r}")
        for record_number, row in enumerate(rows, start=1):
            if len(row) != len(header):
                raise ValueError(
                    f"Record {record_number} has {len(row)} fields; expected {len(header)}"
                )
    except csv.Error as error:
        raise ValueError(f"Malformed CSV: {error}") from error

    frame = pd.read_csv(
        io.BytesIO(content),
        encoding="utf-8-sig",
        keep_default_na=False,
        skip_blank_lines=False,
        dtype={"Compound ID": str, SMILES_COLUMN: str},
    )
    report = audit_dataframe(frame)
    report["source"] = {"path": path.as_posix(), "sha256": hashlib.sha256(content).hexdigest()}
    report["versions"] = {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "rdkit": rdBase.rdkitVersion,
    }
    return report


def main() -> int:
    """Print JSON only to stdout; validation issues produce exit code 1."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path, help="Path to the raw ESOL CSV")
    arguments = parser.parse_args()
    try:
        report = audit_csv(arguments.path)
    except (OSError, ValueError) as error:
        parser.exit(2, f"Audit could not read the dataset: {error}\n")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["passes_validation"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
