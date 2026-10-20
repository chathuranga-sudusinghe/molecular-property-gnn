"""Small synthetic fixtures exercise validation without depending on ignored raw data."""

import csv
import hashlib
import json
import subprocess
import sys

import pandas as pd
import pytest

from src.data.audit import (
    EXPECTED_COLUMNS,
    SMILES_COLUMN,
    TARGET_COLUMN,
    audit_csv,
    audit_dataframe,
    validate_smiles,
)


def make_frame(smiles=("CCO", "O", "CC"), targets=(-3.0, -2.0, -1.0)):
    rows = []
    for position, (structure, target) in enumerate(zip(smiles, targets, strict=True)):
        row = dict.fromkeys(EXPECTED_COLUMNS, 1)
        row.update(
            {"Compound ID": f"compound-{position}", SMILES_COLUMN: structure, TARGET_COLUMN: target}
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=EXPECTED_COLUMNS)


def write_csv(path, frame):
    frame.to_csv(path, index=False)
    return path


@pytest.mark.parametrize("smiles", ["CCO", " c1ccccc1 \t", "[Na+].[Cl-]", "C[C@H](O)C(=O)O"])
def test_valid_smiles(smiles):
    assert validate_smiles(smiles) is None


@pytest.mark.parametrize(
    "smiles",
    [None, float("nan"), 12, "", " \t", "C(C", "[CH5]", "CCO trailing-name", "CCO |garbage|"],
)
def test_invalid_smiles(smiles):
    assert validate_smiles(smiles) is not None


def test_statistics_and_distribution_cover_every_finite_target():
    report = audit_dataframe(make_frame())
    assert report["passes_validation"]
    assert (report["row_count"], report["column_count"]) == (3, 10)
    assert report["columns"] == list(EXPECTED_COLUMNS)
    assert report["dtypes"][TARGET_COLUMN] == "float64"
    assert not any(report["missing_values"].values())
    assert report["target_statistics"] == {
        "count": 3,
        "mean": -2.0,
        "std_sample": 1.0,
        "min": -3.0,
        "q25": -2.5,
        "median": -2.0,
        "q75": -1.5,
        "max": -1.0,
        "skew": 0.0,
    }
    histogram = report["target_histogram"]
    assert sum(histogram["counts"]) == 3
    assert len(histogram["bin_edges"]) == 11
    assert histogram["bin_edges"][0] == -3.0
    assert histogram["bin_edges"][-1] == -1.0


def test_duplicate_definitions_and_no_dataframe_mutation():
    frame = make_frame(("CCO", " CCO ", "CCO"), (-1.0, -1.0, -1.0))
    frame.loc[2] = frame.loc[0]
    original = frame.copy(deep=True)
    report = audit_dataframe(frame)
    pd.testing.assert_frame_equal(frame, original)
    assert report["duplicate_rows"] == {"extra_rows": 1, "affected_rows": 2}
    assert report["duplicate_smiles_exact"]["extra_rows"] == 1
    assert report["duplicate_smiles_exact"]["affected_rows"] == 2
    assert report["duplicate_smiles_trimmed"] == {
        "extra_rows": 2,
        "affected_rows": 3,
        "groups": [{"smiles": "CCO", "count": 3, "record_numbers": [1, 2, 3]}],
    }
    assert report["smiles_surrounding_whitespace_count"] == 1
    assert report["passes_validation"]  # duplicates are findings, not a deletion/failure policy


def test_missing_smiles_not_counted_as_duplicate_keys():
    frame = make_frame((None, " ", "CCO"))
    report = audit_dataframe(frame)
    assert report["missing_values"][SMILES_COLUMN] == 2
    assert len(report["invalid_smiles"]) == 2
    assert report["duplicate_smiles_trimmed"]["extra_rows"] == 0
    assert not report["passes_validation"]


def test_invalid_records_are_retained_and_positioned_independently_of_index():
    frame = make_frame(("CCO", "[CH5]", "CC"), ("bad", "inf", "-2.5"))
    frame.index = [20, 20, 99]
    report = audit_dataframe(frame)
    assert report["row_count"] == len(frame) == 3
    assert report["valid_smiles_count"] == 2
    assert report["invalid_smiles"][0]["record_number"] == 2
    assert report["invalid_smiles"][0]["value"] == "[CH5]"
    assert [item["record_number"] for item in report["invalid_targets"]] == [1, 2]
    assert report["target_statistics"]["count"] == 1
    assert report["target_statistics"]["mean"] == -2.5
    assert report["target_statistics"]["std_sample"] is None
    assert report["target_statistics"]["skew"] is None
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize("target", [None, "", " ", float("inf"), -float("inf"), "NaN"])
def test_missing_or_nonfinite_targets(target):
    report = audit_dataframe(make_frame(("CCO",), (target,)))
    assert len(report["invalid_targets"]) == 1
    assert report["target_statistics"]["count"] == 0
    assert report["target_statistics"]["min"] is None
    assert report["target_histogram"]["counts"] == []
    assert not report["passes_validation"]
    json.dumps(report, allow_nan=False)


def test_auxiliary_column_missingness_is_reported():
    frame = make_frame()
    frame.loc[0, "Compound ID"] = " \t"
    report = audit_dataframe(frame)
    assert report["missing_values"]["Compound ID"] == 1
    assert not report["passes_validation"]


@pytest.mark.parametrize("change", ["missing", "extra", "duplicate", "reordered"])
def test_schema_mismatch_rejected(change):
    frame = make_frame()
    if change == "missing":
        frame = frame.drop(columns=TARGET_COLUMN)
    elif change == "extra":
        frame["unexpected"] = 1
    elif change == "duplicate":
        frame.columns = [*EXPECTED_COLUMNS[:-1], TARGET_COLUMN]
    else:
        frame = frame[list(reversed(frame.columns))]
    with pytest.raises(ValueError, match="Expected columns"):
        audit_dataframe(frame)


def test_csv_bytes_unchanged_and_literal_na_is_preserved(tmp_path):
    frame = make_frame()
    frame.loc[0, "Compound ID"] = "NA"
    path = write_csv(tmp_path / "esol.csv", frame)
    before = path.read_bytes()
    report = audit_csv(path)
    assert path.read_bytes() == before
    assert report["source"]["sha256"] == hashlib.sha256(before).hexdigest()
    assert report["missing_values"]["Compound ID"] == 0
    assert report["passes_validation"]
    assert report["versions"]["rdkit"]


@pytest.mark.parametrize("row", [["too", "short"], ["extra"] * 11, []])
def test_ragged_or_blank_csv_record_rejected(tmp_path, row):
    path = tmp_path / "malformed.csv"
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(EXPECTED_COLUMNS)
        writer.writerow(row)
    with pytest.raises(ValueError, match="Record 1"):
        audit_csv(path)


def test_header_only_dataset_fails_validation_without_nan(tmp_path):
    path = write_csv(tmp_path / "empty.csv", make_frame((), ()))
    report = audit_csv(path)
    assert report["row_count"] == 0
    assert not report["passes_validation"]
    assert report["target_statistics"]["mean"] is None
    json.dumps(report, allow_nan=False)


def test_malformed_quotes_rejected(tmp_path):
    path = tmp_path / "malformed.csv"
    path.write_text(",".join(EXPECTED_COLUMNS) + '\n"unterminated')
    with pytest.raises(ValueError, match="Malformed CSV"):
        audit_csv(path)


def test_nul_byte_rejected_without_truncation(tmp_path):
    path = write_csv(tmp_path / "nul.csv", make_frame(("CCO\x00invalid",), (-1.0,)))
    with pytest.raises(ValueError, match="NUL byte"):
        audit_csv(path)


@pytest.mark.parametrize("content", ["", "smiles,target\nCCO,-1\n"])
def test_empty_file_or_wrong_csv_header_rejected(tmp_path, content):
    path = tmp_path / "wrong.csv"
    path.write_text(content)
    with pytest.raises(ValueError, match="Expected CSV header"):
        audit_csv(path)


@pytest.mark.parametrize("valid, expected_code", [(True, 0), (False, 1)])
def test_cli_reports_json_and_validation_exit_status(tmp_path, valid, expected_code):
    frame = make_frame(("CCO" if valid else "[CH5]",), (-1.0,))
    path = write_csv(tmp_path / "esol.csv", frame)
    result = subprocess.run(
        [sys.executable, "-m", "src.data.audit", str(path)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == expected_code
    assert json.loads(result.stdout)["passes_validation"] is valid


def test_cli_missing_file_exits_with_input_error(tmp_path):
    result = subprocess.run(
        [sys.executable, "-m", "src.data.audit", str(tmp_path / "missing.csv")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 2
    assert "Audit could not read" in result.stderr
