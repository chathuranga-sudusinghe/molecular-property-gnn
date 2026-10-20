"""Synthetic structural diagnostics only; no real data, models, or network access."""

import copy
import csv
import json

import pytest

from src.data.structural_audit import audit_structural_csv, audit_structures, main


def records(*smiles):
    return [
        {"Compound ID": f"synthetic-{position}", "smiles": value}
        for position, value in enumerate(smiles, start=1)
    ]


def test_empty_scaffolds_count_records_and_unique_identities_without_assigning():
    source = records("CCO", "OCC", "CCC", "c1ccccc1", "Cc1ccccc1")
    original = copy.deepcopy(source)
    report = audit_structures(source)
    assert source == original
    assert report["record_count"] == report["valid_record_count"] == 5
    assert report["empty_scaffolds"]["record_count"] == 3
    assert report["empty_scaffolds"]["percentage_of_valid_records"] == 60
    assert report["empty_scaffolds"]["unique_isomeric_identities"] == 2
    assert report["empty_scaffolds"]["percentage_of_unique_isomeric_identities"] == 50
    assert report["nonempty_scaffold_groups"]["sizes_descending"] == [2]
    capacity = report["hypothetical_single_empty_bucket_capacity"]["partitions"]
    assert capacity["train"]["empty_bucket_fits_budget"]
    assert not capacity["test"]["empty_bucket_fits_budget"]
    assert "split" not in report


def test_canonical_and_stereo_free_duplicate_groups_are_separate():
    report = audit_structures(records("CCO", " OCC ", "C[C@H](O)F", "C[C@@H](O)F"))
    assert report["canonical_duplicates"]["group_count"] == 1
    assert report["canonical_duplicates"]["extra_records"] == 1
    assert report["stereo_free_duplicates"]["group_count"] == 2
    assert report["stereo_free_duplicates"]["affected_records"] == 4
    group = report["canonical_duplicates"]["groups"][0]
    assert [record["record_number"] for record in group["records"]] == [1, 2]


def test_typed_scaffolds_do_not_merge_different_ring_elements():
    report = audit_structures(records("c1ccccc1", "Cc1ccccc1", "c1ccncc1"))
    assert report["nonempty_scaffold_groups"]["sizes_descending"] == [2, 1]
    assert report["nonempty_scaffold_groups"]["size_histogram"]["1"]["group_count"] == 1


def test_multicomponent_structures_and_flags_preserve_every_fragment():
    report = audit_structures(records("CCO.[Na+]", "[Na+].[Cl-]", "C=C", "CCC", "CCO"))
    assert report["multicomponent_record_count"] == 2
    assert report["multicomponent_records"][0]["canonical_identity"] == "CCO.[Na+]"
    assert report["multicomponent_records"][1]["fragment_count"] == 2
    categories = report["empty_scaffolds"]["categories"]
    assert categories["contains_formally_charged_atoms"]["record_count"] == 2
    assert categories["carbon_free"]["record_count"] == 1
    assert categories["hydrocarbon_only"]["record_count"] == 2
    assert categories["contains_double_or_triple_bonds"]["record_count"] == 1


def test_atom_maps_ignored_for_keys_and_isotopes_preserved_in_identity():
    report = audit_structures(records("[CH3:1][OH:2]", "CO", "[13CH3]O"))
    assert report["atom_mapped_record_count"] == 1
    assert report["canonical_duplicates"]["groups"][0]["record_count"] == 2
    assert report["unique_valid_isomeric_identities"] == 2


@pytest.mark.parametrize("smiles", ["", "C(C", "[CH5]", "CCO suffix"])
def test_invalid_records_reported_without_assignment_or_silent_exclusion(smiles):
    report = audit_structures(records(smiles, "CCO"))
    assert report["record_count"] == 2
    assert report["valid_record_count"] == 1
    assert report["invalid_records"][0]["record_number"] == 1
    assert report["invalid_records"][0]["smiles"] == smiles


def test_empty_input_is_serializable():
    report = audit_structures([])
    assert report["empty_scaffolds"]["percentage_of_valid_records"] is None
    assert report["hypothetical_single_empty_bucket_capacity"]["partitions"] == {}
    json.dumps(report, allow_nan=False)


def test_unrelated_fields_are_never_accessed():
    class RestrictedRecord(dict):
        def __getitem__(self, key):
            assert key in {"Compound ID", "smiles"}
            return super().__getitem__(key)

    source = RestrictedRecord({"Compound ID": "synthetic", "smiles": "CCO", "unused": object()})
    assert audit_structures([source])["valid_record_count"] == 1


def test_csv_ignores_other_fields_and_preserves_source_bytes(tmp_path):
    path = tmp_path / "synthetic.csv"
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["Compound ID", "smiles", "unused_measurement", "unused_prediction"])
        writer.writerow(["synthetic", "CCO", "not numeric", "never interpreted"])
    original = path.read_bytes()
    first = audit_structural_csv(path)
    assert path.read_bytes() == original
    path.write_text(path.read_text().replace("not numeric", "different value"))
    second = audit_structural_csv(path)
    assert path.read_bytes() != original  # test changes its own synthetic fixture only
    assert first.pop("source")["checksum_unchanged"]
    second.pop("source")
    assert first == second


def test_cli_rejects_writing_inside_input_directory_before_reading(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "sys.argv",
        ["audit", str(tmp_path / "missing.csv"), "--output", str(tmp_path / "report.json")],
    )
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 2
    assert not (tmp_path / "report.json").exists()


def test_cli_refuses_to_overwrite_existing_output(tmp_path, monkeypatch):
    source = tmp_path / "raw" / "missing.csv"
    output = tmp_path / "report.json"
    output.write_text("preserve existing evidence")
    monkeypatch.setattr("sys.argv", ["audit", str(source), "--output", str(output)])
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 2
    assert output.read_text() == "preserve existing evidence"


def test_cli_writes_only_new_structural_report(tmp_path, monkeypatch):
    raw = tmp_path / "raw"
    raw.mkdir()
    source = raw / "synthetic.csv"
    source.write_text("Compound ID,smiles,unused_measurement\nsynthetic,CCO,secret-value\n")
    output = tmp_path / "report.json"
    before = source.read_bytes()
    monkeypatch.setattr("sys.argv", ["audit", str(source), "--output", str(output)])
    assert main() == 0
    assert source.read_bytes() == before
    assert json.loads(output.read_text())["empty_scaffolds"]["record_count"] == 1
    assert "secret-value" not in output.read_text()


@pytest.mark.parametrize(
    "content",
    ["smiles\nCCO\n", "Compound ID,smiles\nx,CCO,extra\n", 'Compound ID,smiles\n"unterminated'],
)
def test_malformed_synthetic_csv_rejected(tmp_path, content):
    source = tmp_path / "malformed.csv"
    source.write_text(content)
    with pytest.raises(ValueError):
        audit_structural_csv(source)
