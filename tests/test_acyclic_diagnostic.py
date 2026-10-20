"""Bounded synthetic tests for label-blind diagnostics; no project datasets used."""

import copy
import hashlib

import pytest
from rdkit import DataStructs

from src.data.acyclic_diagnostic import (
    audit_csv,
    components,
    count_feasible,
    count_tanimoto,
    diagnose_nodes,
    main,
    prepare_nodes,
)


def make_nodes(fingerprints, weights=None):
    if weights is None:
        weights = [1] * len(fingerprints)
    result, position = [], 0
    for index, (counts, weight) in enumerate(zip(fingerprints, weights, strict=True)):
        records = []
        for _ in range(weight):
            position += 1
            records.append(
                {
                    "record_number": position,
                    "Compound ID": f"synthetic-{position}",
                    "smiles": "synthetic-only",
                }
            )
        result.append(
            {
                "identity": str(index),
                "counts": counts,
                "records": records,
                "full_identities": {str(index)},
            }
        )
    return result


def test_count_similarity_matches_rdkit_sparse_count_tanimoto():
    left, right = {1: 2, 2: 1}, {1: 2, 2: 2}
    intersection, union = count_tanimoto(left, right)
    vectors = [DataStructs.UIntSparseIntVect(10), DataStructs.UIntSparseIntVect(10)]
    for vector, counts in zip(vectors, (left, right), strict=True):
        for key, value in counts.items():
            vector[key] = value
    assert (intersection, union) == (3, 4)
    assert intersection / union == DataStructs.TanimotoSimilarity(*vectors)


def test_transitive_components_expose_chaining_and_fixed_thresholds():
    nodes = make_nodes([{1: 2, 2: 1}, {1: 2, 2: 2}, {1: 1, 2: 2}])
    report = diagnose_nodes(nodes)
    assert [item["threshold"] for item in report["threshold_reports"]] == [0.6, 0.7, 0.8]
    assert [item["component_count"] for item in report["threshold_reports"]] == [1, 1, 3]
    for item in report["threshold_reports"][:2]:
        largest = item["largest_five_components"][0]
        assert largest["chaining_below_threshold"]
        assert largest["minimum_pairwise_evidence"]["similarity"] == 0.5
    assert report["threshold_reports"][2]["singleton_component_percentage"] == 100


@pytest.mark.parametrize("intersection, expected", [(700, 1), (699, 2)])
def test_threshold_boundary_is_inclusive_and_exact(intersection, expected):
    report = diagnose_nodes(
        make_nodes([{1: intersection, 2: 1000 - intersection}, {1: intersection}])
    )
    assert report["threshold_reports"][1]["component_count"] == expected


def test_duplicate_records_weight_sizes_but_do_not_duplicate_graph_nodes():
    report = diagnose_nodes(make_nodes([{1: 1}, {2: 1}, {3: 1}], weights=[13, 2, 2]))
    assert report["acyclic_record_count"] == 17
    assert report["conservative_identity_count"] == 3
    item = report["threshold_reports"][0]
    assert item["singleton_component_count"] == 3
    assert item["singleton_component_percentage"] == 100
    assert item["records_in_singleton_components"] == 17
    assert item["largest_component_record_count"] == 13
    assert item["majority_component_over_50_percent"]
    feasibility = item["count_only_feasibility"]
    assert feasibility["rounded_80_10_10_counts"] == [13, 2, 2]
    assert feasibility["exact_rounded_counts_possible"]
    assert "split" not in item["components"][0]


@pytest.mark.parametrize(
    "weights, expected",
    [([8, 1, 1], True), ([4, 4, 1, 1], True), ([9, 1], False), ([7, 2, 1], False), ([10], False)],
)
def test_count_only_feasibility_requires_whole_groups_and_three_nonempty_parts(weights, expected):
    assert count_feasible(weights, [(8, 8), (1, 1), (1, 1)]) is expected


def test_components_cover_every_node_and_keep_every_edge():
    edges = [(0, 1), (1, 2), (3, 4)]
    groups = components(6, edges)
    assert groups == [[0, 1, 2], [3, 4], [5]]
    assert sorted(node for group in groups for node in group) == list(range(6))
    assert all(any(left in group and right in group for group in groups) for left, right in edges)
    assert components(6, list(reversed(edges))) == groups


def test_prepare_only_acyclic_nodes_preserves_source_and_ignores_unrelated_fields():
    class RestrictedRecord(dict):
        def __getitem__(self, key):
            assert key in {"record_number", "Compound ID", "smiles"}
            return super().__getitem__(key)

    source = [
        RestrictedRecord(
            {
                "record_number": i,
                "Compound ID": f"synthetic-{i}",
                "smiles": value,
                "unused": "must not access",
            }
        )
        for i, value in enumerate(["CCO", "OCC", "C[C@H](O)F", "C[C@@H](O)F", "c1ccccc1"], start=1)
    ]
    original = copy.deepcopy(source)
    nodes = prepare_nodes(source)
    assert source == original
    assert len(nodes) == 2
    assert sum(len(node["records"]) for node in nodes) == 4
    assert len(set().union(*(node["full_identities"] for node in nodes))) == 3
    assert nodes == prepare_nodes(list(reversed(source)))


def test_chemical_input_order_does_not_change_graph_evidence():
    source = [
        {"record_number": i, "Compound ID": str(i), "smiles": value}
        for i, value in enumerate(["CCO", "CCC", "CCCO"], start=1)
    ]
    assert diagnose_nodes(prepare_nodes(source)) == diagnose_nodes(prepare_nodes(source[::-1]))


@pytest.mark.parametrize("value", ["", "C(C", "[CH5]", "CCO suffix", "CCO.[Na+]"])
def test_invalid_or_multicomponent_input_requires_attention(value):
    with pytest.raises(ValueError):
        prepare_nodes([{"record_number": 1, "Compound ID": "synthetic", "smiles": value}])


def test_missing_acyclic_population_stops():
    with pytest.raises(ValueError, match="No acyclic"):
        diagnose_nodes([])


def test_csv_label_changes_leave_structural_results_identical(tmp_path):
    source = tmp_path / "synthetic.csv"
    source.write_text(
        "Compound ID,smiles,unused_measurement,unused_prediction\n"
        "first,CCO,never used,also ignored\nsecond,CCC,nothing,ignored\n"
    )
    original = source.read_bytes()
    first = audit_csv(source, hashlib.sha256(original).hexdigest(), 2)
    assert source.read_bytes() == original
    source.write_text(source.read_text().replace("never used", "different contents"))
    second = audit_csv(source, hashlib.sha256(source.read_bytes()).hexdigest(), 2)
    first.pop("source")
    second.pop("source")
    assert first == second


def test_authorized_checksum_and_count_guards(tmp_path):
    source = tmp_path / "synthetic.csv"
    source.write_text("Compound ID,smiles\nfirst,CCO\n")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="checksum"):
        audit_csv(source, "incorrect", 1)
    with pytest.raises(ValueError, match="record count"):
        audit_csv(source, digest, 2)


def test_cli_refuses_existing_output_without_reading_source(tmp_path, monkeypatch):
    output = tmp_path / "report.json"
    output.write_text("existing evidence")
    monkeypatch.setattr(
        "sys.argv",
        [
            "diagnostic",
            str(tmp_path / "raw" / "missing.csv"),
            "--output",
            str(output),
            "--expected-source-sha256",
            "unused",
            "--expected-acyclic-count",
            "2",
        ],
    )
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 2
    assert output.read_text() == "existing evidence"
