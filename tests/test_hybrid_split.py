"""Tiny synthetic fixtures only: no ESOL, models, downloads, or accelerator access."""

import copy
import csv
import hashlib
import io
import json
from dataclasses import replace
from fractions import Fraction

import pytest
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from src.data import hybrid_split as split


def records(*smiles):
    return [
        {"record_number": i, "Compound ID": f"synthetic-{i}", "smiles": value}
        for i, value in enumerate(smiles, 1)
    ]


@pytest.fixture
def population():
    # Three groups per route, weighted 8/1/1; no measured data.
    return split.prepare_structures(
        records(
            *(["CCO"] * 8),
            "[Na+]",
            "Cl",
            *(["c1ccccc1"] * 8),
            "C1CCCCC1",
            "c1ccncc1",
        )
    )


def candidate(structures):
    groups = split.build_structural_groups(structures)
    assignments = {}
    for route in split.ROUTES:
        assignments.update(split.allocate_route([g for g in groups if g.route == route]))
    digest = "0" * 64
    rows = split.manifest_rows(structures, groups, assignments, digest)
    return groups, assignments, rows, digest


def source(digest="0" * 64):
    return {"path": "synthetic.csv", "sha256": digest, "checksum_unchanged": True}


def write_csv(tmp_path, values=None):
    raw = tmp_path / "raw"
    raw.mkdir()
    path = raw / "synthetic.csv"
    values = values or (["CCO"] * 8 + ["[Na+]", "Cl"] + ["c1ccccc1"] * 8 + ["C1CCCCC1", "c1ccncc1"])
    stream = io.StringIO(newline="")
    writer = csv.writer(stream)
    writer.writerow(["Compound ID", "smiles", "unused_target", "unused_prediction"])
    for i, smiles in enumerate(values, 1):
        writer.writerow([f"synthetic-{i}", smiles, "DO NOT INTERPRET", "UNAVAILABLE"])
    path.write_text(stream.getvalue())
    return path


def arguments(path, prefix, digest=None, total=20, acyclic=10):
    return [
        "hybrid_split",
        str(path),
        "--expected-source-sha256",
        digest or hashlib.sha256(path.read_bytes()).hexdigest(),
        "--expected-record-count",
        str(total),
        "--expected-acyclic-count",
        str(acyclic),
        "--output-prefix",
        str(prefix),
    ]


def test_identity_policy_preserves_duplicates_stereo_isotopes_and_source():
    raw = records(
        " CCO ",
        "OCC",
        "C[C@H](O)F",
        "C[C@@H](O)F",
        "CC(O)F",
        "[CH3:1][OH:2]",
        "CO",
        "[13CH3]O",
    )
    original = copy.deepcopy(raw)
    prepared = split.prepare_structures(raw)
    assert raw == original
    assert len(prepared) == 8
    assert prepared[0].full_identity == prepared[1].full_identity
    assert prepared[2].full_identity != prepared[3].full_identity
    assert len({r.conservative_identity for r in prepared[2:5]}) == 1
    assert prepared[5].full_identity == prepared[6].full_identity
    assert prepared[6].conservative_identity != prepared[7].conservative_identity


def test_no_chemical_standardization_or_fragment_removal():
    prepared = split.prepare_structures(records("CC(=O)C", "CC(O)=C", "C[NH3+]", "CN"))
    assert len({r.conservative_identity for r in prepared}) == 4


@pytest.mark.parametrize("value", ["", "C(C", "[CH5]", "CCO suffix", "CCO.[Na+]", "c1ccccc1.Cl"])
def test_invalid_or_disconnected_structures_stop_without_exclusion(value):
    with pytest.raises(ValueError):
        split.prepare_structures(records("CCO", value))


@pytest.mark.parametrize("numbers", [[1, 1], [0, 1], [True, 2]])
def test_invalid_record_identifiers_stop(numbers):
    raw = records("CCO", "CCC")
    for record, number in zip(raw, numbers, strict=True):
        record["record_number"] = number
    with pytest.raises(ValueError, match="Record numbers"):
        split.prepare_structures(raw)


def test_empty_population_stops():
    with pytest.raises(ValueError, match="No source"):
        split.prepare_structures([])


def test_structural_functions_cannot_access_other_source_fields():
    class Restricted(dict):
        def __getitem__(self, key):
            assert key in {"record_number", "Compound ID", "smiles"}
            return super().__getitem__(key)

    raw = [
        Restricted({**r, "unused_target": object(), "unused_prediction": object()})
        for r in records("CCO", "c1ccccc1")
    ]
    assert len(split.prepare_structures(raw)) == 2


def test_typed_scaffold_groups_and_achiral_scaffold_keys():
    prepared = split.prepare_structures(
        records(
            "c1ccccc1",
            "Cc1ccccc1",
            "c1ccncc1",
            "C[C@H]1CCCCO1",
            "C[C@@H]1CCCCO1",
        )
    )
    groups = split.build_structural_groups(prepared)
    assert sorted(group.weight for group in groups) == [1, 2, 2]
    assert prepared[3].scaffold_key == prepared[4].scaffold_key
    assert all(r.route == "nonempty_murcko" for r in prepared)


def test_pinned_sparse_count_configuration_and_multiplicity():
    prepared = split.prepare_structures(records("CCCCCCCC"))
    mol = Chem.MolFromSmiles("CCCCCCCC")
    generator = rdFingerprintGenerator.GetMorganGenerator(**dict(split.FP_PARAMETERS))
    vector = generator.GetSparseCountFingerprint(mol)
    assert prepared[0].counts == dict(vector.GetNonzeroElements())
    assert max(prepared[0].counts.values()) > 1
    assert split.FP_PARAMETERS["radius"] == 2
    assert split.FP_PARAMETERS["includeChirality"] is False
    assert split._configuration()["device"] == "cpu"
    assert split._configuration()["fingerprint"]["output"] == "GetSparseCountFingerprint"


def test_empty_generated_fingerprint_stops(monkeypatch):
    class Empty:
        def GetNonzeroElements(self):
            return {}

    class Generator:
        def GetSparseCountFingerprint(self, mol):
            return Empty()

    monkeypatch.setattr(
        split.rdFingerprintGenerator, "GetMorganGenerator", lambda **kw: Generator()
    )
    with pytest.raises(ValueError, match="fingerprint"):
        split.prepare_structures(records("CCO"))


def artificial(identity, counts, number):
    return split.Structure(number, identity, identity, identity, "", "empty_murcko", counts)


@pytest.mark.parametrize("intersection, expected_groups", [(700, 1), (699, 2)])
def test_inclusive_exact_threshold(intersection, expected_groups):
    prepared = [
        artificial("a", {1: intersection, 2: 1000 - intersection}, 1),
        artificial("b", {1: intersection}, 2),
    ]
    assert len(split.build_structural_groups(prepared)) == expected_groups


def test_count_tanimoto_agrees_with_rdkit_and_rejects_empty_vectors():
    left, right = DataStructs.UIntSparseIntVect(10), DataStructs.UIntSparseIntVect(10)
    left[1], left[2], right[1], right[2] = 2, 1, 2, 2
    i, u = split.count_tanimoto({1: 2, 2: 1}, {1: 2, 2: 2})
    assert Fraction(i, u) == Fraction(3, 4)
    assert i / u == DataStructs.TanimotoSimilarity(left, right)
    with pytest.raises(ValueError):
        split.count_tanimoto({}, {})


def test_chaining_retains_distant_endpoints_and_isolated_nodes():
    prepared = [
        artificial("a", {1: 2, 2: 1}, 1),
        artificial("b", {1: 2, 2: 2}, 2),
        artificial("c", {1: 1, 2: 2}, 3),
        artificial("isolated", {3: 1}, 4),
    ]
    groups = split.build_structural_groups(prepared)
    assert [g.record_numbers for g in groups] == [(1, 2, 3), (4,)]
    diagnostic = split._group_diagnostics(prepared, groups)["largest_five_acyclic_components"][0]
    assert diagnostic["minimum_pairwise_similarity"] == 0.5
    assert diagnostic["chaining_below_threshold"]


def test_identity_routing_conflict_stops():
    first = artificial("same", {1: 1}, 1)
    second = replace(first, record_number=2, route="nonempty_murcko", scaffold_key="C1CC1")
    with pytest.raises(ValueError, match="conflicting"):
        split.build_structural_groups([first, second])


def weighted_groups(weights, keys=None):
    groups, position = [], 0
    for index, weight in enumerate(weights):
        numbers = tuple(range(position + 1, position + weight + 1))
        groups.append(
            split.Group(
                keys[index] if keys else str(index),
                "empty_murcko",
                numbers,
            )
        )
        position += weight
    return groups


def test_exact_allocation_objective():
    counts = (7, 2, 1)
    assert split.allocation_cost(counts, 10) == Fraction(9, 8)
    assert split.allocation_cost((8, 1, 1), 10) == 0


def test_allocation_largest_first_and_partition_ties():
    groups = weighted_groups([8, 1, 1])
    assigned = split.allocate_route(list(reversed(groups)))
    assert {assigned[n] for n in groups[0].record_numbers} == {"train"}
    assert assigned[9] == "validation"
    assert assigned[10] == "test"
    assert split_counts(assigned) == {"train": 8, "validation": 1, "test": 1}


def split_counts(assigned):
    return {name: list(assigned.values()).count(name) for name in split.SPLITS}


def test_group_key_tie_is_complete_key_lexical():
    groups = weighted_groups([8, 1, 1], ["large", "z-key", "a-key"])
    assigned = split.allocate_route(groups)
    assert assigned[10] == "validation"
    assert assigned[9] == "test"


def test_reservation_keeps_every_partition_nonempty():
    # Without reservation the objective would put more than one of these into train.
    groups = weighted_groups([1, 1, 1])
    assert set(split.allocate_route(groups).values()) == set(split.SPLITS)


def test_allocator_preserves_duplicate_record_weights(population):
    groups, assigned, rows, digest = candidate(population)
    assert len(rows) == len(population) == 20
    assert split_counts(assigned) == {"train": 16, "validation": 2, "test": 2}
    assert split.validate_assignments(population, groups, rows, digest)["passed"]


@pytest.mark.parametrize("weights", [[9, 1], [0, 1, 1]])
def test_impossible_or_empty_groups_stop(weights):
    with pytest.raises(ValueError):
        split.allocate_route(weighted_groups(weights))


def test_overlapping_groups_stop():
    groups = weighted_groups([8, 1, 1])
    groups[2] = replace(groups[2], record_numbers=(1,))
    with pytest.raises(ValueError, match="overlap"):
        split.allocate_route(groups)


@pytest.mark.parametrize("train_count, passes", [(78, True), (77, False), (82, True), (83, False)])
def test_tolerance_is_inclusive_exact_and_not_two_percent_relative(train_count, passes):
    prepared = [artificial(str(i), {i: 1}, i) for i in range(1, 101)]
    # Only inspect the train entry; no full split acceptance is asserted here.
    assigned = {i: "train" if i <= train_count else "test" for i in range(1, 101)}
    report = split._population_report(prepared, assigned)
    assert report["splits"]["train"]["within_tolerance"] is passes


def test_iteration_order_does_not_change_groups_or_manifest(population):
    groups, assigned, rows, digest = candidate(population)
    shuffled = list(reversed(population))
    following = split.build_structural_groups(shuffled)
    assert following == groups
    assert split.manifest_rows(shuffled, following, assigned, digest) == rows


@pytest.mark.parametrize(
    "mutation, failed_check",
    [
        ("drop", "record_conservation"),
        ("duplicate", "record_conservation"),
        ("unknown_split", "partition_values"),
        ("sample_id", "sample_ids"),
        ("identity_field", "manifest_fields"),
        ("group_split", "group_integrity"),
    ],
)
def test_corrupted_manifest_is_rejected(population, mutation, failed_check):
    groups, assigned, rows, digest = candidate(population)
    if mutation == "drop":
        rows.pop()
    elif mutation == "duplicate":
        rows.append(dict(rows[0]))
    elif mutation == "unknown_split":
        rows[0]["split"] = "development"
    elif mutation == "sample_id":
        rows[0]["sample_id"] = rows[1]["sample_id"]
    elif mutation == "identity_field":
        rows[0]["full_identity_id"] = "wrong"
    elif mutation == "group_split":
        rows[0]["split"] = "validation"
    report = split.validate_assignments(population, groups, rows, digest)
    assert not report["passed"]
    assert not report["checks"][failed_check]


def test_duplicate_identity_and_scaffold_leakage_are_detected(population):
    groups, assigned, rows, digest = candidate(population)
    for number in (1, 11):
        rows[number - 1]["split"] = "validation"
    checks = split.validate_assignments(population, groups, rows, digest)["checks"]
    assert not checks["full_identity_integrity"]
    assert not checks["conservative_identity_integrity"]
    assert not checks["nonempty_scaffold_integrity"]


def test_conservative_stereo_match_leakage_is_detected():
    prepared = split.prepare_structures(records("C[C@H](O)F", "C[C@@H](O)F"))
    groups = split.build_structural_groups(prepared)
    rows = split.manifest_rows(prepared, groups, {1: "train", 2: "test"}, "0" * 64)
    checks = split.validate_assignments(prepared, groups, rows, "0" * 64)["checks"]
    assert checks["full_identity_integrity"]
    assert not checks["conservative_identity_integrity"]


def test_qualifying_acyclic_edge_leakage_is_recomputed():
    prepared = [artificial("a", {1: 7, 2: 3}, 1), artificial("b", {1: 7}, 2)]
    # Supply deliberately wrong singleton groups: acceptance must not trust the graph.
    groups = [
        split.Group("bad-a", "empty_murcko", (1,)),
        split.Group("bad-b", "empty_murcko", (2,)),
    ]
    rows = split.manifest_rows(prepared, groups, {1: "train", 2: "test"}, "0" * 64)
    report = split.validate_assignments(prepared, groups, rows, "0" * 64)
    assert not report["checks"]["acyclic_edge_integrity"]
    assert len(report["acyclic_cross_partition_edge_violations"]) == 1


def test_missing_route_coverage_and_empty_partition_are_detected(population):
    groups, assigned, rows, digest = candidate(population)
    for row in rows:
        row["split"] = "train"
    checks = split.validate_assignments(population, groups, rows, digest)["checks"]
    assert not checks["all_partitions_nonempty"]
    assert not checks["both_routes_in_all_partitions"]


def test_nearest_neighbor_audit_reports_cross_route_without_enforcing_it():
    prepared = [
        artificial("a", {1: 1}, 1),
        replace(artificial("ring", {2: 1}, 2), route="nonempty_murcko", scaffold_key="C1CC1"),
        artificial("b", {2: 1}, 3),
    ]
    report = split.audit_structural_similarity(prepared, {1: "train", 2: "train", 3: "test"})
    detail = report["heldout_identities"][0]["nearest_neighbors"]
    assert detail["all_training"]["similarity"] == 1
    assert detail["same_route_training"]["similarity"] == 0
    assert detail["cross_route_training"]["similarity"] == 1
    assert "Only empty-to-empty" in report["scope"]
    summary = next(
        x
        for x in report["summaries"]
        if (x["split"], x["route"], x["comparison"]) == ("test", "empty_murcko", "all_training")
    )
    assert summary["at_or_above_0_70"] == 1
    assert sum(item["count"] for item in summary["bins"]) == 1


def test_nearest_neighbor_ties_are_lexical_and_quantiles_have_fixed_definition():
    prepared = [
        artificial("z", {1: 1}, 1),
        artificial("a", {1: 1}, 2),
        artificial("heldout", {1: 1}, 3),
    ]
    report = split.audit_structural_similarity(prepared, {1: "train", 2: "train", 3: "test"})
    neighbor = report["heldout_identities"][0]["nearest_neighbors"]["all_training"]
    assert neighbor["training_identity_id"] == split._key_id("conservative_identity", "a")
    summary = split._similarity_summary([Fraction(0), Fraction(1)])
    assert summary["quantiles"]["0.5"] == 0.5
    assert summary["quantiles"]["0.95"] == 0.95


def test_artifact_schema_hashes_no_labels_and_pending_owner_acceptance(population):
    groups, assigned, rows, digest = candidate(population)
    payloads = split.build_artifacts(population, groups, assigned, source(), {"synthetic": True})
    manifest = list(csv.DictReader(io.StringIO(payloads["manifest"].decode())))
    metadata = json.loads(payloads["metadata"])
    acceptance = json.loads(payloads["acceptance"])
    assert tuple(manifest[0]) == split.COLUMNS
    assert len(manifest) == 20
    assert (
        metadata["artifact_sha256"]["manifest"] == hashlib.sha256(payloads["manifest"]).hexdigest()
    )
    assert (
        metadata["artifact_sha256"]["acceptance"]
        == hashlib.sha256(payloads["acceptance"]).hexdigest()
    )
    assert acceptance["passed"]
    assert acceptance["status"] == "validated_candidate"
    assert acceptance["seal_status"] == "not_sealed"
    assert "pending" in metadata["owner_artifact_acceptance"]
    assert b"DO NOT INTERPRET" not in b"".join(payloads.values())
    assert (
        split.build_artifacts(population, groups, assigned, source(), {"synthetic": True})[
            "manifest"
        ]
        == payloads["manifest"]
    )


def test_hard_failure_produces_report_only(population):
    groups, assigned, rows, digest = candidate(population)
    assigned[1] = "validation"
    payloads = split.build_artifacts(population, groups, assigned, source(), {})
    assert set(payloads) == {"acceptance"}
    report = json.loads(payloads["acceptance"])
    assert not report["hard_checks_passed"]
    assert report["status"] == "rejected_candidate"


def test_hash_collision_stops(monkeypatch, population):
    groups = split.build_structural_groups(population)
    monkeypatch.setattr(split, "_key_id", lambda kind, key: "collision")
    with pytest.raises(ValueError, match="collision"):
        split._key_mappings(population, groups)


def test_tolerance_failure_preserves_rejected_candidate_without_search():
    prepared = split.prepare_structures(
        records(
            *(["CCO"] * 10),
            "[Na+]",
            "Cl",
            *(["c1ccccc1"] * 10),
            "C1CCCCC1",
            "c1ccncc1",
        )
    )
    groups, assigned, rows, digest = candidate(prepared)
    payloads = split.build_artifacts(prepared, groups, assigned, source(), {})
    report = json.loads(payloads["acceptance"])
    assert set(payloads) == {"manifest", "metadata", "acceptance"}
    assert report["hard_checks_passed"]
    assert not report["checks"]["proportion_tolerance"]
    assert report["status"] == "rejected_candidate"
    assert split._configuration()["allocation"]["passes"] == 1


def test_cli_success_only_reads_synthetic_csv_and_preserves_it(tmp_path, monkeypatch):
    path = write_csv(tmp_path)
    before = path.read_bytes()
    prefix = tmp_path / "output"
    monkeypatch.setattr("sys.argv", arguments(path, prefix))
    assert split.main() == 0
    assert path.read_bytes() == before
    report = json.loads((tmp_path / "output_acceptance.json").read_text())
    assert all(report["checks"].values())
    assert report["source"]["checksum_unchanged"]
    metadata = json.loads((tmp_path / "output_metadata.json").read_text())
    assert metadata["execution"]["construction_and_checks_seconds"] >= 0
    assert not list(tmp_path.glob("*.tmp"))


def test_synthetic_target_changes_do_not_change_assignment(tmp_path):
    path = write_csv(tmp_path)
    first_digest = hashlib.sha256(path.read_bytes()).hexdigest()
    first = split.prepare_structures(split.read_structural_records(path, first_digest))
    path.write_text(path.read_text().replace("DO NOT INTERPRET", "DIFFERENT CONTENT"))
    second_digest = hashlib.sha256(path.read_bytes()).hexdigest()
    following = split.prepare_structures(split.read_structural_records(path, second_digest))
    assert first_digest != second_digest
    assert following == first
    assert candidate(first)[1] == candidate(following)[1]


@pytest.mark.parametrize(
    "contents",
    [
        "smiles\nCCO\n",
        "Compound ID,smiles,smiles\nx,CCO,CCO\n",
        "Compound ID,smiles\nx,CCO,extra\n",
        'Compound ID,smiles\n"unterminated',
        "Compound ID,smiles\nx,CCO\x00\n",
    ],
)
def test_malformed_synthetic_csv_stops(tmp_path, contents):
    path = tmp_path / "malformed.csv"
    path.write_text(contents)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError):
        split.read_structural_records(path, digest)


@pytest.mark.parametrize("guard", ["checksum", "total", "acyclic", "version"])
def test_cli_input_guards_publish_nothing(tmp_path, monkeypatch, guard):
    path = write_csv(tmp_path)
    prefix = tmp_path / "output"
    argv = arguments(path, prefix)
    if guard == "checksum":
        argv[argv.index("--expected-source-sha256") + 1] = "incorrect"
    elif guard == "total":
        argv[argv.index("--expected-record-count") + 1] = "21"
    elif guard == "acyclic":
        argv[argv.index("--expected-acyclic-count") + 1] = "9"
    elif guard == "version":
        monkeypatch.setattr(split.rdBase, "rdkitVersion", "unexpected")
    monkeypatch.setattr("sys.argv", argv)
    with pytest.raises(SystemExit) as error:
        split.main()
    assert error.value.code == 2
    assert not list(tmp_path.glob("output*"))


def test_source_change_during_construction_publishes_nothing(tmp_path, monkeypatch):
    path = write_csv(tmp_path)
    original_build = split.build_artifacts

    def changed(*args, **kwargs):
        result = original_build(*args, **kwargs)
        path.write_bytes(path.read_bytes() + b"\n")
        return result

    monkeypatch.setattr(split, "build_artifacts", changed)
    monkeypatch.setattr("sys.argv", arguments(path, tmp_path / "output"))
    with pytest.raises(SystemExit) as error:
        split.main()
    assert error.value.code == 2
    assert not list(tmp_path.glob("output*"))


def test_cli_tolerance_failure_returns_one_and_does_not_retry(tmp_path, monkeypatch):
    path = write_csv(
        tmp_path, (["CCO"] * 10 + ["[Na+]", "Cl"] + ["c1ccccc1"] * 10 + ["C1CCCCC1", "c1ccncc1"])
    )
    original_allocate = split.allocate_route
    calls = []

    def counted(groups):
        calls.append(groups[0].route)
        return original_allocate(groups)

    monkeypatch.setattr(split, "allocate_route", counted)
    monkeypatch.setattr("sys.argv", arguments(path, tmp_path / "output", total=24, acyclic=12))
    assert split.main() == 1
    assert calls == list(split.ROUTES)
    report = json.loads((tmp_path / "output_acceptance.json").read_text())
    assert report["status"] == "rejected_candidate"
    assert report["seal_status"] == "not_sealed"


def test_existing_output_and_raw_location_refused_before_source_read(tmp_path, monkeypatch):
    path = write_csv(tmp_path)
    prefix = tmp_path / "output"
    existing = tmp_path / "output_metadata.json"
    existing.write_text("preserve existing")

    def forbidden(*args):
        raise AssertionError("source should not be read")

    monkeypatch.setattr(split, "read_structural_records", forbidden)
    for chosen in (prefix, path.parent / "output"):
        monkeypatch.setattr("sys.argv", arguments(path, chosen))
        with pytest.raises(SystemExit) as error:
            split.main()
        assert error.value.code == 2
    assert existing.read_text() == "preserve existing"


def test_cli_has_no_methodology_or_sealing_overrides(tmp_path, monkeypatch):
    path = write_csv(tmp_path)
    monkeypatch.setattr("sys.argv", arguments(path, tmp_path / "output") + ["--threshold", "0.8"])
    with pytest.raises(SystemExit) as error:
        split.main()
    assert error.value.code == 2


def test_publication_acceptance_last_and_failure_cleans_created_files(tmp_path, monkeypatch):
    raw = tmp_path / "raw"
    raw.mkdir()
    source_path = raw / "synthetic.csv"
    source_path.write_text("unused synthetic")
    original_link = split.os.link
    order = []

    def failing_link(left, right):
        order.append(right.name)
        if right.name.endswith("_acceptance.json"):
            raise OSError("simulated publication interruption")
        return original_link(left, right)

    monkeypatch.setattr(split.os, "link", failing_link)
    with pytest.raises(OSError, match="interruption"):
        split.write_artifacts(
            tmp_path / "output",
            {"manifest": b"csv", "metadata": b"{}", "acceptance": b"{}"},
            source_path,
        )
    assert order == ["output_manifest.csv", "output_metadata.json", "output_acceptance.json"]
    assert not list(tmp_path.glob("output*"))
    assert not list(tmp_path.glob(".*.tmp"))


def test_publication_race_never_overwrites_other_writer(tmp_path, monkeypatch):
    raw = tmp_path / "raw"
    raw.mkdir()
    path = raw / "synthetic.csv"
    path.write_text("synthetic")
    original_link = split.os.link

    def competing_link(left, right):
        if right.name.endswith("_metadata.json"):
            right.write_bytes(b"other writer")
        return original_link(left, right)

    monkeypatch.setattr(split.os, "link", competing_link)
    with pytest.raises(FileExistsError):
        split.write_artifacts(
            tmp_path / "output",
            {"manifest": b"csv", "metadata": b"{}"},
            path,
        )
    assert (tmp_path / "output_metadata.json").read_bytes() == b"other writer"
    assert not (tmp_path / "output_manifest.csv").exists()


def test_manifest_record_numbers_require_integers(population):
    groups, assigned, rows, digest = candidate(population)
    rows[0]["source_record_number"] = True
    report = split.validate_assignments(population, groups, rows, digest)
    assert not report["checks"]["record_number_types"]
    assert not report["hard_checks_passed"]


def test_unverified_source_cannot_produce_usable_bundle(population):
    groups, assigned, rows, digest = candidate(population)
    unverified = {**source(), "checksum_unchanged": False}
    payloads = split.build_artifacts(population, groups, assigned, unverified, {})
    assert set(payloads) == {"acceptance"}
    report = json.loads(payloads["acceptance"])
    assert not report["checks"]["source_checksum_unchanged"]
    assert not report["passed"]
