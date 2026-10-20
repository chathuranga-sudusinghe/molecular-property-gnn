"""CPU-only, label-blind primary split; real-data execution requires owner approval."""

import argparse
import csv
import hashlib
import io
import json
import os
import platform
import subprocess
import tempfile
import time
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

from src.data.acyclic_diagnostic import FP_PARAMETERS, components, count_tanimoto

PROTOCOL = "esol-hybrid-structure-v1"
RDKIT_VERSION = "2024.09.6"
SPLITS = ("train", "validation", "test")
ROUTES = ("nonempty_murcko", "empty_murcko")
FRACTIONS = (Fraction(4, 5), Fraction(1, 10), Fraction(1, 10))
THRESHOLD = Fraction(7, 10)
TOLERANCE = Fraction(1, 50)
ROOT = Path(__file__).resolve().parents[2]
COLUMNS = (
    "sample_id",
    "source_record_number",
    "compound_id",
    "full_identity_id",
    "conservative_identity_id",
    "structure_route",
    "structural_group_id",
    "split",
)


@dataclass(frozen=True)
class Structure:
    record_number: int
    compound_id: str
    full_identity: str
    conservative_identity: str
    scaffold_key: str
    route: str
    counts: Mapping[int, int]


@dataclass(frozen=True)
class Group:
    key: str
    route: str
    record_numbers: tuple[int, ...]

    @property
    def weight(self) -> int:
        return len(self.record_numbers)


def _json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    )


def _sha(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _key_id(kind: str, key: str) -> str:
    return _sha(_json({"schema": 1, "kind": kind, "key": key}).encode("utf-8"))


def _sample_id(digest: str, number: int) -> str:
    return f"esol:{digest}:{number:06d}"


def read_structural_records(path: Path, expected_sha256: str) -> list[dict[str, Any]]:
    """Project away all other CSV fields before structural functions receive records."""
    content = path.read_bytes()
    if _sha(content) != expected_sha256:
        raise ValueError("Source checksum differs from the approved revision")
    if b"\x00" in content:
        raise ValueError("CSV contains a NUL byte")
    reader = csv.reader(io.StringIO(content.decode("utf-8-sig"), newline=""), strict=True)
    try:
        header = next(reader, [])
        if len(header) != len(set(header)) or not {"Compound ID", "smiles"} <= set(header):
            raise ValueError("Expected unique headers including Compound ID and smiles")
        id_index, smiles_index = header.index("Compound ID"), header.index("smiles")
        records = []
        for number, row in enumerate(reader, 1):
            if len(row) != len(header):
                raise ValueError(f"Unexpected field count at record {number}")
            records.append(
                {
                    "record_number": number,
                    "Compound ID": row[id_index],
                    "smiles": row[smiles_index],
                }
            )
    except csv.Error as error:
        raise ValueError(f"Malformed CSV: {error}") from error
    return records


def prepare_structures(records: Iterable[Mapping[str, Any]]) -> list[Structure]:
    """Preserve records; fingerprint conservative copies for grouping and audit only."""
    parameters = Chem.SmilesParserParams()
    parameters.parseName = False
    parameters.allowCXSMILES = False
    generator = rdFingerprintGenerator.GetMorganGenerator(**dict(FP_PARAMETERS))
    fingerprints: dict[str, dict[int, int]] = {}
    result = []
    numbers = set()
    for source in records:
        number, smiles = source["record_number"], source["smiles"]
        if type(number) is not int or number <= 0 or number in numbers:
            raise ValueError("Record numbers must be unique positive integers")
        numbers.add(number)
        molecule = Chem.MolFromSmiles(smiles.strip(), parameters) if smiles.strip() else None
        if molecule is None or molecule.GetNumAtoms() == 0:
            raise ValueError(f"Invalid molecule at record {number}; no exclusion applied")
        if len(Chem.GetMolFrags(molecule)) != 1:
            raise ValueError(f"Unexpected multicomponent structure at record {number}")
        molecule = Chem.Mol(molecule)
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
        full = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
        scaffold_key = Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=False)
        route = ROUTES[0] if scaffold.GetNumAtoms() else ROUTES[1]
        Chem.RemoveStereochemistry(molecule)
        conservative = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        if conservative not in fingerprints:
            counts = dict(generator.GetSparseCountFingerprint(molecule).GetNonzeroElements())
            if not counts or any(value <= 0 for value in counts.values()):
                raise ValueError(f"Unexpected empty/nonpositive fingerprint at record {number}")
            fingerprints[conservative] = counts
        result.append(
            Structure(
                number,
                source["Compound ID"],
                full,
                conservative,
                scaffold_key,
                route,
                fingerprints[conservative],
            )
        )
    if not result:
        raise ValueError("No source records")
    return sorted(result, key=lambda item: item.record_number)


def _nodes(structures: list[Structure]) -> dict[str, Structure]:
    nodes: dict[str, Structure] = {}
    full_to_conservative: dict[str, str] = {}
    for record in structures:
        old = nodes.setdefault(record.conservative_identity, record)
        if (old.route, old.scaffold_key, old.counts) != (
            record.route,
            record.scaffold_key,
            record.counts,
        ):
            raise ValueError("Conservative identity has conflicting routing/scaffold/fingerprint")
        previous = full_to_conservative.setdefault(
            record.full_identity, record.conservative_identity
        )
        if previous != record.conservative_identity:
            raise ValueError("Full identity has inconsistent conservative identity")
    return {key: nodes[key] for key in sorted(nodes)}


def _qualifies(left: Structure, right: Structure) -> bool:
    intersection, union = count_tanimoto(left.counts, right.counts)
    return 10 * intersection >= 7 * union


def build_structural_groups(structures: list[Structure]) -> list[Group]:
    nodes = _nodes(structures)
    by_identity: dict[str, list[int]] = defaultdict(list)
    by_scaffold: dict[str, list[int]] = defaultdict(list)
    for record in structures:
        by_identity[record.conservative_identity].append(record.record_number)
        if record.route == ROUTES[0]:
            by_scaffold[record.scaffold_key].append(record.record_number)
    groups = [
        Group("murcko:" + key, ROUTES[0], tuple(sorted(members)))
        for key, members in sorted(by_scaffold.items())
    ]
    acyclic = [record for record in nodes.values() if record.route == ROUTES[1]]
    edges = [
        (i, j)
        for i, j in combinations(range(len(acyclic)), 2)
        if _qualifies(acyclic[i], acyclic[j])
    ]
    for members in components(len(acyclic), edges):
        identities = sorted(acyclic[i].conservative_identity for i in members)
        numbers = sorted(number for identity in identities for number in by_identity[identity])
        groups.append(Group("acyclic:" + _json(identities), ROUTES[1], tuple(numbers)))
    return sorted(groups, key=lambda group: (-group.weight, group.key))


def allocation_cost(counts: tuple[int, ...], total: int) -> Fraction:
    return sum(
        (
            (count - fraction * total) ** 2 / (fraction * total)
            for count, fraction in zip(counts, FRACTIONS, strict=True)
        ),
        Fraction(0),
    )


def allocate_route(groups: list[Group]) -> dict[int, str]:
    """One fixed largest-first pass, with exact objective and nonempty reservation."""
    if len(groups) < 3 or len({group.route for group in groups}) != 1:
        raise ValueError("Each route requires at least three indivisible groups")
    if any(not group.weight for group in groups):
        raise ValueError("Empty allocation group")
    ordered = sorted(groups, key=lambda group: (-group.weight, group.key))
    numbers = [number for group in groups for number in group.record_numbers]
    if len(numbers) != len(set(numbers)) or len({group.key for group in groups}) != len(groups):
        raise ValueError("Allocation groups overlap or have duplicate keys")
    total, counts, assigned = len(numbers), [0, 0, 0], {}
    for index, group in enumerate(ordered):
        remaining = len(ordered) - index - 1
        candidates = []
        for split_index in range(3):
            following = counts.copy()
            following[split_index] += group.weight
            if following.count(0) <= remaining:
                candidates.append((allocation_cost(tuple(following), total), split_index))
        if not candidates:
            raise ValueError("Cannot populate all partitions without breaking groups")
        _, chosen = min(candidates)
        counts[chosen] += group.weight
        assigned.update(dict.fromkeys(group.record_numbers, SPLITS[chosen]))
    return assigned


def _group_map(groups: list[Group]) -> dict[int, Group]:
    return {number: group for group in groups for number in group.record_numbers}


def manifest_rows(
    structures: list[Structure], groups: list[Group], assignments: Mapping[int, str], digest: str
) -> list[dict[str, Any]]:
    membership = _group_map(groups)
    return [
        dict(
            zip(
                COLUMNS,
                (
                    _sample_id(digest, record.record_number),
                    record.record_number,
                    record.compound_id,
                    _key_id("full_identity", record.full_identity),
                    _key_id("conservative_identity", record.conservative_identity),
                    record.route,
                    _key_id("structural_group", membership[record.record_number].key),
                    assignments[record.record_number],
                ),
                strict=True,
            )
        )
        for record in sorted(structures, key=lambda item: item.record_number)
    ]


def _population_report(
    structures: list[Structure], assignments: Mapping[int, str]
) -> dict[str, Any]:
    total = len(structures)
    counts = Counter(assignments.get(record.record_number) for record in structures)
    splits = {}
    for name, fraction in zip(SPLITS, FRACTIONS, strict=True):
        actual = Fraction(counts[name], total) if total else None
        deviation = actual - fraction if actual is not None else None
        splits[name] = {
            "record_count": counts[name],
            "percentage": float(100 * actual) if actual is not None else None,
            "deviation_percentage_points": float(100 * deviation)
            if deviation is not None
            else None,
            "within_tolerance": deviation is not None and abs(deviation) <= TOLERANCE,
            "identity_count": len(
                {
                    item.conservative_identity
                    for item in structures
                    if assignments.get(item.record_number) == name
                }
            ),
        }
    return {"record_count": total, "splits": splits}


def validate_assignments(
    structures: list[Structure], groups: list[Group], rows: list[dict[str, Any]], digest: str
) -> dict[str, Any]:
    """Re-enumerate constraints from structures, not the graph used for construction."""
    expected = {record.record_number for record in structures}
    observed = [row.get("source_record_number") for row in rows]
    assigned = {
        row["source_record_number"]: row.get("split")
        for row in rows
        if type(row.get("source_record_number")) is int
    }
    group_numbers = [number for group in groups for number in group.record_numbers]
    membership = _group_map(groups)
    checks = {
        "source_record_numbers": expected == set(range(1, len(structures) + 1)),
        "record_conservation": Counter(observed) == Counter(expected),
        "record_number_types": all(type(number) is int for number in observed),
        "manifest_schema": all(set(row) == set(COLUMNS) for row in rows),
        "partition_values": all(row.get("split") in SPLITS for row in rows),
        "sample_ids": (
            len({row.get("sample_id") for row in rows}) == len(rows)
            and all(
                row.get("sample_id") == _sample_id(digest, row.get("source_record_number"))
                for row in rows
                if type(row.get("source_record_number")) is int
            )
        ),
        "group_conservation": Counter(group_numbers) == Counter(expected),
        "group_keys_unique": len({group.key for group in groups}) == len(groups),
    }
    fields_match = True
    for record in structures:
        matches = [row for row in rows if row.get("source_record_number") == record.record_number]
        group = membership.get(record.record_number)
        if len(matches) != 1 or group is None:
            fields_match = False
            continue
        row = matches[0]
        fields_match &= (
            row.get("compound_id") == record.compound_id
            and row.get("full_identity_id") == _key_id("full_identity", record.full_identity)
            and row.get("conservative_identity_id")
            == _key_id("conservative_identity", record.conservative_identity)
            and row.get("structure_route") == record.route == group.route
            and row.get("structural_group_id") == _key_id("structural_group", group.key)
        )
    checks["manifest_fields"] = bool(fields_match)

    def intact(attribute: str, nonempty_only: bool = False) -> bool:
        partitions: dict[str, set[str | None]] = defaultdict(set)
        for record in structures:
            if not nonempty_only or record.route == ROUTES[0]:
                partitions[getattr(record, attribute)].add(assigned.get(record.record_number))
        return all(len(values) == 1 for values in partitions.values())

    checks["full_identity_integrity"] = intact("full_identity")
    checks["conservative_identity_integrity"] = intact("conservative_identity")
    checks["nonempty_scaffold_integrity"] = intact("scaffold_key", True)
    checks["group_integrity"] = all(
        len({assigned.get(number) for number in group.record_numbers}) == 1 for group in groups
    )
    nodes = list(_nodes(structures).values())
    acyclic = [node for node in nodes if node.route == ROUTES[1]]
    violations = []
    for left, right in combinations(acyclic, 2):
        if _qualifies(left, right) and assigned.get(left.record_number) != assigned.get(
            right.record_number
        ):
            violations.append(
                [
                    _key_id("conservative_identity", left.conservative_identity),
                    _key_id("conservative_identity", right.conservative_identity),
                ]
            )
    checks["acyclic_edge_integrity"] = not violations
    populations = {"all": _population_report(structures, assigned)}
    populations.update(
        {
            route: _population_report(
                [record for record in structures if record.route == route], assigned
            )
            for route in ROUTES
        }
    )
    checks["all_partitions_nonempty"] = all(
        populations["all"]["splits"][split]["record_count"] > 0 for split in SPLITS
    )
    checks["both_routes_in_all_partitions"] = all(
        populations[route]["splits"][split]["record_count"] > 0
        for route in ROUTES
        for split in SPLITS
    )
    checks["proportion_tolerance"] = all(
        population["splits"][split]["within_tolerance"]
        for population in populations.values()
        for split in SPLITS
    )
    return {
        "checks": checks,
        "populations": populations,
        "acyclic_cross_partition_edge_violations": violations,
        "hard_checks_passed": all(
            value for key, value in checks.items() if key != "proportion_tolerance"
        ),
        "passed": all(checks.values()),
    }


def _similarity_summary(values: list[Fraction]) -> dict[str, Any]:
    ordered = sorted(values)

    def quantile(fraction: Fraction) -> float | None:
        if not ordered:
            return None
        position = (len(ordered) - 1) * fraction
        lower = position.numerator // position.denominator
        upper = min(lower + 1, len(ordered) - 1)
        return float(ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower]))

    bounds = (
        Fraction(0),
        Fraction(1, 5),
        Fraction(2, 5),
        Fraction(3, 5),
        THRESHOLD,
        Fraction(4, 5),
        Fraction(1),
    )
    bins = []
    for index, (lower, upper) in enumerate(zip(bounds, bounds[1:])):
        final = index == len(bounds) - 2
        bins.append(
            {
                "lower": float(lower),
                "upper": float(upper),
                "upper_inclusive": final,
                "count": sum(
                    lower <= value and (value <= upper if final else value < upper)
                    for value in ordered
                ),
            }
        )
    return {
        "identity_count": len(values),
        "maximum": float(max(values)) if values else None,
        "quantiles": {
            str(float(q)): quantile(q)
            for q in (
                Fraction(0),
                Fraction(1, 4),
                Fraction(1, 2),
                Fraction(3, 4),
                Fraction(19, 20),
                Fraction(1),
            )
        },
        "at_or_above_0_70": sum(value >= THRESHOLD for value in values),
        "bins": bins,
    }


def audit_structural_similarity(
    structures: list[Structure], assignments: Mapping[int, str]
) -> dict[str, Any]:
    """Stream held-out versus training comparisons; cross-route scores are descriptive."""
    nodes = _nodes(structures)
    training = [node for node in nodes.values() if assignments[node.record_number] == "train"]
    details = []
    summaries: dict[tuple[str, str, str], list[Fraction]] = {
        (split, route, comparison): []
        for split in SPLITS[1:]
        for route in ROUTES
        for comparison in ("all_training", "same_route_training", "cross_route_training")
    }
    for identity, heldout in nodes.items():
        split = assignments[heldout.record_number]
        if split == "train":
            continue
        best: dict[str, tuple[Fraction, str]] = {}
        for neighbor in training:
            similarity = Fraction(*count_tanimoto(heldout.counts, neighbor.counts))
            comparison = (
                "same_route_training" if heldout.route == neighbor.route else "cross_route_training"
            )
            for category in ("all_training", comparison):
                previous = best.get(category)
                if previous is None or similarity > previous[0]:
                    best[category] = (similarity, neighbor.conservative_identity)
        neighbors = {}
        for comparison in ("all_training", "same_route_training", "cross_route_training"):
            match = best.get(comparison)
            neighbors[comparison] = (
                {
                    "similarity": float(match[0]),
                    "training_identity_id": _key_id("conservative_identity", match[1]),
                }
                if match is not None
                else None
            )
            if match is not None:
                summaries[(split, heldout.route, comparison)].append(match[0])
        details.append(
            {
                "identity_id": _key_id("conservative_identity", identity),
                "split": split,
                "route": heldout.route,
                "nearest_neighbors": neighbors,
            }
        )
    return {
        "denominator": "unique conservative held-out identities; not source records",
        "quantiles": "linear interpolation at (n-1)*q; exact fractions before serialization",
        "scope": (
            "Only empty-to-empty separation is threshold-enforced; other scores describe coverage"
        ),
        "summaries": [
            {
                "split": split,
                "route": route,
                "comparison": comparison,
                **_similarity_summary(values),
            }
            for (split, route, comparison), values in sorted(summaries.items())
        ],
        "heldout_identities": details,
    }


def _group_diagnostics(structures: list[Structure], groups: list[Group]) -> dict[str, Any]:
    by_number = {record.record_number: record for record in structures}
    distributions = {}
    for route in ROUTES:
        selected = [group for group in groups if group.route == route]
        distributions[route] = {
            "group_count": len(selected),
            "record_size_histogram": dict(
                sorted(Counter(group.weight for group in selected).items())
            ),
            "identity_size_histogram": dict(
                sorted(
                    Counter(
                        len(
                            {
                                by_number[number].conservative_identity
                                for number in group.record_numbers
                            }
                        )
                        for group in selected
                    ).items()
                )
            ),
        }
    largest = []
    selected = sorted(
        (group for group in groups if group.route == ROUTES[1]),
        key=lambda group: (-group.weight, group.key),
    )[:5]
    for group in selected:
        nodes = list(_nodes([by_number[number] for number in group.record_numbers]).values())
        pairs = [
            (
                Fraction(*count_tanimoto(left.counts, right.counts)),
                left.conservative_identity,
                right.conservative_identity,
            )
            for left, right in combinations(nodes, 2)
        ]
        minimum = min(pairs) if pairs else None
        largest.append(
            {
                "group_id": _key_id("structural_group", group.key),
                "record_count": group.weight,
                "identity_count": len(nodes),
                "minimum_pairwise_similarity": float(minimum[0]) if minimum else None,
                "minimum_pair_identity_ids": [
                    _key_id("conservative_identity", key) for key in minimum[1:]
                ]
                if minimum
                else [],
                "chaining_below_threshold": minimum is not None and minimum[0] < THRESHOLD,
            }
        )
    return {
        "size_distributions": distributions,
        "largest_five_acyclic_components": largest,
        "singleton_definition": "one conservative identity; may have multiple source records",
    }


def _configuration() -> dict[str, Any]:
    return {
        "protocol": PROTOCOL,
        "rdkit": RDKIT_VERSION,
        "columns_used": ["Compound ID", "smiles", "source record position"],
        "parser": {
            "sanitize": True,
            "parseName": False,
            "allowCXSMILES": False,
            "other_defaults": f"RDKit {RDKIT_VERSION}",
        },
        "identity": {
            "full": "canonical isomeric SMILES; atom maps removed on key copies",
            "conservative": "stereochemistry removed; canonical isomeric SMILES retains isotopes",
            "standardization": "none beyond pinned RDKit parsing/sanitization",
        },
        "scaffold": "GetScaffoldForMol; canonical=True,isomericSmiles=False; not generic",
        "fingerprint": {"output": "GetSparseCountFingerprint", "parameters": dict(FP_PARAMETERS)},
        "tanimoto": "sum(min counts)/sum(max counts); edge iff 10*intersection >= 7*union",
        "threshold": "7/10",
        "fractions": ["4/5", "1/10", "1/10"],
        "tolerance": "1/50 absolute fraction, inclusive; global and per route",
        "allocation": {
            "routes": list(ROUTES),
            "order": "decreasing record count, then complete key lexical",
            "objective": "sum((n_s-q_s*N)^2/(q_s*N)); exact rational arithmetic",
            "ties": list(SPLITS),
            "reservation": "remaining groups >= still-empty partitions",
            "passes": 1,
            "seed": None,
        },
        "key_serialization": {
            "schema": 1,
            "encoding": "UTF-8",
            "sort_keys": True,
            "ensure_ascii": False,
            "separators": [",", ":"],
            "digest": "SHA-256",
            "scaffold_key": "murcko: + scaffold SMILES",
            "component_key": "acyclic: + compact JSON of sorted conservative identities",
        },
        "manifest": {
            "columns": list(COLUMNS),
            "encoding": "UTF-8",
            "line_ending": "LF",
            "order": "source_record_number",
        },
        "device": "cpu",
        "learned_preprocessing": "none",
        "model": None,
        "test_policy": {
            "seal": "procedural; public data remain accessible",
            "owner_artifact_acceptance_required": True,
            "refit_train_plus_validation": False,
            "validation": "model/configuration selection and approved calibration decisions",
            "test": "final evaluation only after model and evaluation plan are locked",
        },
    }


def _key_mappings(structures: list[Structure], groups: list[Group]) -> dict[str, Any]:
    mappings = {
        "full_identity": sorted({record.full_identity for record in structures}),
        "conservative_identity": sorted({record.conservative_identity for record in structures}),
        "structural_group": sorted({group.key for group in groups}),
    }
    result = {}
    all_ids = set()
    for kind, keys in mappings.items():
        result[kind] = {}
        for key in keys:
            identifier = _key_id(kind, key)
            if identifier in all_ids:
                raise ValueError("Identity/group digest collision")
            all_ids.add(identifier)
            result[kind][identifier] = key
    return result


def _provenance() -> dict[str, Any]:
    files = [
        Path(__file__),
        Path(__file__).with_name("acyclic_diagnostic.py"),
        ROOT / "pyproject.toml",
        ROOT / "uv.lock",
    ]

    def git(*args: str) -> str | None:
        try:
            result = subprocess.run(
                ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False, timeout=5
            )
        except (OSError, subprocess.TimeoutExpired):
            return None
        return result.stdout.strip() if result.returncode == 0 else None

    return {
        "python": platform.python_version(),
        "rdkit": rdBase.rdkitVersion,
        "system": platform.system(),
        "machine": platform.machine(),
        "device": "cpu",
        "source_hashes": {
            path.relative_to(ROOT).as_posix(): _sha(path.read_bytes()) for path in files
        },
        "git_commit": git("rev-parse", "--verify", "HEAD"),
        "git_status": git("status", "--short"),
        "git_note": "A missing commit is recorded as null; source hashes identify local code",
    }


def build_artifacts(
    structures: list[Structure],
    groups: list[Group],
    assignments: Mapping[int, str],
    source: dict[str, Any],
    execution: dict[str, Any],
) -> dict[str, bytes]:
    mappings = _key_mappings(structures, groups)
    rows = manifest_rows(structures, groups, assignments, source["sha256"])
    report = validate_assignments(structures, groups, rows, source["sha256"])
    report["checks"]["source_checksum_unchanged"] = source.get("checksum_unchanged") is True
    report["hard_checks_passed"] &= report["checks"]["source_checksum_unchanged"]
    report["passed"] &= report["checks"]["source_checksum_unchanged"]
    report["group_diagnostics"] = _group_diagnostics(structures, groups)
    status = "validated_candidate" if report["passed"] else "rejected_candidate"
    report.update(
        {
            "protocol": PROTOCOL,
            "status": status,
            "seal_status": "not_sealed",
            "source": source,
            "target_isolation": {
                "allowed_fields": ["record_number", "Compound ID", "smiles"],
                "construction": "no target/prediction values passed to structural functions",
                "assurance_limit": (
                    "mechanical isolation does not prove absence of human label knowledge"
                ),
            },
        }
    )
    if not report["hard_checks_passed"]:
        return {
            "acceptance": (json.dumps(report, indent=2, allow_nan=False) + "\n").encode("utf-8")
        }
    report["nearest_training_similarity"] = audit_structural_similarity(structures, assignments)
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=COLUMNS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    manifest = stream.getvalue().encode("utf-8")
    report["manifest_sha256"] = _sha(manifest)
    acceptance = (json.dumps(report, indent=2, allow_nan=False) + "\n").encode("utf-8")
    configuration = _configuration()
    metadata = {
        "protocol": PROTOCOL,
        "status": status,
        "owner_artifact_acceptance": "pending at generation; not automatically frozen or sealed",
        "source": source,
        "configuration": configuration,
        "configuration_sha256": _sha(_json(configuration).encode("utf-8")),
        "key_mappings": mappings,
        "groups": [
            {
                "group_id": _key_id("structural_group", group.key),
                "route": group.route,
                "record_count": group.weight,
                "record_numbers": list(group.record_numbers),
            }
            for group in sorted(groups, key=lambda group: (-group.weight, group.key))
        ],
        "execution": execution,
        "artifact_sha256": {"manifest": _sha(manifest), "acceptance": _sha(acceptance)},
        "evidence_preservation": "Historical 0.60/0.70/0.80 diagnostic artifacts are not rewritten",
    }
    return {
        "manifest": manifest,
        "metadata": (json.dumps(metadata, indent=2, allow_nan=False) + "\n").encode("utf-8"),
        "acceptance": acceptance,
    }


def _output_paths(prefix: Path) -> dict[str, Path]:
    return {
        "manifest": prefix.parent / f"{prefix.name}_manifest.csv",
        "metadata": prefix.parent / f"{prefix.name}_metadata.json",
        "acceptance": prefix.parent / f"{prefix.name}_acceptance.json",
    }


def _preflight_outputs(prefix: Path, source: Path) -> dict[str, Path]:
    paths = _output_paths(prefix)
    if not prefix.name or not prefix.parent.is_dir():
        raise ValueError("Output prefix requires an existing parent directory")
    for path in paths.values():
        resolved = path.resolve()
        if resolved.is_relative_to(source.resolve().parent) or resolved.is_relative_to(
            (ROOT / "data" / "raw").resolve()
        ):
            raise ValueError("Outputs must be outside the raw/input data directory")
        if path.exists() or path.is_symlink():
            raise ValueError(f"Output already exists; refusing overwrite: {path}")
    return paths


def write_artifacts(prefix: Path, payloads: dict[str, bytes], source: Path) -> None:
    """Exclusive publication, acceptance last; incomplete bundles are never accepted."""
    paths = _preflight_outputs(prefix, source)
    staged: dict[str, Path] = {}
    created = []
    try:
        for kind, payload in payloads.items():
            with tempfile.NamedTemporaryFile(
                prefix=f".{prefix.name}-", suffix=".tmp", dir=prefix.parent, delete=False
            ) as stream:
                staged[kind] = Path(stream.name)
                stream.write(payload)
                stream.flush()
                os.fsync(stream.fileno())
        for kind in ("manifest", "metadata", "acceptance"):
            if kind in staged:
                os.link(staged[kind], paths[kind])  # Exclusive: never overwrites a racing writer.
                created.append(paths[kind])
    except OSError:
        for path in reversed(created):
            path.unlink()
        raise
    finally:
        for path in staged.values():
            path.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-record-count", type=int, required=True)
    parser.add_argument("--expected-acyclic-count", type=int, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    timestamp = datetime.now(UTC).isoformat()
    try:
        _preflight_outputs(args.output_prefix, args.path)
        if rdBase.rdkitVersion != RDKIT_VERSION:
            raise ValueError(f"This protocol requires RDKit {RDKIT_VERSION}")
        if args.expected_record_count <= 0 or not 0 < args.expected_acyclic_count < (
            args.expected_record_count
        ):
            raise ValueError("Expected counts must describe two nonempty routes")
        raw = read_structural_records(args.path, args.expected_source_sha256)
        if len(raw) != args.expected_record_count:
            raise ValueError("Source record count differs from the authorized population")
        structures = prepare_structures(raw)
        if sum(record.route == ROUTES[1] for record in structures) != args.expected_acyclic_count:
            raise ValueError("Empty-scaffold count differs from the authorized population")
        groups = build_structural_groups(structures)
        assignments = {}
        for route in ROUTES:
            assignments.update(allocate_route([group for group in groups if group.route == route]))
        execution = _provenance()
        execution.update(
            {
                "started_utc": timestamp,
                "command_arguments": {
                    "input": args.path.as_posix(),
                    "output_prefix": args.output_prefix.as_posix(),
                    "expected_source_sha256": args.expected_source_sha256,
                    "expected_record_count": args.expected_record_count,
                    "expected_acyclic_count": args.expected_acyclic_count,
                },
            }
        )
        source = {
            "path": args.path.as_posix(),
            "sha256": args.expected_source_sha256,
            "record_count": len(structures),
            "checksum_unchanged": True,
        }
        payloads = build_artifacts(structures, groups, assignments, source, execution)
        if _sha(args.path.read_bytes()) != args.expected_source_sha256:
            raise ValueError("Source changed during construction; no outputs published")
        finished = datetime.now(UTC).isoformat()
        duration = time.perf_counter() - started
        if "metadata" in payloads:
            metadata = json.loads(payloads["metadata"])
            metadata["execution"].update(
                {
                    "checks_finished_utc": finished,
                    "construction_and_checks_seconds": duration,
                    "duration_scope": "through final source check; excludes publication",
                }
            )
            payloads["metadata"] = (json.dumps(metadata, indent=2, allow_nan=False) + "\n").encode(
                "utf-8"
            )
        write_artifacts(args.output_prefix, payloads, args.path)
        report = json.loads(payloads["acceptance"])
    except (OSError, ValueError, UnicodeError) as error:
        parser.exit(2, f"Hybrid split could not complete: {error}\n")
    print(f"{report['status']}; not sealed; owner artifact acceptance required.")
    for kind, payload in payloads.items():
        print(f"{_output_paths(args.output_prefix)[kind]} SHA-256 {_sha(payload)}")
    if not report["passed"]:
        print("Stopped after the fixed allocation; no alternative orders or thresholds attempted.")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
