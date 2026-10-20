"""Predeclared, label-blind acyclic similarity diagnostics; never assigns sample splits."""

import argparse
import csv
import hashlib
import io
import json
import math
import platform
from collections import Counter
from collections.abc import Iterable, Mapping
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

THRESHOLDS = (Fraction(3, 5), Fraction(7, 10), Fraction(4, 5))
FRACTIONS = (Fraction(4, 5), Fraction(1, 10), Fraction(1, 10))
FP_PARAMETERS = {
    "radius": 2,
    "countSimulation": False,
    "includeChirality": False,
    "useBondTypes": True,
    "onlyNonzeroInvariants": False,
    "includeRingMembership": True,
    "countBounds": None,
    "fpSize": 2048,  # Unused by GetSparseCountFingerprint; no fixed-size folding.
    "atomInvariantsGenerator": None,
    "bondInvariantsGenerator": None,
    "includeRedundantEnvironments": False,
}


def prepare_nodes(records: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Only route all structures; canonicalize/fingerprint zero-atom-scaffold records only."""
    parameters = Chem.SmilesParserParams()
    parameters.parseName = False
    parameters.allowCXSMILES = False
    generator = rdFingerprintGenerator.GetMorganGenerator(**FP_PARAMETERS)
    nodes = {}
    for record in records:
        smiles = record["smiles"]
        molecule = Chem.MolFromSmiles(smiles.strip(), parameters) if smiles.strip() else None
        if molecule is None or molecule.GetNumAtoms() == 0:
            raise ValueError(
                f"Invalid molecule at record {record['record_number']}; no exclusion applied"
            )
        if MurckoScaffold.GetScaffoldForMol(molecule).GetNumAtoms():
            continue
        if len(Chem.GetMolFrags(molecule)) != 1:
            raise ValueError("Unexpected multicomponent acyclic record; owner policy required")
        molecule = Chem.Mol(molecule)
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
        full_identity = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        Chem.RemoveStereochemistry(molecule)
        identity = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        if identity not in nodes:
            counts = dict(generator.GetSparseCountFingerprint(molecule).GetNonzeroElements())
            if not counts or any(count <= 0 for count in counts.values()):
                raise ValueError("Unexpected empty or nonpositive count fingerprint")
            nodes[identity] = {
                "identity": identity,
                "counts": counts,
                "records": [],
                "full_identities": set(),
            }
        nodes[identity]["full_identities"].add(full_identity)
        nodes[identity]["records"].append(
            {key: record[key] for key in ("record_number", "Compound ID", "smiles")}
        )
    for node in nodes.values():
        node["records"].sort(key=lambda record: record["record_number"])
    return [nodes[key] for key in sorted(nodes)]


def count_tanimoto(left: Mapping[int, int], right: Mapping[int, int]) -> tuple[int, int]:
    """Nonnegative sparse-count intersection and union; integer threshold decisions."""
    intersection = sum(min(value, right.get(key, 0)) for key, value in left.items())
    union = sum(left.values()) + sum(right.values()) - intersection
    if union <= 0:
        raise ValueError("Undefined similarity for empty count fingerprints")
    return intersection, union


def components(node_count: int, edges: list[tuple[int, int]]) -> list[list[int]]:
    """Return connected components with stable node ordering; no centroid heuristic."""
    neighbors = [set() for _ in range(node_count)]
    for left, right in edges:
        neighbors[left].add(right)
        neighbors[right].add(left)
    unseen, result = set(range(node_count)), []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        pending, members = [root], []
        while pending:
            node = pending.pop()
            members.append(node)
            for neighbor in sorted(neighbors[node] & unseen):
                unseen.remove(neighbor)
                pending.append(neighbor)
        result.append(sorted(members))
    return result


def count_feasible(weights: list[int], bounds: list[tuple[int, int]]) -> bool:
    """Boolean capacity DP: counts only, no sample IDs, backpointers, or assignments.

    Bounds are train/validation/test count ranges. Every whole group contributes
    its weight exactly once. Both holdouts and training must be nonempty.
    """
    if any(weight <= 0 for weight in weights):
        raise ValueError("Group record counts must be positive")
    if any(upper < max(1, lower) for lower, upper in bounds):
        return False
    (_, _), (_, max_valid), (_, max_test) = bounds
    states = {(0, 0)}
    for weight in weights:
        following = set(states)
        for valid, test in states:
            if valid + weight <= max_valid:
                following.add((valid + weight, test))
            if test + weight <= max_test:
                following.add((valid, test + weight))
        states = following
    total = sum(weights)
    return any(
        all(
            max(1, lower) <= count <= upper
            for count, (lower, upper) in zip(
                (total - valid - test, valid, test), bounds, strict=True
            )
        )
        for valid, test in states
    )


def _pair_evidence(nodes, pair, similarities):
    left, right = pair
    intersection, union = similarities[pair]
    return {
        "similarity": intersection / union,
        "intersection": intersection,
        "union": union,
        "left": {"identity": nodes[left]["identity"], "example": nodes[left]["records"][0]},
        "right": {"identity": nodes[right]["identity"], "example": nodes[right]["records"][0]},
    }


def diagnose_nodes(nodes: list[dict[str, Any]]) -> dict[str, Any]:
    """Evaluate the fixed threshold grid, preserving duplicate record weights."""
    if not nodes:
        raise ValueError("No acyclic identities to diagnose")
    similarities = {
        pair: count_tanimoto(nodes[pair[0]]["counts"], nodes[pair[1]]["counts"])
        for pair in combinations(range(len(nodes)), 2)
    }
    record_count = sum(len(node["records"]) for node in nodes)
    quotas = [record_count * fraction for fraction in FRACTIONS]
    rounded = [math.floor(quota) for quota in quotas]
    order = sorted(range(3), key=lambda index: (-(quotas[index] - rounded[index]), index))
    for index in order[: record_count - sum(rounded)]:
        rounded[index] += 1
    relaxed = [
        (
            math.ceil(record_count * (fraction - Fraction(1, 50))),
            math.floor(record_count * (fraction + Fraction(1, 50))),
        )
        for fraction in FRACTIONS
    ]
    reports = []
    for threshold in THRESHOLDS:
        edges = [
            pair
            for pair, (intersection, union) in similarities.items()
            if threshold.denominator * intersection >= threshold.numerator * union
        ]
        groups = components(len(nodes), edges)
        groups.sort(
            key=lambda members: (
                -sum(len(nodes[i]["records"]) for i in members),
                tuple(nodes[i]["identity"] for i in members),
            )
        )
        summaries = []
        for members in groups:
            keys = [nodes[i]["identity"] for i in members]
            summaries.append(
                {
                    "component_id": hashlib.sha256(json.dumps(keys).encode()).hexdigest(),
                    "identity_count": len(members),
                    "record_count": sum(len(nodes[i]["records"]) for i in members),
                    "record_numbers": sorted(
                        record["record_number"] for i in members for record in nodes[i]["records"]
                    ),
                }
            )
        weights = [summary["record_count"] for summary in summaries]
        largest_details = []
        for members, summary in zip(groups[:5], summaries[:5], strict=True):
            internal_pairs = list(combinations(members, 2))
            minimum_pair = (
                min(internal_pairs, key=lambda pair: (Fraction(*similarities[pair]), pair))
                if internal_pairs
                else None
            )
            evidence = _pair_evidence(nodes, minimum_pair, similarities) if minimum_pair else None
            largest_details.append(
                {
                    **summary,
                    "minimum_pairwise_evidence": evidence,
                    "chaining_below_threshold": bool(
                        evidence
                        and Fraction(evidence["intersection"], evidence["union"]) < threshold
                    ),
                }
            )
        near_pairs = sorted(
            edges, key=lambda pair: (Fraction(*similarities[pair]) - threshold, pair)
        )[:5]
        near_evidence = []
        for pair in near_pairs:
            evidence = _pair_evidence(nodes, pair, similarities)
            evidence["gap_above_threshold"] = evidence["similarity"] - float(threshold)
            evidence["within_0_05_above_threshold"] = Fraction(
                *similarities[pair]
            ) - threshold <= Fraction(1, 20)
            near_evidence.append(evidence)
        singleton_groups = sum(len(members) == 1 for members in groups)
        singleton_records = sum(
            summary["record_count"] for summary in summaries if summary["identity_count"] == 1
        )
        reports.append(
            {
                "threshold": float(threshold),
                "identity_count": len(nodes),
                "edge_count": len(edges),
                "component_count": len(groups),
                "singleton_component_count": singleton_groups,
                "singleton_component_percentage": 100 * singleton_groups / len(groups),
                "records_in_singleton_components": singleton_records,
                "record_size_histogram": dict(sorted(Counter(weights).items())),
                "identity_size_histogram": dict(sorted(Counter(map(len, groups)).items())),
                "largest_component_record_count": weights[0],
                "largest_component_percentage": 100 * weights[0] / record_count,
                "components_above_nominal_10_percent_budget": sum(
                    10 * weight > record_count for weight in weights
                ),
                "majority_component_over_50_percent": 2 * weights[0] > record_count,
                "largest_component_over_80_percent": 5 * weights[0] > 4 * record_count,
                "near_threshold_connected_pairs": near_evidence,
                "largest_five_components": largest_details,
                "components": summaries,
                "count_only_feasibility": {
                    "rounded_80_10_10_counts": rounded,
                    "exact_rounded_counts_possible": count_feasible(
                        weights, [(n, n) for n in rounded]
                    ),
                    "illustrative_two_percentage_point_ranges": relaxed,
                    "within_illustrative_ranges_possible": count_feasible(weights, relaxed),
                    "note": "Boolean group-size capacity checks only; no assignments retained. "
                    "Two-percentage-point ranges are diagnostic, not an accepted policy. "
                    "This checks the acyclic subset, not joint full-dataset allocation.",
                },
            }
        )
    return {
        "acyclic_record_count": record_count,
        "full_canonical_identity_count": len(
            set().union(*(node["full_identities"] for node in nodes))
        ),
        "conservative_identity_count": len(nodes),
        "pair_comparison_count": len(similarities),
        "nominal_record_budgets": [float(quota) for quota in quotas],
        "threshold_reports": reports,
        "definitions": {
            "singleton": "one conservative identity; may retain more than one original record",
            "singleton_percentage_denominator": "number of connected components",
            "giant_flag": "descriptive majority flag: over 50% of acyclic records",
            "minimum_similarity": "distinct identities in each of five largest components",
            "near_pairs": "five closest connected pairs at/above threshold; gaps reported",
        },
    }


def audit_csv(path: Path, expected_sha256: str, expected_count: int) -> dict[str, Any]:
    content = path.read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    if digest != expected_sha256:
        raise ValueError("Input checksum differs from the owner-authorized dataset revision")
    if b"\x00" in content:
        raise ValueError("NUL byte in CSV")
    reader = csv.reader(io.StringIO(content.decode("utf-8-sig"), newline=""), strict=True)
    try:
        header = next(reader, [])
        if len(header) != len(set(header)) or not {"Compound ID", "smiles"} <= set(header):
            raise ValueError("Expected unique headers including Compound ID and smiles")
        id_index, smiles_index = header.index("Compound ID"), header.index("smiles")
        records = []
        for number, row in enumerate(reader, start=1):
            if len(row) != len(header):
                raise ValueError(f"Unexpected field count at record {number}")
            records.append(
                {"record_number": number, "Compound ID": row[id_index], "smiles": row[smiles_index]}
            )
    except csv.Error as error:
        raise ValueError(f"Malformed CSV: {error}") from error
    nodes = prepare_nodes(records)
    if sum(len(node["records"]) for node in nodes) != expected_count:
        raise ValueError("Empty-scaffold record count differs from the authorized population")
    report = diagnose_nodes(nodes)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError("Source changed during diagnostic; no report saved")
    report["source"] = {"path": path.as_posix(), "sha256": digest, "checksum_unchanged": True}
    report["method"] = {
        "status": "diagnostic only; no method, threshold, or sample split frozen",
        "columns_used": ["Compound ID", "smiles"],
        "fingerprint_output": "GetSparseCountFingerprint",
        "parameters": FP_PARAMETERS,
        "similarity": "sum(min(counts))/sum(max(counts)); exact rational edge comparisons",
        "grouping": "connected components on conservative identities; original records retained",
        "thresholds": [float(threshold) for threshold in THRESHOLDS],
        "routing": "parse structures for Murcko routing; fingerprint only zero-atom scaffolds",
        "preprocessing": (
            "approved parsing/sanitization; atom maps and stereo ignored in key copies; "
            "no desalting, neutralization, tautomer normalization, or record removal"
        ),
    }
    report["versions"] = {
        "python": platform.python_version(),
        "rdkit": rdBase.rdkitVersion,
        "diagnostic_code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-acyclic-count", type=int, required=True)
    args = parser.parse_args()
    try:
        if rdBase.rdkitVersion != "2024.09.6":
            raise ValueError("This diagnostic requires the approved RDKit 2024.09.6")
        if args.output.resolve().is_relative_to(args.path.resolve().parent):
            raise ValueError("Output must be outside the source data directory")
        if args.output.exists() or not args.output.parent.is_dir():
            raise ValueError("Output must be new, with an existing parent directory")
        report = audit_csv(args.path, args.expected_source_sha256, args.expected_acyclic_count)
        payload = json.dumps(report, indent=2, allow_nan=False) + "\n"
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(payload)
    except (OSError, ValueError) as error:
        parser.exit(2, f"Acyclic diagnostic could not complete: {error}\n")
    print(f"Diagnostic saved to {args.output}; no sample split assignments created.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
