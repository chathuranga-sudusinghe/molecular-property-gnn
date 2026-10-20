"""Label-blind structural diagnostics; no allocation or empty-scaffold policy is applied."""

import argparse
import csv
import hashlib
import io
import json
import platform
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold


def _reference(record: dict[str, Any]) -> dict[str, Any]:
    return {key: record[key] for key in ("record_number", "compound_id", "smiles")}


def _duplicates(records: list[dict[str, Any]], key: str) -> dict[str, Any]:
    groups = defaultdict(list)
    for record in records:
        groups[record[key]].append(record)
    repeated = sorted(
        ((identity, members) for identity, members in groups.items() if len(members) > 1),
        key=lambda item: (-len(item[1]), item[0]),
    )
    return {
        "group_count": len(repeated),
        "extra_records": sum(len(members) - 1 for _, members in repeated),
        "affected_records": sum(len(members) for _, members in repeated),
        "groups": [
            {
                "identity": identity,
                "record_count": len(members),
                "records": [_reference(member) for member in members],
            }
            for identity, members in repeated
        ],
    }


def audit_structures(records: Iterable[Mapping[str, str]]) -> dict[str, Any]:
    """Access only Compound ID and smiles; preserve every observation and its position.

    Canonical identities use sanitized, as-supplied graphs. Atom maps are removed
    from key copies only; no fragments, charges, tautomers, or records are normalized.
    Empty scaffolds are counted, not assigned a final grouping policy or split.
    """
    parameters = Chem.SmilesParserParams()
    parameters.parseName = False
    parameters.allowCXSMILES = False
    valid, invalid, disconnected, empty = [], [], [], []
    categories = defaultdict(list)
    mapped_count = 0
    total = 0
    for total, source in enumerate(records, start=1):
        smiles = source["smiles"]
        reference = {"record_number": total, "compound_id": source["Compound ID"], "smiles": smiles}
        molecule = Chem.MolFromSmiles(smiles.strip(), parameters) if smiles.strip() else None
        if molecule is None or molecule.GetNumAtoms() == 0:
            invalid.append({**reference, "reason": "empty, unparsable, or unsanitizable molecule"})
            continue
        molecule = Chem.Mol(molecule)
        mapped_count += int(any(atom.GetAtomMapNum() for atom in molecule.GetAtoms()))
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
        identity = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        stereo_free = Chem.Mol(molecule)
        Chem.RemoveStereochemistry(stereo_free)
        conservative = Chem.MolToSmiles(stereo_free, canonical=True, isomericSmiles=True)
        scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
        scaffold_key = Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=False)
        fragment_count = len(Chem.GetMolFrags(molecule))
        record = {
            **reference,
            "canonical_identity": identity,
            "stereo_free_identity": conservative,
            "scaffold_key": scaffold_key,
            "fragment_count": fragment_count,
        }
        valid.append(record)
        if fragment_count > 1:
            disconnected.append(
                {
                    **_reference(record),
                    "canonical_identity": identity,
                    "fragment_count": fragment_count,
                }
            )
        if scaffold.GetNumAtoms() == 0:
            empty.append(record)
            elements = {atom.GetAtomicNum() for atom in molecule.GetAtoms()}
            if 6 not in elements:
                categories["carbon_free"].append(record)
            elif elements <= {1, 6}:
                categories["hydrocarbon_only"].append(record)
            else:
                categories["carbon_and_heteroatoms"].append(record)
            if any(atom.GetFormalCharge() != 0 for atom in molecule.GetAtoms()):
                categories["contains_formally_charged_atoms"].append(record)
            if fragment_count > 1:
                categories["disconnected"].append(record)
            if any(
                bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)
                for bond in molecule.GetBonds()
            ):
                categories["contains_double_or_triple_bonds"].append(record)

    scaffold_sizes = Counter(record["scaffold_key"] for record in valid if record["scaffold_key"])
    sizes = sorted(scaffold_sizes.values(), reverse=True)
    valid_count, empty_count = len(valid), len(empty)
    unique_count = len({record["canonical_identity"] for record in valid})
    unique_empty = len({record["canonical_identity"] for record in empty})
    histogram = {}
    for lower, upper, label in [
        (1, 1, "1"),
        (2, 5, "2-5"),
        (6, 10, "6-10"),
        (11, 25, "11-25"),
        (26, 50, "26-50"),
        (51, None, "51+"),
    ]:
        selected = [size for size in sizes if size >= lower and (upper is None or size <= upper)]
        histogram[label] = {"group_count": len(selected), "record_count": sum(selected)}
    return {
        "record_count": total,
        "valid_record_count": valid_count,
        "unique_valid_isomeric_identities": unique_count,
        "invalid_records": invalid,
        "atom_mapped_record_count": mapped_count,
        "canonical_duplicates": _duplicates(valid, "canonical_identity"),
        "stereo_free_duplicates": _duplicates(valid, "stereo_free_identity"),
        "empty_scaffolds": {
            "record_count": empty_count,
            "percentage_of_valid_records": 100 * empty_count / valid_count if valid_count else None,
            "unique_isomeric_identities": unique_empty,
            "percentage_of_unique_isomeric_identities": (
                100 * unique_empty / unique_count if unique_count else None
            ),
            "categories_are_nonexclusive": True,
            "categories": {
                name: {
                    "record_count": len(members),
                    "examples": [_reference(record) for record in members[:3]],
                }
                for name, members in sorted(categories.items())
            },
        },
        "nonempty_scaffold_groups": {
            "group_count": len(sizes),
            "record_count": sum(sizes),
            "sizes_descending": sizes,
            "size_histogram": histogram,
            "largest_groups": [
                {"scaffold_key": key, "record_count": count}
                for key, count in sorted(
                    scaffold_sizes.items(), key=lambda item: (-item[1], item[0])
                )[:10]
            ],
        },
        "multicomponent_records": disconnected,
        "multicomponent_record_count": len(disconnected),
        "hypothetical_single_empty_bucket_capacity": {
            "note": "Count-only scenario, not a split or accepted empty-scaffold policy. "
            "Denominator is valid records, not a finalized eligibility population.",
            "partitions": {
                name: {
                    "ideal_record_budget": valid_count * fraction,
                    "empty_bucket_fits_budget": empty_count <= valid_count * fraction,
                    "budget_remaining_after_bucket": valid_count * fraction - empty_count,
                }
                for name, fraction in (("train", 0.8), ("validation", 0.1), ("test", 0.1))
            }
            if valid_count
            else {},
        },
        "policy": {
            "status": "structural audit only; empty-scaffold policy remains open",
            "columns_used": ["Compound ID", "smiles"],
            "scaffold": "RDKit typed Murcko; canonical non-isomeric scaffold SMILES",
            "identity": "canonical isomeric SMILES; stereo-free key retains isotopes",
            "parser": "default sanitization; parseName=False; allowCXSMILES=False",
            "atom_maps": "ignored in key copies only",
            "chemistry_standardization": "none beyond RDKit parsing/sanitization",
            "allocation": "none; no assignments or seed/split search",
        },
    }


def audit_structural_csv(path: Path) -> dict[str, Any]:
    """Read one snapshot; discard all fields except ID/SMILES before structural analysis."""
    content = path.read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    if b"\x00" in content:
        raise ValueError("CSV contains a NUL byte")
    reader = csv.reader(io.StringIO(content.decode("utf-8-sig"), newline=""), strict=True)
    try:
        header = next(reader, [])
        if len(header) != len(set(header)) or not {"Compound ID", "smiles"} <= set(header):
            raise ValueError("CSV needs unique headers including Compound ID and smiles")
        identifier_index, smiles_index = header.index("Compound ID"), header.index("smiles")
        records = []
        for number, row in enumerate(reader, start=1):
            if len(row) != len(header):
                raise ValueError(f"Record {number} has an unexpected field count")
            records.append({"Compound ID": row[identifier_index], "smiles": row[smiles_index]})
    except csv.Error as error:
        raise ValueError(f"Malformed CSV: {error}") from error
    report = audit_structures(records)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError("Source file changed during audit; report not saved")
    report["source"] = {"path": path.as_posix(), "sha256": digest, "checksum_unchanged": True}
    report["versions"] = {"python": platform.python_version(), "rdkit": rdBase.rdkitVersion}
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.output.resolve().is_relative_to(args.path.resolve().parent):
            raise ValueError("Report must be outside the input data directory")
        if args.output.exists():
            raise ValueError("Output already exists; refusing to overwrite it")
        if not args.output.parent.is_dir():
            raise ValueError("Output directory must already exist")
        report = audit_structural_csv(args.path)
        payload = json.dumps(report, indent=2, allow_nan=False) + "\n"
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(payload)
    except (OSError, ValueError) as error:
        parser.exit(2, f"Structural audit could not complete: {error}\n")
    print(f"Structural report saved to {args.output}; no split created.")
    return 1 if report["invalid_records"] or not report["valid_record_count"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
