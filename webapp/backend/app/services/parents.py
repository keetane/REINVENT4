from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

from app.errors import AppError


@dataclass(frozen=True, slots=True)
class ValidatedParent:
    input_smiles: str
    canonical_smiles: str
    molblock: str
    atom_count: int
    molecular_formula: str


def validate_parent_smiles(
    smiles: str, *, max_atoms: int, model_key: str = "mol2mol_high"
) -> ValidatedParent:
    value = smiles.strip()
    if not value:
        raise _invalid("Enter a parent SMILES before continuing.", "empty")
    if ">" in value:
        raise _invalid(
            "Reaction SMILES are not supported. Enter one parent molecule.",
            "reaction_not_supported",
        )
    if model_key == "linkinvent":
        return _validate_linkinvent_input(value, max_atoms=max_atoms)
    if "|" in value:
        raise _invalid(
            "The '|' separator is only supported for LinkInvent inputs.",
            "unexpected_fragment_separator",
        )
    with rdBase.BlockLogs():
        molecule = Chem.MolFromSmiles(value)
    if molecule is None:
        raise _invalid(
            "The parent SMILES could not be parsed. Check the syntax and valence.",
            "parse_failed",
        )
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise _invalid(
            "Multiple fragments are not supported. Enter one connected parent molecule.",
            "multiple_fragments",
        )
    if model_key == "libinvent":
        _validate_attachment_points(
            molecule,
            source="LibInvent scaffold",
            exact_count=None,
        )
    atom_count = molecule.GetNumAtoms()
    if atom_count == 0:
        raise _invalid("The parent structure is empty.", "empty_structure")
    if atom_count > max_atoms:
        raise _invalid(
            f"The parent molecule exceeds the {max_atoms}-atom limit.",
            "atom_limit_exceeded",
        )

    return ValidatedParent(
        input_smiles=value,
        canonical_smiles=Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True),
        molblock=Chem.MolToMolBlock(molecule),
        atom_count=atom_count,
        molecular_formula=rdMolDescriptors.CalcMolFormula(molecule),
    )


def _validate_linkinvent_input(value: str, *, max_atoms: int) -> ValidatedParent:
    if "." in value:
        raise _invalid(
            "Separate LinkInvent warheads with '|' instead of '.'.",
            "linkinvent_pipe_required",
        )
    fragments = [fragment.strip() for fragment in value.split("|")]
    if len(fragments) != 2 or any(not fragment for fragment in fragments):
        raise _invalid(
            "LinkInvent requires exactly two warheads separated by '|'.",
            "linkinvent_warhead_count",
        )

    molecules: list[Chem.Mol] = []
    canonical_fragments: list[str] = []
    for index, fragment in enumerate(fragments, start=1):
        with rdBase.BlockLogs():
            molecule = Chem.MolFromSmiles(fragment)
        if molecule is None:
            raise _invalid(
                f"LinkInvent warhead {index} could not be parsed.",
                "parse_failed",
            )
        _validate_attachment_points(
            molecule,
            source=f"LinkInvent warhead {index}",
            exact_count=1,
        )
        molecules.append(molecule)
        canonical_fragments.append(
            Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
        )

    atom_count = sum(molecule.GetNumAtoms() for molecule in molecules)
    if atom_count > max_atoms:
        raise _invalid(
            f"The LinkInvent input exceeds the {max_atoms}-atom limit.",
            "atom_limit_exceeded",
        )

    combined = Chem.MolFromSmiles(".".join(canonical_fragments))
    assert combined is not None
    return ValidatedParent(
        input_smiles=value,
        canonical_smiles="|".join(canonical_fragments),
        molblock=Chem.MolToMolBlock(combined),
        atom_count=atom_count,
        molecular_formula=rdMolDescriptors.CalcMolFormula(combined),
    )


def _validate_attachment_points(
    molecule: Chem.Mol, *, source: str, exact_count: int | None
) -> None:
    attachment_points = [
        atom for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 0
    ]
    if exact_count is None and not attachment_points:
        raise _invalid(
            f"{source} requires at least one '*' attachment point.",
            "attachment_point_required",
        )
    if exact_count is not None and len(attachment_points) != exact_count:
        raise _invalid(
            f"{source} requires exactly {exact_count} '*' attachment point.",
            "attachment_point_count",
        )
    if any(atom.GetDegree() != 1 for atom in attachment_points):
        raise _invalid(
            f"Every '*' attachment point in {source} must have exactly one bond.",
            "invalid_attachment_point",
        )


def _invalid(message: str, reason: str) -> AppError:
    return AppError(
        status_code=422,
        code="invalid_parent_smiles",
        message=message,
        context={"reason": reason},
    )
