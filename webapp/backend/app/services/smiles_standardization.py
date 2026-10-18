"""Shared, versioned SMILES standardization for TL, QSAR and RL.

Generated molecules are never rewritten in their source output.  Callers use the
``standardized_smiles`` value only for model input/scoring and retain
``original_smiles`` for traceability.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Iterable

import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize


STANDARDIZATION_VERSION = "1.0.0"
STANDARDIZATION_CONFIG: dict[str, Any] = {
    "cleanup": True,
    "normalize": True,
    "fragment_policy": "largest_heavy_atom_then_mw_then_canonical",
    "neutralize": True,
    "retain_permanent_charge": True,
    "minimum_molecular_weight_exclusive": 200.0,
    "canonical_isomeric_smiles": True,
    "reject_cxsmiles_pipe": True,
    "reject_reaction_smiles": True,
}


@dataclass(frozen=True, slots=True)
class StandardizationResult:
    original_smiles: str
    standardized_smiles: str | None
    accepted: bool
    reason: str | None
    molecular_weight: float | None
    formal_charge: int | None
    fragment_extracted: bool = False
    neutralized: bool = False
    charged_structure_retained: bool = False
    history: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["history"] = list(self.history)
        return value


def standardization_metadata() -> dict[str, Any]:
    return {
        "version": STANDARDIZATION_VERSION,
        "rdkit_version": rdkit.__version__,
        "config": dict(STANDARDIZATION_CONFIG),
    }


def assert_standardization_compatible(metadata: dict[str, Any] | None) -> None:
    stored = (metadata or {}).get("smiles_standardization")
    if not isinstance(stored, dict):
        raise ValueError("model_smiles_standardization_metadata_missing")
    if stored.get("version") != STANDARDIZATION_VERSION:
        raise ValueError("model_smiles_standardization_version_mismatch")
    if stored.get("config") != STANDARDIZATION_CONFIG:
        raise ValueError("model_smiles_standardization_config_mismatch")


def _fragment_key(molecule: Chem.Mol) -> tuple[int, float, str]:
    return (
        molecule.GetNumHeavyAtoms(),
        float(Descriptors.MolWt(molecule)),
        Chem.MolToSmiles(molecule, isomericSmiles=True, canonical=True),
    )


def standardize_smiles(value: object) -> StandardizationResult:
    original = str(value or "").strip()
    if not original:
        return StandardizationResult(original, None, False, "empty_smiles", None, None)
    # A pipe can start CXSMILES annotations.  The application accepts normal
    # molecular SMILES here; LinkInvent attachment syntax is handled elsewhere.
    if "|" in original:
        return StandardizationResult(original, None, False, "unexpected_cxsmiles_separator", None, None)
    if ">" in original:
        return StandardizationResult(original, None, False, "reaction_smiles_not_supported", None, None)

    try:
        molecule = Chem.MolFromSmiles(original)
        if molecule is None:
            return StandardizationResult(original, None, False, "invalid_smiles", None, None)
        molecule = rdMolStandardize.Cleanup(molecule)
        fragments = list(Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True))
        if not fragments:
            return StandardizationResult(original, None, False, "no_molecular_fragment", None, None)
        molecule = max(fragments, key=_fragment_key)
        fragment_extracted = len(fragments) > 1
        before = Chem.MolToSmiles(molecule, isomericSmiles=True, canonical=True)
        before_charge = int(Chem.GetFormalCharge(molecule))
        molecule = rdMolStandardize.Uncharger(canonicalOrder=True).uncharge(molecule)
        after = Chem.MolToSmiles(molecule, isomericSmiles=True, canonical=True)
        formal_charge = int(Chem.GetFormalCharge(molecule))
        charged_structure_retained = any(atom.GetFormalCharge() != 0 for atom in molecule.GetAtoms())
        neutralized = before_charge != formal_charge or before != after
        molecular_weight = float(Descriptors.MolWt(molecule))
        history: list[str] = ["parsed", "cleanup_normalize"]
        if fragment_extracted:
            history.append("largest_fragment_selected")
        if neutralized:
            history.append("uncharged")
        if charged_structure_retained:
            history.append("permanent_charge_retained")
        history.append("canonical_isomeric_smiles")
        if molecular_weight <= float(STANDARDIZATION_CONFIG["minimum_molecular_weight_exclusive"]):
            history.append("excluded_low_molecular_weight")
            return StandardizationResult(
                original, after, False, "molecular_weight_not_above_200", molecular_weight,
                formal_charge, fragment_extracted, neutralized, charged_structure_retained, tuple(history)
            )
        return StandardizationResult(
            original, after, True, None, molecular_weight, formal_charge,
            fragment_extracted, neutralized, charged_structure_retained, tuple(history)
        )
    except Exception:
        return StandardizationResult(original, None, False, "standardization_failed", None, None)


def standardization_counts(results: Iterable[StandardizationResult]) -> dict[str, int]:
    rows = list(results)
    return {
        "input_count": len(rows),
        "standardized_count": sum(row.standardized_smiles is not None for row in rows),
        "fragment_extracted_count": sum(row.fragment_extracted for row in rows),
        "neutralized_count": sum(row.neutralized for row in rows),
        "charged_structure_retained_count": sum(row.charged_structure_retained for row in rows),
        "low_molecular_weight_excluded_count": sum(row.reason == "molecular_weight_not_above_200" for row in rows),
        "invalid_count": sum(row.reason not in {None, "molecular_weight_not_above_200"} for row in rows),
        "accepted_before_deduplication_count": sum(row.accepted for row in rows),
    }
