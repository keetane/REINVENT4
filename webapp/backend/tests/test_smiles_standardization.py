from __future__ import annotations

import pytest

from app.services.smiles_standardization import (
    STANDARDIZATION_CONFIG,
    assert_standardization_compatible,
    standardization_metadata,
    standardize_smiles,
)


IBUPROFEN = "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"


def test_salt_uses_largest_fragment_and_neutralizes_deterministically():
    result = standardize_smiles(f"[Na+].CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)[O-]")

    assert result.accepted
    assert result.fragment_extracted
    assert result.neutralized
    assert result.formal_charge == 0
    assert result.standardized_smiles == standardize_smiles(IBUPROFEN).standardized_smiles


@pytest.mark.parametrize(
    "counterion",
    ["CS(=O)(=O)[O-]", "Cc1ccc(S(=O)(=O)[O-])cc1"],
)
def test_protonated_tertiary_amine_salts_are_desalted_and_neutralized(counterion):
    result = standardize_smiles(f"C[NH+](C)CCOC(c1ccccc1)c1ccccc1.{counterion}")
    assert result.accepted
    assert result.fragment_extracted
    assert result.neutralized
    assert result.formal_charge == 0


def test_protonated_pyridinium_is_neutralized_when_possible():
    result = standardize_smiles("CCCCCCCCCc1cc[nH+]cc1")
    assert result.accepted
    assert result.neutralized
    assert result.formal_charge == 0


def test_permanent_charge_is_retained():
    result = standardize_smiles("CCCCCCCCCCCC[N+](C)(C)C")

    assert result.accepted
    assert result.formal_charge == 1
    assert result.charged_structure_retained


def test_n_alkyl_pyridinium_and_zwitterion_keep_unavoidable_charge():
    pyridinium = standardize_smiles("CCCCCCCCCC[n+]1ccccc1")
    zwitterion = standardize_smiles("C[N+](C)(C)CCCCCCCCCC(=O)[O-]")
    assert pyridinium.accepted and pyridinium.formal_charge == 1
    assert pyridinium.charged_structure_retained
    assert zwitterion.accepted
    assert zwitterion.formal_charge == 0
    assert zwitterion.charged_structure_retained


@pytest.mark.parametrize("value", ["", "not-smiles", "CCO |atomProp:0.foo.bar|", "CCO>>CC=O"])
def test_unsafe_or_invalid_input_is_rejected(value):
    result = standardize_smiles(value)
    assert not result.accepted
    assert result.reason


def test_molecular_weight_cutoff_is_strictly_above_200():
    below = standardize_smiles("CCCCCCCCCCCCCC")
    above = standardize_smiles("CCCCCCCCCCCCCCC")

    assert below.molecular_weight <= 200
    assert not below.accepted
    assert below.reason == "molecular_weight_not_above_200"
    assert above.molecular_weight > 200
    assert above.accepted


def test_exact_molecular_weight_200_is_rejected(monkeypatch):
    from app.services import smiles_standardization as module

    monkeypatch.setattr(module.Descriptors, "MolWt", lambda molecule: 200.0)
    result = module.standardize_smiles(IBUPROFEN)
    assert result.molecular_weight == 200.0
    assert not result.accepted


def test_metadata_requires_exact_version_and_config():
    metadata = {"smiles_standardization": standardization_metadata()}
    assert_standardization_compatible(metadata)
    incompatible = {"smiles_standardization": {**metadata["smiles_standardization"], "config": {**STANDARDIZATION_CONFIG, "neutralize": False}}}
    with pytest.raises(ValueError, match="config_mismatch"):
        assert_standardization_compatible(incompatible)
