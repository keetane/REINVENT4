"""Zero the final RL reward for invalid or low-molecular-weight molecules."""

from __future__ import annotations

__all__ = ["SmilesStandardizationEligibility"]

from typing import List

import numpy as np
from pydantic.dataclasses import dataclass

from app.services.smiles_standardization import standardize_smiles

from .add_tag import add_tag
from .component_results import ComponentResults


@add_tag("__parameters")
@dataclass
class Parameters:
    standardization_version: List[str]


@add_tag("__component", "penalty")
class SmilesStandardizationEligibility:
    def __init__(self, params: Parameters):
        self.number_of_endpoints = len(params.standardization_version)

    def __call__(self, smilies: List[str]) -> ComponentResults:
        values = np.asarray([1.0 if standardize_smiles(value).accepted else 0.0 for value in smilies])
        return ComponentResults([values.copy() for _ in range(self.number_of_endpoints)])
