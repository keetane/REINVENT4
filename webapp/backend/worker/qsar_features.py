from __future__ import annotations

import math

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, MACCSkeys, rdFingerprintGenerator
from sklearn.base import BaseEstimator, TransformerMixin


class SmilesFeaturizer(BaseEstimator, TransformerMixin):
    """Trusted, persisted sklearn transformer for canonical SMILES."""

    def __init__(
        self,
        fingerprint: str = "morgan",
        feature_mode: str = "fingerprint",
        morgan_radius: int = 2,
        bit_length: int = 2048,
        use_chirality: bool = True,
    ):
        self.fingerprint = fingerprint
        self.feature_mode = feature_mode
        self.morgan_radius = morgan_radius
        self.bit_length = bit_length
        self.use_chirality = use_chirality

    def fit(self, X, y=None):
        return self

    def _fingerprint(self, molecule: Chem.Mol) -> np.ndarray:
        if self.fingerprint == "morgan":
            generator = rdFingerprintGenerator.GetMorganGenerator(
                radius=self.morgan_radius,
                fpSize=self.bit_length,
                includeChirality=self.use_chirality,
            )
            fingerprint = generator.GetFingerprint(molecule)
        elif self.fingerprint == "rdkit":
            fingerprint = Chem.RDKFingerprint(molecule, fpSize=self.bit_length)
        elif self.fingerprint == "maccs":
            fingerprint = MACCSkeys.GenMACCSKeys(molecule)
        else:
            raise ValueError(f"unsupported fingerprint: {self.fingerprint}")
        array = np.zeros((fingerprint.GetNumBits(),), dtype=np.float64)
        DataStructs.ConvertToNumpyArray(fingerprint, array)
        return array

    @staticmethod
    def _descriptors(molecule: Chem.Mol) -> np.ndarray:
        values = []
        for _, function in Descriptors._descList:
            try:
                value = float(function(molecule))
            except Exception:
                value = math.nan
            values.append(value if math.isfinite(value) else math.nan)
        return np.asarray(values, dtype=np.float64)

    def transform(self, X):
        rows = []
        for smiles in X:
            molecule = Chem.MolFromSmiles(str(smiles))
            if molecule is None:
                raise ValueError("invalid SMILES reached the trusted QSAR pipeline")
            parts = []
            if self.feature_mode in {"fingerprint", "both"}:
                parts.append(self._fingerprint(molecule))
            if self.feature_mode in {"descriptors", "both"}:
                parts.append(self._descriptors(molecule))
            if not parts:
                raise ValueError("no QSAR features were selected")
            rows.append(np.concatenate(parts))
        return np.asarray(rows, dtype=np.float64)


def ad_fingerprint(smiles: str):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return generator.GetFingerprint(molecule)
