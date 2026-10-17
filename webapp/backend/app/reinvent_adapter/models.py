from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ModelSpec:
    key: str
    mode: str
    relative_path: str


MODEL_SPECS: dict[str, ModelSpec] = {
    spec.key: spec
    for spec in (
        ModelSpec("mol2mol_high", "mol2mol", "priors/mol2mol_high_similarity.prior"),
        ModelSpec("mol2mol_medium", "mol2mol", "priors/mol2mol_medium_similarity.prior"),
        ModelSpec("mol2mol_mmp", "mol2mol", "priors/mol2mol_mmp.prior"),
        ModelSpec("mol2mol_scaffold", "mol2mol", "priors/mol2mol_scaffold.prior"),
        ModelSpec(
            "mol2mol_scaffold_generic",
            "mol2mol",
            "priors/mol2mol_scaffold_generic.prior",
        ),
        ModelSpec("mol2mol_similarity", "mol2mol", "priors/mol2mol_similarity.prior"),
        ModelSpec("reinvent", "reinvent", "priors/reinvent.prior"),
        ModelSpec("libinvent", "libinvent", "priors/libinvent.prior"),
        ModelSpec("linkinvent", "linkinvent", "priors/linkinvent.prior"),
    )
}

