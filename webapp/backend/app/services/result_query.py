from __future__ import annotations

import csv
import io
import math
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import rdkit
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Crippen, Descriptors, Draw, Lipinski, QED, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Contrib.SA_Score import sascorer

from app.services.results import MAX_RESULT_BYTES, ResultValidationError


CORE_DESCRIPTORS = ("MW", "RotB", "HBD", "HBA", "LogP", "TPSA")
ADDED_DESCRIPTORS = (
    "QED",
    "SAScore",
    "FractionCSP3",
    "RingCount",
    "AromaticRingCount",
    "AliphaticRingCount",
    "HeavyAtomCount",
    "FormalCharge",
    "LipinskiViolations",
    "MolMR",
    "ParentSimilarity",
    "ScaffoldSimilarity",
)
REQUIRED_DESCRIPTORS = CORE_DESCRIPTORS + ADDED_DESCRIPTORS
OPTIONAL_DESCRIPTORS = (
    "NLL",
    "CtAmides",
    "CtRings",
    "CtAromaticRings",
    "n_member_ring",
    "Tanimoto_Similarity_to_Parent",
)
ROTATABLE_ALIASES = ("RotB", "RoB", "CtRotB", "CtRotBonds", "Rotatable_Bonds")
ALL_DESCRIPTORS = REQUIRED_DESCRIPTORS + OPTIONAL_DESCRIPTORS
MAX_ROWS = 100_000
MAX_SCATTER_POINTS = 500


@dataclass(frozen=True, slots=True)
class DescriptorInfo:
    name: str
    source: str
    aliases: tuple[str, ...]
    minimum: float | None
    maximum: float | None
    missing_count: int


@dataclass(frozen=True, slots=True)
class ResultRow:
    row_number: int
    original: dict[str, str]
    smiles: str
    canonical_smiles: str
    molecule: Chem.Mol
    values: dict[str, float | None]


@dataclass(frozen=True, slots=True)
class ResultDataset:
    columns: tuple[str, ...]
    rows: tuple[ResultRow, ...]
    total_rows: int
    invalid_smiles: int
    descriptors: tuple[DescriptorInfo, ...]
    warnings: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ResultFilters:
    ranges: dict[str, tuple[float | None, float | None]]
    sort_by: str | None = None
    sort_direction: str = "asc"


def load_result_dataset(path: Path, parent_smiles: str | None = None) -> ResultDataset:
    if path.is_symlink() or not path.is_file():
        raise ResultValidationError("result CSV is unavailable")
    size = path.stat().st_size
    if size <= 0 or size > MAX_RESULT_BYTES:
        raise ResultValidationError("result CSV size is invalid")

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "SMILES" not in reader.fieldnames:
            raise ResultValidationError("result CSV must contain a SMILES column")
        raw_rows: list[dict[str, str]] = []
        for index, row in enumerate(reader, start=1):
            if index > MAX_ROWS:
                raise ResultValidationError("result CSV exceeds the row limit")
            raw_rows.append({key: value or "" for key, value in row.items() if key})

    warnings: list[str] = []
    rotatable_source = next((name for name in ROTATABLE_ALIASES if name in reader.fieldnames), None)
    present_rotatable = [name for name in ROTATABLE_ALIASES if name in reader.fieldnames]
    if len(present_rotatable) > 1:
        conflicts = sum(
            1
            for row in raw_rows
            if len({_finite(row.get(name, "")) for name in present_rotatable if _finite(row.get(name, "")) is not None}) > 1
        )
        if conflicts:
            warnings.append(
                f"RotB aliases conflict in {conflicts} rows; {rotatable_source} takes precedence."
            )

    reference = _reference_molecule(parent_smiles)
    reference_scaffold = _scaffold(reference) if reference is not None else None
    rows: list[ResultRow] = []
    invalid = 0
    sources: dict[str, str] = {}
    for row_number, raw in enumerate(raw_rows, start=1):
        smiles = raw.get("SMILES", "").strip()
        with rdBase.BlockLogs():
            molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            invalid += 1
            continue
        computed = _computed_descriptors(molecule, reference, reference_scaffold)
        values: dict[str, float | None] = {}
        for name in REQUIRED_DESCRIPTORS:
            source = rotatable_source if name == "RotB" else (name if name in raw else None)
            existing = _finite(raw.get(source, "")) if source else None
            values[name] = existing if existing is not None else computed[name]
            if source:
                sources.setdefault(name, source)
            elif name == "SAScore":
                sources.setdefault(name, "RDKit Contrib SA_Score")
            elif name in {"ParentSimilarity", "ScaffoldSimilarity"}:
                sources.setdefault(name, "RDKit Morgan/Tanimoto")
            else:
                sources.setdefault(name, f"RDKit {rdkit.__version__}")
        for name in OPTIONAL_DESCRIPTORS:
            if name in raw:
                values[name] = _finite(raw[name])
                sources.setdefault(name, name)
        rows.append(
            ResultRow(
                row_number=row_number,
                original=raw,
                smiles=smiles,
                canonical_smiles=Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True),
                molecule=molecule,
                values=values,
            )
        )

    descriptor_names = list(REQUIRED_DESCRIPTORS)
    descriptor_names.extend(name for name in OPTIONAL_DESCRIPTORS if name in reader.fieldnames)
    descriptors = tuple(
        DescriptorInfo(
            name=name,
            source=sources.get(name, name),
            aliases=ROTATABLE_ALIASES if name == "RotB" else (),
            minimum=_minimum(row.values.get(name) for row in rows),
            maximum=_maximum(row.values.get(name) for row in rows),
            missing_count=sum(row.values.get(name) is None for row in rows),
        )
        for name in descriptor_names
    )
    return ResultDataset(
        columns=tuple(reader.fieldnames),
        rows=tuple(rows),
        total_rows=len(raw_rows),
        invalid_smiles=invalid,
        descriptors=descriptors,
        warnings=tuple(warnings),
    )


def query_rows(dataset: ResultDataset, filters: ResultFilters) -> list[ResultRow]:
    allowed = {item.name for item in dataset.descriptors}
    for name, (minimum, maximum) in filters.ranges.items():
        if name not in allowed:
            raise ValueError(f"unsupported filter descriptor: {name}")
        if minimum is not None and maximum is not None and minimum > maximum:
            raise ValueError(f"minimum exceeds maximum for {name}")
    if filters.sort_by is not None and filters.sort_by not in allowed:
        raise ValueError("unsupported sort descriptor")
    if filters.sort_direction not in {"asc", "desc"}:
        raise ValueError("unsupported sort direction")

    rows = [
        row
        for row in dataset.rows
        if all(_within(row.values.get(name), minimum, maximum) for name, (minimum, maximum) in filters.ranges.items())
    ]
    if filters.sort_by:
        name = filters.sort_by
        rows.sort(
            key=lambda row: (
                row.values.get(name) is None,
                row.values.get(name) if row.values.get(name) is not None else 0.0,
            ),
            reverse=filters.sort_direction == "desc",
        )
    return rows


def molecule_svg(molecule: Chem.Mol) -> str:
    drawer = Draw.MolDraw2DSVG(260, 180)
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def statistics_for(dataset: ResultDataset, rows: list[ResultRow]) -> tuple[dict[str, dict[str, Any]], dict[str, list[dict[str, Any]]]]:
    summaries: dict[str, dict[str, Any]] = {}
    histograms: dict[str, list[dict[str, Any]]] = {}
    for descriptor in dataset.descriptors:
        values = [row.values[descriptor.name] for row in rows if row.values.get(descriptor.name) is not None]
        numeric = [float(value) for value in values]
        summaries[descriptor.name] = _summary(numeric, len(rows) - len(numeric))
        if descriptor.name in REQUIRED_DESCRIPTORS:
            histograms[descriptor.name] = _histogram(numeric)
    return summaries, histograms


def scatter_for(
    dataset: ResultDataset,
    rows: list[ResultRow],
    x_descriptor: str,
    y_descriptor: str,
    limit: int = MAX_SCATTER_POINTS,
) -> tuple[list[dict[str, float | int]], bool]:
    allowed = {item.name for item in dataset.descriptors}
    if x_descriptor not in allowed or y_descriptor not in allowed:
        raise ValueError("unsupported scatter descriptor")
    candidates = [
        row
        for row in rows
        if row.values.get(x_descriptor) is not None and row.values.get(y_descriptor) is not None
    ]
    sampled = len(candidates) > limit
    if sampled:
        # Deterministic, evenly-spaced sampling keeps the same query reproducible.
        candidates = [candidates[index * len(candidates) // limit] for index in range(limit)]
    return [
        {
            "row_number": row.row_number,
            "x": float(row.values[x_descriptor]),
            "y": float(row.values[y_descriptor]),
        }
        for row in candidates
    ], sampled


def export_csv(
    columns: Iterable[str],
    rows: Iterable[ResultRow],
    descriptor_names: Iterable[str] = (),
) -> str:
    output = io.StringIO(newline="")
    fieldnames = list(columns)
    fieldnames.extend(name for name in descriptor_names if name not in fieldnames)
    writer = csv.DictWriter(output, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        values = {key: _safe_csv_cell(value) for key, value in row.original.items()}
        values.update({name: value for name, value in row.values.items() if name not in values})
        writer.writerow(values)
    return output.getvalue()


def export_sdf(rows: Iterable[ResultRow]) -> str:
    output = io.StringIO()
    writer = Chem.SDWriter(output)
    for row in rows:
        molecule = Chem.Mol(row.molecule)
        molecule.SetProp("SMILES", row.canonical_smiles)
        for name, value in row.values.items():
            if value is not None:
                molecule.SetProp(name, str(value))
        writer.write(molecule)
    writer.flush()
    writer.close()
    return output.getvalue()


def _computed_descriptors(
    molecule: Chem.Mol,
    reference: Chem.Mol | None = None,
    reference_scaffold: Chem.Mol | None = None,
) -> dict[str, float | None]:
    mw = float(Descriptors.MolWt(molecule))
    logp = float(Crippen.MolLogP(molecule))
    hbd = float(Lipinski.NumHDonors(molecule))
    hba = float(Lipinski.NumHAcceptors(molecule))
    scaffold = _scaffold(molecule)
    return {
        "MW": mw,
        "RotB": float(Lipinski.NumRotatableBonds(molecule)),
        "HBD": hbd,
        "HBA": hba,
        "LogP": logp,
        "TPSA": float(rdMolDescriptors.CalcTPSA(molecule)),
        "QED": float(QED.qed(molecule)),
        "SAScore": float(sascorer.calculateScore(molecule)),
        "FractionCSP3": float(rdMolDescriptors.CalcFractionCSP3(molecule)),
        "RingCount": float(rdMolDescriptors.CalcNumRings(molecule)),
        "AromaticRingCount": float(rdMolDescriptors.CalcNumAromaticRings(molecule)),
        "AliphaticRingCount": float(rdMolDescriptors.CalcNumAliphaticRings(molecule)),
        "HeavyAtomCount": float(molecule.GetNumHeavyAtoms()),
        "FormalCharge": float(Chem.GetFormalCharge(molecule)),
        "LipinskiViolations": float(sum((mw > 500, logp > 5, hbd > 5, hba > 10))),
        "MolMR": float(Crippen.MolMR(molecule)),
        "ParentSimilarity": _similarity(molecule, reference),
        "ScaffoldSimilarity": _similarity(scaffold, reference_scaffold),
    }


def _reference_molecule(smiles: str | None) -> Chem.Mol | None:
    if not smiles:
        return None
    with rdBase.BlockLogs():
        molecule = Chem.MolFromSmiles(smiles.replace("|", "."))
    if molecule is None:
        return None
    editable = Chem.RWMol(molecule)
    for atom_index in reversed([atom.GetIdx() for atom in editable.GetAtoms() if atom.GetAtomicNum() == 0]):
        editable.RemoveAtom(atom_index)
    cleaned = editable.GetMol()
    try:
        Chem.SanitizeMol(cleaned)
    except ValueError:
        return None
    return cleaned if cleaned.GetNumAtoms() else None


def _scaffold(molecule: Chem.Mol | None) -> Chem.Mol | None:
    if molecule is None:
        return None
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    return scaffold if scaffold.GetNumAtoms() else None


def _similarity(left: Chem.Mol | None, right: Chem.Mol | None) -> float | None:
    if left is None or right is None:
        return None
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return float(DataStructs.TanimotoSimilarity(
        generator.GetFingerprint(left),
        generator.GetFingerprint(right),
    ))


def _finite(value: str | float | None) -> float | None:
    try:
        parsed = float(value) if value not in {None, ""} else None
    except (TypeError, ValueError):
        return None
    return parsed if parsed is not None and math.isfinite(parsed) else None


def _minimum(values: Iterable[float | None]) -> float | None:
    present = [value for value in values if value is not None]
    return min(present) if present else None


def _maximum(values: Iterable[float | None]) -> float | None:
    present = [value for value in values if value is not None]
    return max(present) if present else None


def _within(value: float | None, minimum: float | None, maximum: float | None) -> bool:
    if minimum is None and maximum is None:
        return True
    if value is None:
        return False
    return (minimum is None or value >= minimum) and (maximum is None or value <= maximum)


def _summary(values: list[float], missing: int) -> dict[str, Any]:
    if not values:
        return {"count": 0, "missing_count": missing}
    ordered = sorted(values)
    quartiles = statistics.quantiles(ordered, n=4, method="inclusive") if len(ordered) > 1 else [ordered[0]] * 3
    return {
        "count": len(ordered),
        "missing_count": missing,
        "minimum": ordered[0],
        "maximum": ordered[-1],
        "mean": statistics.fmean(ordered),
        "median": statistics.median(ordered),
        "standard_deviation": statistics.pstdev(ordered),
        "q1": quartiles[0],
        "q3": quartiles[2],
    }


def _histogram(values: list[float], bins: int = 10) -> list[dict[str, Any]]:
    if not values:
        return []
    minimum, maximum = min(values), max(values)
    if minimum == maximum:
        return [{"start": minimum, "end": maximum, "count": len(values)}]
    width = (maximum - minimum) / bins
    counts = [0] * bins
    for value in values:
        index = min(bins - 1, int((value - minimum) / width))
        counts[index] += 1
    return [
        {"start": minimum + index * width, "end": minimum + (index + 1) * width, "count": count}
        for index, count in enumerate(counts)
    ]


def _safe_csv_cell(value: str) -> str:
    if value.startswith(("=", "+", "-", "@")):
        try:
            float(value)
        except ValueError:
            return "'" + value
    return value
