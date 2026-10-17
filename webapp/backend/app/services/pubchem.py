from __future__ import annotations

from dataclasses import dataclass
from urllib.parse import quote

import httpx

from app.errors import AppError


PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


@dataclass(frozen=True, slots=True)
class PubChemCompound:
    cid: int
    molecular_name: str
    isomeric_smiles: str


def lookup_pubchem_compound(name: str) -> PubChemCompound:
    normalized_name = name.strip()
    if not normalized_name:
        raise AppError(
            status_code=422,
            code="invalid_molecular_name",
            message="Molecular Name を入力してください。",
        )

    encoded_name = quote(normalized_name, safe="")
    url = (
        f"{PUBCHEM_BASE_URL}/compound/name/{encoded_name}"
        "/property/Title,IsomericSMILES/JSON"
    )
    try:
        response = httpx.get(
            url,
            headers={"Accept": "application/json", "User-Agent": "Molecular-Generator/0.1"},
            timeout=7.5,
        )
    except httpx.HTTPError as error:
        raise AppError(
            status_code=502,
            code="pubchem_unavailable",
            message="PubChem に接続できませんでした。時間をおいて再試行してください。",
        ) from error

    if response.status_code == 404:
        raise AppError(
            status_code=404,
            code="pubchem_not_found",
            message=f"PubChem に「{normalized_name}」が見つかりませんでした。",
        )
    if response.status_code >= 400:
        raise AppError(
            status_code=502,
            code="pubchem_unavailable",
            message="PubChem から分子情報を取得できませんでした。",
            context={"pubchem_status": response.status_code},
        )

    try:
        properties = response.json()["PropertyTable"]["Properties"]
        first = properties[0]
        cid = int(first["CID"])
        smiles = (
            first.get("IsomericSMILES")
            or first.get("SMILES")
            or first.get("CanonicalSMILES")
        )
        title = str(first.get("Title") or normalized_name)
    except (KeyError, IndexError, TypeError, ValueError) as error:
        raise AppError(
            status_code=502,
            code="invalid_pubchem_response",
            message="PubChem の応答からSMILESを取得できませんでした。",
        ) from error

    if not isinstance(smiles, str) or not smiles.strip():
        raise AppError(
            status_code=502,
            code="invalid_pubchem_response",
            message="PubChem の応答にSMILESが含まれていません。",
        )

    return PubChemCompound(
        cid=cid,
        molecular_name=title,
        isomeric_smiles=smiles.strip(),
    )
