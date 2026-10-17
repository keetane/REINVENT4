from typing import Annotated

from fastapi import APIRouter, Query, Request

from app.schemas.parents import (
    ParentValidationRequest,
    ParentValidationResponse,
    PubChemLookupResponse,
)
from app.services.parents import validate_parent_smiles
from app.services.pubchem import lookup_pubchem_compound


router = APIRouter(prefix="/parents", tags=["parents"])


@router.get("/pubchem", response_model=PubChemLookupResponse)
def lookup_pubchem(
    name: Annotated[str, Query(min_length=1, max_length=200)],
) -> PubChemLookupResponse:
    compound = lookup_pubchem_compound(name)
    return PubChemLookupResponse(
        cid=compound.cid,
        molecular_name=compound.molecular_name,
        isomeric_smiles=compound.isomeric_smiles,
    )


@router.post("/validate", response_model=ParentValidationResponse)
def validate_parent(
    request: Request, payload: ParentValidationRequest
) -> ParentValidationResponse:
    parent = validate_parent_smiles(
        payload.smiles,
        max_atoms=request.app.state.settings.parent_max_atoms,
        model_key=payload.model_key,
    )
    return ParentValidationResponse(
        input_smiles=parent.input_smiles,
        canonical_smiles=parent.canonical_smiles,
        molblock=parent.molblock,
        atom_count=parent.atom_count,
        molecular_formula=parent.molecular_formula,
    )
