from pydantic import BaseModel, Field


class ParentValidationRequest(BaseModel):
    smiles: str = Field(min_length=1, max_length=10_000)
    model_key: str = Field(default="mol2mol_high", min_length=1, max_length=100)


class ParentValidationResponse(BaseModel):
    input_smiles: str
    canonical_smiles: str
    molblock: str
    atom_count: int = Field(ge=1)
    molecular_formula: str


class PubChemLookupResponse(BaseModel):
    cid: int = Field(ge=1)
    molecular_name: str
    isomeric_smiles: str
