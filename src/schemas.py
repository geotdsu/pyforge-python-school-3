from pydantic import BaseModel
from typing import Optional, List


class Molecule(BaseModel):
    molecule_id: int
    smiles: str


class UpdateMolecule(BaseModel):
    smiles: str


class SubstructureQuery(BaseModel):
    substructure: str


class SearchResponse(BaseModel):
    message: str
    molecules: Optional[List[Molecule]] = None
