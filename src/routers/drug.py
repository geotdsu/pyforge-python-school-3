from fastapi import status, HTTPException, Depends, APIRouter
from .. import models, schemas
from ..database import get_db
from sqlalchemy.orm import Session
from typing import List, Optional
from rdkit.Chem import AllChem
from rdkit import Chem


router = APIRouter(tags=['Drugs'])


@router.get("/drugs", response_model=List[schemas.DrugResponse])
def get_drugs(db: Session = Depends(get_db), 
              limit: int = 10, skip: int = 0, search: Optional[str] = ""):

    drugs = db.query(models.Drug).filter(models.Drug.name.contains(search)).limit(limit).offset(skip).all()
    return drugs
    

@router.post("/drugs", status_code=status.HTTP_201_CREATED, response_model=schemas.DrugResponse)
def create_drugs(drug: schemas.DrugAdd, db: Session = Depends(get_db)):
     
    new_drug = models.Drug(**drug.model_dump())
    db.add(new_drug)
    db.commit()
    db.refresh(new_drug)

    return new_drug



@router.get("/drugs/{id}", response_model=schemas.DrugResponse)
def get_drug(id: int, db: Session = Depends(get_db)):

    drug = db.query(models.Drug).filter(models.Drug.id == id).first()

    if not drug:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                            detail=f"drug with id: {id} was not found")

    return drug


@router.delete("/drugs/{id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_drug(id: int, db: Session = Depends(get_db)):

    drug = db.query(models.Drug).filter(models.Drug.id == id) 

    if drug.first() == None:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                            detail=f"drug with id: {id} was not found")
    
    drug.delete(synchronize_session=False)
    db.commit()


@router.put("/drugs/{id}", response_model=schemas.DrugResponse)
def update_drug(id: int, updated_drug: schemas.DrugAdd, db: Session = Depends(get_db)):

    drug_query = db.query(models.Drug).filter(models.Drug.id == id)

    drug = drug_query.first()

    if drug == None:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                            detail=f"post with id: {id} does not exist")
    
    
    drug_query.update(updated_drug.model_dump(), synchronize_session=False)
    db.commit()
    return drug_query.first()


@router.get("/substructure_search", response_model=List[schemas.DrugResponse])
def search_drugs_by_substructure(substructure: str, db: Session = Depends(get_db)):
    query_mol = Chem.MolFromSmiles(substructure)
    if query_mol is None:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid substructure query")

    drugs = db.query(models.Drug).all()
    
    matching_drugs = []

    for drug in drugs:
        mol = Chem.MolFromSmiles(drug.smiles)
        if mol and mol.HasSubstructMatch(query_mol):
            matching_drugs.append(drug)

    return matching_drugs