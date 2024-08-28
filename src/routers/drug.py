from fastapi import Query, status, HTTPException, Depends, APIRouter
from .. import models, schemas
from ..database import get_db
from sqlalchemy.orm import Session
from typing import List, Optional
from rdkit import Chem
from ..logger import logger


router = APIRouter(tags=['Drugs'])


def drugs_iterator(db: Session, limit: int, skip: int):
    offset = skip
    total_fetched = 0

    while total_fetched < limit:
        drugs = db.query(models.Drug).offset(offset).limit(
            min(limit - total_fetched, 100)).all()

        if not drugs:
            break

        for drug in drugs:
            yield drug
            total_fetched += 1

            if total_fetched >= limit:
                return

        offset += len(drugs)


def substructure_search_iterator(db: Session, query_mol: Chem.Mol, limit: int):
    offset = 0
    results_count = 0

    while results_count < limit:
        drugs = db.query(models.Drug).offset(offset).limit(100).all()

        if not drugs:
            break

        for drug in drugs:
            mol = Chem.MolFromSmiles(drug.smiles)
            if mol and mol.HasSubstructMatch(query_mol):
                yield drug
                results_count += 1
                if results_count >= limit:
                    return

        offset += len(drugs)


@router.get("/drugs", response_model=List[schemas.DrugResponse])
def get_drugs(
    db: Session = Depends(get_db),
    limit: int = Query(10, le=100),
    skip: int = Query(0, ge=0),
    search: Optional[str] = ""
):
    logger.info(
        "Fetching drugs with limit: %s, skip: %s, search query: %s", limit, skip, search)

    drugs_iter = drugs_iterator(db, limit, skip)
    drugs = list(drugs_iter)

    return drugs


@router.post("/drugs", status_code=status.HTTP_201_CREATED,
             response_model=schemas.DrugResponse)
def create_drugs(drug: schemas.DrugAdd, db: Session = Depends(get_db)):
    logger.info(f'Creating new drug with name="{drug.name}"')

    new_drug = models.Drug(**drug.model_dump())
    db.add(new_drug)
    db.commit()
    db.refresh(new_drug)

    return new_drug


@router.get("/drugs/{id}", response_model=schemas.DrugResponse)
def get_drug(id: int, db: Session = Depends(get_db)):
    logger.info(f'Retrieving drug with id={id}')
    drug = db.query(models.Drug).filter(models.Drug.id == id).first()

    if not drug:
        logger.error(f'Drug with id={id} not found')
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                            detail=f"drug with id: {id} was not found")

    return drug


@router.delete("/drugs/{id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_drug(id: int, db: Session = Depends(get_db)):
    logger.info(f'Deleting drug with id={id}')

    drug = db.query(models.Drug).filter(models.Drug.id == id)

    if drug.first() is None:
        logger.error(f'Drug with id={id} not found for deletion')
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                            detail=f"drug with id: {id} was not found")

    drug.delete(synchronize_session=False)
    db.commit()


@router.put("/drugs/{id}", response_model=schemas.DrugResponse)
def update_drug(
        id: int,
        updated_drug: schemas.DrugAdd,
        db: Session = Depends(get_db)):
    logger.info(f'Updating drug with id={id}')

    drug_query = db.query(models.Drug).filter(models.Drug.id == id)

    drug = drug_query.first()

    if drug is None:
        logger.error(f'Drug with id={id} not found for update')
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                            detail=f"post with id: {id} does not exist")

    drug_query.update(updated_drug.model_dump(), synchronize_session=False)
    db.commit()
    return drug_query.first()


@router.get('/substructure_search', response_model=List[schemas.DrugResponse])
def get_substructure_match(substructure: str = Query(..., description="SMILES structure to search for"),
                           db: Session = Depends(get_db),
                           limit: int = Query(100, le=100)):
    logger.info(
        f'Substructure search for substructure="{substructure}" with limit={limit}')

    substructure = substructure.strip()

    query_mol = Chem.MolFromSmiles(substructure)
    if not query_mol:
        logger.error(
            f'Invalid molecule for substructure search: {substructure}')
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f'{status.HTTP_400_BAD_REQUEST} BAD REQUEST - not a valid molecule'
        )

    sub_matches = []
    for drug in substructure_search_iterator(db, query_mol, limit):
        sub_matches.append(drug)
        if len(sub_matches) >= limit:
            break

    return sub_matches
