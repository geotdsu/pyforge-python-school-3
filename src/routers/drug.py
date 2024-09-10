from fastapi import Query, status, HTTPException, Depends, APIRouter
from .. import models, schemas
from fastapi_cache.decorator import cache
from ..database import get_db
from sqlalchemy.orm import Session
from typing import List
from rdkit import Chem
from ..logger import logger
from celery.result import AsyncResult
from src.tasks import substructure_search_task

router = APIRouter(tags=['Drugs'])


@router.get("/substructure_search", response_model=dict)
@cache(expire=30)
def initiate_substructure_search(substructure: str = Query(..., description="SMILES structure to search for"),
                                 limit: int = Query(100, le=100)):
    logger.info(
        f"Initiating substructure search for: {substructure}, limit: {limit}")

    # Start Celery task
    task = substructure_search_task.delay(substructure, limit)

    return {"task_id": task.id, "status": task.status}

# Fetch Results of Substructure Search by Task ID


@router.get("/substructure_search/results/{task_id}", response_model=dict)
def get_substructure_search_results(task_id: str):
    logger.info(f"Fetching results for task_id: {task_id}")

    task_result = AsyncResult(task_id)

    if task_result.state == 'PENDING':
        return {"task_id": task_id, "status": "Task is still processing"}
    elif task_result.state == 'SUCCESS':
        drug_ids = task_result.result
        return {
            "task_id": task_id,
            "status": "Task completed",
            "result": drug_ids  # Return the list of matching drug IDs
        }
    else:
        return {"task_id": task_id, "status": task_result.state}


@router.get("/drugs", response_model=List[schemas.DrugResponse])
@cache(expire=10)
def get_drugs(db: Session = Depends(get_db),
              limit: int = Query(10, le=1000),
              skip: int = Query(0, ge=0)):
    drugs = db.query(models.Drug).offset(skip).limit(limit).all()
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
