from src.celery_worker import celery
from src.models import Drug
from src.database import SessionLocal
from rdkit import Chem

@celery.task
def substructure_search_task(substructure: str, limit: int):
    query_mol = Chem.MolFromSmiles(substructure)
    if not query_mol:
        raise ValueError("Invalid SMILES structure")
    
    db = SessionLocal()
    results = []
    
    offset = 0
    results_count = 0
    
    while results_count < limit:
        drugs = db.query(Drug).offset(offset).limit(100).all()
        if not drugs:
            break
        
        for drug in drugs:
            mol = Chem.MolFromSmiles(drug.smiles)
            if mol and mol.HasSubstructMatch(query_mol):
                results.append(drug)
                results_count += 1
                if results_count >= limit:
                    break
        
        offset += len(drugs)
    
    db.close()
    return [drug.id for drug in results]
