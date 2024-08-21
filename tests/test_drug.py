import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from src.main import app
from src.database import get_db, Base
from src import models

# Set up a test database
DATABASE_URL = "sqlite:///tests/test.db"

engine = create_engine(DATABASE_URL)
TestingSessionLocal = sessionmaker(
    autocommit=False,
    autoflush=False,
    bind=engine)


# Override the get_db function to use the test database
def override_get_db():
    try:
        db = TestingSessionLocal()
        yield db
    finally:
        db.close()


app.dependency_overrides[get_db] = override_get_db

client = TestClient(app)


@pytest.fixture(scope="module")
def setup_db():
    # Create tables
    Base.metadata.create_all(bind=engine)
    yield
    # Drop tables after tests are done
    Base.metadata.drop_all(bind=engine)


@pytest.fixture(autouse=True)
def clear_database():
    db = TestingSessionLocal()
    try:
        db.query(models.Drug).delete()
        db.commit()
    finally:
        db.close()


# Test the GET /drugs endpoint
def test_get_drugs(setup_db, clear_database):
    response = client.get("/drugs")
    assert response.status_code == 200
    assert response.json() == []


# Test the POST /drugs endpoint
def test_create_drug(setup_db, clear_database):
    drug_data = {
        "name": "TestDrug",
        "smiles": "CCO"
    }
    response = client.post("/drugs", json=drug_data)
    assert response.status_code == 201
    assert response.json()["name"] == "TestDrug"


# Test the GET /drugs/{id} endpoint
def test_get_drug(setup_db, clear_database):
    response = client.post("/drugs", json={"name": "TestDrug2", "smiles": "C"})
    drug_id = response.json()["id"]

    response = client.get(f"/drugs/{drug_id}")
    assert response.status_code == 200
    assert response.json()["name"] == "TestDrug2"


# Test the DELETE /drugs/{id} endpoint
def test_delete_drug(setup_db, clear_database):
    response = client.post(
        "/drugs",
        json={
            "name": "TestDrug3",
            "smiles": "CC"})
    drug_id = response.json()["id"]

    response = client.delete(f"/drugs/{drug_id}")
    assert response.status_code == 204

    response = client.get(f"/drugs/{drug_id}")
    assert response.status_code == 404


# Test the PUT /drugs/{id} endpoint
def test_update_drug(setup_db, clear_database):
    response = client.post(
        "/drugs",
        json={
            "name": "TestDrug4",
            "smiles": "CO"})
    drug_id = response.json()["id"]

    updated_drug_data = {
        "name": "UpdatedDrug",
        "smiles": "COO"
    }
    response = client.put(f"/drugs/{drug_id}", json=updated_drug_data)
    assert response.status_code == 200
    assert response.json()["name"] == "UpdatedDrug"


# Test the GET /substructure_search endpoint
def test_search_drugs_by_substructure(setup_db, clear_database):
    # Add some test data
    client.post("/drugs", json={"name": "Drug1", "smiles": "CCO"})
    client.post("/drugs", json={"name": "Drug2", "smiles": "CCC"})

    response = client.get("/substructure_search?substructure=CCO")
    assert response.status_code == 200
    results = response.json()
    assert len(results) == 1
    assert results[0]["name"] == "Drug1"


# Test creating a drug with empty smiles (empty 'smiles')
def test_create_drug_empty_smiles(setup_db, clear_database):
    drug_data = {"name": "DrugWithEmptySMILES", "smiles": ""}
    response = client.post("/drugs", json=drug_data)
    assert response.status_code == 422
    assert "detail" in response.json()


# Test creating a drug with invalid data (missing 'name')
def test_create_drug_invalid_data_missing_name(setup_db, clear_database):
    drug_data = {
        "smiles": "CCO"
    }
    response = client.post("/drugs", json=drug_data)
    assert response.status_code == 422
    assert "detail" in response.json()


# Test pagination with a limit on the number of retrieved drugs
def test_get_drugs_pagination(setup_db, clear_database):
    for i in range(15):
        client.post("/drugs", json={"name": f"Drug{i}", "smiles": "CCO"})

    response = client.get("/drugs?limit=10")
    assert response.status_code == 200
    assert len(response.json()) == 10


# Test creating a drug with invalid data (invalid SMILES)
def test_create_drug_invalid_data_invalid_smiles(setup_db, clear_database):
    drug_data = {
        "name": "TestDrugInvalidSMILES",
        "smiles": "invalid_smiles"
    }
    response = client.post("/drugs", json=drug_data)
    assert response.status_code == 422
    assert "detail" in response.json()


# Test retrieving a drug with a non-existent ID
def test_get_drug_non_existent(setup_db, clear_database):
    response = client.get("/drugs/99999")
    assert response.status_code == 404
    assert "detail" in response.json()


# Test deleting a drug with a non-existent ID
def test_delete_drug_non_existent(setup_db, clear_database):
    response = client.delete("/drugs/99999")
    assert response.status_code == 404
    assert "detail" in response.json()


# Test updating a drug with a non-existent ID
def test_update_drug_non_existent(setup_db, clear_database):
    updated_drug_data = {
        "name": "NonExistentDrug",
        "smiles": "COO"
    }
    response = client.put("/drugs/99999", json=updated_drug_data)
    assert response.status_code == 404
    assert "detail" in response.json()


# Test updating a drug with invalid data (missing 'smiles')
def test_update_drug_invalid_data_missing_smiles(setup_db, clear_database):
    response = client.post(
        "/drugs",
        json={
            "name": "TestDrug5",
            "smiles": "CO"})
    drug_id = response.json()["id"]

    updated_drug_data = {
        "name": "UpdatedDrugMissingSmiles"
    }
    response = client.put(f"/drugs/{drug_id}", json=updated_drug_data)
    assert response.status_code == 422
    assert "detail" in response.json()


# Test substructure search with no matching drugs
def test_search_drugs_by_substructure_no_match(setup_db, clear_database):
    # Add some test data
    client.post("/drugs", json={"name": "Drug3", "smiles": "CCO"})
    client.post("/drugs", json={"name": "Drug4", "smiles": "CCC"})

    response = client.get("/substructure_search?substructure=N")
    assert response.status_code == 200
    results = response.json()
    assert len(results) == 0
