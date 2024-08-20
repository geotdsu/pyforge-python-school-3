import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from src.main import app
from src.database import get_db, Base
from src import models

# Set up a test database
DATABASE_URL = "sqlite:///./test.db"

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

# Fixture to create the database tables


@pytest.fixture(scope="module")
def create_tables():
    Base.metadata.create_all(bind=engine)
    yield
    # Do nothing here; drop happens in another fixture

# Fixture to drop the database tables


@pytest.fixture(scope="module")
def drop_tables():
    yield
    Base.metadata.drop_all(bind=engine)

# Fixture to clear the database before each test


@pytest.fixture(autouse=True)
def clear_database():
    db = TestingSessionLocal()
    try:
        db.query(models.Drug).delete()
        db.commit()
    finally:
        db.close()


# Test the GET /drugs endpoint
def test_get_drugs(create_tables, drop_tables):
    response = client.get("/drugs")
    assert response.status_code == 200
    assert isinstance(response.json(), list)

# Test the POST /drugs endpoint


def test_create_drug(create_tables, drop_tables):
    drug_data = {
        "name": "TestDrug",
        "smiles": "CCO"
    }
    response = client.post("/drugs", json=drug_data)
    assert response.status_code == 201
    assert response.json()["name"] == "TestDrug"

# Test the GET /drugs/{id} endpoint


def test_get_drug(create_tables, drop_tables):
    response = client.post("/drugs", json={"name": "TestDrug2", "smiles": "C"})
    drug_id = response.json()["id"]

    response = client.get(f"/drugs/{drug_id}")
    assert response.status_code == 200
    assert response.json()["name"] == "TestDrug2"

# Test the DELETE /drugs/{id} endpoint


def test_delete_drug(create_tables, drop_tables):
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


def test_update_drug(create_tables, drop_tables):
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


def test_search_drugs_by_substructure(create_tables, drop_tables):
    # Add some test data
    client.post("/drugs", json={"name": "Drug1", "smiles": "CCO"})
    client.post("/drugs", json={"name": "Drug2", "smiles": "CCC"})

    response = client.get("/substructure_search?substructure=CCO")
    assert response.status_code == 200
    results = response.json()
    print("Results:", results)  # Debugging output
    assert len(results) == 1
    assert results[0]["name"] == "Drug1"


# Test creating a drug with invalid data (missing 'name')
def test_create_drug_invalid_data_missing_name(create_tables, drop_tables):
    drug_data = {
        "smiles": "CCO"
    }
    response = client.post("/drugs", json=drug_data)
    assert response.status_code == 422
    assert "detail" in response.json()

# Test creating a drug with invalid data (invalid SMILES)


def test_create_drug_invalid_data_invalid_smiles(create_tables, drop_tables):
    drug_data = {
        "name": "TestDrugInvalidSMILES",
        "smiles": "invalid_smiles"
    }
    response = client.post("/drugs", json=drug_data)
    assert response.status_code == 422
    assert "detail" in response.json()

# Test retrieving a drug with a non-existent ID


def test_get_drug_non_existent(create_tables, drop_tables):
    response = client.get("/drugs/99999")
    assert response.status_code == 404
    assert "detail" in response.json()

# Test deleting a drug with a non-existent ID


def test_delete_drug_non_existent(create_tables, drop_tables):
    response = client.delete("/drugs/99999")
    assert response.status_code == 404
    assert "detail" in response.json()

# Test updating a drug with a non-existent ID


def test_update_drug_non_existent(create_tables, drop_tables):
    updated_drug_data = {
        "name": "NonExistentDrug",
        "smiles": "COO"
    }
    response = client.put("/drugs/99999", json=updated_drug_data)
    assert response.status_code == 404
    assert "detail" in response.json()

# Test updating a drug with invalid data (missing 'smiles')


def test_update_drug_invalid_data_missing_smiles(create_tables, drop_tables):
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


def test_search_drugs_by_substructure_no_match(create_tables, drop_tables):
    # Add some test data
    client.post("/drugs", json={"name": "Drug3", "smiles": "CCO"})
    client.post("/drugs", json={"name": "Drug4", "smiles": "CCC"})

    response = client.get("/substructure_search?substructure=N")
    assert response.status_code == 200
    results = response.json()
    assert len(results) == 0
