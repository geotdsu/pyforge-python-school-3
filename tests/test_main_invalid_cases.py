import pytest
from fastapi.testclient import TestClient
from src.main import app, molecules
from tests.test_main_valid_cases import reset_molecules
import os


client = TestClient(app)


### Testing add_molecule Function

@pytest.mark.parametrize("molecule_id,smiles_string,expected_status,expected_detail", [
    (5, "Saba", 400, "Invalid SMILES string"),                # Invalid SMILES syntax
    (6, "CCO@", 400, "Invalid SMILES string"),                # Invalid characters
    (7, "CCO123", 400, "Invalid SMILES string"),             # Invalid SMILES format
    (8, "", 400, "SMILES string cannot be empty")           # Empty string
])
def test_add_invalid_smiles_syntax(molecule_id, smiles_string, expected_status, expected_detail):
    response = client.post("/add", json={"molecule_id": molecule_id, "smiles": smiles_string})
    assert response.status_code == expected_status
    assert response.json() == {"detail": expected_detail}


def test_add_non_string_smiles():
    response = client.post("/add", json={"molecule_id": 8, "smiles": 12345})
    assert response.status_code == 422
    assert "detail" in response.json()


def test_add_duplicate_molecule_id():
    client.post("/add", json={"molecule_id": 1, "smiles": "CCO"})
    response = client.post("/add", json={"molecule_id": 1, "smiles": "CCO"})
    assert response.status_code == 400
    assert response.json() == {"detail": "Molecule with the ID 1 already exists"}


def test_add_missing_molecule_id():
    response = client.post("/add", json={"smiles": "CCO"})
    assert response.status_code == 422
    assert "detail" in response.json()  


def test_add_missing_smiles():
    response = client.post("/add", json={"molecule_id": 10})
    assert response.status_code == 422
    assert "detail" in response.json()  



### Testing get_molecule Function

@pytest.mark.parametrize("molecule_id", [
    999,   # ID not in the list
    -1,    # Negative ID
    0      # Zero ID
])
def test_get_non_existing_molecule(molecule_id):
    response = client.get(f"/molecule/{molecule_id}")
    assert response.status_code == 404  


def test_get_non_integer_molecule_id():
    response = client.get("/molecule/Testosterone")
    assert response.status_code == 422
    assert "detail" in response.json()



### Testing update_molecule Function

@pytest.mark.parametrize("molecule_id,smiles,expected_status,expected_detail", [
    (1, "Protein", 400, "Invalid SMILES string"),              # Invalid SMILES syntax
    (2, "CCO@", 400, "Invalid SMILES string"),              # Invalid characters
    (3, "CCO123", 400, "Invalid SMILES string"),            # Invalid SMILES format
    (4, "", 400, "SMILES string cannot be empty"),          # Empty string
    (999, "CCO", 404, "Molecule with the ID 999 is not found")  # Non-existing molecule ID
])
def test_update_invalid_smiles(molecule_id, smiles, expected_status, expected_detail):
    response = client.put(f"/molecule/{molecule_id}", json={"smiles": smiles})
    assert response.status_code == expected_status
    assert response.json() == {"detail": expected_detail}    


def test_update_non_string_smiles():
    response = client.put("/molecule/1", json={"smiles": 200315})
    assert response.status_code == 422
    assert "detail" in response.json()


def test_update_missing_smiles():
    response = client.put("/molecule/1", json={})
    assert response.status_code == 422
    assert "detail" in response.json()
    


## Testing delete_molecule Function

@pytest.mark.parametrize("molecule_id,expected_status,expected_detail", [
    (999, 404, "Molecule with the ID 999 is not found"),
    (-1, 404, "Molecule with the ID -1 is not found"),
    (0, 404, "Molecule with the ID 0 is not found"),
])
def test_delete_non_existing_molecule(molecule_id, expected_status, expected_detail):
    response = client.delete(f"/molecule/{molecule_id}")
    assert response.status_code == expected_status
    assert response.json() == {"detail": expected_detail}


def test_delete_non_integer_id():
    response = client.delete("/molecule/Testosterone")
    assert response.status_code == 422
    assert "detail" in response.json()



## Test get_molecules Function

def test_get_molecules_empty(reset_molecules):
    response = client.get("/molecules")
    assert response.status_code == 200
    assert response.json() == []


def test_get_molecules_invalid_query(reset_molecules):
    # Adding valid molecules
    client.post("/add", json={"molecule_id": 1, "smiles": "CCO"})
    client.post("/add", json={"molecule_id": 2, "smiles": "C[C@H](O)[C@@H](N)C(=O)O"})
    
    response = client.get("/molecules")
    assert response.status_code == 200
    assert response.json() == [
        {"molecule_id": 1, "smiles": "CCO"},
        {"molecule_id": 2, "smiles": "C[C@H](O)[C@@H](N)C(=O)O"}
    ]



## Testing Substructue Search Function

def test_substructure_search_empty_substructure():
    response = client.post("/search", json={"substructure": ""})
    assert response.status_code == 400
    assert response.json() == {"detail": "Substructure query cannot be empty."}


def test_substructure_search_invalid_substructure():
    response = client.post("/search", json={"substructure": "Saba"})
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid substructure SMILES."}
    


## Testing upload_file Function

@pytest.mark.parametrize("file_path, expected_status, expected_detail", [
    ("real_csv.txt", 400, "Invalid file format. Only CSV files are supported"), #
    ("empty.csv", 400, "Uploaded file is empty"),
    ("missing_smiles_column.csv", 400, "CSV file must contain 'SMILES', 'smiles', or 'Smiles' header."),
    ("missing_smiles_row.csv", 400, "Missing 'SMILES' or 'smiles' or 'Smiles' value in row 3"),
    ("incorrect_smiles_value.csv", 400, 'Invalid SMILES string in row 3: Saba')
])
def test_upload_file(file_path, expected_status, expected_detail):
    # This function will test if a CSV file has incorrect SMILES in a row

    with open(file_path, 'rb') as file:
        response = client.post("/uploadfile", files={"file": file})

    assert response.status_code == expected_status
    assert response.json() == {"detail": expected_detail}


def test_upload_file_no_file_provided():
    response = client.post("/uploadfile", data={})

    assert response.status_code == 422
    assert "detail" in response.json()


@pytest.mark.xfail
def test_upload_file_large_file_size():
    # This function will test if the CSV is too large or not and the limit is 100 MB
    # This test will fail because I don't have a CSV file that is over 100MB :) 
    # if you do, put it in tests
    file_path = 'your_big_file.csv'


    with open(file_path, 'rb') as file:
        response = client.post("/uploadfile", files={"file": file})

    assert response.status_code == 400
    assert response.json() == {"detail": "File is too large. Maximum size is 100 MB"}



# Testing get_server Function

def test_get_server_invalid():
    os.environ["SERVER_ID"] = ""
    response = client.get("/server")
    assert response.status_code == 200
    assert response.json() == {"server_id": ""}

    del os.environ["SERVER_ID"]
    response = client.get("/server")
    assert response.status_code == 200
    assert response.json() == {"server_id": "1"}