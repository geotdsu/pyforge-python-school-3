import pytest
from fastapi.testclient import TestClient
from src.main import app, molecules
import os

client = TestClient(app)


@pytest.fixture
def reset_molecules():
    # Reset the in-memory storage before each test
    molecules.clear()
    yield



# Testing the root of the application
def test_root():
    response = client.get("/")
    assert response.status_code == 200
    assert response.json() == {'message': 'Molecule Fast API application'}



### Testing add_molecule Function

@pytest.mark.parametrize("molecule_id,smiles,expected_status", [
    (5, "CCO", 201),                          # Simple SMILES
    (6, "C1CCCCC1", 201),                     # Complex SMILES (Cyclohexane)
    (7, "CC(C)C1=CC2=C(C=C1)C(=O)C=C2", 201), # SMILES with rings and branches
    (8, "C[C@H](O)[C@@H](N)C(=O)O", 201)      # SMILES with stereo chemistry
])
def test_add_valid_simple_smiles(molecule_id, smiles, expected_status):
    response = client.post("/add", json={"molecule_id": molecule_id, "smiles": smiles})
    assert response.status_code == expected_status
    assert response.json() == {"molecule_id": molecule_id, "smiles": smiles}



### Testing get_molecule Function

@pytest.mark.parametrize("molecule_id,expected_status,expected_response", [
    (1, 200, {"molecule_id": 1, "smiles": "CCO"}),
    (2, 200, {"molecule_id": 2, "smiles": "c1ccccc1"}),
    (3, 200, {"molecule_id": 3, "smiles": "CC(=O)O"})
])
def test_get_molecule(molecule_id, expected_status, expected_response):
    response = client.get(f"/molecule/{molecule_id}")
    assert response.status_code == expected_status
    assert response.json() == expected_response
    


### Testing get_molecule Function
@pytest.mark.parametrize("molecule_id,smiles,expected_status", [
    (1, "C1=CC=CC=C1", 200),                      # Update with valid benzene SMILES
    (2, "CCO", 200),                              # Update with simple ethanol SMILES
    (3, "CC(C)C1=CC2=C(C=C1)C(=O)C=C2", 200),     # Update with complex SMILES
    (4, "C[C@H](O)[C@@H](N)C(=O)O", 200)          # Update with stereochemistry SMILES
])
def test_update_valid_smiles(molecule_id, smiles, expected_status):
    response = client.put(f"/molecule/{molecule_id}", json={"smiles": smiles})
    assert response.status_code == expected_status
    assert response.json() == {"molecule_id": molecule_id, "smiles": smiles}
    assert molecules[molecule_id]["smiles"] == smiles



## Testing delete_molecule Function

def test_delete_existing_molecule():
    assert 1 in molecules
    response = client.delete("/molecule/1")
    assert response.status_code == 204
    assert 1 not in molecules



## Testing get_molecules Function

def test_get_molecules(reset_molecules):
    # notice that with reset_molecules now our In-memory storage is empty

    # adding few molecules before printing them out
    client.post("/add", json={"molecule_id": 1, "smiles": "CCO"})
    client.post("/add", json={"molecule_id": 2, "smiles": "C[C@H](O)[C@@H](N)C(=O)O"})

    response = client.get("/molecules")
    assert response.status_code == 200
    assert response.json() == [
        {"molecule_id": 1, "smiles": "CCO"},
        {"molecule_id": 2, "smiles": "C[C@H](O)[C@@H](N)C(=O)O"}
    ]



## Testing Substructue Search Function

def test_substructure_search(reset_molecules):
    # adding few molecules to the database
    client.post("/add", json={"molecule_id": 1, "smiles": "c1ccccc1"})
    client.post("/add", json={"molecule_id": 2, "smiles": "O=[N+]([O-])c1nccn1COCC=CCO"})
    client.post("/add", json={"molecule_id": 3, "smiles": "OCCN1CCN(CCCOn2nnc3ccccc32)CC1"})


    response = client.post("/search", json={"substructure": "CCO"})
    assert response.status_code == 200
    assert response.json() == {
        "message": "Matching molecules found",
        "molecules": [
            {"molecule_id": 2, "smiles": "O=[N+]([O-])c1nccn1COCC=CCO"},
            {"molecule_id": 3, "smiles": "OCCN1CCN(CCCOn2nnc3ccccc32)CC1"}
        ]
    }
    

def test_substructure_search_no_matches(reset_molecules):
    # adding few molecules to the database
    client.post("/add", json={"molecule_id": 1, "smiles": "S=C=NCCCCCCCCc1ccccc1"})
    client.post("/add", json={"molecule_id": 2, "smiles": "S=c1[nH]nc(Cn2ccc3ccccc32)n1-c1ccccc1"})
    client.post("/add", json={"molecule_id": 3, "smiles": "OCCN1CCN(CCCOn2nnc3ccccc32)CC1"})

    response = client.post("/search", json={"substructure": "O=c1cc(CCl)occ1O"})
    response.status_code == 200
    assert response.json() == {
        "message": "No matching molecules found",
        "molecules": None
    }



## Testing upload_file Function

def test_upload_file():
    # This csv file 
    file_path = 'valid_smiles_dataset.csv'


    with open(file_path, 'rb') as file:
        response = client.post("/uploadfile", files={"file": file})

    assert response.status_code == 200
    assert response.json() == {"message": "Successfully uploaded 3 molecules"}



# Testing get_server Function

def test_get_server():
    os.environ["SERVER_ID"] = "SERVER-123"  
    response = client.get("/server")
    assert response.status_code == 200
    assert response.json() == {"server_id": "SERVER-123"} 

    del os.environ["SERVER_ID"]  
    response = client.get("/server")
    assert response.status_code == 200
    assert response.json() == {"server_id": "1"}

