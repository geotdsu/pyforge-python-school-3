from typing import List
from fastapi import FastAPI, HTTPException, status, UploadFile, File
from rdkit import Chem
from . import schemas
import csv
from os import getenv


tags_metadata = [
    {"name": "Root", "description": "Root endpoint."},
    {"name": "Molecule Management", "description": "Manage molecules: add, update, retrieve, and delete."},
    {"name": "Substructure Search", "description": "Search for molecules containing a substructure."},
    {"name": "Upload a file", "description": "Upload a CSV file with molecule SMILES."}
]

app = FastAPI(openapi_tags=tags_metadata)



# In-memory storage for molecules
molecules = {
    1: {"smiles": "CCO"},
    2: {"smiles": "c1ccccc1"},
    3: {"smiles": "CC(=O)O"},
    4: {"smiles": "CC(=O)Oc1ccccc1C(=O)O"}
}


@app.get("/", status_code=status.HTTP_200_OK, tags=["Root"])
def root():
    """Root endpoint for the application"""
    return {"message": "Molecule Fast API application"}


@app.get("/server", tags=["Root"])
def get_server():
    """Check which server is responding"""
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/add", status_code=status.HTTP_201_CREATED, response_model=schemas.Molecule, tags=["Molecule Management"])
def add_molecule(molecule: schemas.Molecule):
    """Add a new molecule to the database

    **Request Body:**
    - `molecule_id` (int): Unique identifier for the molecule.
    - `smiles` (str): SMILES representation of the molecule.

    **Response:**
    - Returns the added molecule with the same `molecule_id` and `smiles`.

    **Exceptions:**
    - `400 Bad Request`: If the molecule ID already exists.
    - `400 Bad Request`: If the SMILES string is empty or invalid.
    """

    # Check if the molecule ID already exists
    if molecule.molecule_id in molecules:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Molecule with the ID {molecule.molecule_id} already exists")
    
    # Check if SMILES string is provided and valid
    if not molecule.smiles:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="SMILES string cannot be empty")
    
    molecule_object = Chem.MolFromSmiles(molecule.smiles)
    if molecule_object is None:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid SMILES string")
    
    # Store the new molecule
    molecules[molecule.molecule_id] = {"smiles": molecule.smiles}
    return molecule

@app.get("/molecule/{molecule_id}", response_model=schemas.Molecule, tags=["Molecule Management"])
def get_molecule(molecule_id: int):
    """Retrieve a molecule by its ID

    **Path Parameter:**
    - `molecule_id` (int): Unique identifier for the molecule.

    **Response:**
    - Returns the molecule with the specified `molecule_id`.

    **Exceptions:**
    - `404 Not Found`: If no molecule with the given `molecule_id` is found.
    """
    
    # Check if the molecule ID exists
    if molecule_id not in molecules:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Molecule with the ID {molecule_id} is not found")
    
    # Return the molecule details
    return {"molecule_id": molecule_id, "smiles": molecules[molecule_id]["smiles"]}


@app.put("/molecule/{molecule_id}", response_model=schemas.Molecule, tags=["Molecule Management"])
def update_molecule(molecule_id: int, molecule: schemas.UpdateMolecule):
    """Update a molecule by its ID

    **Path Parameter:**
    - `molecule_id` (int): Unique identifier for the molecule.

    **Request Body:**
    - `smiles` (str): New SMILES representation of the molecule.

    **Response:**
    - Returns the updated molecule with the same `molecule_id` and new `smiles`.

    **Exceptions:**
    - `400 Bad Request`: If the SMILES string is empty or invalid.
    - `404 Not Found`: If no molecule with the given `molecule_id` is found.
    """

    # Check if the molecule ID exists
    if molecule_id not in molecules:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Molecule with the ID {molecule_id} is not found")
    
    # Validate the new SMILES string
    if not molecule.smiles:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="SMILES string cannot be empty")
    
    molecule_object = Chem.MolFromSmiles(molecule.smiles)
    if molecule_object is None:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid SMILES string")
    
    # Update the molecule information
    molecules[molecule_id]["smiles"] = molecule.smiles
    return {"molecule_id": molecule_id, "smiles": molecule.smiles}


@app.delete("/molecule/{molecule_id}", status_code=status.HTTP_204_NO_CONTENT, tags=["Molecule Management"])
def delete_molecule(molecule_id: int):
    """Delete a molecule by its ID

    **Path Parameter:**
    - `molecule_id` (int): Unique identifier for the molecule.

    **Exceptions:**
    - `404 Not Found`: If no molecule with the given `molecule_id` is found.
    """

    # Check if the molecule ID exists
    if molecule_id not in molecules:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Molecule with the ID {molecule_id} is not found")
    
    # Remove the molecule
    del molecules[molecule_id]
    return


@app.get("/molecules", response_model=List[schemas.Molecule] ,tags=["Molecule Management"])
def get_molecules():
    """List all molecules currently stored"""
    return [schemas.Molecule(molecule_id=id, smiles=info["smiles"]) for id, info in molecules.items()]


@app.post("/search", status_code=status.HTTP_200_OK, response_model=schemas.SearchResponse, tags=["Substructure Search"])
def substructure_search(query: schemas.SubstructureQuery):
    """Perform substructure search using a query molecule

    **Request Body:**
    - `substructure` (str): The SMILES string representing the substructure to search for within existing molecules.

    **Response:**
    - Returns a message indicating the result of the search and, if applicable, the list of matching molecules.

    **Exceptions:**
    - `400 Bad Request`: If the `substructure` query is empty or invalid.
    """


    # Validate the substructure query
    if not query.substructure:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Substructure query cannot be empty.")
    
    substructure_molecule = Chem.MolFromSmiles(query.substructure)
    if substructure_molecule is None:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid substructure SMILES.")
    
    # Find molecules matching the substructure
    matching_molecules = [
        schemas.Molecule(molecule_id=id, smiles=info["smiles"])
        for id, info in molecules.items()
        if Chem.MolFromSmiles(info["smiles"]) and Chem.MolFromSmiles(info["smiles"]).HasSubstructMatch(substructure_molecule)
    ]
    
    if not matching_molecules:
        return {"message": "No matching molecules found", "molecules": None}
    
    return {"message": "Matching molecules found", "molecules": matching_molecules}


@app.post("/uploadfile", status_code=status.HTTP_200_OK, tags=["Upload a file"])
async def upload_file(file: UploadFile = File(...)):
    """Upload a CSV file containing molecule data
    
    **Request Body:**
    - `file` (UploadFile): A CSV file containing molecule SMILES strings. The CSV file must have a column named 'SMILES' or 'smiles'.

    **Response:**
    - Returns a message indicating the number of successfully uploaded molecules and any errors encountered.

    **Exceptions:**
    - `400 Bad Request`: If the file format is not CSV.
    - `400 Bad Request`: If the file is empty or exceeds the size limit.
    - `400 Bad Request`: If the CSV file does not contain a 'SMILES' or 'smiles' column.
    - `400 Bad Request`: If there are invalid or missing SMILES strings in the CSV file.
    - `500 Internal Server Error`: For any unexpected errors.

    **Details:**
    - The CSV file should contain a column named 'SMILES' or 'smiles' with molecule SMILES strings. Rows with invalid or missing SMILES strings will be skipped, and a summary of errors will be provided in the response.

    **Example Dataset:**
    - For an example CSV file with SMILES data, you can download the dataset from [this Kaggle link](https://www.kaggle.com/datasets/yanmaksi/big-molecules-smiles-dataset).
    """
    
    # Maximum allowed file size (100 MB)
    MAX_FILE_SIZE = 100 * 1024 * 1024
    
    
    # Validate the file format
    if file.content_type != "text/csv":
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid file format. Only CSV files are supported") 
    
    # Check if a file was selected
    if file.filename == "":
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail="No file selected for upload") 
    
    # Read the file content
    contents = await file.read()

    # Check if the file is empty
    if not contents:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Uploaded file is empty")   
    
    # Check if the file exceeds the size limit
    if len(contents) > MAX_FILE_SIZE:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="File is too large. Maximum size is 100 MB")    
    
    try:
        # Decode the contents and create a CSV reader
        reader = csv.DictReader(contents.decode().splitlines())

        # Ensure the CSV contains a 'SMILES', 'smiles', or 'Smiles' header
        valid_headers = {"SMILES", "smiles", "Smiles"}
        if not any(header in reader.fieldnames for header in valid_headers):
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="CSV file must contain 'SMILES', 'smiles', or 'Smiles' header.")
        
        count = 0

        # Process each row in the CSV file
        for index, row in enumerate(reader, start=len(molecules) + 1):
            # Retrieve the SMILES value from the row
            smiles_value = row.get("SMILES") or row.get("smiles") or row.get("Smiles")

            # Validate the SMILES value
            if not smiles_value:
                raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Missing 'SMILES' or 'smiles' or 'Smiles' value in row {index}")
            
            # Validate the SMILES string using RDKit
            molecule_object = Chem.MolFromSmiles(smiles_value)
            if molecule_object is None:
                raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Invalid SMILES string in row {index}: {smiles_value}")
            
            # Store the valid molecule
            molecules[index] = {"smiles": smiles_value}
            count += 1

        # Return a success message with the number of uploaded molecules
        return {"message": f"Successfully uploaded {count} molecules"}
    
    except HTTPException:
        # Re-raise HTTP exceptions to return proper status codes and messages
        raise
    except Exception as e:
        # Handle any unexpected errors
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Unexpected error: {str(e)}")

