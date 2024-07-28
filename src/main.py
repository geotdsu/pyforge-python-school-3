from typing import List, Optional
from fastapi import FastAPI, HTTPException, status, UploadFile, File
from pydantic import BaseModel
from rdkit import Chem
import csv
from os import getenv


tags_metadata = [
    {
        "name": "Root",
        "description": "Root endpoint.",
    },
    {
        "name": "Molecule Management",
        "description": "Manage molecules: add, update, retrieve, and delete.",
    },
    {
        "name": "Substructure Search",
        "description": "Search for molecules containing a substructure.",
    },
    {
        "name": "Upload a file",
        "description": "Upload a CSV file with molecule SMILES.",
    }
]

app = FastAPI(openapi_tags=tags_metadata)


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


molecules = [
    {"molecule_id": 1, "smiles": "CCO"},
    {"molecule_id": 2, "smiles": "c1ccccc1"},
    {"molecule_id": 3, "smiles": "CC(=O)O"},
    {"molecule_id": 4, "smiles": "CC(=O)Oc1ccccc1C(=O)O"}
]


@app.get("/", tags=["Root"])
def root():
    """Root endpoint for the application"""
    return {"message": "Molecule Fast API application"}


@app.get("/server", tags=["Root"])
def get_server():
    """Check which server is responding"""
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/add", status_code=status.HTTP_201_CREATED, response_model=Molecule, tags=["Molecule Management"])
def add_molecule(molecule: Molecule):
    """Add a new molecule to the database
    
    **Request Body:**
    - `molecule_id` (int): Unique identifier for the molecule.
    - `smiles` (str): SMILES representation of the molecule.
    
    **Response:**
    - Returns the added molecule with the same `molecule_id` and `smiles`.

    **Example Response:**
    ```json
    {
        "molecule_id": 5,
        "smiles": "CCN"
    }
    ```

    *Exceptions:**
    - `400 Bad Request`: If the molecule ID already exists.
    - `400 Bad Request`: If the SMILES string is empty.
    - `400 Bad Request`: If the SMILES string is invalid.
    - `500 Internal Server Error`: For any unexpected errors.
    """

    try:
        for mol in molecules:
            if mol["molecule_id"] == molecule.molecule_id:
                raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Molecule with the ID {molecule.molecule_id} already exists")
        
        if not molecule.smiles:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="SMILES string cannot be empty")
        
        molecule_object = Chem.MolFromSmiles(molecule.smiles)
        if molecule_object is None:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid SMILES string")
        
        molecules.append(molecule.model_dump())
        return molecule
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Unexpected error: {str(e)}")


@app.get("/molecule/{molecule_id}", response_model=Molecule, tags=["Molecule Management"])
def get_molecule(molecule_id: int):
    """Retrieve a molecule by its ID
    
    **Path Parameter:**
    - `molecule_id` (int): Unique identifier for the molecule.
    
    **Response:**
    - Returns the molecule with the specified `molecule_id`.
    
    **Example Response:**
    ```json
    {
        "molecule_id": 2,
        "smiles": "c1ccccc1"
    }
    ```
    
    **Exceptions:**
    - `400 Bad Request`: If `molecule_id` is not a valid integer.
    - `404 Not Found`: If no molecule with the given `molecule_id` is found.
    - `500 Internal Server Error`: For any unexpected errors.
    """
    
    try:
        if not isinstance(molecule_id, int):
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid molecule ID type. It must be an integer.")

        for mol in molecules:
            if mol["molecule_id"] == molecule_id:
                return mol
        
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Molecule with the ID {molecule_id} is not found")
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Unexpected error: {str(e)}")


@app.put("/molecule/{molecule_id}", response_model=Molecule, tags=["Molecule Management"])
def update_molecule(molecule_id: int, molecule: UpdateMolecule):
    """Update a molecule by its ID
    
    **Path Parameter:**
    - `molecule_id` (int): Unique identifier for the molecule.
    
    **Request Body:**
    - `smiles` (str): New SMILES representation of the molecule.
    
    **Response:**
    - Returns the updated molecule with the same `molecule_id` and new `smiles`.
    
    **Example Response:**
    ```json
    {
        "molecule_id": 1,
        "smiles": "CCN"
    }
    ```

    **Exceptions:**
    - `400 Bad Request`: If `molecule_id` is not a valid integer.
    - `400 Bad Request`: If the SMILES string is invalid.
    - `404 Not Found`: If no molecule with the given `molecule_id` is found.
    - `500 Internal Server Error`: For any unexpected errors.
    """

    try:
        if not isinstance(molecule_id, int):
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid molecule ID type. It must be an integer.")
        
        
        for index, mol in enumerate(molecules):
            if mol["molecule_id"] == molecule_id:
                if not molecule.smiles:
                    raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="SMILES string cannot be empty.")
                
                molecule_object = Chem.MolFromSmiles(molecule.smiles)
                if molecule_object is None:
                    raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid SMILES string.")
                
                molecules[index]["smiles"] = molecule.smiles
                return {"molecule_id": molecule_id, "smiles": molecule.smiles}
        
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Molecule with ID {molecule_id} is not found")
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Unexpected error: {str(e)}")


@app.delete("/molecule/{molecule_id}", status_code=status.HTTP_204_NO_CONTENT, tags=["Molecule Management"])
def delete_molecule(molecule_id: int):
    """Delete a molecule by its ID

    **Path Parameter:**
    - `molecule_id` (int): Unique identifier for the molecule.
    
    **Exceptions:**
    - `400 Bad Request`: If `molecule_id` is not a valid integer.
    - `404 Not Found`: If no molecule with the given `molecule_id` is found.
    - `500 Internal Server Error`: For any unexpected errors.
    """

    try:
        if not isinstance(molecule_id, int):
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid molecule ID type. It must be an integer.")
        
        for index, mol in enumerate(molecules):
            if mol["molecule_id"] == molecule_id:
                molecules.pop(index)
                return
        
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Molecule with ID {molecule_id} is not found")
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Unexpected error: {str(e)}")


@app.get("/molecules", response_model=List[Molecule] ,tags=["Molecule Management"])
def get_molecules():
    """List all molecules currently stored"""
    return molecules


@app.post("/search", status_code=status.HTTP_200_OK, response_model=SearchResponse, tags=["Substructure Search"])
def substructure_search(query: SubstructureQuery):
    """Perform substructure search using a query molecule
    
    **Request Body:**
    - `substructure` (str): The SMILES string representing the substructure to search for within existing molecules.
    
    **Example Response (Matches Found):**
    ```json
    {
        "message": "Matching molecules found",
        "molecules": [
            {
                "id": 1,
                "name": "Aspirin",
                "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
            },
            {
                "id": 2,
                "name": "Ibuprofen",
                "smiles": "CC(C)CC1=CC2=C(C=C1)C(=O)O2"
            }
        ]
    }
    ```

    **Example Response (No Matches Found):**
    ```json
    {
        "message": "No matching molecules found",
        "molecules": null
    }
    ```

    **Exceptions:**
    - `400 Bad Request`: If the `substructure` query is empty.
    - `400 Bad Request`: If the `substructure` SMILES string is invalid.
    - `500 Internal Server Error`: For any unexpected errors during the search process.
    """

    if not query.substructure:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Substructure query cannot be empty.")
    
    try:
        substructure_molecule = Chem.MolFromSmiles(query.substructure)
        if substructure_molecule is None:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid substructure SMILES.")
        
        matching_molecules = []
        for mol in molecules:
                molecule = Chem.MolFromSmiles(mol["smiles"])
                if molecule and molecule.HasSubstructMatch(substructure_molecule):
                    matching_molecules.append(mol)
        
        if not matching_molecules:
            return {"message": "No matching molecules found", "molecules": None}
        
        return {"message": "Matching molecules found", "molecules": matching_molecules}
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Unexpected error: {str(e)}")


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
    
    # File size limit (10 MB)
    MAX_FILE_SIZE = 10 * 1024 * 1024
    
    if file.content_type != "text/csv":
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid file format. Only CSV files are supported.")
    
    if file.filename == "":
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="No file selected for upload.")
    
    contents = await file.read()
    if not contents:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Uploaded file is empty.")
    
    if len(contents) > MAX_FILE_SIZE:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="File is too large. Maximum size is 10 MB.")
    
    try:
        reader = csv.DictReader(contents.decode().splitlines())
        if "SMILES" not in reader.fieldnames and "smiles" not in reader.fieldnames:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="CSV file must contain 'SMILES' or 'smiles' header.")
        
        new_molecules = []
        errors = []
        for index, row in enumerate(reader, start=len(molecules) + 1):
            smiles_value = row.get("SMILES") or row.get("smiles")
            if not smiles_value:
                errors.append(f"Missing 'SMILES' or 'smiles' value in row {index}.")
                continue
            
            try:
                molecule_object = Chem.MolFromSmiles(smiles_value)
                if molecule_object is None:
                    errors.append(f"Invalid SMILES string in row {index}: {smiles_value}.")
                    continue
            except Exception as e:
                errors.append(f"Error processing SMILES in row {index}: {str(e)}.")
                continue
            
            molecule = {
                "molecule_id": index,
                "smiles": smiles_value
            }
            new_molecules.append(molecule)
        
        molecules.extend(new_molecules)
        if errors:
            return {"message": f"Successfully uploaded {len(new_molecules)} molecules, with errors in the following rows: {errors}"}
        return {"message": f"Successfully uploaded {len(new_molecules)} molecules"}
    
    except csv.Error as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Error processing CSV file: {str(e)}")
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Unexpected error: {str(e)}")