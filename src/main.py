from fastapi import FastAPI
from .routers import drug


app = FastAPI()

app.include_router(drug.router) 
