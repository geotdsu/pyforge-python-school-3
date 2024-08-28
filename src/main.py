from fastapi import FastAPI
from .routers import drug
from .models import Base
from .database import engine
from .logger import logger


Base.metadata.create_all(bind=engine)

app = FastAPI()

app.include_router(drug.router)

logger.info("Application started")
