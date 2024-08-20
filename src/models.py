from sqlalchemy import Column, String, Integer, TIMESTAMP, func
from sqlalchemy.sql.expression import text

from .database import Base



class Drug(Base):
    __tablename__ = "drugs"

    id = Column(Integer, primary_key=True, nullable=False)
    name = Column(String, nullable=False)
    smiles = Column(String, nullable=False)
    created_at = Column(TIMESTAMP, server_default=func.now(), nullable=False)
    updated_at = Column(TIMESTAMP, server_default=func.now(), onupdate=func.now(), nullable=False)
    
