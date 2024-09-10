import datetime
from typing import Annotated
from sqlalchemy import create_engine, func
from sqlalchemy.orm import sessionmaker, DeclarativeBase, mapped_column
from src.config import get_db_url

DATABASE_URL = get_db_url()

engine = create_engine(DATABASE_URL)

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


int_pk = Annotated[int, mapped_column(primary_key=True)]
created_at = Annotated[datetime.datetime,
                       mapped_column(server_default=func.now())]
updated_at = Annotated[
    datetime.datetime, mapped_column(
        server_default=func.now(), onupdate=datetime.datetime.now)
]
str_uniq = Annotated[str, mapped_column(unique=True, nullable=False)]
str_null_true = Annotated[str, mapped_column(nullable=True)]


class Base(DeclarativeBase):
    pass


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
