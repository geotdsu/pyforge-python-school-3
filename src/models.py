from sqlalchemy.orm import Mapped
from sqlalchemy import Column, DateTime, func
from src.database import Base, str_uniq, int_pk

class Drug(Base):
    __tablename__ = "drugs"

    id: Mapped[int_pk]
    name: Mapped[str_uniq]
    smiles: Mapped[str_uniq]
    created_at: Mapped[DateTime] = Column(DateTime, server_default=func.now())
    updated_at: Mapped[DateTime] = Column(DateTime, server_default=func.now(), onupdate=func.now())

    def __str__(self):
        return (
            f"{self.__class__.__name__}(id={self.id}, "
            f"name={self.name!r},"
            f"smiles={self.smiles!r},"
            f"created_at={self.created_at!r},"
            f"updated_at={self.updated_at!r})"
        )

    def __repr__(self):
        return str(self)
