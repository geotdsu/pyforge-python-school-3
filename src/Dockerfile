FROM continuumio/miniconda3

WORKDIR /app

COPY requirements.txt .

RUN conda install conda-forge::rdkit && \
    conda clean -afy

RUN pip install --no-cache-dir -r requirements.txt

COPY . /app

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
