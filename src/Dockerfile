# Use Miniconda base image
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy requirements file
COPY requirements.txt .

# Install RDKit and clean conda cache
RUN conda install conda-forge::rdkit && \
    conda clean -afy

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy project files
COPY . /app

# Expose port 8000
EXPOSE 8000

# Run the app with uvicorn
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]