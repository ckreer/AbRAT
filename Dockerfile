# Basis-Image mit Miniconda verwenden und Plattform festlegen (optional für Apple Silicon)
FROM continuumio/miniconda3
# --platform=linux/amd64 continuumio/miniconda3

# Define Workdir in Container
WORKDIR /app

# Install Linux-packages
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    unzip \
    perl \
    tar \
    && rm -rf /var/lib/apt/lists/*

# Create and activate conda environemnt
COPY environment.yml .
RUN conda env create -f environment.yml
ENV PATH /opt/conda/envs/abrat_env/bin:$PATH

# Set ENV variables
ENV IGDATA=/opt/conda/envs/abrat_env/share/igblast
ENV BLASTDB=/app/data/database/blastdb:/app/data/database/igblastdb

# TODO: Remove before production!
ENV PYTHONDONTWRITEBYTECODE=1

# Copy setup.py und MANIFEST.in
COPY setup.py MANIFEST.in /app/

# Copy den gesamten package-Ordner
COPY abrat /app/abrat

# Projekt als Paket installieren (damit "import abrat" klappt)
# TODO: -e für production entfernen!
RUN pip install -e . --no-cache-dir

# Update igblast db with most recent data
RUN python3 /app/abrat/core/build_igblast_db.py --output /app/data/database/igblastdb

# Streamlit Port freigeben
EXPOSE 8501

# Check
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1

# Define Entrypoint
ENTRYPOINT ["streamlit", "run", "abrat/gui/abrat_app.py"]

# Standardarguments
CMD ["--server.port=8501", "--server.address=0.0.0.0"]
# ["python", "testfile.py"]