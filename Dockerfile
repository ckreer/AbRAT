# Base-Image miniconda
FROM continuumio/miniconda3

LABEL org.opencontainers.image.source="https://github.com/ckreer/AbRAT"
LABEL org.opencontainers.image.description="AbRAT – The Antibody Repertoire Analysis Toolkit"
LABEL org.opencontainers.image.licenses="GPL-3.0"
ARG VERSION=v1.0.1
LABEL org.opencontainers.image.version="${VERSION}"
LABEL org.opencontainers.image.documentation="https://abrat.readthedocs.io"
LABEL org.opencontainers.image.authors="Christoph Kreer"

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
ENV PATH=/opt/conda/envs/abrat_env/bin:$PATH
RUN conda clean -afy

# Set ENV variables
ENV IGDATA=/opt/conda/envs/abrat_env/share/igblast
ENV BLASTDB=/app/data/database/blastdb:/app/data/database/igblastdb

# Copy setup.py and MANIFEST.in
COPY setup.py MANIFEST.in /app/

# streamlit config
COPY .streamlit /root/.streamlit
ENV HOME=/root

# Copy package
COPY abrat /app/abrat

# Install abrat package
RUN pip install . --no-cache-dir

# Update igblast db with most recent data
RUN python3 /app/abrat/core/build_igblast_db.py --output /app/data/database/igblastdb

# Streamlit Port
EXPOSE 8501

# Check
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health || exit 1

# Define Entrypoint
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "abrat_env", "streamlit", "run", "abrat/gui/abrat_app.py"]

# Standardarguments
CMD ["--server.port=8501", "--server.address=0.0.0.0"]