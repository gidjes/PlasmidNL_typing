#!/usr/bin/env bash
set -euo pipefail

echo "Creating conda environments..."

# Main workflow environment
conda env create -f environment.yml

# CARD / RGI
conda create -y -n CARD \
    -c conda-forge \
    -c bioconda \
    rgi=6.0.5

# MOB-suite
conda create -y -n mob_suite \
    -c conda-forge \
    -c bioconda \
    mob_suite=3.1.9

# MGE cluster
conda create -y -n mge_cluster \
    -c conda-forge \
    -c bioconda \
    mge-cluster=1.1.0

echo "All environments have been created!"

# Databases
echo "Downloading databases"

# AMRFinderPlus database
AMRFINDER_DB_VERSION="2026-03-24.1"
AMRFINDER_DB_URL="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/4.2/${AMRFINDER_DB_VERSION}"
AMRFINDER_DB_DIR="$(conda run -n plasmidnl_typing bash -c 'echo $CONDA_PREFIX')/share/amrfinderplus/data"

echo "Installing AMRFinderPlus database"

mkdir -p "${AMRFINDER_DB_DIR}/${AMRFINDER_DB_VERSION}"

wget -q -r -np -nH --cut-dirs=8 \
    --reject="index.html*" \
    "${AMRFINDER_DB_URL}/" \
    -P "${AMRFINDER_DB_DIR}/${AMRFINDER_DB_VERSION}"

ln -sfn "${AMRFINDER_DB_VERSION}" "${AMRFINDER_DB_DIR}/latest"

# Check
conda run -n plasmidnl_typing amrfinder -V

# Resfinder database
echo "Installing ResFinder database"
cd data/databases
git clone --branch resfinder-4.7.2 --depth 1 \
    https://bitbucket.org/genomicepidemiology/resfinder_db/

# Resfinder database
echo "Installing PlasmidFinder database"
git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git
cd plasmidfinder_db
python3 INSTALL.py 
cd ../../../

# Done!
echo "Setup complete."