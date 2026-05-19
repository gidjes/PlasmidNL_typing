#!/usr/bin/env bash
set -euo pipefail

echo "Creating conda environments..."

# Main workflow environment
conda env create -f environment.yml

# CARD / RGI
conda create -y -n CARD \
    -c conda-forge \
    -c bioconda \
    rgi

# MOB-suite
conda create -y -n mob_suite \
    -c conda-forge \
    -c bioconda \
    mob_suite

# MGE cluster
conda create -y -n mge_cluster \
    -c conda-forge \
    -c bioconda \
    mge-cluster

echo "All environments have been created!"

# Databases
echo "Downloading databases"
conda run -n plasmidnl_typing download-db.sh
conda run -n plasmidnl_typing amrfinder -u
cd data/databases
git clone https://bitbucket.org/genomicepidemiology/resfinder_db/
cd ../../

# Install package into environment
conda run -n plasmidnl_typing pip install -e .

echo "Setup complete."