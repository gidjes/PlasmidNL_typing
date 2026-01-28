#!/usr/bin/env bash
set -euo pipefail

echo "Creating conda environments..."

# environments
conda env create -f plasmidnl_typing.yml
conda env create -n CARD rgi
conda env create -n mob_suite mob_suite
conda env create -n mge_cluster mge_cluster

echo "All environments have been created!"

echo "Downloading databases"
conda run -n plasmidnl_typing download-db.sh
conda run -n plasmidnl_typing amrfinder -u
cd data/databases
git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
cd ../../

conda activate plasmidnl_typing
pip install -e .