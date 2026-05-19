#!/bin/bash

# Read out arguments
plasmid_name=$1
in_dir=$2
use_custom=$3

# Check if output directories exists
mkdir -p output/${1}/mge_cluster/

# Set up mge input file
echo ${in_dir}/${1}.fasta > output/${1}/mge_cluster/${1}_path.txt 

# Determine scheme directory
if [ "$use_custom" = "1" ]; then
    model_folder="data/mge_model_custom"
else
    model_folder="data/mge_model_CPO"
fi

# Run mge-cluster
mge_cluster --existing --input output/${1}/mge_cluster/${1}_path.txt  --model_folder "$model_folder" --model_prefix mge-cluster --prefix ${1} --outdir output/${1}/mge_cluster/

echo mge_cluster --existing --input output/${1}/mge_cluster/${1}_path.txt  --model_folder "$model_folder" --model_prefix mge-cluster --prefix ${1} --outdir output/${1}/mge_cluster/


# If mge-cluster fails (e.g. no matching unitigs in fasta) generate NA output file
if [ $? -ne 0 ]; then
    if [[ ! -f "$output_file" ]]; then
        echo "Error: Output file $output_file does not exist."
        exit 1
    fi
    if awk -F'\t' '{ if ($2 != 0) exit 1 }' "$output_file"; then
        echo "tsne1D,tsne2D,Standard_Cluster,Sample_Name" > output/${1}/mge_cluster/${1}_prediction.csv
        echo "-,-,-,${1}" >> output/${1}/mge_cluster/${1}_prediction.csv
    else
        echo "Error: Not all entries in the second column are 0."
        exit 1
    fi
fi
