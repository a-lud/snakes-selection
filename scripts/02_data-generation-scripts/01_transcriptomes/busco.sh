#!/usr/bin/env bash

# ------------------------------------------------------------------------------------------------ #
# Directories
PROJ='/home/a1645424/al/00_phd/genome-data/transcriptomes'
PROT='/home/a1645424/al/00_phd/genome-data/transcriptomes/completeORF'
OUT=${PROJ}/busco

mkdir -p ${OUT}

# ------------------------------------------------------------------------------------------------ #
# BUSCO assessment
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

conda activate BIOTOOLS_env

for sample in ${PROT}/*.pep; do

    BN=$(basename ${sample%.*})

    # Aipysurus
    busco \
        -i ${sample} \
        -o ${BN}_tetrapoda \
        -m prot \
        -l tetrapoda_odb10 \
        --cpu 20 \
        --out_path ${OUT}

    busco \
        -i ${sample} \
        -o ${BN}_metazoa \
        -m prot \
        -l metazoa_odb10 \
        --cpu 20 \
        --out_path ${OUT}
done

conda deactivate