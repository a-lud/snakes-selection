#!/usr/bin/env bash

# ------------------------------------------------------------------------------------------------ #
# Conda parameters
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate KRAKEN2_env

# ------------------------------------------------------------------------------------------------ #
# Directories
DIR='/home/a1645424/al/00_phd/orthologue-analysis/oma-16-sample/hogs-gain-loss'
DIRS=$(find ${DIR} -type d -name 'ancestral-node-*' | tr '\n' ' ')
KRAKENDB='/home/a1645424/al/databases/k2_standard_20210517'

for ANC in ${DIRS}; do
    HOGS=$(cat ${ANC}/blast-gained-samples-threshold-no-hit.txt | tr '\n' ' ')
    HOGSEQS="${ANC}/samples/gained_threshold_fasta"
    OUT="${ANC}/kraken2-check"

    mkdir -p ${OUT}

    # -------------------------------------------------------------------------------------------- #
    # Concatenate HOG fasta files into master one
    if [[ ! -f ${ANC}/concatentated-gained-samples-threshold-no-hits.fa || ! -s ${ANC}/concatentated-gained-samples-threshold-no-hits.fa ]]; then
        printf "Concatentating FASTA files\n"

        # Append each HOG to concatenated file
        for H in ${HOGS}; do
            cat "${HOGSEQS}/${H}.fa" >>${ANC}/concatentated-gained-samples-threshold-no-hits.fa
        done
    fi

    # -------------------------------------------------------------------------------------------- #
    # Search against Kraken2
    kraken2 \
        --db ${KRAKENDB} \
        --threads 10 \
        --unclassified-out ${OUT}/unclassified.out \
        --classified-out ${OUT}/classified.out \
        --output ${OUT}/output.txt \
        --report ${OUT}/report.out \
        --use-names \
        ${ANC}/concatentated-gained-samples-threshold-no-hits.fa
done

conda deactivate
