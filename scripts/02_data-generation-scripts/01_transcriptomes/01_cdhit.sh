#!/usr/bin/env bash

GENOMEDATA='/home/a1645424/al/00_phd/genome-data/transcriptomes'
CDHITOUT="${GENOMEDATA}/cdhit"

mkdir -p ${CDHITOUT}

CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

conda activate TRINITY_NEW

for transcripts in "${GENOMEDATA}"/*fasta.gz; do

    BN=$(basename "${transcripts%%.*}")

    if [[ ! -f ${CDHITOUT}/${BN}.cdhit.fasta ]];then
        cd-hit-est \
            -o ${CDHITOUT}/${BN}.cdhit.fasta \
            -M 2000M \
            -c 0.95 \
            -i "${transcripts}" \
            -p 1 \
            -d 0 \
            -b 3 \
            -T 20
    fi
done

conda deactivate

