#!/usr/bin/env bash

## Location of where the pipeline is installed
DIR="/home/a1645424/al/alastair-phd/orthologue-analysis/oma-16-sample"
MSA="${DIR}/msa/nine-sample/nucleotide_msa"
TREE="${DIR}/trees/nine-sample-relax.nwk"
OUT="${DIR}/selection/nine-sample"

mkdir -p ${OUT}

CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

conda activate /home/a1645424/al/nf-condaEnvs/nf-workflows-ebc7570a00a9f01eadb4bbef4dae4640

for FILE in ${MSA}/*.fasta; do
    BN=$(basename ${FILE%_*})
    printf "Current File: %s\n\n" "${BN}"
    hyphy relax \
        --alignment ${FILE} \
        --tree ${TREE} \
        --model minimal \
        --test Foreground \
        --reference Background \
        --output "${OUT}/RELAX-${BN}-nine-sample-relax.json" >"${OUT}/${BN}-stderr.txt"
    printf "\n"
done

conda deactivate
