#!/usr/bin/env bash

COMPLETEORF="/Users/alastair.ludington/Documents/phd/00_papers/snakes-selection-hogs/data/oma-16-samples/oma-cds"

for fasta in ${COMPLETEORF}/*.cds; do
    BN=$(basename ${fasta%.*})
    if [[ ! -f ${COMPLETEORF}/${BN}.gc.stats ]]; then
        printf "GC-Content: %s\n" ${BN}
        seqkit fx2tab \
            --length \
            --name \
            --gc \
            --gc-skew \
            --threads 2 \
            --out-file ${COMPLETEORF}/${BN}.gc.stats \
            ${fasta}
    fi
done
