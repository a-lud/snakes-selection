#!/usr/bin/env bash

DIR='/Users/alastair.ludington/Documents/phd/00_papers/snakes-selection-hogs/data/HOG-data/hogs-gain-loss'
DIRS=$(find ${DIR} -type d -name 'ancestral-node-*' | tr '\n' ' ')

for ANC in ${DIRS}; do
    HOGS=$(cat ${ANC}/blast-gained-samples-threshold-no-hit.txt | tr '\n' ' ')
    HOGSEQS="${ANC}/samples/gained_threshold_fasta"
    OUT=${ANC}/distance-matrix
    mkdir -p ${OUT}

    # Remove the file if it exists already
    if [[ -f ${ANC}/concatentated-gained-samples-threshold-no-hits.fa ]]; then
        rm ${ANC}/concatentated-gained-samples-threshold-no-hits.fa
    fi

    # Re-create with seqs
    for H in ${HOGS}; do
        cat "${HOGSEQS}/${H}.fa" >>${ANC}/concatentated-gained-samples-threshold-no-hits.fa
    done

    clustalo \
        -i ${ANC}/concatentated-gained-samples-threshold-no-hits.fa \
        --percent-id \
        --distmat-out=${OUT}/gained-samples-threshold-no-hits-pid.txt \
        --out ${OUT}/gained-samples-threshold-no-hits-pid.clw \
        --full \
        --verbose \
        --force
done
