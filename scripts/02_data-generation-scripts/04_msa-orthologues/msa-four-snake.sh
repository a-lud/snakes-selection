#!/usr/bin/env bash

## Required Phoenix modules for the pipeline to run
module load arch/haswell
module load Java/11.0.6

## Location of where the pipeline is installed
PIPE="/home/a1645424/hpcfs/software/nf-pipelines"
DIR="/home/a1645424/hpcfs/orthologue-analysis/oma-16-samples"

## MSA pipeline
nextflow run ${PIPE}/main.nf \
    -profile slurm,conda \
    --pipeline msa \
    -work-dir ${DIR}/work-four-sample \
    --outdir ${DIR}/msa/four-sample \
    --email a1645424@adelaide.edu.au \
    --files_dir ${DIR}/complete-oma-groups/four-sample/protein \
    --files_ext '*.fa' \
    --files_type fasta \
    --aligner mafft \
    --pep2nuc \
    --nucleotide_dir ${DIR}/complete-oma-groups/four-sample/nucleotide \
    --nucleotide_ext '*.fa' \
    --partition test \
    -resume
