#!/usr/bin/env bash

## Required Phoenix modules for the pipeline to run
module load arch/haswell
module load Java/11.0.6

## Location of where the pipeline is installed
PIPE="/home/a1645424/hpcfs/software/nf-pipelines"
DIR="/home/a1645424/hpcfs/orthologue-analysis/oma-16-samples"
MSA="${DIR}/msa/nine-sample/nucleotide_msa"

nextflow run ${PIPE}/main.nf \
    -profile conda,slurm \
    --pipeline hyphy \
    -work-dir ${DIR}/work-nine-sample-hyphy \
    --files_dir ${MSA} \
    --files_ext '*.fasta' \
    --files_type fasta \
    --outdir ${DIR}/selection/nine-sample \
    --email alastair.ludington@adelaide.edu.au \
    --partition skylake \
    --method busted,absrel,meme \
    --tree ${DIR}/trees/nine-sample-fg.nwk \
    --absrel_optional '--branches Foreground' \
    --meme_optional '--branches Foreground' \
    -resume
