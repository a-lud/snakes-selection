#!/usr/bin/env bash

# CDHITDIR="/home/a1645424/al/00_phd/genome-data/transcriptomes/cdhit"
TRANSDECODEROUT="/home/a1645424/al/00_phd/genome-data/transcriptomes/transdecoder"
COMPLETEORF="/home/a1645424/al/00_phd/genome-data/transcriptomes/completeORF"

mkdir -p ${COMPLETEORF}

SAMPLES=$(find ${TRANSDECODEROUT} -type f -name '*.transdecoder.cds' )

for fasta in ${SAMPLES}; do
    BN=$(basename "${fasta%%.*}")

    if [[ ! -f ${COMPLETEORF}/${BN}.completeORF.cds ]]; then
        # Complete CDS
        grep "type:complete" ${fasta} | sed 's/^>//g' > ${TRANSDECODEROUT}/${BN}/complete.txt

        # Subset sequences
        seqtk subseq ${fasta} ${TRANSDECODEROUT}/${BN}/complete.txt > ${COMPLETEORF}/${BN}.completeORF.cds
        seqtk subseq ${fasta%.*}.pep ${TRANSDECODEROUT}/${BN}/complete.txt > ${COMPLETEORF}/${BN}.completeORF.pep
    fi

done