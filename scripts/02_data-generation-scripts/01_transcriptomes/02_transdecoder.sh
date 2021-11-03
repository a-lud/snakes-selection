#!/usr/bin/env bash

CDHITDATA="/home/a1645424/al/00_phd/genome-data/transcriptomes/cdhit"
TRANSDECODEROUT="/home/a1645424/al/00_phd/genome-data/transcriptomes/transdecoder"
DATABASE="/home/a1645424/al/databases"

mkdir -p "${TRANSDECODEROUT}"

source activate TRINITY_NEW

for transcripts in "${CDHITDATA}"/*.fasta; do

    BN=$(basename "${transcripts%%.*}")

    if [[ ! -d ${TRANSDECODEROUT}/${BN} ]]; then

        if [[ ! -f "${TRANSDECODEROUT}/${BN}/longest_orfs.cds" ]]; then
            # Create gene-trans-map
            get_Trinity_gene_to_trans_map.pl \
                "${transcripts} >${CDHITDATA}/${BN}.gene_trans_map"

            # Extract longest open reading frames
            TransDecoder.LongOrfs \
                -t ${transcripts} \
                --gene_trans_map "${CDHITDATA}/${BN}.gene_trans_map" \
                --output_dir ${TRANSDECODEROUT}/${BN}
        fi

        # BLAST
        if [[ ! -f ${TRANSDECODEROUT}/${BN}/${BN}.outfmt6 ]]; then
            echo "BLAST: ${BN}\n"
            blastp \
                -query "${TRANSDECODEROUT}/${BN}/longest_orfs.pep" \
                -db "${DATABASE}/uniprot_sprot.fasta" \
                -max_target_seqs 1 \
                -outfmt 6 \
                -evalue 1e-5 \
                -num_threads 40 >${TRANSDECODEROUT}/${BN}/${BN}.outfmt6
        fi

        # PFAM
        if [[ ! -f ${TRANSDECODEROUT}/${BN}/${BN}.domtblout ]]; then
            echo "HMMER: ${BN}\n"
            hmmscan \
                --cpu 40 \
                --domtblout "${TRANSDECODEROUT}/${BN}/${BN}.domtblout" \
                "${DATABASE}/Pfam-A.hmm" \
                "${TRANSDECODEROUT}/${BN}/longest_orfs.pep"
        fi

        # Predict coding regions
        if [[ ! -f "${TRANSDECODEROUT}/${BN}/${BN}.cdhit.fasta.transdecoder.cds" ]]; then
            TransDecoder.Predict \
                -t ${transcripts} \
                --retain_pfam_hits "${TRANSDECODEROUT}/${BN}/${BN}.domtblout" \
                --retain_blastp_hits "${TRANSDECODEROUT}/${BN}/${BN}.outfmt6" \
                --single_best_only \
                --output_dir ${TRANSDECODEROUT}/${BN}

            cp "${BN}.cdhit.fasta.transdecoder*" ${TRANSDECODEROUT}/${BN}
        fi

    fi
done

conda deactivate
