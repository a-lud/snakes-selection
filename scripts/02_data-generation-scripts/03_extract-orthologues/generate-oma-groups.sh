#!/usr/bin/env bash

BASE='/Users/alastair.ludington/Documents/phd/00_papers/snakes-selection'
OMA="${BASE}/data/oma-16-samples/Output"
NUC="${BASE}/data/oma-cds"
SCRIPT="${BASE}/scripts_paper/util"
OUT="${OMA}/../complete-oma-groups"

# -------------------------------------------------------- #
# Generate orthologue combinations
# -------------------------------------------------------- #

# Four samples - genomes
${SCRIPT}/getCompleteOrthologues.py \
    -oma ${OMA} \
    -s 'aipysurusLaevis notechisScutatus pseudonajaTextilis najaNaja' \
    -n ${NUC} \
    -e '.cds' \
    -o "${OUT}/four-sample"

# Nine samples - genomes + transcriptomes (original)
${SCRIPT}/getCompleteOrthologues.py \
    -oma ${OMA} \
    -s 'Nscutatus_L1_eye Ptextilis_L1_eye aipysurusLaevis aipysurusLaevis_ACAGTG_ACTTGA_CGATGT 
    aipysurusLaevis_KLS0468 aipysurusLaevis_L1_eye najaNaja notechisScutatus pseudonajaTextilis' \
    -n ${NUC} \
    -e '.cds' \
    -o "${OUT}/nine-sample"

# Sixteen samples - all data
${SCRIPT}/getCompleteOrthologues.py \
    -oma ${OMA} \
    -n ${NUC} \
    -e '.cds' \
    -o "${OUT}/sixteen-sample"
