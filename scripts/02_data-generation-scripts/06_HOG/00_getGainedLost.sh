#!/usr/bin/env bash

# --------------------------------------------------------------------------- #
# Directories
UTIL='/Users/alastair.ludington/Documents/phd/00_papers/snakes-selection-hogs/scripts_paper/util'
OMA='/Users/alastair.ludington/Documents/phd/00_papers/snakes-selection-hogs/data/oma-16-samples/Output'
TREE='/Users/alastair.ludington/Documents/phd/00_papers/snakes-selection-hogs/data/trees'
OUT='/Users/alastair.ludington/Documents/phd/00_papers/snakes-selection-hogs/data/HOG-data'

# --------------------------------------------------------------------------- #
# Get HOGs at different levels

if [[ ! -f "${OUT}/hogs-gain-loss/nodes.txt" ]]; then
    ${UTIL}/pyham-getHogs.py \
        -pn \
        -tr "${TREE}/sixteen-sample.nwk" \
        "${OMA}" \
        "${OUT}/hogs-gain-loss"
fi

${UTIL}/pyham-getHogs.py \
    -n "${OUT}/hogs-gain-loss/nodes.txt" \
    -tr "${TREE}/sixteen-sample.nwk" \
    "${OMA}" \
    "${OUT}/hogs-gain-loss"
