#!/usr/bin/env bash
#SBATCH --job-name=subset-rna-bam
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -a 1-11
#SBATCH --ntasks-per-core=1
#SBATCH --time=2:00:00
#SBATCH --mem=4GB
#SBATCH -o /home/a1645424/hpcfs/orthologue-analysis/oma-16-samples/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

# --------------------------------------------------------------------------- #
# Directories
BAMS="/home/a1645424/hpcfs/orthologue-analysis/oma-16-samples/rna-alignment"
DIR="/home/a1645424/hpcfs/orthologue-analysis/oma-16-samples/hog-gain-loss"

# --------------------------------------------------------------------------- #
# Pipeline
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate BIOTOOLS_env

# Get BAM file
BAM=$(find ${BAMS} -type f -name '*sorted.bam' | tr '\n' ' ' | cut -d ' ' -f ${SLURM_ARRAY_TASK_ID})
BAMBN=$(basename "${BAM%%.*}")
printf "BAM File: %s\n" ${BAMBN}

DIRS=$(find ${DIR} -type d -name 'ancestral-node-*' | tr '\n' ' ')

for ANC in ${DIRS}; do
    HEADER_FILES=$(find "${ANC}/hog-headers-by-snake" -type f -name '*.txt')
    OUT="${ANC}/rna-subset"
    mkdir -p ${OUT}

    for FILE in ${HEADER_FILES}; do
        FILEBN=$(basename ${FILE} .txt)
        if [[ ${FILEBN} == "${BAMBN}" ]]; then
            GENES=$(cat ${FILE} | tr '\n' ' ')
            samtools view -b ${BAM} ${GENES} |
                samtools sort \
                    -l 9 \
                    -o ${OUT}/${BAMBN}.sorted.bam \
                    -T ${OUT}/${BAMBN}
            samtools index ${OUT}/${BAMBN}.sorted.bam
            samtools depth -a -Q 20 -o ${OUT}/${BAMBN}.depth ${OUT}/${BAMBN}.sorted.bam
        fi
    done
done

[main_samview] region "TRINITY_DN12305_c0_g1_i1.p2:1-360" specifies an invalid region or unknown reference. Continue anyway.
[main_samview] region "TRINITY_DN18603_c0_g1_i1.p1:1-253" spec
