#!/usr/bin/env bash
#SBATCH --job-name=validate-HOGS
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -a 1-5
#SBATCH --ntasks-per-core=1
#SBATCH --time=72:00:00
#SBATCH --mem=8GB
#SBATCH -o /home/a1645424/hpcfs/orthologue-analysis/oma-16-samples/slurm/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

# --------------------------------------------------------------------------- #
# Script:
#
# BLAST the HOGs belonging to different sample-groups to their respective
# counter samples. A gained hog in a group should not have homology to any
# sequences in the BLAST-DB. If it does, this indicates that the HOG has not
# been gained since the ancestral-grouping-node (as predicted by PyHam).
#
## --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Directories
PROJ="/home/a1645424/hpcfs/orthologue-analysis/oma-16-samples"
REFDIR="${PROJ}/data/protein"

# --------------------------------------------------------------------------- #
# Array - dataset to work on
DIR=$(ls -d ${PROJ}/hog-gain-loss/* | tr '\n' ' ' | cut -d ' ' -f ${SLURM_ARRAY_TASK_ID})

OUT_GAIN=${DIR}/blast-gained-samples
mkdir -p ${OUT_GAIN}

OUT_GAIN_THRESHOLD=${DIR}/blast-gained-samples-threshold
mkdir -p ${OUT_GAIN_THRESHOLD}

OUT_LOST=${DIR}/blast-lost-samples
mkdir -p ${OUT_LOST}

# --------------------------------------------------------------------------- #
# Pipeline
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate BLAST_env

# References to check gain/lost HOGs to - should be no presence in either set
# (respectively)
REF_GAINS=$(cat ${DIR}/refs-validate-gains.txt)
REF_LOST=$(cat ${DIR}/refs-validate-lost.txt)

# Gained sequences
for REF in ${REF_GAINS}; do
    find "${DIR}/samples/gained_fasta" -type f -name '*.fa' | parallel \
        -j ${SLURM_CPUS_PER_TASK} \
        --joblog "${OUT_GAIN}/parallel-gain.log" \
        "blastp -query {} -db ${REFDIR}/${REF} -out ${OUT_GAIN}/{/.}-${REF}.outfmt6 -evalue 1e-10 -outfmt \"6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident qcovs qcovhsp mismatch gaps\" -num_threads 1"
done

# Gained threshold
for REF in ${REF_GAINS}; do
    find "${DIR}/samples/gained_threshold_fasta" -type f -name '*.fa' | parallel \
        -j ${SLURM_CPUS_PER_TASK} \
        --joblog "${OUT_GAIN_THRESHOLD}/parallel-gain-threshold.log" \
        "blastp -query {} -db ${REFDIR}/${REF} -out ${OUT_GAIN_THRESHOLD}/{/.}-${REF}.outfmt6 -evalue 1e-10 -outfmt \"6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident qcovs qcovhsp mismatch gaps\" -num_threads 1"
done

# Lost sequences
for REF in ${REF_LOST}; do
    find "${DIR}/samples/lost_fasta" -type f -name '*.fa' | parallel \
        -j ${SLURM_CPUS_PER_TASK} \
        --joblog "${OUT_LOST}/parallel-lost.log" \
        "blastp -query {} -db ${REFDIR}/${REF} -out ${OUT_LOST}/{/.}-${REF}.outfmt6 -evalue 1e-10 -outfmt \"6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident qcovs qcovhsp mismatch gaps\" -num_threads 1"
done

conda deactivate
