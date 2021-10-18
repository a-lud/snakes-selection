#!/usr/bin/env bash
#SBATCH --job-name=mitchochondria-check
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -a 1-5
#SBATCH --ntasks-per-core=1
#SBATCH --time=06:00:00
#SBATCH --mem=4GB
#SBATCH -o /home/a1645424/hpcfs/orthologue-analysis/oma-16-samples/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

# --------------------------------------------------------------------------- #
# Conda parameters
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate BIOTOOLS_env

# --------------------------------------------------------------------------- #
# Directories
DIR='/home/a1645424/hpcfs/orthologue-analysis/oma-16-samples'
REF='/home/a1645424/hpcfs/database/mitochondrion.1-2021-05-20'

# --------------------------------------------------------------------------- #
# Get array + remaining directories
ANC=$(find ${DIR}/hog-gain-loss -type d -name 'ancestral-node-*' | tr '\n' ' ' | cut -d ' ' -f ${SLURM_ARRAY_TASK_ID})
HOGDIR="${ANC}/samples/gained_threshold_fasta"
OUT="${ANC}/mitochondria-check"

mkdir -p ${OUT}

# --------------------------------------------------------------------------- #
# Create BLASTDB
if [[ ! -f ${REF}.phr ]]; then
    makeblastdb \
        -dbtype prot \
        -in ${REF}.protein.faa \
        -parse_seqids \
        -out ${REF} \
        -logfile ${REF}.makeblastdb.log
fi

# --------------------------------------------------------------------------- #
# HOGs to search against Hydrophis snakesx
HOGIDS=$(cat ${ANC}/blast-gained-samples-threshold-no-hit.txt | tr '\n' ' ')

for hog in ${HOGIDS}; do
    printf "Query: %s\n" ${hog}

    blastp \
        -query ${HOGDIR}/${hog}.fa \
        -db ${REF} \
        -out ${OUT}/${hog}.outfmt6 \
        -outfmt "6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident qcovs qcovhsp mismatch gaps" \
        -num_threads ${SLURM_CPUS_PER_TASK}

    blastp \
        -query ${HOGDIR}/${hog}.fa \
        -db ${REF} \
        -out ${OUT}/${hog}-blosum45.outfmt6 \
        -outfmt "6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore lengthpident qcovs qcovhsp mismatch gaps" \
        -num_threads ${SLURM_CPUS_PER_TASK} \
        -matrix BLOSUM45
done

conda deactivate
