#!/usr/bin/env bash
#SBATCH --job-name=blast-nr
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH -a 1-5
#SBATCH --ntasks-per-core=1
#SBATCH --time=72:00:00
#SBATCH --mem=40GB
#SBATCH -o /home/a1645424/hpcfs/orthologue-analysis/oma-16-samples/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

# --------------------------------------------------------------------------- #
# Conda parameters
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate BLAST_env

# --------------------------------------------------------------------------- #
# Directories
DIR='/home/a1645424/hpcfs/orthologue-analysis/oma-16-samples'
REF='/hpcfs/archive/biorefs/blast_db/nr'

# --------------------------------------------------------------------------- #
# Get array + remaining directories
ANC=$(find ${DIR}/hog-gain-loss -type d -name 'ancestral-node-*' | tr '\n' ' ' | cut -d ' ' -f ${SLURM_ARRAY_TASK_ID})
HOGDIR="${ANC}/samples/gained_threshold_fasta"
OUT="${ANC}/hog-homology/blast-nr"

mkdir -p ${OUT}

# --------------------------------------------------------------------------- #
# HOGs to search against NR database
HOGIDS=$(cat ${ANC}/blast-gained-samples-threshold-no-hit.txt | tr '\n' ' ')

for hog in ${HOGIDS}; do
    printf "Query: %s\n" ${hog}

    blastp \
        -query ${HOGDIR}/${hog}.fa \
        -db ${REF} \
        -out ${OUT}/${hog}.outfmt6 \
        -outfmt "6 qseqid sseqid qlen slen qstart qend sstart send qseq sseq evalue bitscore length pident qcovs qcovhsp mismatch gaps sacc sgi stitle" \
        -num_threads ${SLURM_CPUS_PER_TASK}
done

conda deactivate
