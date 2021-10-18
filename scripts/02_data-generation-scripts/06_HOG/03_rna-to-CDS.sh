#!/usr/bin/env bash
#SBATCH --job-name=validation-rna
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -a 1-16
#SBATCH --ntasks-per-core=1
#SBATCH --time=24:00:00
#SBATCH --mem=16GB
#SBATCH -o /home/a1645424/hpcfs/orthologue-analysis/oma-16-samples/slurm/%x_%a_%A_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alastair.ludington@adelaide.edu.au

# --------------------------------------------------------------------------- #
# Directories
REF='/home/a1645424/hpcfs/orthologue-analysis/oma-16-samples/data/cds'
READS='/home/a1645424/hpcfs/genome-data/raw/rna'
OUT="/home/a1645424/hpcfs/orthologue-analysis/oma-16-samples/rna-alignment"

mkdir -p ${OUT}

# --------------------------------------------------------------------------- #
# Pipeline
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

conda activate BIOTOOLS_env

# Get reference file
REFSEQ=$(find ${REF} -type f -name '*.cds' | tr '\n' ' ' | cut -d ' ' -f ${SLURM_ARRAY_TASK_ID})
REFBN=$(basename ${REFSEQ} .cds)

# Read files for each reference
case ${REFBN} in
aipysurusLaevis_ACAGTG_ACTTGA_CGATGT)
    R1=${READS}/${REFBN}_R1.fastq.gz
    R2=${READS}/${REFBN}_R2.fastq.gz
    ;;
aipysurusLaevis_KLS0468)
    R1=${READS}/${REFBN}_R1.fastq.gz
    R2=${READS}/${REFBN}_R2.fastq.gz
    ;;
aipysurusLaevis_L1_eye)
    R1=${READS}/${REFBN}_R1.fastq.gz
    R2=${READS}/${REFBN}_R2.fastq.gz
    ;;
aipysurusMosaicus)
    R1=${READS}/${REFBN}_R1.fastq.gz
    R2=${READS}/${REFBN}_R2.fastq.gz
    ;;
Atenuis_CAGATC_TTAGGC_GCCAAT)
    R1=${READS}/${REFBN}_R1.fastq.gz
    R2=${READS}/${REFBN}_R2.fastq.gz
    ;;
hydrelapsDarwiniensis)
    R1=${READS}/${REFBN}_R1.fastq.gz
    R2=${READS}/${REFBN}_R2.fastq.gz
    ;;
hydrophisKingii)
    R1=${READS}/${REFBN}_R1.fastq.gz
    R2=${READS}/${REFBN}_R2.fastq.gz
    ;;
hydrophisMajor)
    R1=${READS}/${REFBN}_R1.fastq.gz
    R2=${READS}/${REFBN}_R2.fastq.gz
    ;;
hydrophisMajor_KLS0460)
    R1=${READS}/${REFBN}_R1.fastq.gz
    R2=${READS}/${REFBN}_R2.fastq.gz
    ;;
hydrophisPeronii)
    R1=${READS}/${REFBN}_R1.fastq.gz
    R2=${READS}/${REFBN}_R2.fastq.gz
    ;;
aipysurusLaevis)
    cat ${READS}/aipysurusLaevis_ACAGTG_ACTTGA_CGATGT_R1.fastq.gz ${READS}/aipysurusLaevis_KLS0468_R1.fastq.gz ${READS}/aipysurusLaevis_L1_R1.fastq.gz >${READS}/reads_R1.fastq.gz
    cat ${READS}/aipysurusLaevis_ACAGTG_ACTTGA_CGATGT_R2.fastq.gz ${READS}/aipysurusLaevis_KLS0468_R2.fastq.gz ${READS}/aipysurusLaevis_L1_R2.fastq.gz >${READS}/reads_R2.fastq.gz
    R1=${READS}/reads_R1.fastq.gz
    R2=${READS}/reads_R2.fastq.gz
    ;;
*)
    printf "Unused file: %s\n" ${REFBN}
    exit 0
    ;;
esac

printf "REF: %s\nBASENAME: %s\nR1: %s\nR2: %s\n\n" ${REFSEQ} ${REFBN} ${R1} ${R2}

# Create BWA index
if [[ ! -f ${REFSEQ}.0123 ]]; then
    printf "Creating index: %s\n\n" ${REFBN}
    bwa-mem2 index ${REFSEQ}
fi

# Align data
bwa-mem2 mem \
    -t ${SLURM_CPUS_PER_TASK} \
    -M \
    ${REFSEQ} \
    ${R1} \
    ${R2} |
    samtools sort -T ${OUT}/${REFBN} |
    samtools view -F 4 -q 20 -o ${OUT}/${REFBN}.bam

conda deactivate
