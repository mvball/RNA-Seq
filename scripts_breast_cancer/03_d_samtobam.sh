#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=samtools_conversion_bam
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/samtools_conversion_bam%E.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/samtools_conversion_bam%E.err
#SBATCH --partition=pibu_el8
#SBATCH --array=1-12
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL

# Setting the constants for directories and required files
WORKDIR="/data/users/${USER}/rnaseq_course"
LOGDIR="$WORKDIR/log"
SAMPLELIST="$WORKDIR/FASTQ_files/samplelist/sample_list.txt"
MAPDIR=$WORKDIR/MAPPING
OUTDIR=$WORKDIR/MAPPING_BAMS

# Create the directory for error and output files if not present
mkdir -p $LOGDIR
mkdir -p $OUTDIR

# Extract sample name and file paths for this task
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)

# Log the current task
echo "Starting task for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "Processing sample: $SAMPLE"

# Check if the SAM file exists before proceeding
SAM_FILE="$MAPDIR/$SAMPLE.sam"
if [[ ! -f "$SAM_FILE" ]]; then
    echo "Error: SAM file $SAM_FILE does not exist!" >&2
    exit 1
fi

# Log the SAM to BAM conversion process
echo "Converting $SAM_FILE to BAM format..."
apptainer exec --bind $MAPDIR /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
    samtools view -hbS $SAM_FILE > $OUTDIR/$SAMPLE.bam

# Confirm successful completion
if [[ $? -eq 0 ]]; then
    echo "Successfully converted $SAM_FILE to $OUTDIR/$SAMPLE.bam"
else
    echo "Error occurred during conversion of $SAM_FILE to BAM" >&2
    exit 1
fi

echo "Task for $SAMPLE completed."
