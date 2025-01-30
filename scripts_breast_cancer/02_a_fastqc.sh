#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=02:10:00
#SBATCH --mem=2g
#SBATCH --cpus-per-task=2
#SBATCH --job-name=fastqc_attempt1
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/fastqc_attempt1%J.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/fastqc_attempt1%J.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL

#Note: SBATCH array above allows you to process 12 samples in parallel, which is efficient.

# define variables
MAINDIR="/data/users/${USER}/rnaseq_course"

WORKDIR="/data/users/${USER}/rnaseq_course/FASTQ_files"

OUTDIR="$MAINDIR/QC_RESULTS"

# Create the directory that we specified above if they don't exist already
mkdir -p $OUTDIR

SAMPLELIST="$WORKDIR/samplelist/sample_list.txt"

SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# Load the FastQC module
module load FastQC/0.11.9-Java-11

# Run FastQC on the reads
fastqc "$READ1" "$READ2" -o "$OUTDIR"

# Log completion
echo "FastQC completed for $SAMPLE. Results saved to $OUTDIR."


