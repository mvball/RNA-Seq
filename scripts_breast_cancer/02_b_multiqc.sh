#!/bin/bash
#SBATCH --job-name=multiQC1
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/multiQC1%A.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/multiQC1%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL

# Define directories
MAINDIR="/data/users/${USER}/rnaseq_course"
WORKDIR="$MAINDIR/QC_RESULTS"  # Directory containing FastQC results
OUTDIR="$MAINDIR/MULTIQC_RESULTS"  # Directory for MultiQC results

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Run MultiQC, use bind so that the host system can be accessible inside the container
apptainer exec --bind "$MAINDIR":"$MAINDIR" /containers/apptainer/multiqc-1.19.sif \
    multiqc -o "$OUTDIR" "$WORKDIR"

