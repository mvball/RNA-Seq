#!/bin/bash
#SBATCH --job-name=featureCounts_multiQC
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/featureCounts_multiQC%A.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/featureCounts_multiQC%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL

# Define directories
MAINDIR="/data/users/${USER}/rnaseq_course"
WORKDIR="$MAINDIR/FEATURECOUNTS"  # Directory containing .txt mapping results
OUTDIR="$MAINDIR/FEATURECOUNTS"  # Directory for MultiQC results

# Ensure output directory exists
mkdir -p "$OUTDIR"

# Check if input directory exists
if [ ! -d "$WORKDIR" ]; then
    echo "Error: $WORKDIR does not exist." >&2
    exit 1
fi

# Run MultiQC
apptainer exec --bind "$MAINDIR":"$MAINDIR" /containers/apptainer/multiqc-1.19.sif \
    multiqc -o "$OUTDIR" "$WORKDIR"