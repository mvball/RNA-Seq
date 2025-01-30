#!/bin/bash
#SBATCH --job-name=multiQC_corrected
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/multiQC_corrected%D.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/multiQC_corrected%D.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --time=10:00:00
#SBATCH --mem=24G
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL

# Define directories
MAINDIR="/data/users/${USER}/rnaseq_course"
WORKDIR="$MAINDIR/FASTP_RESULTS"  # Directory containing cleaned FastQ files
OUTDIRFASTQ="$MAINDIR/corrected_FASTQC_RESULTS"  # Directory for FastQC results
OUTDIRMULTI="$MAINDIR/corrected_MULTIQC_RESULTS"  # Directory for MultiQC results
LOGDIR="$WORKDIR/log"  # Log directory for FastQC and MultiQC

# Create directories if they do not exist
mkdir -p "$OUTDIRFASTQ" "$OUTDIRMULTI" "$LOGDIR"

# Check the contents of the cleaned FastQC root directory
echo "Checking cleaned FastQ root directory: $WORKDIR"
ls -lh "$WORKDIR"

# Load FastQC module
module load FastQC/0.11.9-Java-11

# Run FastQC on cleaned FastQ files
echo "Running FastQC on cleaned data..."
fastqc -o "$OUTDIRFASTQ" "$WORKDIR"/*_clean_R1.fastq.gz "$WORKDIR"/*_clean_R2.fastq.gz || { echo "FastQC failed. Exiting."; exit 1; }

# Run MultiQC on all FastQC results
echo "Running MultiQC on FastQC results..."
apptainer exec --bind "$OUTDIRFASTQ":"$OUTDIRFASTQ" \
               --bind "$OUTDIRMULTI":"$OUTDIRMULTI" \
               /containers/apptainer/multiqc-1.19.sif \
               multiqc -o "$OUTDIRMULTI" "$OUTDIRFASTQ" || { echo "MultiQC failed. Exiting."; exit 1; }

# Check if the MultiQC run was successful
echo "MultiQC run complete. Checking output directory for results..."
ls -lh "$OUTDIRMULTI"
