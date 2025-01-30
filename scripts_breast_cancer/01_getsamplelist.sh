#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500M
#SBATCH --time=00:05:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=sample_list
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/samplelist_%J.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/samplelist_%J.err

# Set constants for directories
WORKDIR="/data/users/mschlaepfer/rnaseq_course"
READSDIR="/data/courses/rnaseq_course/breastcancer_de/reads"
OUTDIR="$WORKDIR/FASTQ_files"
LOGDIR="$WORKDIR/log"
SAMPLE_LIST="$OUTDIR/sample_list.txt"  # Changed to a file path

# Create necessary directories if they don't exist
mkdir -p "$LOGDIR"
mkdir -p "$OUTDIR"

# Create symbolic links for each FASTQ file
for file in "$READSDIR"/*.fastq.gz; do
    FILENAME=$(basename "$file")  # Fixed variable name
    ln -s "$file" "$OUTDIR/$FILENAME"
done

# Generate the sample list
for file in "$OUTDIR"/*_*1.fastq.gz; do
    PREFIX="${file%_*.fastq.gz}"  # Get the prefix of the file name
    SAMPLE=$(basename "$PREFIX") # Extract the sample name
    echo -e "${SAMPLE}\t${file}\t${file%1.fastq.gz}2.fastq.gz"
done > "$SAMPLE_LIST"

echo "Sample list written to: $SAMPLE_LIST"



#Example:
#If the input directory contains the files:
    #/path/to/FASTQ_FOLDER/sample1_1.fastq.gz
    #/path/to/FASTQ_FOLDER/sample1_2.fastq.gz
    #/path/to/FASTQ_FOLDER/sample2_1.fastq.gz
    #/path/to/FASTQ_FOLDER/sample2_2.fastq.gz

#The output will be:
    #sample1    /path/to/FASTQ_FOLDER/sample1_1.fastq.gz    /path/to/FASTQ_FOLDER/sample1_2.fastq.gz
    #sample2    /path/to/FASTQ_FOLDER/sample2_1.fastq.gz    /path/to/FASTQ_FOLDER/sample2_2.fastq.gz