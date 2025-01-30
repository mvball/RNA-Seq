#!/bin/bash
#SBATCH --time=14:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=hisat_mapping
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/hisat_mapping%Z.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/hisat_mapping%Z.err
#SBATCH --partition=pibu_el8
#SBATCH --array=1-12
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL

# Print task ID for debugging
echo "Task ID: $SLURM_ARRAY_TASK_ID"

# Define paths
WORKDIR="/data/users/${USER}/rnaseq_course"
INDEXDIR="$WORKDIR/index_hisat"
OUTDIR="$WORKDIR/MAPPING"
BASE_FASTQ_DIR="$WORKDIR/FASTP_RESULTS"
SAMPLELIST="$WORKDIR/FASTQ_files/samplelist/sample_list.txt"
CONTAINER_IMAGE="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# Create output directory if it doesn't exist
mkdir -p $OUTDIR

# Get the sample name and file paths
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)
R1_FILE="$BASE_FASTQ_DIR/${SAMPLE}_clean_R1.fastq.gz"
R2_FILE="$BASE_FASTQ_DIR/${SAMPLE}_clean_R2.fastq.gz"

# Print sample details for debugging
echo "Sample name: $SAMPLE"
echo "R1 file: $R1_FILE"
echo "R2 file: $R2_FILE"

# Check if the FastQ files exist
if [[ ! -f "$R1_FILE" ]]; then
    echo "Error: $R1_FILE does not exist!" >&2
    exit 1
fi
if [[ ! -f "$R2_FILE" ]]; then
    echo "Error: $R2_FILE does not exist!" >&2
    exit 1
fi

# Map the reads
apptainer exec --bind /data:/data $CONTAINER_IMAGE \
    hisat2 -x $INDEXDIR/genome_index \
           -1 $R1_FILE \
           -2 $R2_FILE \
           -S $OUTDIR/${SAMPLE}.sam \
           --threads 4 \
           --summary-file $OUTDIR/${SAMPLE}_mapping_summary.txt




#apptainer exec: Runs a command within a container image. The container image here is /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif.
# --bind /data: Maps the /data directory on your host system into the container's filesystem. Ensures that the container can access /data while running.
# hisat2: The primary tool being executed. It performs the alignment of RNA-Seq reads to a reference genome.

# -x $INDEXDIR/genome_index: Specifies the path to the HISAT2 index built from the reference genome. This is the reference the reads will be aligned to.
# -1 $READ1: Specifies the file containing the first reads (R1) from paired-end sequencing.
# -2 $READ2: Specifies the file containing the second reads (R2) from paired-end sequencing.
# -S $OUTDIR/$SAMPLE.sam: Specifies the output SAM file where the aligned reads will be written.
# --threads 4: Allocates 4 threads for parallel processing, speeding up the alignment process.
# --rna-strandness : Indicates the strand-specificity of the RNA-Seq reads: We have unstranded data so this is omitted, but can be RF or RF
# --summary-file $OUTDIR/${SAMPLE}mapping_summary.txt: Writes a summary of the mapping results (e.g., percentage of reads mapped) to a specified file.

#For our data, The library preparation protocol did not preserve information on the transcribed strand (i.e. unstranded). The libraries were sequenced on an Illumina HiSeq 2000 in paired-end mode.
