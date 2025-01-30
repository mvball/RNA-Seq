#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --mem=144G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=featurecount
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/featurecounts%A.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/featurecounts%A.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-12

# Define directories and files
WORKDIR="/data/users/${USER}/rnaseq_course"
LOGDIR="$WORKDIR/log/featurecounts"
MAPPINGDIR="$WORKDIR/MAPPING_BAMS"
OUTDIR="$WORKDIR/FEATURECOUNTS"
ANNOTATIONFILE="$WORKDIR/reference_genome/Homo_sapiens.GRCh38.113.gtf.gz"
SAMPLELIST="$WORKDIR/FASTQ_files/samplelist/sample_list.txt"

# Create output and log directories if not present
mkdir -p "$LOGDIR"
mkdir -p "$OUTDIR"

# Unzip the annotation file (if not already unzipped)
if [[ ! -f "$ANNOTATIONFILE" ]]; then
     gunzip -c "$WORKDIR/reference_genome/Homo_sapiens.GRCh38.113.gtf.gz" > "$ANNOTATIONFILE"
fi

# Print the SLURM_ARRAY_TASK_ID to the log for debugging
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID" >> "$LOGDIR/featurecounts_${SLURM_ARRAY_TASK_ID}.log"

# Get sample name from samplelist.txt
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)

# Log the extracted sample name
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID, SAMPLE: $SAMPLE" >> "$LOGDIR/featurecounts_${SLURM_ARRAY_TASK_ID}.log"

echo "Using file: $MAPPINGDIR/${SAMPLE}_sorted.bam"

# Check if the BAM file exists
if [[ ! -f "$MAPPINGDIR/${SAMPLE}_sorted.bam" ]]; then
    echo "ERROR: The file $MAPPINGDIR/${SAMPLE}_sorted.bam does not exist!" >> "$LOGDIR/featurecounts_${SLURM_ARRAY_TASK_ID}.log"
    exit 1
fi

# Log when starting
echo "Starting featureCounts for sample: $SAMPLE at $(date)" >> "$LOGDIR/featurecounts_${SLURM_ARRAY_TASK_ID}.log"

# Run featureCounts for the specific BAM file
apptainer exec --bind "$WORKDIR:$WORKDIR" \
    /containers/apptainer/subread_2.0.1--hed695b0_0.sif featureCounts \
    -T 4 \
    -p \
    -s 0 \
    -Q 10 \
    -t exon \
    -g gene_id \
    -a "$ANNOTATIONFILE" \
    -o "$OUTDIR/gene_count_${SAMPLE}.txt" \
    "$MAPPINGDIR/${SAMPLE}_sorted.bam"

# Logging when finished
echo "Finished featureCounts for sample: $SAMPLE at $(date)" >> "$LOGDIR/featurecounts_${SLURM_ARRAY_TASK_ID}.log"


#Key featureCounts Parameters:
# -T 4: Specifies the use of 4 threads (CPU cores) for processing.
# -p: Indicates paired-end sequencing, meaning reads come in pairs (e.g., R1 and R2).
# -s 0: Unstranded.
# -Q 10: Filters reads with a mapping quality below 10, ensuring only confidently mapped reads are counted.
# -t exon: Specifies the feature type to count as exons (defined in the GTF file).
# -g gene_id : Groups counts by the gene_id attribute in the GTF file.
# -a $ANNOTATIONFILE : Specifies the annotation file (GTF) used to map reads to features.
# -o $OUTDIR/gene_count.txt : Writes the output (gene count matrix) to the specified file.