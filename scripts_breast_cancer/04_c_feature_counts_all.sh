#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --mem=144G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=featurecount
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/featurecounts%B.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/featurecounts%B.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-12

#Setting the constant for the directory and the required files
WORKDIR="/data/users/${USER}/rnaseq_course"
LOGDIR="$WORKDIR/log"
SAMPLELIST="$WORKDIR/FASTQ_files/samplelist/sample_list.txt"
MAPPINGDIR="$WORKDIR/MAPPING_BAMS"
OUTDIR="$WORKDIR/FEATURECOUNTS"
ANNOTATIONFILE="$WORKDIR/reference_genome/Homo_sapiens.GRCh38.113.gtf.gz"


#Create the directory for the error and output file if not present
mkdir -p $LOGDIR

mkdir -p $OUTDIR

#take the sample name, path to the read1 and read2 line by line 
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# Log when starting
echo "Starting featureCounts for sample: $SAMPLE at $(date)" >> "$LOGDIR/featurecounts_${SLURM_ARRAY_TASK_ID}.log"

apptainer exec --bind $WORKDIR /containers/apptainer/subread_2.0.1--hed695b0_0.sif featureCounts -T4 -p -s0 -Q10 -t exon -g gene_id -a $ANNOTATIONFILE -o "$OUTDIR/gene_count.txt" $MAPPINGDIR/*sorted.bam

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