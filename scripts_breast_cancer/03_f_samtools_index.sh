#!/bin/bash
#SBATCH --time=02:30:00
#SBATCH --mem=8G        
#SBATCH --cpus-per-task=1
#SBATCH --job-name=index_samtools
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/index_samtools%D.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/index_samtools%D.err
#SBATCH --partition=pibu_el8
#SBATCH --array=1-12     
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL

# Setting the constant for the directories and required files
WORKDIR="/data/users/${USER}/rnaseq_course"
LOGDIR="$WORKDIR/log"
SAMPLELIST="$WORKDIR/FASTQ_files/samplelist/sample_list.txt"
OUTDIR="$WORKDIR/MAPPING_BAMS"

# Create the directory for the error and output file if not present
mkdir -p $LOGDIR

# Take the sample name, path to read1 and read2 line by line 
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)
READ1=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST)
READ2=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST)

# Debug and give an update
echo "Starting task for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "Processing sample: $SAMPLE for indexing" >> $LOGDIR/${SLURM_ARRAY_TASK_ID}.log

# Check if the index file already exists
if [ ! -f "$OUTDIR/${SAMPLE}_sorted.bam.bai" ]; then
  # Indexing the bam file
  echo "Indexing BAM file for sample: $SAMPLE" >> $LOGDIR/${SLURM_ARRAY_TASK_ID}.log
  
  apptainer exec --bind $WORKDIR /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
    samtools index $OUTDIR/${SAMPLE}_sorted.bam
  
  echo "Finished indexing sample: $SAMPLE" >> $LOGDIR/${SLURM_ARRAY_TASK_ID}.log
else
  echo "Index file for $SAMPLE already exists. Skipping indexing." >> $LOGDIR/${SLURM_ARRAY_TASK_ID}.log
fi
