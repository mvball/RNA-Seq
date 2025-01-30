#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --mem=512G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=samtools_sort
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/samtools_sort%A_%a.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/samtools_sort%A_%a.err
#SBATCH --partition=pibu_el8
#SBATCH --array=1-12
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL

# Setting constants for directories and required files
WORKDIR="/data/users/${USER}/rnaseq_course"
LOGDIR="$WORKDIR/log"
SAMPLELIST="$WORKDIR/FASTQ_files/samplelist/sample_list.txt"
OUTDIR=$WORKDIR/MAPPING_BAMS
TEMP_DIR="/scratch/${USER}/temp"  # Using scratch directory for temporary files

# Create necessary directories
mkdir -p $LOGDIR
mkdir -p $TEMP_DIR

# Take sample name, path to read1, and read2 line by line
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)
READ1=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST)
READ2=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST)

#Log when starting
echo "Starting samtools sort for sample: $SAMPLE" >> $LOGDIR/$SLURM_ARRAY_TASK_ID.out
  
# Sorting the BAM file with temporary files stored in scratch directory
apptainer exec --bind $WORKDIR /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
    samtools sort -m 32G -@ 8 -o $OUTDIR/${SAMPLE}_sorted.bam -T $TEMP_DIR $OUTDIR/$SAMPLE.bam
  
# Logging when finished
echo "Finished sorting sample: $SAMPLE" >> $LOGDIR/$SLURM_ARRAY_TASK_ID.out


#apptainer exec: This command runs a program inside a container. exec stands for "execute" and is the primary command to launch processes inside containers.
#--bind $WORKDIR: This binds (mounts) a local directory ($WORKDIR) to the container's filesystem. This allows the container to access files and directories on the host system.
#samtools is a widely used bioinformatics tool for manipulating SAM/BAM/CRAM files. The sort command is used to sort BAM files, which is an essential step in many sequencing workflows (like alignment or variant calling).

