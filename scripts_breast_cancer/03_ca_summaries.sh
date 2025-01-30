#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem=100M
#SBATCH --cpus-per-task=1
#SBATCH --job-name=merging_summaries
#SBATCH --job-name=hisat_mapping
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/merging_summaries%C.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/merging_summaries%C.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL

#Setting the constant for the directories and required files
WORKDIR="/data/users/${USER}/rnaseq_course"
MAPPINGDIR="$WORKDIR/MAPPING"
LOGDIR="$WORKDIR/log"
OUTFILE="$MAPPINGDIR/all_summary_mapping.txt"

#Create the directory for the error and output file if not present
mkdir -p $LOGDIR

echo "Alignement rate summary for all samples" > $OUTFILE
for FILE in $(ls $MAPPINGDIR/*mapping_summary.txt);
do
    SAMPLE=$(basename "$FILE" | sed 's/mapping_summary.txt//')
    echo "####################################################################################################################################################################" >> $OUTFILE
    echo $SAMPLE >> $OUTFILE
    cat $FILE >> $OUTFILE
    echo "####################################################################################################################################################################" >> $OUTFILE
done