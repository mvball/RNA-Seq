#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --partition=pibu_el8
#SBATCH --mem=10M
#SBATCH --cpus-per-task=1
#SBATCH --job-name=get_reference
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/get_reference%A.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/get_reference%A.err
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL

#Setting the constant for the directories and required files
WORKDIR="/data/users/${USER}/rnaseq_course"
REFGENDIR="$WORKDIR/reference_genome"
LOGDIR="$WORKDIR/log"
REFGENOMEFILE="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
ANNOTATIONFILE="Homo_sapiens.GRCh38.113.gtf.gz"

#Create the directory for the error and output file if not present
mkdir -p $LOGDIR

#Create the directory where the reference genome and the annotation will be downloaded
mkdir -p $REFGENDIR

#enter the created folder for the reference genome and download the fa and gtf file from ensembl there
cd $REFGENDIR
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/$REFGENOMEFILE
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/$ANNOTATIONFILE


#Perform the checksum for the files to be sure there was no error/mixup of files
echo "Checksum for fasta file"
sum $REFGENDIR/$REFGENOMEFILE
echo "Checksum for gtf file"
sum $REFGENDIR/$ANNOTATIONFILE

#unzip the reference genome for later use
gunzip $REFGENDIR/$REFGENOMEFILE

#unzip the annotation file for later use
gunzip $REFGENDIR/$ANNOTATIONFILE