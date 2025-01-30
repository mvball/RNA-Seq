#!/usr/bin/env bash
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=hisat_index
#SBATCH --output=/data/users/mschlaepfer/rnaseq_course/log/hisat_index%A.out
#SBATCH --error=/data/users/mschlaepfer/rnaseq_course/log/hisat_index%A.err
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=marina.schlaepferrubio@students.unibe.ch
#SBATCH --mail-type=END,FAIL

# Set the constants for directories and files
WORKDIR="/data/users/${USER}/rnaseq_course"
REFGENDIR="$WORKDIR/reference_genome"
LOGDIR="$WORKDIR/log"
INDEXDIR="$WORKDIR/index_hisat"
COMPRESSED_GENOMEFILE="$REFGENDIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
DECOMPRESSED_GENOMEFILE="$REFGENDIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Create directories if they don't exist
mkdir -p "$LOGDIR"
mkdir -p "$INDEXDIR"

# Check if the decompressed genome file already exists from a previous step
if [[ -f "$DECOMPRESSED_GENOMEFILE" ]]; then
    echo "Decompressed genome file already exists: $DECOMPRESSED_GENOMEFILE"
else
    echo "Decompressing genome file..."
    gunzip -k "$COMPRESSED_GENOMEFILE" || { echo "Error decompressing genome file. Exiting."; exit 1; }
fi

# Build the HISAT2 index
echo "Building HISAT2 index..."
apptainer exec --bind "$WORKDIR" \
    /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
    hisat2-build "$DECOMPRESSED_GENOMEFILE" "$INDEXDIR/genome_index" || { echo "Error building HISAT2 index. Exiting."; exit 1; }

echo "HISAT2 index creation completed successfully."
