#!/bin/bash

set -euo pipefail

# Usage check
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <species_name> <reference_name> <base_output_dir> [--keep-temp] [--no-cache] [--mapq <value>] [--save-alignment]"
    exit 1
fi

# Required arguments
SPECIES=$1
REFERENCE=$2
BASE_DIR=$3

# Optional flags
KEEP_TEMP=false
NO_CACHE=false
SAVE_ALIGNMENT=false
MAPQ=""

i=4
while [ $i -le $# ]; do
    arg="${!i}"
    if [[ "$arg" == "--keep-temp" ]]; then
        KEEP_TEMP=true
    elif [[ "$arg" == "--no-cache" ]]; then
        NO_CACHE=true
    elif [[ "$arg" == "--save-alignment" ]]; then
        SAVE_ALIGNMENT=true
    elif [[ "$arg" == "--mapq" ]]; then
        next=$((i+1))
        MAPQ="${!next}"
        i=$((i+1))  # Skip next arg in loop
    fi
    i=$((i+1))
done

# Paths
SPECIES_FASTA="${BASE_DIR}/${SPECIES}/${SPECIES}.fasta"
REFERENCE_FASTA="${BASE_DIR}/${REFERENCE}/${REFERENCE}.fasta"
FASTQ="${BASE_DIR}/${SPECIES}/${SPECIES}_reads_to_${REFERENCE}.fastq"
RAW_BAM="${BASE_DIR}/${SPECIES}_to_${REFERENCE}_raw.bam"
BAM="${BASE_DIR}/${SPECIES}_to_${REFERENCE}.bam"
LOG="${BAM%.bam}.log"

# Skip if BAM already exists (unless --no-cache is set)
if [[ -f "$BAM" && "$NO_CACHE" == false ]]; then
    echo "✅ Final BAM already exists: $BAM"
    exit 0
fi

# Step 1: Generate FASTQ
echo "🔧 Generating FASTQ from $SPECIES_FASTA → $FASTQ"
python create_fastq_fragments.py "$SPECIES_FASTA" --output "$FASTQ" --force

# Step 2: Set filtering command
FILTER_CMD="python filter_sam.py -"
if [[ -n "$MAPQ" ]]; then
    FILTER_CMD+=" --mapq $MAPQ"
fi

# Step 3: Align, filter, and sort (with optional saving of raw BAM)
CORES=$(nproc)
if [[ "$SAVE_ALIGNMENT" == true ]]; then
    echo "📡 Saving raw BAM before filtering..."
    bwa mem "$REFERENCE_FASTA" "$FASTQ" \
      | samtools sort -@ "$CORES" -o "$RAW_BAM"
    samtools index "$RAW_BAM"
    echo "✅ Raw BAM saved: $RAW_BAM"

    echo "🔬 Filtering $RAW_BAM and writing final BAM..."
    samtools view -h "$RAW_BAM" \
      | eval "$FILTER_CMD" 2> "$LOG" \
      | samtools sort -@ "$CORES" -o "$BAM"
else
    echo "🔄 Aligning, filtering, and sorting in stream..."
    bwa mem "$REFERENCE_FASTA" "$FASTQ" \
      | eval "$FILTER_CMD" 2> "$LOG" \
      | samtools sort -@ "$CORES" -o "$BAM"
fi

# Step 4: Index BAM
echo "📌 Indexing final BAM..."
samtools index "$BAM"

# Step 5: Cleanup
if ! $KEEP_TEMP; then
    rm -f "$FASTQ"
    echo "🧹 Removed intermediate FASTQ: $FASTQ"
else
    echo "🗃️ Keeping intermediate FASTQ (--keep-temp enabled)"
fi

echo "✅ Finished: $BAM"
echo "📄 Filter stats written to: $LOG"
if [[ "$SAVE_ALIGNMENT" == true ]]; then
    echo "📁 Unfiltered alignment saved: $RAW_BAM"
fi
