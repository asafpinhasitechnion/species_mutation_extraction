#!/bin/bash

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <species_name> <base_output_dir>"
    exit 1
fi

# Input arguments
SPECIES_NAME=$1
BASE_DIR=$2
FASTA_FILE="${BASE_DIR}/${SPECIES_NAME}/${SPECIES_NAME}.fasta"

# Check tool availability
command -v bwa >/dev/null 2>&1 || { echo >&2 "❌ bwa not found in PATH."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "❌ samtools not found in PATH."; exit 1; }

# Check that FASTA exists
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "❌ FASTA file not found: $FASTA_FILE"
    exit 1
fi

# Clean old BWA indices
rm -f "${FASTA_FILE}".{amb,ann,bwt,pac,sa}

# Run BWA index
echo "🔧 Indexing with BWA..."
bwa index "$FASTA_FILE" || { echo "❌ BWA indexing failed."; exit 1; }

# Run Samtools faidx
echo "🔧 Indexing with Samtools..."
samtools faidx "$FASTA_FILE" || { echo "❌ Samtools indexing failed."; exit 1; }

echo "✅ Indexing complete for $SPECIES_NAME in $BASE_DIR"
