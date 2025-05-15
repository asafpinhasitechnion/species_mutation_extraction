#!/bin/bash

set -euo pipefail

# === Usage ===
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <reference_name> <taxa1_name> <taxa2_name> <base_output_dir> [--no-cache]"
    exit 1
fi

# === Input arguments ===
REFERENCE=$1
TAXA1=$2
TAXA2=$3
BASE_OUTPUT_DIR=$4

# === Optional flags ===
NO_CACHE=false
if [[ "${5:-}" == "--no-cache" ]]; then
    NO_CACHE=true
fi

# === Paths ===
REF_FASTA="${BASE_OUTPUT_DIR}/${REFERENCE}/${REFERENCE}.fasta"
BAM1="${BASE_OUTPUT_DIR}/${TAXA1}_to_${REFERENCE}.bam"
BAM2="${BASE_OUTPUT_DIR}/${TAXA2}_to_${REFERENCE}.bam"


PILEUP_FILE="${BASE_OUTPUT_DIR}/${REFERENCE}__${TAXA1}__${TAXA2}.pileup.gz"

# === Check existence of input files ===
for FILE in "$REF_FASTA" "$BAM1" "$BAM2"; do
    if [[ ! -f "$FILE" ]]; then
        echo "❌ Required file not found: $FILE"
        exit 1
    fi
done

# === Skip if exists ===
if [[ -f "$PILEUP_FILE" && "$NO_CACHE" == false ]]; then
    echo "✅ Pileup already exists: $PILEUP_FILE"
    exit 0
fi

# === Generate pileup ===
echo "🔬 Generating pileup:"
echo "    Reference: $REFERENCE"
echo "    Taxa: $TAXA1, $TAXA2"
echo "    Output: $PILEUP_FILE"

samtools mpileup -f "$REF_FASTA" -B -d 100 "$BAM1" "$BAM2" | gzip > "$PILEUP_FILE"

echo "✅ Pileup written to: $PILEUP_FILE"
