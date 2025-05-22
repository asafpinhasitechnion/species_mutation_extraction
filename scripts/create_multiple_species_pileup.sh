#!/bin/bash

set -euo pipefail

# === Usage ===
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <reference_name> <base_output_dir> <taxa1> [<taxa2> ...] [--no-cache]"
    exit 1
fi

# === Required arguments ===
REFERENCE=$1
BASE_OUTPUT_DIR=$2
shift 2

# === Handle taxa and optional args ===
NO_CACHE=false
TAXA=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-cache)
            NO_CACHE=true
            shift
            ;;
        *)
            TAXA+=("$1")
            shift
            ;;
    esac
done

if [ "${#TAXA[@]}" -lt 3 ]; then
    echo "❌ At least three taxa must be provided."
    exit 1
fi

# === Paths ===
REF_FASTA="${BASE_OUTPUT_DIR}/${REFERENCE}/${REFERENCE}.fasta"
BAM_FOLDER="${BASE_OUTPUT_DIR}/BAMs"

joined=$(printf "__%s" "${TAXA[@]}")
joined=${joined:2}  # Remove the leading __
PILEUP_FILE="${BASE_OUTPUT_DIR}/${REFERENCE}__${joined}.pileup.gz"

# === Check input reference ===
if [[ ! -f "$REF_FASTA" ]]; then
    echo "❌ Reference FASTA not found: $REF_FASTA"
    exit 1
fi

# === Build list of BAM files ===
BAM_FILES=()
for TAXON in "${TAXA[@]}"; do
    BAM="${BAM_FOLDER}/${TAXON}_to_${REFERENCE}.bam"
    if [[ ! -f "$BAM" ]]; then
        echo "❌ Required BAM not found: $BAM"
        exit 1
    fi
    BAM_FILES+=("$BAM")
done

# === Skip if exists and caching is allowed ===
if [[ -f "$PILEUP_FILE" && "$NO_CACHE" == false ]]; then
    echo "✅ Pileup already exists: $PILEUP_FILE"
    exit 0
fi

# === Run mpileup ===
echo "🔬 Generating pileup:"
echo "    Reference: $REFERENCE"
echo "    Taxa: ${TAXA[*]}"
echo "    Output: $PILEUP_FILE"

samtools mpileup -f "$REF_FASTA" -B -d 100 "${BAM_FILES[@]}" | gzip > "$PILEUP_FILE"

echo "✅ Pileup written to: $PILEUP_FILE"
