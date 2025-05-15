#!/bin/bash

set -euo pipefail

# === Usage check ===
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <species_name> <reference_name> <base_output_dir> [--no-cache] [--mapq <value>] [--remove-temp]"
    exit 1
fi

# === Required arguments ===
SPECIES=$1
REFERENCE=$2
BASE_DIR=$3
shift 3  # Remove positional args

# === Optional flags ===
NO_CACHE=false
REMOVE_TEMP=false
MAPQ=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-cache)
            NO_CACHE=true
            shift
            ;;
        --remove-temp)
            REMOVE_TEMP=true
            shift
            ;;
        --mapq)
            if [[ -n "${2:-}" && "$2" != --* ]]; then
                MAPQ="$2"
                shift 2
            else
                echo "❌ --mapq requires a value"
                exit 1
            fi
            ;;
        *)
            echo "❌ Unknown argument: $1"
            exit 1
            ;;
    esac
done

# === Paths ===
SPECIES_FASTA="${BASE_DIR}/${SPECIES}/${SPECIES}.fasta"
REFERENCE_FASTA="${BASE_DIR}/${REFERENCE}/${REFERENCE}.fasta"
BAM_DIR="${BASE_DIR}/BAMs"
FASTQ="${BASE_DIR}/${SPECIES}/${SPECIES}_reads_to_${REFERENCE}.fastq"
RAW_BAM="${BAM_DIR}/${SPECIES}_to_${REFERENCE}_raw.bam"
BAM="${BAM_DIR}/${SPECIES}_to_${REFERENCE}.bam"
LOG="${BAM%.bam}.log"

# === Input checks ===
if [[ ! -f "$SPECIES_FASTA" ]]; then
    echo "❌ Species FASTA not found: $SPECIES_FASTA"
    exit 1
fi

if [[ ! -f "$REFERENCE_FASTA" ]]; then
    echo "❌ Reference FASTA not found: $REFERENCE_FASTA"
    exit 1
fi

# === Final BAM caching ===
if [[ -f "$BAM" && "$NO_CACHE" == false ]]; then
    echo "✅ Final BAM already exists: $BAM"
    exit 0
fi

mkdir -p "$BAM_DIR"

# === FASTQ caching ===
if [[ -f "$FASTQ" && "$NO_CACHE" == false ]]; then
    echo "📂 Reusing existing FASTQ file: $FASTQ"
else
    echo "🔧 Generating FASTQ from $SPECIES_FASTA → $FASTQ"
    python create_fastq_fragments.py "$SPECIES_FASTA" --output "$FASTQ" --force
fi

# === Prepare filter command ===
FILTER_CMD="python filter_sam.py -"
if [[ -n "$MAPQ" ]]; then
    FILTER_CMD+=" --mapq $MAPQ"
fi

CORES=$(nproc)

# === Alignment and filtering logic ===
if [[ "$REMOVE_TEMP" == false ]]; then
    # Use cached raw BAM if exists
    if [[ -f "$RAW_BAM" && "$NO_CACHE" == false ]]; then
        echo "📡 Using cached raw BAM: $RAW_BAM"
    else
        echo "📡 Aligning and saving raw BAM..."
        bwa mem "$REFERENCE_FASTA" "$FASTQ" \
            | samtools sort -@ "$CORES" -o "$RAW_BAM"
        samtools index "$RAW_BAM"
        echo "✅ Raw BAM saved: $RAW_BAM"
    fi

    echo "🔬 Filtering $RAW_BAM to produce final BAM..."
    samtools view -h "$RAW_BAM" \
        | eval "$FILTER_CMD" 2> "$LOG" \
        | samtools sort -@ "$CORES" -o "$BAM"

else
    echo "🔄 Aligning, filtering, and sorting in stream..."
    bwa mem "$REFERENCE_FASTA" "$FASTQ" \
        | eval "$FILTER_CMD" 2> "$LOG" \
        | samtools sort -@ "$CORES" -o "$BAM"
    if [[ "$REMOVE_TEMP" == true && -f "$FASTQ" ]]; then
        rm "$FASTQ"
    fi

fi

# === Index final BAM ===
echo "📌 Indexing final BAM..."
samtools index "$BAM"

# === Done ===
echo "✅ Finished: $BAM"
echo "📄 Filter stats written to: $LOG"
[[ "$REMOVE_TEMP" == false ]] && echo "📁 Unfiltered alignment saved: $RAW_BAM"
