#!/bin/bash

set -euo pipefail

# === ARGUMENT PARSING ===
if [ "$#" -lt 6 ]; then
    echo "Usage: $0 <taxa1_name> <taxa1_accession> <taxa2_name> <taxa2_accession> <outgroup_name> <outgroup_accession> [optional args...]"
    echo "Optional args: --keep-temp --no-cache --mapq <value> --save-alignment --no-full-mutations"
    exit 1
fi

# === CORE ARGUMENTS ===
T1_NAME=$1
T1_ACC=$2
T2_NAME=$3
T2_ACC=$4
OUT_NAME=$5
OUT_ACC=$6
shift 6  # Shift to access optional args

# === BASE DIRECTORY SETUP ===
RUN_ID="${T1_NAME}__${T2_NAME}__${OUT_NAME}"
BASE_OUTPUT_DIR="../Output/${RUN_ID}"
mkdir -p "$BASE_OUTPUT_DIR"

echo "📁 Run ID: $RUN_ID"
echo "📂 Base output directory: $BASE_OUTPUT_DIR"

# === GENOME DOWNLOADS ===
echo "⬇️ Downloading genomes..."
bash download_genome.sh "$T1_NAME" "$T1_ACC" "$BASE_OUTPUT_DIR"
bash download_genome.sh "$T2_NAME" "$T2_ACC" "$BASE_OUTPUT_DIR"
bash download_genome.sh "$OUT_NAME" "$OUT_ACC" "$BASE_OUTPUT_DIR"

# === REFERENCE INDEXING ===
echo "🧬 Indexing outgroup genome: $OUT_NAME"
bash index_reference_genome.sh "$OUT_NAME" "$BASE_OUTPUT_DIR"

# === ALIGNMENT + FILTERING ===
echo "🔗 Aligning $T1_NAME to $OUT_NAME"
bash align_and_filter.sh "$T1_NAME" "$OUT_NAME" "$BASE_OUTPUT_DIR" "$@"

echo "🔗 Aligning $T2_NAME to $OUT_NAME"
bash align_and_filter.sh "$T2_NAME" "$OUT_NAME" "$BASE_OUTPUT_DIR" "$@"

echo "✅ Alignment and filtering complete for $RUN_ID"

# === PILEUP CREATION ===
echo "📊 Creating pileup..."
bash create_pileup.sh "$OUT_NAME" "$T1_NAME" "$T2_NAME" "$BASE_OUTPUT_DIR" "$@"


# Optional argument for full output CSV
FULL_OUTPUT_ARG=""
for ((i=1; i<=$#; i++)); do
    if [[ "${!i}" == "--full-output-dir" ]]; then
        j=$((i+1))
        FULL_OUTPUT_ARG="--full-output-dir ${!j}"
        break
    fi
done

# === MUTATION EXTRACTION ===
echo "🧪 Extracting mutations..."
python extract_mutations.py "$OUT_NAME" "$T1_NAME" "$T2_NAME" \
  --pileup-dir "$BASE_OUTPUT_DIR" \
  --output-dir "$BASE_OUTPUT_DIR/Mutations" \
  $FULL_OUTPUT_ARG

# === TRIPLET CONTEXT EXTRACTION ===
echo "🧬 Extracting triplet contexts..."
python extract_triplets.py "$OUT_NAME" "$T1_NAME" "$T2_NAME" \
  --pileup-dir "$BASE_OUTPUT_DIR" \
  --output-dir "$BASE_OUTPUT_DIR/Triplets"

