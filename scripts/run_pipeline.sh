#!/bin/bash

set -euo pipefail

# === ARGUMENT PARSING ===
if [ "$#" -lt 6 ]; then
    echo "Usage: $0 <taxa1_name> <taxa1_accession> <taxa2_name> <taxa2_accession> <outgroup_name> <outgroup_accession> [optional args...]"
    echo "Optional args: --keep-temp --no-cache --mapq <value> --save-alignment --no-full-mutations --full-output-dir <subdir>"
    exit 1
fi

# === CORE ARGUMENTS ===
T1_NAME=$1
T1_ACC=$2
T2_NAME=$3
T2_ACC=$4
OUT_NAME=$5
OUT_ACC=$6
shift 6

# === BASE DIRECTORY SETUP ===
RUN_ID="${T1_NAME}__${T2_NAME}__${OUT_NAME}"
BASE_OUTPUT_DIR="../Output/${RUN_ID}"
mkdir -p "$BASE_OUTPUT_DIR"

echo "📁 Run ID: $RUN_ID"
echo "📂 Base output directory: $BASE_OUTPUT_DIR"

# === COLLECT OPTIONAL ARGUMENTS ===
ALIGN_FILTER_ARGS=()
PILEUP_ARGS=()
MUTATION_ARGS=()
TRIPLET_ARGS=()
FULL_OUTPUT_SUBDIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --keep-temp|--save-alignment)
            ALIGN_FILTER_ARGS+=("$1")
            PILEUP_ARGS+=("$1")
            shift
            ;;
        --no-cache)
            ALIGN_FILTER_ARGS+=("$1")
            PILEUP_ARGS+=("$1")
            MUTATION_ARGS+=("$1")
            TRIPLET_ARGS+=("$1")
            shift
            ;;
        --mapq)
            ALIGN_FILTER_ARGS+=("$1" "$2")
            PILEUP_ARGS+=("$1" "$2")
            shift 2
            ;;
        --no-full-mutations)
            MUTATION_ARGS+=("$1")
            shift
            ;;
        --full-output-dir)
            shift
            FULL_OUTPUT_SUBDIR="$1"
            MUTATION_ARGS+=(--full-output-dir "$BASE_OUTPUT_DIR/$FULL_OUTPUT_SUBDIR")
            shift
            ;;
        *)
            echo "❗ Unknown argument: $1"
            exit 1
            ;;
    esac
done



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
bash align_and_filter.sh "$T1_NAME" "$OUT_NAME" "$BASE_OUTPUT_DIR" "${ALIGN_FILTER_ARGS[@]}"

echo "🔗 Aligning $T2_NAME to $OUT_NAME"
bash align_and_filter.sh "$T2_NAME" "$OUT_NAME" "$BASE_OUTPUT_DIR" "${ALIGN_FILTER_ARGS[@]}"

echo "✅ Alignment and filtering complete for $RUN_ID"

# === PILEUP CREATION ===
echo "📊 Creating pileup..."
bash create_pileup.sh "$OUT_NAME" "$T1_NAME" "$T2_NAME" "$BASE_OUTPUT_DIR" "${PILEUP_ARGS[@]}"

# === MUTATION EXTRACTION ===
echo "🧪 Extracting mutations..."
python extract_mutations.py "$OUT_NAME" "$T1_NAME" "$T2_NAME" \
  --pileup-dir "$BASE_OUTPUT_DIR" \
  --output-dir "$BASE_OUTPUT_DIR/Mutations" \
  "${MUTATION_ARGS[@]}"

# === TRIPLET CONTEXT EXTRACTION ===
echo "🧬 Extracting triplet contexts..."
python extract_triplets.py "$OUT_NAME" "$T1_NAME" "$T2_NAME" \
  --pileup-dir "$BASE_OUTPUT_DIR" \
  --output-dir "$BASE_OUTPUT_DIR/Triplets" \
  "${TRIPLET_ARGS[@]}"

# === MUTATION NORMALIZATION ===
echo "📐 Normalizing mutation spectra..."
python normalize_extracted_mutations.py --input-dir "$BASE_OUTPUT_DIR"

# === PLOTTING MUTATION AND TRIPLET SPECTRA ===
echo "🖼️ Plotting mutation and triplet spectra..."
python plot_spectra.py --input-dir "$BASE_OUTPUT_DIR/Tables"

# === INTERVAL EXTRACTION ===
echo "📏 Extracting genomic intervals..."
INTERVAL_DIR="$BASE_OUTPUT_DIR/Intervals"
mkdir -p "$INTERVAL_DIR"

python bam_to_intervals.py "$BASE_OUTPUT_DIR/${T1_NAME}_to_${OUT_NAME}.bam" "$INTERVAL_DIR"
python bam_to_intervals.py "$BASE_OUTPUT_DIR/${T2_NAME}_to_${OUT_NAME}.bam" "$INTERVAL_DIR"



# === GENOMIC DENSITY PLOTS ===
echo "📊 Plotting read coverage and mutation distributions..."
CHROMOSOME="NW_021628937.1"
FAI_FILE="$BASE_OUTPUT_DIR/$OUT_NAME/${OUT_NAME}.fasta.fai"
PLOT_DIR="$BASE_OUTPUT_DIR/Plots"
mkdir -p "$PLOT_DIR"

# Plot coverage only
python plot_genomic_coverage.py \
  --output-dir "$BASE_OUTPUT_DIR" \
  --chromosome "$CHROMOSOME" \
  --bin-size 100000 \
  --slide 100000

# Plot mutation density
python plot_genomic_mutation_distribution.py \
  --mutation-dir "$BASE_OUTPUT_DIR/Mutations" \
  --fai-file "$FAI_FILE" \
  --chromosome "$CHROMOSOME" \
  --output-dir "$PLOT_DIR" \
  --bin-size 100000 \
  --coverage-dir "$BASE_OUTPUT_DIR/Intervals"

# Plot CpG C>T mutations specifically
python plot_genomic_mutation_distribution.py \
  --mutation-dir "$BASE_OUTPUT_DIR/Mutations" \
  --fai-file "$FAI_FILE" \
  --chromosome "$CHROMOSOME" \
  --output-dir "$PLOT_DIR" \
  --bin-size 100000 \
  --coverage-dir "$BASE_OUTPUT_DIR/Intervals" \
  --mutation_category "[ACGT]\\[C>[T]\\]G"

# Plot all CpG mutations
python plot_genomic_mutation_distribution.py \
  --mutation-dir "$BASE_OUTPUT_DIR/Mutations" \
  --fai-file "$FAI_FILE" \
  --chromosome "$CHROMOSOME" \
  --output-dir "$PLOT_DIR" \
  --bin-size 100000 \
  --coverage-dir "$BASE_OUTPUT_DIR/Intervals" \
  --mutation_category "[ACGT]\\[C>[ACGT]\\]G"
