#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <genome_name> <genome_accession> <output_directory>"
    exit 1
fi

# Assign input arguments to variables
GENOME_NAME=$1
ACCESSION=$2
BASE_OUTPUT_DIR=$3
OUTPUT_DIR="${BASE_OUTPUT_DIR}/${GENOME_NAME}"
TEMP_DIR="${BASE_OUTPUT_DIR}/temp_${GENOME_NAME}"

# Step 1: Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Step 2: Create a temporary directory for unzipping
mkdir -p "$TEMP_DIR"

# Step 3: Download the genome data using NCBI datasets CLI
echo "Downloading genome data for accession: $ACCESSION"
datasets download genome accession "$ACCESSION" --filename "${TEMP_DIR}/${GENOME_NAME}.zip"

# Check if the download was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to download genome data for accession $ACCESSION"
    rm -rf "$TEMP_DIR"
    exit 1
fi

# Step 4: Unzip the downloaded file into the temporary folder
echo "Unzipping genome data to $TEMP_DIR..."
unzip -q "${TEMP_DIR}/${GENOME_NAME}.zip" -d "$TEMP_DIR"

# Step 5: Locate the FASTA file within the unzipped directory
FASTA_FILE=$(find "$TEMP_DIR" -name "*.fna" -print -quit)

# Check if a FASTA file was found
if [ -z "$FASTA_FILE" ]; then
    echo "Error: No FASTA file found in the downloaded data."
    rm -rf "$TEMP_DIR"
    exit 1
fi

# Step 6: Move the FASTA file to the final output directory
mv "$FASTA_FILE" "${OUTPUT_DIR}/${GENOME_NAME}.fasta"

# Step 7: Remove the temporary directory and zip file
echo "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

# Step 8: Notify the user of success
echo "FASTA file for $GENOME_NAME has been extracted and saved in: ${OUTPUT_DIR}/${GENOME_NAME}.fasta"
