#!/bin/bash

# Usage: bash organize_fastq.sh <fastq_dir> <sample_csv> <output_dir>
# Example: bash organize_fastq.sh ./fastq_files samples.csv ./matched_samples

FASTQ_DIR="${1}"
SAMPLE_CSV="${2}"
OUTPUT_DIR="${3:-./matched_samples}"

# --- Validate inputs ---
if [[ -z "$FASTQ_DIR" || -z "$SAMPLE_CSV" ]]; then
    echo "Usage: bash organize_fastq.sh <fastq_dir> <sample_csv> [output_dir]"
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_DIR"
    exit 1
fi

if [[ ! -f "$SAMPLE_CSV" ]]; then
    echo "ERROR: Sample CSV not found: $SAMPLE_CSV"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "FASTQ directory : $FASTQ_DIR"
echo "Sample CSV      : $SAMPLE_CSV"
echo "Output directory: $OUTPUT_DIR"
echo "----------------------------------------"

matched=0
missing=0

# Read each sample ID from the CSV (handles single-column or multi-column CSVs)
# Skips blank lines and strips carriage returns / surrounding whitespace
while IFS=',' read -r sample_id _rest; do
    # Clean up the sample ID
    sample_id=$(echo "$sample_id" | tr -d '\r' | xargs)

    # Skip empty lines or header-like lines if needed
    [[ -z "$sample_id" ]] && continue

    r1="${FASTQ_DIR}/${sample_id}_R1.fq.gz"
    r2="${FASTQ_DIR}/${sample_id}_R2.fq.gz"

    found=0

    if [[ -f "$r1" ]]; then
        cp "$r1" "$OUTPUT_DIR/"
        echo "  [COPIED] $(basename "$r1")"
        ((found++))
    fi

    if [[ -f "$r2" ]]; then
        cp "$r2" "$OUTPUT_DIR/"
        echo "  [COPIED] $(basename "$r2")"
        ((found++))
    fi

    if [[ $found -gt 0 ]]; then
        ((matched++))
    else
        echo "  [MISSING] No files found for sample: $sample_id"
        ((missing++))
    fi

done < "$SAMPLE_CSV"

echo "----------------------------------------"
echo "Done. Matched: $matched sample(s) | Missing: $missing sample(s)"
echo "Files copied to: $OUTPUT_DIR"
