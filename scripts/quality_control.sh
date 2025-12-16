#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Script: rnaseq_qc_trimming.sh
# Purpose: Perform RNA-seq quality control and trimming
#          - FastQC on raw reads
#          - MultiQC summary
#          - Trimmomatic trimming (paired-end)
#          - FastQC on trimmed reads
#
# Usage:
#   bash rnaseq_qc_trimming.sh \
#     <raw_reads_dir> \
#     <output_dir> \
#     <trimmomatic_jar> \
#     <adapter_fasta>
#
# Example:
#   bash rnaseq_qc_trimming.sh \
#     data/reads \
#     qc_results \
#     tools/Trimmomatic/trimmomatic.jar \
#     tools/Trimmomatic/adapters/TruSeq2-PE.fa
# ============================================================

SECONDS=0

# ----------------------------
# Argument check
# ----------------------------
if [ "$#" -ne 4 ]; then
  echo "Usage: bash $0 <raw_reads_dir> <output_dir> <trimmomatic_jar> <adapter_fasta>"
  exit 1
fi

RAW_READS_DIR="$1"
OUT_DIR="$2"
TRIMMOMATIC_JAR="$3"
ADAPTERS="$4"

# ----------------------------
# Directory structure
# ----------------------------
RAW_FASTQC="${OUT_DIR}/fastqc_raw"
RAW_MULTIQC="${OUT_DIR}/multiqc_raw"
TRIMMED_READS="${OUT_DIR}/trimmed_reads"
TRIM_FASTQC="${OUT_DIR}/fastqc_trimmed"

mkdir -p "$RAW_FASTQC" "$RAW_MULTIQC" "$TRIMMED_READS" "$TRIM_FASTQC"

# ============================================================
# STEP 1: FastQC on raw reads
# ============================================================
echo "Running FastQC on raw reads..."
fastqc "$RAW_READS_DIR"/*.fastq.gz -o "$RAW_FASTQC"

# ============================================================
# STEP 2: MultiQC on raw FastQC results
# ============================================================
echo "Running MultiQC on raw FastQC outputs..."
multiqc "$RAW_FASTQC" -o "$RAW_MULTIQC"

# ============================================================
# STEP 3: Trimmomatic (paired-end)
# ============================================================
echo "Running Trimmomatic..."

for R1 in "$RAW_READS_DIR"/*_R1_*.fastq.gz; do
  R2="${R1/_R1_/_R2_}"

  if [ ! -f "$R2" ]; then
    echo "WARNING: Missing R2 for $R1. Skipping..."
    continue
  fi

  SAMPLE=$(basename "$R1" _R1_*.fastq.gz)

  echo "Processing sample: $SAMPLE"

  java -jar "$TRIMMOMATIC_JAR" PE -phred33 \
    "$R1" "$R2" \
    "$TRIMMED_READS/${SAMPLE}_R1_paired.fastq.gz" \
    "$TRIMMED_READS/${SAMPLE}_R1_unpaired.fastq.gz" \
    "$TRIMMED_READS/${SAMPLE}_R2_paired.fastq.gz" \
    "$TRIMMED_READS/${SAMPLE}_R2_unpaired.fastq.gz" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done

# ============================================================
# STEP 4: FastQC on trimmed reads
# ============================================================
echo "Running FastQC on trimmed reads..."
fastqc "$TRIMMED_READS"/*_paired.fastq.gz -o "$TRIM_FASTQC"

echo "RNA-seq QC and trimming completed in ${SECONDS} seconds."
