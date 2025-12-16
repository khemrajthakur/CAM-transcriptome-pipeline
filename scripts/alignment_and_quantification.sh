#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Script: rnaseq_alignment_quantification.sh
# Purpose:
#   - Download HISAT2 genome index (GRCh38)
#   - Align trimmed paired-end reads using HISAT2
#   - Convert SAM → sorted BAM
#   - Quantify gene expression using featureCounts
#
# Usage:
#   bash rnaseq_alignment_quantification.sh \
#     <trimmed_reads_dir> \
#     <output_dir> \
#     <threads>
#
# Example:
#   bash rnaseq_alignment_quantification.sh \
#     results/trimmed_reads \
#     results/alignment \
#     8
# ============================================================

SECONDS=0

# ----------------------------
# Argument check
# ----------------------------
if [ "$#" -ne 3 ]; then
  echo "Usage: bash $0 <trimmed_reads_dir> <output_dir> <threads>"
  exit 1
fi

TRIMMED_READS="$1"
OUT_DIR="$2"
THREADS="$3"

# ----------------------------
# Directory structure
# ----------------------------
GENOME_DIR="${OUT_DIR}/genome"
ALIGN_DIR="${OUT_DIR}/hisat2"
BAM_DIR="${OUT_DIR}/bam"
COUNT_DIR="${OUT_DIR}/counts"
GTF_DIR="${OUT_DIR}/annotation"

mkdir -p "$GENOME_DIR" "$ALIGN_DIR" "$BAM_DIR" "$COUNT_DIR" "$GTF_DIR"

# ============================================================
# STEP 1: Download HISAT2 genome index (GRCh38)
# ============================================================
echo "Downloading HISAT2 genome index (GRCh38)..."

cd "$GENOME_DIR"
wget -q https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xzf grch38_genome.tar.gz
rm grch38_genome.tar.gz

GENOME_INDEX="${GENOME_DIR}/grch38/genome"

# ============================================================
# STEP 2: HISAT2 alignment
# ============================================================
echo "Starting HISAT2 alignment..."

for R1 in "$TRIMMED_READS"/*_R1_paired.fastq.gz; do
  R2="${R1/_R1_paired/_R2_paired}"

  if [ ! -f "$R2" ]; then
    echo "WARNING: Missing R2 for $R1. Skipping..."
    continue
  fi

  SAMPLE=$(basename "$R1" _R1_paired.fastq.gz)
  echo "Aligning sample: $SAMPLE"

  hisat2 -p "$THREADS" \
    -x "$GENOME_INDEX" \
    -1 "$R1" -2 "$R2" \
    --rna-strandness R \
    -S "$ALIGN_DIR/${SAMPLE}.sam"
done

# ============================================================
# STEP 3: Convert SAM → sorted BAM
# ============================================================
echo "Converting SAM to sorted BAM..."

for SAM in "$ALIGN_DIR"/*.sam; do
  SAMPLE=$(basename "$SAM" .sam)

  samtools view -@ "$THREADS" -bS "$SAM" | \
  samtools sort -@ "$THREADS" -o "$BAM_DIR/${SAMPLE}.sorted.bam"

  samtools index "$BAM_DIR/${SAMPLE}.sorted.bam"
  rm "$SAM"
done

# ============================================================
# STEP 4: Download annotation (GTF)
# ============================================================
echo "Downloading GTF annotation..."

cd "$GTF_DIR"
wget -q http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
gunzip Homo_sapiens.GRCh38.106.gtf.gz

GTF_FILE="${GTF_DIR}/Homo_sapiens.GRCh38.106.gtf"

# ============================================================
# STEP 5: Quantification with featureCounts
# ============================================================
echo "Running featureCounts..."

featureCounts -T "$THREADS" \
  -p -s 2 \
  -a "$GTF_FILE" \
  -o "$COUNT_DIR/gene_counts.txt" \
  "$BAM_DIR"/*.sorted.bam

# ============================================================
# Final report
# ============================================================
echo "Pipeline completed successfully!"
echo "Total elapsed time: ${SECONDS} seconds"
