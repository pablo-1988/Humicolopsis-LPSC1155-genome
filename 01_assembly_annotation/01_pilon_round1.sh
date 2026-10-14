#!/bin/bash
set -euo pipefail

GENOME="/Users/pablo/Documents/Proyecto-Hongo/SPAdes/contigs_filtered.fasta"
READ1="/Users/pablo/Documents/Proyecto-Hongo/cleaned-hongo/sub_set2/hongo_sub2_1.fq.gz"
READ2="/Users/pablo/Documents/Proyecto-Hongo/cleaned-hongo/sub_set2/hongo_sub2_2.fq.gz"
OUTPUT_DIR="/Users/pablo/Documents/Proyecto-Hongo/Pilon"
OUTPUT_NAME="hongo_pilon"
THREADS=10

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

SORTED_BAM="aln.primary.q20.sorted.bam"

# Index referencia (una vez)
if [ ! -f "${GENOME}.bwt" ]; then
  bwa index "$GENOME"
fi

# BAM (recomendado: primary + MAPQ>=20)
if [ ! -f "$SORTED_BAM" ]; then
  bwa mem -t "$THREADS" "$GENOME" "$READ1" "$READ2" \
    | samtools view -@ "$THREADS" -b -q 20 -F 0x900 \
    | samtools sort -@ "$THREADS" -o "$SORTED_BAM"
  samtools index "$SORTED_BAM"
fi

# Heap Java (ajusta: 10g, 12g, 14g…)
export JAVA_TOOL_OPTIONS="-Xmx12g -XX:+UseG1GC -XX:MaxGCPauseMillis=200"

# Pilon con chunks pequeños (menos RAM)
pilon \
  --genome "$GENOME" \
  --frags "$SORTED_BAM" \
  --output "$OUTPUT_NAME" \
  --outdir "$OUTPUT_DIR" \
  --vcf --changes --tracks \
  --fix snps,indels \
  --chunksize 10000
