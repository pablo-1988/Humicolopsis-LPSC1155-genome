#!/bin/bash
set -euo pipefail

# ====== INPUTS ======
GENOME="/Users/pablo/Documents/Proyecto-Hongo/Pilon/hongo_assembly.pilon1.fasta"   # genoma pulido ronda 1
READ1="/Users/pablo/Documents/Proyecto-Hongo/cleaned-hongo/sub_set2/hongo_sub2_1.fq.gz"
READ2="/Users/pablo/Documents/Proyecto-Hongo/cleaned-hongo/sub_set2/hongo_sub2_2.fq.gz"

OUTDIR="/Users/pablo/Documents/Proyecto-Hongo/Pilon"
PREFIX="hongo_pilon2"
THREADS=10
JAVA_MEM="10g"   # sube a 16g si tu Mac aguanta

# ====== WORKDIR ======
mkdir -p "$OUTDIR"
cd "$OUTDIR"

BAM="${PREFIX}.q20.sorted.bam"

# ====== INDEX + ALIGN ======
if [ ! -f "${GENOME}.bwt" ]; then
  bwa index "$GENOME"
fi

bwa mem -t "$THREADS" "$GENOME" "$READ1" "$READ2" \
  | samtools view -@ 4 -b -q 20 \
  | samtools sort -@ 4 -o "$BAM"

samtools index "$BAM"

# QC rápido
samtools flagstat "$BAM" > "${PREFIX}.flagstat.txt"

# ====== PILON ======
export _JAVA_OPTIONS="-Xmx${JAVA_MEM}"

pilon \
  --genome "$GENOME" \
  --frags "$BAM" \
  --output "$PREFIX" \
  --outdir "$OUTDIR" \
  --changes --vcf \
  --fix snps,indels

