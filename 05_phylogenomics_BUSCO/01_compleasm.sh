#!/usr/bin/env bash
# A1 paso 1: BUSCO leotiomycetes_odb10 (3,156 loci) con compleasm (miniprot + hmmsearch) sobre cada ensamble.
# Uniforme para genomas anotados y no anotados. Directorio de trabajo local (no Drive).
# Nota: odb12 (4,305 loci) resultó inviable en Mac local (>35 h por genoma); odb10 + binarios arm64 nativos.
# Ejecutar con: caffeinate -i bash 01_compleasm.sh
set -euo pipefail
W=~/hongo_comparativa
C=/Users/pablo/.local/share/mamba/envs/comparativa/bin/compleasm        # compleasm (python, osx-64)
A=/Users/pablo/.local/share/mamba/envs/comparativa_arm/bin   # miniprot/hmmsearch arm64 nativos
mkdir -p $W/compleasm
cd $W
for f in fna/*.fna; do
  code=$(basename "$f" .fna)
  out=compleasm/$code
  [ -f "$out/summary.txt" ] && { echo "skip $code"; continue; }
  echo "[$(date '+%d %H:%M')] $code"
  $C run -a "$f" -o "$out" -l leotiomycetes --odb odb10 -L mb_lib -t 11 \
     --miniprot_execute_path $A/miniprot --hmmsearch_execute_path $A/hmmsearch > "$out.log" 2>&1
done
echo "[$(date '+%d %H:%M')] fin"
