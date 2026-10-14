#!/usr/bin/env bash
# A7.2a — Arma el set de proteomas de A7 (A7/prot/<código>.faa, encabezados ">código|id")
#   7 proteomas nuevos de Funannotate + 17 que se mantienen desde A3/prot/ (ya con isoforma más larga por gen).
# Control: completitud BUSCO de cada proteoma (compleasm protein, leotiomycetes_odb10) -> resultados/busco_proteomas.tsv
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
mkdir -p "$A7/prot" "$D/resultados"

for code in $MANTENER; do cp "$W/A3/prot/$code.faa" "$A7/prot/"; done
for code in $(echo "$GENOMAS_7" | awk '{print $1}'); do
  f=$(ls "$A7"/funannotate/$code/predict/predict_results/*.proteins.fa)
  awk -v c="$code" '/^>/{print ">" c "|" substr($1,2); next}{gsub(/\*/,""); print}' "$f" > "$A7/prot/$code.faa"
done
for f in "$A7"/prot/*.faa; do printf "%s\t%s\n" "$(basename "$f" .faa)" "$(grep -c '>' "$f")"; done \
  | sort > "$D/resultados/n_proteinas.tsv"
cat "$D/resultados/n_proteinas.tsv"

# BUSCO de proteomas
mkdir -p "$A7/busco_prot"
conda run -n "$ENV_BUSCO" compleasm download leotiomycetes -L "$A7/busco_prot/mb_lib"
printf "codigo\tresumen\n" > "$D/resultados/busco_proteomas.tsv"
for f in "$A7"/prot/*.faa; do
  c=$(basename "$f" .faa); o=$A7/busco_prot/$c
  [ -s "$o/summary.txt" ] || conda run -n "$ENV_BUSCO" compleasm protein -p "$f" -l leotiomycetes -L "$A7/busco_prot/mb_lib" \
      -t "$THREADS" -o "$o" > /dev/null 2>&1
  printf "%s\t%s\n" "$c" "$(grep -E '^S:|^D:|^F:|^M:' "$o/summary.txt" | tr '\n' ' ')" >> "$D/resultados/busco_proteomas.tsv"
done
cat "$D/resultados/busco_proteomas.tsv"
