#!/usr/bin/env bash
# A7.2b — OrthoFinder sobre los 24 proteomas de A7 (mismos parámetros que A3: diamond, MSA con MAFFT, FastTree).
# Salida: A7/prot/OrthoFinder/Results_A7/Orthogroups/Orthogroups.GeneCount.tsv
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
cd "$A7"
[ -d prot/OrthoFinder/Results_A7 ] || conda run -n "$ENV" orthofinder -f prot -n A7 -t "$THREADS" -a 8 \
    -S diamond -M msa -A mafft -T fasttree > orthofinder.log 2>&1
cp prot/OrthoFinder/Results_A7/Orthogroups/Orthogroups.GeneCount.tsv "$D/resultados/"
echo "ortogrupos: $(($(wc -l < "$D/resultados/Orthogroups.GeneCount.tsv") - 1))"
