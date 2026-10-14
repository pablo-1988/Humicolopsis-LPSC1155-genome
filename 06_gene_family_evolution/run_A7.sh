#!/usr/bin/env bash
# A7 completo en el PC de 256 GB / 70 núcleos. Cada paso se salta si su salida ya existe (se puede relanzar).
# Uso:  nohup bash run_A7.sh > ~/hongo_comparativa/A7_run.log 2>&1 &
set -euo pipefail
H=$(dirname "$(readlink -f "$0")")
source "$H/00_config.sh"
mkdir -p "$A7"
echo "[$(date)] 1/5 Funannotate (7 genomas, ~6-8 h)";        bash "$H/01_reanotar_funannotate.sh"
echo "[$(date)] 2/5 Proteomas + BUSCO";                       bash "$H/02_proteomas.sh"
echo "[$(date)] 3/5 OrthoFinder y eggNOG-mapper en paralelo"
bash "$H/03_orthofinder.sh" & P1=$!
bash "$H/06_eggnog.sh"      & P2=$!
wait $P1; wait $P2
echo "[$(date)] 4/5 CAFE5";                                   bash "$H/05_cafe5.sh"
echo "[$(date)] 5/5 Enriquecimiento"
python3 "$H/07_enriquecimiento.py" "$A7/cafe/final" "$A7/prot/OrthoFinder/Results_A7/Orthogroups/Orthogroups.tsv" \
    "$A7/eggnog" "$FUNANNOTATE_DB_HOST/go.obo" "$D/resultados/enriquecimiento"
echo "[$(date)] A7 terminado. Resultados en $D/resultados/"
