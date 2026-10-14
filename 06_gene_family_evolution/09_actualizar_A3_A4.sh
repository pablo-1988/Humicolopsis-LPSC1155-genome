#!/usr/bin/env bash
# A7.6b — Rehace los conteos del repertorio de P (A3 + A4) con los proteomas de A7 (7 genomas reanotados con Funannotate)
# sin tocar los resultados originales: todo sale en carpetas *_A7prot (PROT_SET=_A7prot).
#   1. A3: dominios Pfam (hmmsearch --cut_ga) + referencias (DIAMOND) sobre A7/prot; reutiliza el OrthoFinder de A7
#   2. A3: familias por ortogrupo          -> A3_repertorio_fosforo_comparativo/resultados_A7prot/
#   3. A4: árboles génicos (20 familias en paralelo) -> A4_curacion_candidatos/arboles_A7prot/
#   4. A4: clados ortólogos + tabla final  -> A4_curacion_candidatos/resultados_A7prot/
# Uso: bash 09_actualizar_A3_A4.sh
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
export PROT_SET=_A7prot
BIN=/home/ubuntu/miniforge3/envs/comparativa/bin; PY=$BIN/python; export A4_BIN=$BIN A4_THREADS=3
A3D=$PROJ/Paper/analisis/A3_repertorio_fosforo_comparativo; A4D=$PROJ/Paper/analisis/A4_curacion_candidatos
X=$W/A3$PROT_SET; mkdir -p "$X" "$W/A4$PROT_SET"
[ -e "$X/prot" ] || ln -s "$A7/prot" "$X/prot"          # mismos 24 proteomas (">código|id") + OrthoFinder Results_A7

echo "[$(date +%T)] 1/4 dominios Pfam y referencias"
cat "$A7"/prot/*.faa > "$X/all_prot.faa"
[ -s "$X/pfam_familias.domtbl" ] || $BIN/hmmsearch --cut_ga --cpu "$THREADS" --noali --domtblout "$X/pfam_familias.domtbl" \
    "$A3D/hmm/familias.hmm" "$X/all_prot.faa" > /dev/null
$BIN/diamond makedb --in "$X/all_prot.faa" -d "$X/all_prot" --quiet
$BIN/diamond blastp -q "$A3D/refs/referencias.faa" -d "$X/all_prot" --very-sensitive -k 500 -e 1e-5 -p "$THREADS" \
  --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore -o "$X/refs_vs_all.tsv" --quiet

echo "[$(date +%T)] 2/4 familias por ortogrupo (A3)"
$PY "$A3D/04_repertorio_familias.py" > "$X/04_repertorio.log"
tail -1 "$X/04_repertorio.log" | head -1; grep -m1 OrthoFinder "$X/04_repertorio.log"

echo "[$(date +%T)] 3/4 árboles génicos (A4), 20 familias en paralelo"
FAMS=$($PY - <<EOF
import ast, re
s = open("$A4D/01_arboles_familia.py").read()
print(" ".join(ast.literal_eval(re.search(r"FAM = (\{.*?\n\})", s, re.S).group(1)).keys()))
EOF
)
echo $FAMS | tr ' ' '\n' | xargs -P 20 -I{} sh -c "cd '$A4D' && $PY 01_arboles_familia.py {} > '$W/A4$PROT_SET/{}.log' 2>&1 || echo 'FALLÓ {}'"

echo "[$(date +%T)] 4/4 clados ortólogos + tabla final (A4)"
( cd "$A4D" && $PY 02_clados_ortologos.py && $PY 03_tabla_final.py ) | tee "$X/A4_tabla.log"
echo "[$(date +%T)] listo: $A4D/resultados$PROT_SET/tabla_final_repertorio.tsv"
