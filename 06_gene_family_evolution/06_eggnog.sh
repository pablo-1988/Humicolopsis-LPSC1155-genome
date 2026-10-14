#!/usr/bin/env bash
# A7.5a — Anotación funcional uniforme de los 24 proteomas con eggNOG-mapper v2 (GO, KEGG, Pfam).
# Una sola corrida sobre todos los proteomas concatenados; --dbmem carga la base (~45 GB) en RAM.
# Se puede correr en paralelo con 03_orthofinder.sh (ambos usan solo A7/prot/).
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
mkdir -p "$EGGNOG_DB" "$A7/eggnog"
# download_eggnog_data.py (v2.1.13) puede "terminar bien" dejando archivos vacíos si falla el DNS del servidor.
# Se descarga directo de eggnog5.embl.de y se verifica que la base quede completa antes de seguir.
if [ ! -s "$EGGNOG_DB/eggnog.db" ] || [ ! -s "$EGGNOG_DB/eggnog_proteins.dmnd" ]; then
  U=http://eggnog5.embl.de/download/emapperdb-5.0.2
  ( cd "$EGGNOG_DB"
    # el servidor corta descargas largas: se reanuda (-C -) hasta completar cada archivo
    for f in eggnog.db.gz eggnog.taxa.tar.gz eggnog_proteins.dmnd.gz; do
      for i in $(seq 1 200); do curl -sS -L -C - -o "$f" "$U/$f" && break; sleep 10; done
    done
    gunzip -f eggnog.db.gz; tar -zxf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz; gunzip -f eggnog_proteins.dmnd.gz )
fi
[ -s "$EGGNOG_DB/eggnog.taxa.db" ] || { echo "ERROR: base eggNOG incompleta (eggnog.taxa.db)"; exit 1; }
for f in eggnog.db eggnog_proteins.dmnd; do
  [ "$(stat -c %s "$EGGNOG_DB/$f" 2>/dev/null || echo 0)" -gt 100000000 ] || { echo "ERROR: base eggNOG incompleta ($f)"; exit 1; }
done
cat "$A7"/prot/*.faa > "$A7/eggnog/todos.faa"
[ -s "$A7/eggnog/todos.emapper.annotations" ] || conda run -n "$ENV" emapper.py -i "$A7/eggnog/todos.faa" \
    -o todos --output_dir "$A7/eggnog" --data_dir "$EGGNOG_DB" -m diamond --dbmem --cpu "$THREADS" --temp_dir "$A7/eggnog" --override \
    > "$A7/eggnog/emapper.log" 2>&1
echo "proteínas anotadas: $(grep -vc '^#' "$A7/eggnog/todos.emapper.annotations")"
