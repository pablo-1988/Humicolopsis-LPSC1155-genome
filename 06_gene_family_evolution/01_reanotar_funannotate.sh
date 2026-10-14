#!/usr/bin/env bash
# A7.1 — Reanotación de los 7 genomas sin anotación con el MISMO protocolo que LPSC 1155
# (Proyecto-Hongo/protocolo_anotacion_humicolopsis_M4.txt):
#   renombrar scaffolds (scaffold_N) -> funannotate mask -m tantan -> funannotate predict
#   (--busco_seed_species neurospora_crassa --busco_db sordariomycetes_odb10 --optimize_augustus)
# Uso:  bash 01_reanotar_funannotate.sh            # los 7 en paralelo (FA_JOBS x FA_CPUS)
#       bash 01_reanotar_funannotate.sh Allacyti   # uno solo
# Requiere docker (usuario en el grupo docker) o, si docker no está, Apptainer sin root (misma imagen, ver 00_config.sh),
# y `datasets` de NCBI si faltan los genomas.
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
mkdir -p "$A7/funannotate" "$W/fna"

runtime() { if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then echo docker; else echo apptainer; fi; }

uno() {
  local code=$1 acc sp strain
  read -r _ acc sp strain < <(echo "$GENOMAS_7" | awk -v c="$code" '$1==c')
  local out=$A7/funannotate/$code
  mkdir -p "$out"; cd "$out"
  if [ -s predict/predict_results/*.proteins.fa ] 2>/dev/null; then echo "$code: ya anotado"; return; fi

  # 1. genoma (si falta, descargar de NCBI)
  if [ ! -s "$W/fna/$code.fna" ]; then
    datasets download genome accession "$acc" --include genome --filename "$code.zip"
    unzip -o -q "$code.zip" -d "${code}_ncbi"
    cat "${code}_ncbi"/ncbi_dataset/data/"$acc"/*_genomic.fna > "$W/fna/$code.fna"
    rm -rf "$code.zip" "${code}_ncbi"
  fi
  # 2. renombrar scaffolds (funannotate/AUGUSTUS fallan con nombres > 16 caracteres); se guarda la equivalencia
  awk -v m=renamed_map.tsv '/^>/{n="scaffold_" ++i; print n "\t" substr($1,2) > m; print ">" n; next}{print}' \
      "$W/fna/$code.fna" > renamed.fasta

  local dk
  if [ "$(runtime)" = docker ]; then
    dk="docker run --rm -u $(id -u):$(id -g) -v $out:/data -v $FUNANNOTATE_DB_HOST:/db -w /data -e FUNANNOTATE_DB=/db $FA_IMAGE"
  else   # Apptainer: corre como el usuario y la imagen es de solo lectura. --optimize_augustus escribe los parámetros
         # entrenados en $AUGUSTUS_CONFIG_PATH (/usr/share/augustus/config), así que cada genoma usa su propia copia
         # escribible de esa carpeta montada en la misma ruta (las 7 corridas en paralelo no se pisan).
    [ -d "$out/augustus_config" ] || "$APPTAINER" exec "$FA_SIF" cp -r /usr/share/augustus/config "$out/augustus_config"
    # Apptainer monta el /tmp del host (tmpfs con ~1 M de inodos): el alineamiento de proteínas (p2g_*) de 7 genomas en
    # paralelo lo agota. Cada genoma usa su propio /tmp en disco.
    mkdir -p "$out/tmp"
    # optimize_augustus.pl usa #!/usr/bin/perl, que no tiene Parallel::ForkManager -> corre secuencial (1 núcleo; también
    # pasó con LPSC). Se monta una copia idéntica salvo el intérprete (/venv/bin/perl, que sí lo tiene): mismas
    # particiones k-fold y mismo resultado, pero en paralelo (--cpus).
    local fix=$A7/augustus_scripts_fix/optimize_augustus.pl
    if [ ! -s "$fix" ]; then
      mkdir -p "$(dirname "$fix")"
      sed '1s|^#!.*perl.*$|#!/venv/bin/perl|' "$FA_SIF/usr/share/augustus/scripts/optimize_augustus.pl" > "$fix"; chmod +x "$fix"
    fi
    dk="$APPTAINER exec --cleanenv --writable-tmpfs --bind $out:/data --bind $out/tmp:/tmp --bind $FUNANNOTATE_DB_HOST:/db --bind $out/augustus_config:/usr/share/augustus/config --bind $A7/augustus_scripts_fix/optimize_augustus.pl:/usr/share/augustus/scripts/optimize_augustus.pl --pwd /data --env FUNANNOTATE_DB=/db --env TMPDIR=/tmp $FA_SIF"
  fi
  # 3. enmascarado de repeticiones (tantan, igual que LPSC)
  [ -s masked.fasta ] || $dk funannotate mask -i renamed.fasta -o masked.fasta -m tantan --cpus "$FA_CPUS" > mask.log 2>&1
  # 4. predicción génica
  $dk funannotate predict -i masked.fasta -o predict --species "${sp//_/ }" --strain "$strain" \
      --busco_seed_species neurospora_crassa --busco_db sordariomycetes_odb10 \
      --cpus "$FA_CPUS" --optimize_augustus > predict.log 2>&1
  echo "$code: listo ($(grep -c '>' predict/predict_results/*.proteins.fa) proteínas)"
}

if [ $# -ge 1 ]; then
  uno "$1"
else
  if [ "$(runtime)" = docker ]; then docker pull "$FA_IMAGE"
  elif [ ! -d "$FA_SIF" ]; then
    "$APPTAINER" pull "$A7/funannotate_v1.8.15.sif" "docker://$FA_IMAGE"
    "$APPTAINER" build --sandbox "$FA_SIF" "$A7/funannotate_v1.8.15.sif"
  fi
  echo "Contenedor: $(runtime)"
  echo "$GENOMAS_7" | awk '{print $1}' | xargs -P "$FA_JOBS" -I{} bash "$(readlink -f "$0")" {}
fi
