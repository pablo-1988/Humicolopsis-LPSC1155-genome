#!/usr/bin/env bash
# A5 — Secretoma (herramientas sin licencia; SignalP 6 / DeepTMHMM / NetGPI requieren licencia DTU).
#   1. DeepSig (péptidos señal, modelo eucarionte) sobre los 24 proteomas de A7 (7 reanotados con Funannotate + 17 originales)
#   2. TMbed (ProtT5; hélices/barriles transmembrana y péptido señal) SOLO sobre las proteínas con péptido señal de DeepSig
#      (en CPU tarda ~2–3 s por proteína: correrlo sobre los ~290.000 proteínas no es viable; la GPU Quadro P400 de 2 GB no alcanza)
# Entornos: `deepsig` (python 3.8, tensorflow 2.2, keras 2.4.3, numpy 1.18, deepsig-biocomp 0.9 + modelos de
#            github.com/BolognaBiocomp/deepsig en DEEPSIG_ROOT) y `tmbed` (python 3.10, torch CPU, TMbed 1.0.2 desde GitHub,
#            transformers 4.49).
# Uso: nohup bash 01_predicciones.sh > ~/hongo_comparativa/A5_run.log 2>&1 &
set -euo pipefail
W=$HOME/hongo_comparativa; A=$W/A5; PROT=$W/A7/prot
export DEEPSIG_ROOT=$A/deepsig_src
DS=$HOME/miniforge3/envs/deepsig/bin/deepsig; TM=$HOME/miniforge3/envs/tmbed/bin/tmbed
mkdir -p "$A/deepsig" "$A/tmbed"

echo "[$(date +%T)] 1/2 DeepSig en 24 proteomas (paralelo)"
ls "$PROT"/*.faa | xargs -P 24 -I{} sh -c 'c=$(basename {} .faa); o='"$A"'/deepsig/$c.gff3;
  [ -s "$o" ] || OMP_NUM_THREADS=3 TF_NUM_INTRAOP_THREADS=3 '"$DS"' -f {} -o "$o" -k euk -m gff3 > '"$A"'/deepsig/$c.log 2>&1 || echo "FALLÓ DeepSig $c"'
cut -f1 "$A"/deepsig/*.gff3 | sort -u | wc -l | xargs echo "proteínas evaluadas:"
grep -h "Signal peptide" "$A"/deepsig/*.gff3 | cut -f1 | sort -u > "$A/sp_ids.txt"
echo "con péptido señal: $(wc -l < "$A/sp_ids.txt")"

echo "[$(date +%T)] 2/2 TMbed sobre proteínas con péptido señal (4 procesos x 16 hilos)"
[ -s "$A/tmbed/chunk00.faa" ] || python3 - "$A/sp_ids.txt" "$PROT" "$A/tmbed" <<'EOF'
import sys, glob
ids = set(open(sys.argv[1]).read().split()); out = sys.argv[3]
seqs = []
for f in sorted(glob.glob(sys.argv[2] + "/*.faa")):
    n = None
    for l in open(f):
        if l.startswith(">"): n = l[1:].split()[0]; keep = n in ids; s = []
        if keep:
            if l.startswith(">"): seqs.append([n, ""])
            else: seqs[-1][1] += l.strip()
seqs.sort(key=lambda x: len(x[1]))                       # lotes por largo = menos relleno
# 40 bloques: TMbed solo escribe al terminar cada archivo, así un corte pierde a lo más un bloque por proceso
for k in range(40):
    with open(f"{out}/chunk{k:02d}.faa", "w") as fh:
        for n, s in seqs[k::40]: fh.write(f">{n}\n{s[:3000]}\n")   # TMbed/ProtT5: proteínas muy largas truncadas a 3.000 aa
print("secuencias a TMbed:", len(seqs))
EOF
# 4 procesos x 16 hilos; cada bloque se salta si ya tiene predicción (se puede relanzar tras un corte)
ls "$A"/tmbed/chunk[0-9][0-9].faa | xargs -P 4 -I{} sh -c 'p=${1%.faa}.pred; [ -s "$p" ] || OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 '"$TM"' predict -f "$1" -p "$p.tmp" --no-use-gpu --out-format 0 > "${1%.faa}.log" 2>&1 && mv "$p.tmp" "$p"' _ {}
echo "[$(date +%T)] listo: $(cat "$A"/tmbed/chunk*.pred | grep -c '>') proteínas con predicción TMbed"
