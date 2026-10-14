#!/usr/bin/env bash
# A7.6 — Robustez de CAFE5 y del enriquecimiento.
#   (a) Más réplicas del modelo gamma (k = 3) sobre el conjunto principal (22 genomas), con el error fijo del mejor
#       modelo base: las 3 réplicas de 05_cafe5.sh no convergieron bien (λ = 0,445–0,475).
#   (b) Anotación uniforme: CAFE5 + enriquecimiento solo con los 8 genomas anotados con Funannotate (mismo protocolo):
#       LPSC 1155 + los 7 reanotados en A7.1. Quita el sesgo de método de anotación (los proteomas de Funannotate tienen
#       78–88 % de BUSCO completos frente a 91–99 % de las anotaciones de NCBI/JGI).
# Uso: bash 08_robustez.sh     (requiere haber corrido run_A7.sh)
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
C=$A7/cafe; R=$A7/robustez; mkdir -p "$R" "$D/resultados/robustez"
CAFE="conda run -n $ENV cafe5"
ERR=$(grep -H "Final Likelihood" "$C"/base_e_[0-9]/Base_results.txt | sort -t: -k3,3g | head -1 | cut -d: -f1 | xargs dirname)/Base_error_model.txt
echo "error fijo: $ERR"

# ---- (a) réplicas gamma 4–8 en paralelo (14 núcleos cada una) ----
for r in 4 5 6 7 8; do
  [ -s "$C/gamma_e_$r/Gamma_results.txt" ] || $CAFE -i "$C/cafe_input.tsv" -t "$C/ultra.nwk" -e"$ERR" -k 3 -c 14 \
      -o "$C/gamma_e_$r" > "$C/gamma_e_$r.log" 2>&1 &
done
wait
grep -H -E "Final Likelihood|Lambda|Alpha" "$C"/gamma_e_*/Gamma_results.txt | tee "$D/resultados/robustez/cafe_gamma_replicas.txt"

# ---- (b) solo genomas de Funannotate ----
FUN="HumcLPSC HumcCH11 Allacyti Allalyco Xenopyri Leotlubr Cyclminu Chaeraph"
Rscript -e "suppressMessages(library(ape)); t <- read.tree('$C/ultra.nwk');
  t <- keep.tip(t, strsplit('$FUN', ' ')[[1]]); stopifnot(is.ultrametric(t, tol = 1e-5)); write.tree(t, '$R/ultra_funannotate.nwk')"
python3 - "$C/cafe_input.tsv" "$R/cafe_input_funannotate.tsv" "$FUN" <<'EOF'
import sys, csv
inp, out, keep = sys.argv[1], sys.argv[2], sys.argv[3].split()
rows = list(csv.reader(open(inp), delimiter="\t")); h = rows[0]
idx = [h.index(k) for k in keep]
n = 0
with open(out, "w") as f:
    f.write("Desc\tFamily ID\t" + "\t".join(keep) + "\n")
    for r in rows[1:]:
        c = [int(r[i]) for i in idx]
        if sum(x > 0 for x in c) < 2 or max(c) >= 100: continue   # mismos filtros que 05a
        f.write(f"{r[0]}\t{r[1]}\t" + "\t".join(map(str, c)) + "\n"); n += 1
print("familias (solo Funannotate):", n)
EOF
for r in 1 2 3; do
  [ -s "$R/base_e_$r/Base_results.txt" ] || $CAFE -i "$R/cafe_input_funannotate.tsv" -t "$R/ultra_funannotate.nwk" -e -c 20 \
      -o "$R/base_e_$r" > "$R/base_e_$r.log" 2>&1 &
done
wait
ERRF=$(grep -H "Final Likelihood" "$R"/base_e_[0-9]/Base_results.txt | sort -t: -k3,3g | head -1 | cut -d: -f1 | xargs dirname)/Base_error_model.txt
for r in 1 2 3; do
  [ -s "$R/gamma_e_$r/Gamma_results.txt" ] || $CAFE -i "$R/cafe_input_funannotate.tsv" -t "$R/ultra_funannotate.nwk" \
      -e"$ERRF" -k 3 -c 20 -o "$R/gamma_e_$r" > "$R/gamma_e_$r.log" 2>&1 || true &
done
wait
grep -H -E "Final Likelihood|Lambda|Epsilon|Alpha" "$R"/base_e_*/Base_results.txt "$R"/gamma_e_*/Gamma_results.txt 2>/dev/null \
  | tee "$D/resultados/robustez/cafe_funannotate_modelos.txt"
# mejor modelo por AIC (base: λ+ε; gamma: λ+α con ε fijo), igual que 05_cafe5.sh
BEST=$(python3 - "$R" <<'EOF'
import sys, re, glob
R = sys.argv[1]; best = None
for m, pref in (("base", "Base"), ("gamma", "Gamma")):
    for d in glob.glob(f"{R}/{m}_e_[0-9]"):
        f = f"{d}/{pref}_results.txt"
        try: t = open(f).read()
        except OSError: continue
        nll = float(re.search(r"\(-lnL\):\s*([0-9.eE+-]+)", t).group(1)); aic = 4 + 2 * nll
        if best is None or aic < best[0]: best = (aic, d)
print(best[1])
EOF
)
echo "modelo elegido (solo Funannotate): $BEST" | tee -a "$D/resultados/robustez/cafe_funannotate_modelos.txt"
ln -sfn "$BEST" "$R/final"
python3 "$D/07_enriquecimiento.py" "$R/final" "$A7/prot/OrthoFinder/Results_A7/Orthogroups/Orthogroups.tsv" \
    "$A7/eggnog" "$FUNANNOTATE_DB_HOST/go.obo" "$D/resultados/robustez/enriquecimiento_funannotate"
cp "$R/ultra_funannotate.nwk" "$R"/final/*_family_results.txt "$R"/final/*_change.tab "$D/resultados/robustez/" 2>/dev/null || true
echo "Robustez lista: $D/resultados/robustez/"
