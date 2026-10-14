#!/usr/bin/env bash
# A7.4 — CAFE5: expansión/contracción de familias sobre el árbol ultramétrico de A1.
#   1. árbol ultramétrico (04_arbol_ultrametrico.R)
#   2. entrada filtrada (05a_preparar_cafe.py): principal (< 100 copias) y familias grandes aparte
#   3. modelo base con error estimado (-e), 3 réplicas para ver si λ converge  -> modelo principal
#      modelo gamma (-k 3) con el error del mejor base FIJO (CAFE5 no estima error con gamma). Es opcional:
#      en la prueba del Mac falló ("Failed to initialize any reasonable values") por familias con gran
#      diferencia max-min; si falla, se sigue solo con el base.
#   4. elección por AIC (base: λ + ε = 2 parámetros; gamma: λ + α = 2, con ε fijo); familias grandes con λ fijo
# Uso: bash 05_cafe5.sh [GeneCount.tsv]   (por defecto el de 03_orthofinder.sh)
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/00_config.sh"
GC=${1:-$D/resultados/Orthogroups.GeneCount.tsv}
C=$A7/cafe; mkdir -p "$C"
CAFE="conda run -n $ENV cafe5"

Rscript "$D/04_arbol_ultrametrico.R" "$PROJ/Paper/analisis/A1_filogenomica_ortologos/resultados_particionado/part_tree.treefile" "$C/ultra.nwk"
python3 "$D/05a_preparar_cafe.py" "$GC" "$C" "$NO_CAFE"

for r in 1 2 3; do
  [ -s "$C/base_e_$r/Base_results.txt" ] || $CAFE -i "$C/cafe_input.tsv" -t "$C/ultra.nwk" -e -c "$THREADS" -o "$C/base_e_$r" > "$C/base_e_$r.log" 2>&1
done
ERR=$(grep -H "Final Likelihood" "$C"/base_e_[0-9]/Base_results.txt | sort -t: -k3,3g | head -1 | cut -d: -f1 | xargs dirname)/Base_error_model.txt
for r in 1 2 3; do
  [ -s "$C/gamma_e_$r/Gamma_results.txt" ] || $CAFE -i "$C/cafe_input.tsv" -t "$C/ultra.nwk" -e"$ERR" -k 3 -c "$THREADS" -o "$C/gamma_e_$r" > "$C/gamma_e_$r.log" 2>&1 \
    || echo "gamma réplica $r falló (ver $C/gamma_e_$r.log); se usa el modelo base"
done
grep -H -E "Model|Lambda|Epsilon|Alpha" "$C"/base_e_*/Base_results.txt "$C"/gamma_e_*/Gamma_results.txt 2>/dev/null | tee "$D/resultados/cafe_modelos.txt"

# Elección por AIC, réplica de mejor verosimilitud de cada modelo
python3 - "$C" > "$C/elegido.txt" <<'EOF'
import sys, re, glob
C = sys.argv[1]; best = {}
for m, pref, k in (("base", "Base", 2), ("gamma", "Gamma", 2)):
    for d in glob.glob(f"{C}/{m}_e_[0-9]"):
        if not glob.glob(f"{d}/{pref}_results.txt"): continue
        t = open(f"{d}/{pref}_results.txt").read()
        nll = float(re.search(r"\(-lnL\):\s*([0-9.eE+-]+)", t).group(1))   # "Model Base Final Likelihood (-lnL): 145487"
        lam = re.search(r"^Lambda:\s*([0-9.eE+-]+)", t, re.M).group(1)
        if m not in best or nll < best[m][1]: best[m] = (d, nll, lam, pref, 2 * k + 2 * nll)
d, nll, lam, pref, aic = min(best.values(), key=lambda b: b[4])
print(d, pref, lam)
for m, b in best.items(): print(f"# {m}: -lnL={b[1]:.2f} AIC={b[4]:.2f} lambda={b[2]}", file=sys.stderr)
EOF
read -r BEST PREF LAMBDA < "$C/elegido.txt"
echo "Modelo elegido: $BEST (λ = $LAMBDA)" | tee -a "$D/resultados/cafe_modelos.txt"
ln -sfn "$BEST" "$C/final"

# Familias grandes (>= 100 copias en alguna especie) con λ fijo. En la prueba del Mac (10 familias) CAFE5 no
# pudo inicializarlas; si vuelve a fallar, se reportan como excluidas (lista en resultados/cafe_large.tsv).
cp "$C/cafe_large.tsv" "$D/resultados/"
if [ "$(wc -l < "$C/cafe_large.tsv")" -gt 1 ]; then
  [ -s "$C/large/Base_results.txt" ] || $CAFE -i "$C/cafe_large.tsv" -t "$C/ultra.nwk" -l "$LAMBDA" -c "$THREADS" -o "$C/large" > "$C/large.log" 2>&1
fi
cp "$C/ultra.nwk" "$D/resultados/"
cp "$C"/final/*_family_results.txt "$C"/final/*_change.tab "$C"/final/*_asr.tre "$D/resultados/" 2>/dev/null || true
