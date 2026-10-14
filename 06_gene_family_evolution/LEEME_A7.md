# A7 — Expansión de familias génicas (CAFE5) y enriquecimiento funcional

Objetivo y justificación: `Paper/00_plan/PLAN_ANALISIS.md` (sección A7). Decisión del 2026-10-02: se reanotan **solo los 7 genomas sin anotación** (los más cercanos a *Humicolopsis*: CH-1.1, *Allantophomopsis* ×2, *Xenosphaeropsis*, *Leotia*, *Cyclaneusma*, *Chaetomella*). Los demás conservan su anotación original.

## Al llegar al PC de 256 GB (70 núcleos)

1. **Datos:** seguir `00_plan/traspaso_nuevo_PC/LEEME_TRASPASO.md`, sección 0. Se necesitan `~/hongo_comparativa/A3/prot/` (proteomas de las 17 especies que se mantienen) y el proyecto en `~/proyectos/Proyecto-Hongo`. Si faltan los genomas de los 7 (`~/hongo_comparativa/fna/`), el paso 1 los descarga de NCBI.
2. **Revisar las rutas** en `00_config.sh` (`PROJ`, `W`, `EGGNOG_DB`). `nproc` debe dar 70.
3. **Entorno:** `mamba env create -f entorno_a7.yml`. Comprobar `conda run -n a7 cafe5 --help`. Si `cafe` no instala, compilar CAFE5 desde github.com/hahnlab/CAFExp con cmake y dejar `cafe5` en el PATH del entorno.
4. **Docker:** `docker run --rm hello-world` sin sudo (el usuario debe estar en el grupo `docker`).
5. **Lanzar:**
   ```bash
   conda activate a7     # Rscript (ape), datasets y unzip se usan desde el PATH
   cd ~/proyectos/Proyecto-Hongo/Paper/analisis/A7_expansion_enriquecimiento
   nohup bash run_A7.sh > ~/hongo_comparativa/A7_run.log 2>&1 &
   ```
   Cada paso se salta si su salida ya existe, así que se puede relanzar tras un corte.

## Pasos y tiempos estimados (70 núcleos)

| Script | Qué hace | Tiempo |
|---|---|---|
| `01_reanotar_funannotate.sh` | Protocolo de LPSC (`protocolo_anotacion_humicolopsis_M4.txt`): renombrar scaffolds → `funannotate mask -m tantan` → `funannotate predict` (BUSCO sordariomycetes_odb10, semilla *N. crassa*, `--optimize_augustus`). Los 7 en paralelo, 10 núcleos cada uno | ~6–8 h (*Leotia*, 86 Mb, es el más lento) |
| `02_proteomas.sh` | Arma `A7/prot/` (7 nuevos + 17 de A3) con encabezados `código\|id`, y la completitud BUSCO de cada proteoma | ~1 h |
| `03_orthofinder.sh` | OrthoFinder con los mismos parámetros que A3 | ~2–3 h |
| `06_eggnog.sh` | eggNOG-mapper en todos los proteomas (`--dbmem`); descarga la base (~50 GB) | ~3–5 h (en paralelo con OrthoFinder) |
| `05_cafe5.sh` | Árbol ultramétrico (`04_arbol_ultrametrico.R`), filtros (`05a_preparar_cafe.py`), CAFE5 base con error estimado (3 réplicas; modelo principal) y gamma k=3 con ese error fijo (opcional), elección por AIC, familias grandes con λ fijo | ~1–2 h |
| `07_enriquecimiento.py` | Familias expandidas o contraídas en el tronco de *Humicolopsis*, en la rama terminal de LPSC y en el tronco de Phacidiales; Fisher + BH sobre GO, rutas KEGG, KO y Pfam | minutos |

**Total: ~12–16 h**, casi todo sin supervisión.

## Controles al terminar
- `resultados/n_proteinas.tsv`: los 7 reanotados deberían quedar en ~9.000–13.000 proteínas (antes ~7.200–8.200 con miniprot).
- `resultados/busco_proteomas.tsv`: completitud de los proteomas nuevos comparable a la de LPSC.
- `resultados/cafe_modelos.txt`: λ parecido entre réplicas; si no, correr más réplicas.
- `resultados/enriquecimiento/resumen_ramas.tsv` y `enriquecimiento.tsv`: buscar términos de P (transporte de fosfato, fosfatasas, fitasas, polifosfato, oxalato y ácidos orgánicos).
- Después, rehacer los conteos de A3/A4 con `A7/prot/` y actualizar la Fig. 2 si cambian.

## Prueba en el Mac (2026-10-02)
Los scripts 04, 05a, CAFE5 y 07 se probaron con los proteomas actuales (miniprot) para validar el código. **Esos resultados no son válidos para el paper.** Lo que se aprendió:
- **Árbol ultramétrico:** `chronos` con modelo `relaxed` y λ = 0,1 converge y no deja ramas de longitud ~0. El modelo `correlated` dejaba la rama basal de Leotiomycetes en 0.
- **CAFE5 base con `-e`:** funcionó (~1 min con 10 núcleos; λ = 0,57, ε = 0,013).
- **Gamma:** CAFE5 no permite `-e` (estimar el error) junto con `-k`. Con el error fijo tampoco inicializó ("Failed to initialize any reasonable values") por familias con gran diferencia max-min. Por eso es opcional.
- **Familias de ≥ 100 copias** (10 en la prueba): tampoco inicializan con λ fijo, así que se reportan como excluidas.
- **`07_enriquecimiento.py`** funcionó con la anotación eggNOG de LPSC. Como adelanto, solo para orientar: en el tronco de *Humicolopsis* aparecieron enriquecidos los transportadores de membrana y una fosfatasa ácida (KEGG K01078). Hay que confirmarlo con los proteomas reanotados.
