#!/usr/bin/env bash
# A1 (b) fase final: búsqueda de árbol con el esquema MFP+MERGE ya elegido (part.best_scheme.nex: 2346 particiones con modelo),
# UFBoot 1000, luego (c) gCF/sCF. Equivale a la continuación de 03_concatenado_particionado.sh paso (b)+(c).
set -uo pipefail
export PATH=/home/ubuntu/miniforge3/envs/comparativa/bin:$PATH
IQ=iqtree3; T=70
cd ~/hongo_comparativa/A1_phylogenomics/particionado
log(){ echo "[$(date '+%F %T')] $*"; }
log "Árbol con esquema final, -T $T"
$IQ -s ../concat.faa -p part.best_scheme.nex -B 1000 -T $T --prefix part_tree --quiet > ../../A1_arbol.log 2>&1 \
  || { log "falló el árbol (ver part_tree.log)"; exit 1; }
log "Árbol listo"; grep -E "Log-likelihood of the tree" part_tree.iqtree
log "(c) gCF/sCF"
$IQ -t part_tree.treefile --gcf ../genetrees.nwk -s ../concat.faa --scf 100 -T $T --prefix part_concord -redo --quiet > ../../A1_c.log 2>&1 \
  && log "(c) listo: part_concord.cf.tree" || log "(c) falló, ver A1_c.log"
