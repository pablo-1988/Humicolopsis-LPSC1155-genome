#!/usr/bin/env bash
# A1 paso 3 (PC con >= 128 GB RAM): concatenado PARTICIONADO por gen de los 2,742 loci BUSCO (23 taxa).
# Responde a la limitación de la corrida en Mac (modelo único Q.YEAST+F+I+R5 elegido en 300 loci por falta de RAM).
#
#   (a) ModelFinder sobre la matriz completa, modelo único (verifica la elección hecha en la submuestra)
#   (b) Particionado por gen: ModelFinder por partición + fusión de particiones (MFP+MERGE, rcluster 10)
#   (c) Factores de concordancia gCF/sCF sobre el árbol particionado
#
# Requisitos: IQ-TREE >= 2.2 (probado con 3.1.3). Uso:
#   bash 03_concatenado_particionado.sh <carpeta_A1_phylogenomics> <hilos>
# La carpeta debe contener concat.faa, partitions.txt y genetrees.nwk (vienen en hongo_comparativa_esencial.tar.gz).
set -euo pipefail
D=${1:?carpeta con concat.faa}; T=${2:-32}
IQ=${IQTREE:-iqtree3}
cd "$D"
mkdir -p particionado && cd particionado

# (a) modelo único con ModelFinder sobre la matriz completa (conjunto de modelos de proteína habituales)
[ -f mf_full.iqtree ] || $IQ -s ../concat.faa -m MF -mset LG,WAG,JTT,Q.pfam,Q.yeast,Q.insect -T "$T" \
    --prefix mf_full -redo --quiet
grep "Best-fit model according to BIC" mf_full.iqtree

# (b) particionado por gen con fusión de particiones (relaxed clustering 10%)
[ -f part.contree ] || $IQ -s ../concat.faa -p ../partitions.txt -m MFP+MERGE -rcluster 10 \
    -mset LG,WAG,JTT,Q.pfam,Q.yeast,Q.insect -B 1000 -T "$T" --prefix part -redo --quiet
# (IQ-TREE puede fallar con -o + MFP+MERGE; el árbol se enraíza en Neurcras en la figura)

# (c) concordancia
$IQ -t part.treefile --gcf ../genetrees.nwk -s ../concat.faa --scf 100 -T "$T" --prefix part_concord -redo --quiet
echo "Listo: $(date)"; grep -E "Best-fit|Log-likelihood of the tree" part.iqtree
