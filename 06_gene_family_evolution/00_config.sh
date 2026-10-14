#!/usr/bin/env bash
# A7 — configuración común. Ajustar estas rutas al llegar al PC de 256 GB (Linux, 70 núcleos).
# Todos los scripts hacen: source "$(dirname "$0")/00_config.sh"

# Carpeta del proyecto (copia de Drive; ver 00_plan/traspaso_nuevo_PC/LEEME_TRASPASO.md)
PROJ=${PROJ:-$HOME/Documentos/Projects/Project-Hongo}   # PC Linux 256 GB (antes: $HOME/proyectos/Proyecto-Hongo)
# Datos de trabajo (fuera de Drive): faa/, fna/, A3/prot/ ...
W=${W:-$HOME/hongo_comparativa}
# Salidas pesadas de A7
A7=$W/A7
# Scripts y resultados livianos de A7 (esta carpeta)
D=$PROJ/Paper/analisis/A7_expansion_enriquecimiento

THREADS=${THREADS:-70}
# Funannotate: 7 genomas en paralelo x 10 núcleos = 70
FA_JOBS=7
FA_CPUS=10

# Base de datos de Funannotate usada para LPSC 1155 (sordariomycetes_odb10, UniProt, Pfam, dbCAN, MEROPS)
FUNANNOTATE_DB_HOST=$PROJ/humicolopsis_annot/funannotate_db
FA_IMAGE=nextgenusfs/funannotate:v1.8.15
# Sin docker (PC Linux sin sudo): Apptainer instalado sin root en el env `apptainer` corre la misma imagen como .sif
APPTAINER=${APPTAINER:-$HOME/miniforge3/envs/apptainer/bin/apptainer}
# Sin squashfuse, Apptainer descomprimiría la .sif (4 GB) en cada llamada: se usa una sandbox (carpeta) construida una vez
FA_SIF=${FA_SIF:-$W/A7/funannotate_v1.8.15_sandbox}

# Base de datos de eggNOG-mapper (~50 GB)
EGGNOG_DB=${EGGNOG_DB:-$HOME/db/eggnog}

# Entorno conda/mamba con: orthofinder diamond mafft fasttree cafe eggnog-mapper compleasm r-ape python>=3.9 biopython scipy pandas
ENV=${ENV:-a7}
# compleasm no resuelve dentro de a7 (conflicto dendropy/sepp): se usa el del env de A1 (misma versión 0.2.9)
ENV_BUSCO=${ENV_BUSCO:-comparativa}

# Genomas a reanotar (proteomas actuales de miniprot): código  accesión  especie  cepa
GENOMAS_7="HumcCH11 GCA_050941815.1 Humicolopsis_cephalosporioides CH-1.1
Allacyti GCA_024752435.1 Allantophomopsis_cytisporea ATCC66955
Allalyco GCA_024752465.1 Allantophomopsis_lycopodina ATCC66958
Xenopyri GCA_022609375.1 Xenosphaeropsis_pyriputrescens CBS115176
Leotlubr GCA_965240495.1 Leotia_lubrica gyLeoLubr2
Cyclminu GCA_037179335.1 Cyclaneusma_minus NZFS725
Chaeraph GCA_040369945.1 Chaetomella_raphigera JCM9995"

# Proteomas que se mantienen (anotación original o de este estudio) y entran desde A3/prot/
MANTENER="HumcLPSC Clausp Pseudest Oidimaiu Hyalhepa Phiasuba Cadosp Sclescle Botrcine Glarlozo Hymefrax Mollscop Chloaeru Erysneca Neurcras Sacccere Aspenige"
# Referencias funcionales: entran a OrthoFinder y eggNOG, NO a CAFE5
NO_CAFE="Sacccere Aspenige"
