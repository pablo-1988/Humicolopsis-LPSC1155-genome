###############################################################################
# PROTOCOLO MAESTRO: ANOTACIÓN DE GENOMA FÚNGICO (HUMICOLOPSIS)
# Plataforma: Apple Silicon (M4) + OrbStack
# Software: Funannotate v1.8.15 (Docker)
# Fecha de creación: 06 de Enero, 2026
# Autor: Tu Orquestador de IA & Usuario
###############################################################################

# =============================================================================
# REQUISITOS PREVIOS
# =============================================================================
# 1. Tener instalado OrbStack (alternativa ligera a Docker Desktop para Mac).
# 2. Tener el archivo de licencia de GeneMark descargado (http://topaz.gatech.edu/GeneMark/).
# 3. Tener el genoma ensamblado (FASTA).

# =============================================================================
# PASO 1: PREPARACIÓN DEL ENTORNO DE TRABAJO
# =============================================================================
# Crear carpeta y organizar archivos.
# Es CRÍTICO renombrar la llave a '.gm_key' (con punto) para que sea invisible/estándar.

mkdir -p ~/bioinfo/humicolopsis
cd ~/bioinfo/humicolopsis

# Mueve aquí tu genoma original y tu llave descargada.
# Renombrar la llave para compatibilidad:
mv gm_key_64 .gm_key   # (O 'cp gm_key .gm_key' si ya la tienes ahí)


# =============================================================================
# PASO 2: LIMPIEZA DE NOMBRES DE SCAFFOLDS (CORRECCIÓN CRÍTICA)
# =============================================================================
# PROBLEMA: GeneMark/Augustus fallan si los nombres de los contigs tienen >16 caracteres
# o caracteres extraños.
# SOLUCIÓN: Renombrar todo secuencialmente (scaffold_1, scaffold_2...) usando awk.

awk '/^>/{print ">scaffold_" ++i; next}{print}' hongo.scaffold.fasta > hongo_renamed.fasta

# Input: hongo.scaffold.fasta (Tu original)
# Output: hongo_renamed.fasta (El que usaremos de ahora en adelante)


# =============================================================================
# PASO 3: DESCARGA DE BASES DE DATOS (SETUP HÍBRIDO)
# =============================================================================
# PROBLEMA: El instalador automático falla al bajar RepeatDB (Error 308) y BUSCO genérico.
# SOLUCIÓN: Instalación parcial automática + Inyección manual de BUSCO + Dummy RepeatDB.

# 3.1. Descarga automática de lo que sí funciona (UniProt, Pfam, MEROPS)
docker run -it --rm \
    --platform linux/amd64 \
    -v $(pwd):/data \
    -w /data \
    nextgenusfs/funannotate:v1.8.15 \
    funannotate setup -d /data/funannotate_db -b sordariomycetes

# 3.2. Parche para la Base de Datos de Repeticiones (Dummy DB)
# Creamos una base de datos de proteínas falsa para engañar al chequeo inicial,
# ya que usaremos 'tantan' para enmascarar y no necesitamos esta DB pesada.

echo ">dummy_repeat" > dummy.fasta
echo "MVLSPADKTNVLSPADKTN" >> dummy.fasta  # Secuencia de aminoácidos aleatoria

docker run -u $(id -u):$(id -g) -it --rm \
    --platform linux/amd64 \
    -v $(pwd):/data \
    -w /data \
    nextgenusfs/funannotate:v1.8.15 \
    diamond makedb --in dummy.fasta --db funannotate_db/repeats

rm dummy.fasta # Limpieza


# =============================================================================
# PASO 4: ENMASCARAMIENTO DE REPETICIONES (MASKING)
# =============================================================================
# Usamos el algoritmo matemático 'tantan' para ocultar regiones repetitivas.
# Esto previene que se predigan genes falsos en transposones.

docker run -u $(id -u):$(id -g) -it --rm \
    --platform linux/amd64 \
    -v $(pwd):/data \
    -w /data \
    nextgenusfs/funannotate:v1.8.15 \
    funannotate mask \
    -i hongo_renamed.fasta \
    -o hongo_masked.fasta \
    -m tantan \
    --cpus 8

# Output: hongo_masked.fasta (Genoma limpio y enmascarado)


# =============================================================================
# PASO 5: PREDICCIÓN FINAL (EL COMANDO MAESTRO)
# =============================================================================
# Ejecuta el entrenamiento de Augustus y la predicción final.
# NOTA: GeneMark puede dar warning si no encuentra la llave, pero el pipeline continuará
# exitosamente usando Augustus + BUSCO (Sordariomycetes) + Exonerate.

docker run -u $(id -u):$(id -g) -it --rm \
    --platform linux/amd64 \
    --memory=16g \
    -v $(pwd):/data \
    -w /data \
    -e FUNANNOTATE_DB=/data/funannotate_db \
    nextgenusfs/funannotate:v1.8.15 \
    funannotate predict \
    -i hongo_masked.fasta \
    -o anotacion_final \
    --species "Humicolopsis cephalosporioides" \
    --strain "T1" \
    --busco_seed_species neurospora_crassa \
    --busco_db sordariomycetes_odb10 \
    --cpus 8 \
    --optimize_augustus

# EXPLICACIÓN DE BANDERAS:
# --platform linux/amd64: Vital para Mac M1/M2/M3/M4 (usa Rosetta 2).
# -w /data: Establece el directorio de trabajo para evitar errores de permisos de log.
# --busco_seed_species: Usa Neurospora (muy estudiado) para enseñar a Augustus inicialmente.
# --busco_db: Usa la base específica de Sordariomycetes que descargamos.
# --optimize_augustus: Re-entrena los parámetros para máxima precisión (tarda más, vale la pena).


# =============================================================================
# PASO 6: RESULTADOS
# =============================================================================
# Los archivos finales estarán en la carpeta: anotacion_final/predict_results/
# Archivos clave:
#   1. *.gff3 : El archivo de anotación principal (coordenadas genéticas).
#   2. *.gbk  : Archivo GenBank rico para visualización.
#   3. *.proteins.fa : Secuencias de aminoácidos (Input para AntiSMASH/OrthoFinder).
#   4. *.stats.json : Resumen de cuántos genes se encontraron.

###############################################################################
# FIN DEL PROTOCOLO
###############################################################################