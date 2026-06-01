#!/usr/bin/env bash

# 06_star_pe.sh
#
# Author: Alondra Márquez
# Date: 2026-05-30
#
# Description:
#   Alinea lecturas paired-end de RNA-seq usando STAR.
#   El script recorre carpetas SRR dentro de un directorio de FASTQ,
#   identifica los archivos R1 y R2 de cada muestra, y genera un archivo
#   BAM ordenado por coordenadas para cada corrida.
#
# Input:
#   1. Directorio del índice STAR.
#   2. Directorio con las carpetas SRR que contienen los FASTQ paired-end.
#   3. Directorio de salida para los alineamientos.
#   4. Directorio de logs.
#   5. Número de threads para STAR.
#
# Output:
#   Para cada muestra se genera:
#     - BAM ordenado por coordenadas.
#     - Log estándar de STAR.
#     - Log de error de STAR.
#     - Archivo con tiempo y uso de recursos.
#
# Use:
#   nohup bash 03_scripts/06_star_pe.sh \
#     02_data/reference/mm39.gencode.M36.star \
#     02_data/clean_fastq \
#     04_results/06_star_pe \
#     06_logs/star_pe \
#     6 \
#     > 06_logs/star_pe/nohup_star_pe.log 2>&1 &

set -euo pipefail

# Argumentos de entrada
index_dir="$1"
srr_dir="$2"
output_dir="$3"
log_dir="$4"
threads="$5"

# Crear carpetas generales si no existen
mkdir -p "$output_dir"
mkdir -p "$log_dir"


echo "Alineamiento paired-end con STAR"

# Recorrer todas las carpetas que empiezan con SRR
for dir in "$srr_dir"/SRR*
do
    # Extraer nombre de la carpeta, que corresponde al ID de la corrida
    srr=$(basename "$dir")

    # Construir rutas de los archivos paired-end
    r1="$dir/${srr}_1.clean.fastq.gz"
    r2="$dir/${srr}_2.clean.fastq.gz"

    # Crear carpeta de salida para esta muestra
    mkdir -p "$output_dir/$srr"

    echo "Procesando $srr"
    echo "R1: $r1"
    echo "R2: $r2"

    # Parametros
    # /usr/bin/time -v:registra tiempo de ejecución y uso de recursos.
    # --runThreadN:número de hilos usados por STAR.
    # --genomeDir:directorio donde se encuentra el índice STAR.
    # --readFilesIn: archivos FASTQ paired-end de entrada.
    # --readFilesCommand zcat:permite leer archivos FASTQ comprimidos en formato .gz.
    # --outFileNamePrefix: prefijo de los archivos de salida para cada muestra.
    # --outSAMtype BAM SortedByCoordinate:genera directamente un archivo BAM ordenado por coordenadas

    /usr/bin/time -v -o "$log_dir/${srr}.time.txt" \
    STAR \
        --runThreadN "$threads" \
        --genomeDir "$index_dir" \
        --readFilesIn "$r1" "$r2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$output_dir/$srr/${srr}_" \
        --outSAMtype BAM SortedByCoordinate \
        > "$log_dir/${srr}.stdout.log" \
        2> "$log_dir/${srr}.stderr.log"

    echo "$srr terminado"
    echo ""
done
echo "Alineamiento finalizado."