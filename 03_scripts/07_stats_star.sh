#!/usr/bin/env bash

# 07_star_pe_stats.sh
#
# Author: Alondra Márquez
# Date: 2026-05-30
#
# Description:
#   Extrae estadísticas generales de los alineamientos paired-end con STAR.
#   Usa los archivos Log.final.out y genera una tabla por muestra
#   y una tabla resumen para el reporte.
#
# Use:
#   bash 03_scripts/07_stats_star.sh \
#     04_results/06_star_pe \
#     04_results/07_star_pe_stats

set -euo pipefail

# Argumentos de entrada
star_dir="$1"
out_dir="$2"

# Crear carpeta de salida si no existe
mkdir -p "$out_dir"

# Archivos de salida
per_sample="$out_dir/star_pe_stats.tsv"
summary="$out_dir/star_pe_summary.tsv"

# Header de la tabla por muestra
echo -e "sample\tinput_reads\tunique_reads\tunique_percent\tmulti_reads\tmulti_percent" > "$per_sample"

# Recorrer los archivos de estadísticas finales de STAR
for log in "$star_dir"/SRR*/*Log.final.out
do
    # Extraer el ID de la corrida
    sample=$(basename "$log" "_Log.final.out")

    # Extraer estadísticas principales del archivo Log.final.out
    input_reads=$(awk -F '|' '/Number of input reads/ {gsub(/[ \t]/,"",$2); print $2}' "$log")
    unique_reads=$(awk -F '|' '/Uniquely mapped reads number/ {gsub(/[ \t]/,"",$2); print $2}' "$log")
    unique_pct=$(awk -F '|' '/Uniquely mapped reads %/ {gsub(/[ \t%]/,"",$2); print $2}' "$log")
    multi_loci=$(awk -F '|' '/Number of reads mapped to multiple loci/ {gsub(/[ \t]/,"",$2); print $2}' "$log")
    multi_pct=$(awk -F '|' '/% of reads mapped to multiple loci/ {gsub(/[ \t%]/,"",$2); print $2}' "$log")

    # Agregar los resultados de la muestra a la tabla
    echo -e "$sample\t$input_reads\t$unique_reads\t$unique_pct\t$multi_loci\t$multi_pct" >> "$per_sample"
done

# Crear tabla resumen
echo -e "Metrica\tResultado" > "$summary"

# Calcular número de muestras y archivos FASTQ usados
n_samples=$(awk -F '\t' 'NR > 1 {n++} END {print n}' "$per_sample")
n_fastq=$((n_samples * 2))

# Calcular promedios generales de alineamiento
avg_unique=$(awk -F '\t' 'NR > 1 {sum += $4; n++} END {printf "%.2f", sum/n}' "$per_sample")
avg_multi=$(awk -F '\t' 'NR > 1 {sum += $6; n++} END {printf "%.2f", sum/n}' "$per_sample")

# Escribir tabla resumen para el reporte
echo -e "Número de muestras alineadas\t$n_samples" >> "$summary"
echo -e "Tipo de datos\tPaired-end" >> "$summary"
echo -e "Archivos FASTQ limpios usados\t$n_fastq" >> "$summary"
echo -e "Porcentaje promedio de alineamiento único\t${avg_unique}%" >> "$summary"
echo -e "Porcentaje promedio de lecturas multimapeadas\t${avg_multi}%" >> "$summary"

echo "Tablas generadas:" >&2
echo "$per_sample" >&2
echo "$summary" >&2
