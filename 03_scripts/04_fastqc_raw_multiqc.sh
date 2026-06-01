#!/usr/bin/env bash

# 04_fastqc_raw_multiqc.sh
#
# Author: Alondra Márquez
# Date: 2026-05-29
#
# Description:
#   Ejecuta FastQC sobre los FASTQ crudos descargados desde ENA
#   y genera un reporte MultiQC inicial.
#
# Use:
#   nohup bash 03_scripts/04_fastqc_raw_multiqc.sh \
#     01_metadata/runs_to_download.txt \
#     02_data/raw_fastq \
#     04_results \
#     6 \
#     > 06_logs/qc/fastqc_raw_multiqc.log 2>&1 &

set -euo pipefail

runs_file="$1"
raw_dir="$2"
results_dir="$3"
threads="${4:-6}"

fastqc_raw_dir="${results_dir}/01_fastqc_raw"
multiqc_raw_dir="${results_dir}/02_multiqc_raw"

mkdir -p "$fastqc_raw_dir"
mkdir -p "$multiqc_raw_dir"
mkdir -p 06_logs/qc

if [ ! -f "$runs_file" ]; then
    echo "Error no se encontro el archivo de corridas: $runs_file"
    exit 1
fi

if [ ! -d "$raw_dir" ]; then
    echo "Error no se encontro el directorio de FASTQ crudos: $raw_dir"
    exit 1
fi

echo "FastQC de datos crudos"
echo "Runs file: $runs_file"


while IFS= read -r run
do
    run=$(echo "$run" | tr -d '\r')
    [ -z "$run" ] && continue

    echo "Procesando $run"

    r1="${raw_dir}/${run}/${run}_1.fastq.gz"
    r2="${raw_dir}/${run}/${run}_2.fastq.gz"

    if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
        echo "[ERROR] No se encontraron ambos FASTQ paired-end para $run"
        continue
    fi

    fastqc \
        -t "$threads" \
        -o "$fastqc_raw_dir" \
        "$r1" "$r2"

    echo "$run terminado"
    echo ""

done < "$runs_file"


echo "Reporte MultiQC"

multiqc \
    "$fastqc_raw_dir" \
    -o "$multiqc_raw_dir" \
    -n "multiqc_raw_GSE115354.html"

echo "FastQC crudo + MultiQC finalizado."
