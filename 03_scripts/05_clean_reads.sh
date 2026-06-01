#!/usr/bin/env bash

# 05_fastp_clean_multiqc.sh
#
# Author: Alondra Márquez
# Date: 2026-05-30
#
# Description:
#   Realiza una limpieza conservadora de lecturas paired-end con fastp
#   y evalúa los FASTQ limpios con FastQC y MultiQC.
#
#   Esta limpieza no aplica recortes fijos de extremos ni eliminación
#   automática de duplicados. Se enfoca en detección de adaptadores,
#   filtrado moderado por calidad, límite de bases N y eliminación de
#   lecturas demasiado cortas después del procesamiento.
#
# Use:
#   nohup bash 03_scripts/05_clean_reads.sh \
#     01_metadata/runs_to_download.txt \
#     02_data/raw_fastq \
#     02_data/clean_fastq \
#     04_results \
#     6 \
#     > 06_logs/qc/fastp_clean_multiqc.log 2>&1 &

set -euo pipefail

runs_file="$1"
raw_dir="$2"
clean_dir="$3"
results_dir="$4"
threads="${5:-6}"

fastp_report_dir="${results_dir}/03_fastp_reports"
fastqc_clean_dir="${results_dir}/04_fastqc_clean"
multiqc_clean_dir="${results_dir}/05_multiqc_clean"

mkdir -p "$clean_dir"
mkdir -p "$fastp_report_dir"
mkdir -p "$fastqc_clean_dir"
mkdir -p "$multiqc_clean_dir"
mkdir -p 06_logs/qc

if [ ! -f "$runs_file" ]; then
    echo "Error: no se encontro el archivo de corridas: $runs_file"
    exit 1
fi

if [ ! -d "$raw_dir" ]; then
    echo "Error: no se encontro el directorio de FASTQ crudos: $raw_dir"
    exit 1
fi


echo "Limpieza con fastp"

while IFS= read -r run
do
    run=$(echo "$run" | tr -d '\r')
    [ -z "$run" ] && continue

    echo "Procesando fastp para $run"

    r1="${raw_dir}/${run}/${run}_1.fastq.gz"
    r2="${raw_dir}/${run}/${run}_2.fastq.gz"

    sample_clean_dir="${clean_dir}/${run}"
    mkdir -p "$sample_clean_dir"

    clean_r1="${sample_clean_dir}/${run}_1.clean.fastq.gz"
    clean_r2="${sample_clean_dir}/${run}_2.clean.fastq.gz"

    html_report="${fastp_report_dir}/${run}.fastp.html"
    json_report="${fastp_report_dir}/${run}.fastp.json"


    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "$clean_r1" \
        -O "$clean_r2" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --n_base_limit 5 \
        --length_required 35 \
        --thread "$threads" \
        --html "$html_report" \
        --json "$json_report"

    echo "$run fastp terminado"
    echo ""

done < "$runs_file"

echo "FastQC de datos limpios"

while IFS= read -r run
do
    run=$(echo "$run" | tr -d '\r')
    [ -z "$run" ] && continue

    echo "Procesando FastQC limpio para $run"

    clean_r1="${clean_dir}/${run}/${run}_1.clean.fastq.gz"
    clean_r2="${clean_dir}/${run}/${run}_2.clean.fastq.gz"


    fastqc \
        -t "$threads" \
        -o "$fastqc_clean_dir" \
        "$clean_r1" "$clean_r2"

    echo "$run FastQC limpio terminado"
    echo ""

done < "$runs_file"


echo "Reporte MultiQC de datos limpios"

multiqc \
    "$fastqc_clean_dir" \
    "$fastp_report_dir" \
    -o "$multiqc_clean_dir" \
    -n "multiqc_clean_GSE115354.html"

echo "Proceso finalizado."