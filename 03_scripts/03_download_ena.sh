#!/usr/bin/env bash

# 03_download_ena.sh
#
# uso:
#
# nohup bash 03_scripts/03_download_ena.sh \
#   metadata/runs_to_download.txt \
#   01_data/raw_fastq \
#   3 \
#   > logs/download/download_ena.log 2>&1 &

set -u

# paso de argumentos en linea de comandos
txt_file="$1"
outdir="$2"
max_jobs="${3:-2}"

# creacion de la carpeta de salida
mkdir -p "$outdir"

if [ ! -f "$txt_file" ]; then
    echo "Error: no se encontró el archivo $txt_file"
    exit 1
fi

download_run() {
    run="$1"
    outdir="$2"

    run=$(echo "$run" | tr -d '\r')

    # se revisa que la variable no esté vacía
    [ -z "$run" ] && return 0

    echo ""
    echo "========================================"
    echo "Procesando $run"
    echo "========================================"

    # Creacion del directorio para cada ID
    run_dir="$outdir/$run"
    mkdir -p "$run_dir"

    # Consulta al API de ENA para obtener los archivos FASTQ asociados al run
    urls=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,fastq_ftp&format=tsv" \
        | awk -F '\t' 'NR==2 {print $2}' \
        | tr -d '\r')

    # Si no se encuentra nada en ENA para el ID, se imprime mensaje y se pasa al siguiente
    [ -z "$urls" ] && echo "No se encontraron FASTQ para $run" && return 0

    # Se separan las URLs en un array
    IFS=';' read -r -a fastq_array <<< "$urls"

    echo "Se encontraron ${#fastq_array[@]} archivo(s) para $run"

    for url in "${fastq_array[@]}"
    do
        filename=$(basename "$url")
        outfile="$run_dir/$filename"

        # Si el archivo ya existe y está bien comprimido, no se descarga otra vez
        if [ -s "$outfile" ] && gzip -t "$outfile" 2>/dev/null; then
            echo "Archivo ya existente y válido: $outfile"
            continue
        fi

        # Descarga del archivo
        # -c: reanuda si una descarga quedó incompleta
        # -O: guarda con el nombre esperado dentro de la carpeta de la accesión
        echo "Descargando https://$url"

        wget -c -O "$outfile" "https://$url"

        # Verificación básica de integridad del gzip descargado
        if gzip -t "$outfile"; then
            echo "Archivo verificado: $outfile"
        else
            echo "[ERROR] Archivo corrupto o incompleto: $outfile"
        fi
    done
}

# ciclo para lectura y procesamiento de cada ID del archivo txt
while IFS= read -r run
do
    run=$(echo "$run" | tr -d '\r')
    [ -z "$run" ] && continue

    # Ejecuta la descarga de cada run en segundo plano
    download_run "$run" "$outdir" &

    # Limita el número de descargas simultáneas
    while [ "$(jobs -rp | wc -l)" -ge "$max_jobs" ]
    do
        wait -n
    done

done < "$txt_file"

# Espera a que terminen las descargas restantes
wait

echo ""
echo "Descarga finalizada."
