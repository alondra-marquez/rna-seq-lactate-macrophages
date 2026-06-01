
"""
02_download.py

Author: Alondra Márquez
Date: 2026-05-29

Description:
    Descarga corridas SRR con prefetch, las convierte a FASTQ con fasterq-dump
    y comprime los archivos FASTQ usando pigz.

Use:
    python3 03_scripts/02_download.py \
      --txt metadata/runs_to_download.txt \
      --outdir 01_data/raw_fastq \
      --threads 6
"""

import subprocess
import argparse
import sys
import gzip
from pathlib import Path


def check_file(file_path, min_size=1000):
    """
    Verifica que un archivo exista y no esté vacío.
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"No se encontró el archivo {file_path}")

    if file_path.stat().st_size < min_size:
        raise ValueError(f"Archivo muy pequeño o incompleto: {file_path}")

    return True


def verify_fastq(file_path):
    """
    Comprueba que un archivo sea un FASTQ válido.
    Acepta:
    - .fastq
    - .fastq.gz
    """
    file_path = Path(file_path)

    if not (
        str(file_path).endswith(".fastq")
        or str(file_path).endswith(".fastq.gz")
    ):
        raise ValueError(f"{file_path} no tiene una extensión FASTQ válida")

    check_file(file_path)

    try:
        if str(file_path).endswith(".fastq.gz"):
            with gzip.open(file_path, "rt") as f:
                header = f.readline()
        else:
            with open(file_path, "r") as f:
                header = f.readline()

        if not header.startswith("@"):
            raise ValueError(f"{file_path} no parece ser un FASTQ válido")

        return True

    except OSError as e:
        raise RuntimeError(f"Error al abrir/verificar {file_path}: {e}") from e


def compress_with_pigz(file_path, threads=6):
    """
    Comprime un archivo FASTQ usando pigz.
    Si el archivo .gz ya existe y es válido, no repite la compresión.
    """
    file_path = Path(file_path)
    gz_path = Path(str(file_path) + ".gz")

    if gz_path.exists():
        verify_fastq(gz_path)
        return gz_path

    if not file_path.exists():
        raise FileNotFoundError(f"No se encontró el archivo para comprimir: {file_path}")

    subprocess.run(
        ["pigz", "-p", str(threads), str(file_path)],
        check=True
    )

    verify_fastq(gz_path)

    return gz_path


def find_sra_file(srr_dir, srr_id):
    """
    Busca el archivo .sra generado por prefetch.
    """
    srr_dir = Path(srr_dir)

    possible_files = [
        srr_dir / srr_id / f"{srr_id}.sra",
        srr_dir / f"{srr_id}.sra"
    ]

    for file_path in possible_files:
        if file_path.exists():
            return file_path

    sra_files = list(srr_dir.rglob("*.sra"))

    if sra_files:
        return sra_files[0]

    raise FileNotFoundError(f"No se encontró archivo .sra para {srr_id}")


def prefetch_fasterqdump(srr_id, fastq_dir, retries=2, threads=6):
    """
    Descarga un SRR con prefetch, lo convierte a FASTQ con fasterq-dump
    y comprime los FASTQ usando pigz.

    Returns:
        (file_r1, file_r2) si los datos son paired-end.
        (file_single, None) si los datos son single-end.
        (None, None) si falla la descarga o conversión.
    """
    srr_id = srr_id.strip()
    fastq_dir = Path(fastq_dir)

    srr_dir = fastq_dir / srr_id
    srr_dir.mkdir(parents=True, exist_ok=True)

    file_r1 = srr_dir / f"{srr_id}_1.fastq"
    file_r2 = srr_dir / f"{srr_id}_2.fastq"
    file_single = srr_dir / f"{srr_id}.fastq"

    file_r1_gz = srr_dir / f"{srr_id}_1.fastq.gz"
    file_r2_gz = srr_dir / f"{srr_id}_2.fastq.gz"
    file_single_gz = srr_dir / f"{srr_id}.fastq.gz"

    # Si ya existen archivos comprimidos, no repetir descarga
    try:
        if verify_fastq(file_r1_gz) and verify_fastq(file_r2_gz):
            return str(file_r1_gz), str(file_r2_gz)
    except Exception:
        pass

    try:
        if verify_fastq(file_single_gz):
            return str(file_single_gz), None
    except Exception:
        pass

    # Si existen archivos no comprimidos, solo comprimir
    try:
        if verify_fastq(file_r1) and verify_fastq(file_r2):
            r1_gz = compress_with_pigz(file_r1, threads)
            r2_gz = compress_with_pigz(file_r2, threads)
            return str(r1_gz), str(r2_gz)
    except Exception:
        pass

    try:
        if verify_fastq(file_single):
            single_gz = compress_with_pigz(file_single, threads)
            return str(single_gz), None
    except Exception:
        pass

    for intento in range(1, retries + 1):
        try:
            print(f"\n[{srr_id}] Intento {intento}/{retries}")

            # Descargar accession SRA
            subprocess.run(
                ["prefetch", srr_id, "-O", str(srr_dir)],
                check=True
            )

            # Localizar archivo .sra descargado
            sra_file = find_sra_file(srr_dir, srr_id)

            # Convertir a FASTQ
            subprocess.run(
                [
                    "fasterq-dump",
                    str(sra_file),
                    "--split-files",
                    "--threads",
                    str(threads),
                    "-O",
                    str(srr_dir)
                ],
                check=True
            )

            # Verificar y comprimir paired-end
            if verify_fastq(file_r1) and verify_fastq(file_r2):
                r1_gz = compress_with_pigz(file_r1, threads)
                r2_gz = compress_with_pigz(file_r2, threads)
                return str(r1_gz), str(r2_gz)

            # Verificar y comprimir single-end, por si alguna corrida no fuera paired-end
            if verify_fastq(file_single):
                single_gz = compress_with_pigz(file_single, threads)
                return str(single_gz), None

            raise RuntimeError(f"No se encontraron FASTQ válidos para {srr_id}")

        except subprocess.CalledProcessError as e:
            print(
                f"[ERROR] Falló la descarga/conversión de {srr_id}: {e}",
                file=sys.stderr
            )

            if intento == retries:
                return None, None

        except Exception as e:
            print(f"[ERROR] Problema con {srr_id}: {e}", file=sys.stderr)

            if intento == retries:
                return None, None

    return None, None


def main():
    parser = argparse.ArgumentParser(
        description="Descarga SRR con prefetch, convierte a FASTQ y comprime con pigz."
    )

    parser.add_argument(
        "--txt",
        required=True,
        help="Ruta al archivo .txt que contiene los IDs SRR."
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Directorio de salida donde se guardarán los archivos descargados."
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=6,
        help="Número de hilos para fasterq-dump y pigz. Default: 6."
    )

    args = parser.parse_args()

    txt_path = Path(args.txt)
    outdir = Path(args.outdir)

    if not txt_path.exists():
        sys.exit(f"Error: no se encontró el archivo {txt_path}")

    outdir.mkdir(parents=True, exist_ok=True)

    with open(txt_path, "r") as file:
        for line in file:
            srr_id = line.strip()

            if not srr_id:
                continue

            r1, r2 = prefetch_fasterqdump(
                srr_id=srr_id,
                fastq_dir=outdir,
                threads=args.threads
            )

            if r1 and r2:
                print(f"{srr_id}: paired-end descargado y comprimido correctamente")
                print(f"  R1: {r1}")
                print(f"  R2: {r2}")

            elif r1:
                print(f"{srr_id}: single-end descargado y comprimido correctamente")
                print(f"  Archivo: {r1}")

            else:
                print(f"{srr_id}: descarga o conversión falló", file=sys.stderr)


if __name__ == "__main__":
    main()
