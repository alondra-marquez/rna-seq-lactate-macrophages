
"""
Author: Alondra Márquez
Date: 2026-05-29
Description:
    Filtra la SraRunTable de GSE115354 para conservar únicamente las muestras de RNA-seq.
    Genera una metadata simple y una lista de corridas SRR para descargar.

Use:
    python3 03_scripts/01_filter_metadata.py \
      --csv metadata/SraRunTable_original.csv \
      --outdir metadata/
"""

import sys
import re
import pandas as pd
import argparse
from pathlib import Path


def clean_col_name(df):
    """
    Mejora la sintaxis del nombre de las columnas para evitar errores en el filtrado.
    """
    df = df.copy()
    new_cols = []

    for col in df.columns:
        s = str(col).strip().lower()

        if s.startswith("characteristics:"):
            s = s.replace("characteristics:", "", 1).strip()

        # Reemplaza espacios, diagonales, paréntesis y otros símbolos por "_"
        s = re.sub(r"[^a-z0-9]+", "_", s)

        # Evita dobles guiones bajos y guiones bajos al inicio/final
        s = re.sub(r"_+", "_", s).strip("_")

        new_cols.append(s)

    df.columns = new_cols
    return df


def clean_values(df):
    """
    Limpia espacios y caracteres escapados en columnas de texto.
    """
    df = df.copy()

    for col in df.columns:
        df[col] = (
            df[col]
            .astype(str)
            .str.strip()
            .str.replace("\\,", ",", regex=False)
        )

    return df


def filter_samples(df):
    """
    Conserva únicamente las muestras de RNA-seq.
    """
    if "assay_type" not in df.columns:
        sys.exit("Error: no se encontró la columna 'assay_type'.")

    filtered = df.loc[
        df["assay_type"].str.lower() == "rna-seq"
    ].copy()


    return filtered


def write_csv(df, outdir):
    """
    Escribe una metadata simple con las columnas más importantes.
    """
    output_file = outdir / "samples_metadata.tsv"

    selected_cols = [
        "sample_name",
        "run",
        "experiment",
        "biosample",
        "bioproject",
        "sra_study",
        "organism",
        "source_name",
        "cell_type",
        "strain_background",
        "assay_type",
        "librarylayout",
        "libraryselection",
        "librarysource",
        "instrument",
        "time_after_lps_ifng_treatment",
        "treatment"
    ]

    # Conserva solo las columnas que sí existan en el archivo
    selected_cols = [col for col in selected_cols if col in df.columns]

    selected = df[selected_cols].copy()

    selected = selected.rename(columns={
        "sample_name": "geo_accession",
        "run": "run_accession",
        "librarylayout": "library_layout",
        "libraryselection": "library_selection",
        "librarysource": "library_source",
        "time_after_lps_ifng_treatment": "time"
    })

    output_file.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(output_file, sep="\t", index=False)


def write_txt(df, outdir):
    """
    Escribe un archivo .txt con los SRR de las corridas seleccionadas.
    """
    output_file = outdir / "runs_to_download.txt"

    if "run" not in df.columns:
        sys.exit("Error: no se encontró la columna 'run'.")

    runs = df["run"].dropna().astype(str).str.strip()
    runs.to_csv(output_file, index=False, header=False)



def main():
    parser = argparse.ArgumentParser(
        description="Filtra la SraRunTable de GSE115354 para seleccionar muestras RNA-seq."
    )

    parser.add_argument(
        "--csv",
        required=True,
        help="Ruta al archivo CSV de metadata descargado desde SRA Run Selector."
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Directorio de salida donde se guardarán los archivos generados."
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        df = pd.read_csv(Path(args.csv), dtype=str)
    except FileNotFoundError:
        sys.exit("Error: incorrect path, file was not found.")

    df = clean_col_name(df)

    df = clean_values(df)

    filtered = filter_samples(df)

    write_csv(filtered, outdir)
    write_txt(filtered, outdir)


if __name__ == "__main__":
    main()
