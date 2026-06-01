# 17_filter_go_results.R
#
# Author: Alondra Márquez
# Date: 2026-05-31
#
# Description:
#   Filtra resultados de enriquecimiento GO Biological Process obtenidos
#   en ShinyGO.
#
#   El script procesa las tablas descargadas para genes UP y DOWN de los
#   contrastes asociados con ácido láctico.
#
#   Solo considera GO Biological Process.
#
# Input:
#   - 04_results/12_go_bp_results/raw/
#
# Output:
#   - 04_results/12_go_bp_results/filtered/
#   - 04_results/12_go_bp_results/go_bp_terms_for_report_combined.csv
#
# Use:
#   Rscript 03_scripts/17_filter_go_results.R

# Rutas de entrada y salida
raw_dir <- "04_results/12_go_bp_results/raw"
filtered_dir <- "04_results/12_go_bp_results/filtered"
out_combined <- "04_results/12_go_bp_results/go_bp_terms_for_report_combined.csv"

# Crear carpeta de salida
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

# Parámetros de filtrado
FDR_cutoff <- 0.05
min_genes <- 5
top_n <- 10

# Contrastes principales
lactate_contrasts <- c(
  "Lac_8h_vs_M1_8h",
  "Lac_16h_vs_M1_16h",
  "Lac_24h_vs_M1_24h"
)

# Direcciones de cambio
directions <- c("UP", "DOWN")

# Columnas esperadas en las tablas de ShinyGO
expected_cols <- c(
  "Enrichment FDR",
  "nGenes",
  "Pathway Genes",
  "Fold Enrichment",
  "Pathway",
  "URL",
  "Genes"
)

# Función para filtrar una tabla GO BP
filter_go_bp_table <- function(input_file, output_file, contrast, direction) {
  # Verificar que exista el archivo
  if (!file.exists(input_file)) {
    stop(paste("No se encontró el archivo:", input_file))
  }

  # Cargar tabla de ShinyGO
  go <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Revisar columnas necesarias
  missing_cols <- setdiff(expected_cols, colnames(go))

  if (length(missing_cols) > 0) {
    stop(paste(
      "Faltan columnas en la tabla:",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Convertir columnas numéricas
  go$`Enrichment FDR` <- as.numeric(go$`Enrichment FDR`)
  go$nGenes <- as.numeric(go$nGenes)
  go$`Pathway Genes` <- as.numeric(go$`Pathway Genes`)
  go$`Fold Enrichment` <- as.numeric(go$`Fold Enrichment`)

  # Extraer GO ID
  go$GO_ID <- ifelse(
    grepl("^GO:[0-9]+", go$Pathway),
    sub("^(GO:[0-9]+).*", "\\1", go$Pathway),
    NA
  )

  # Extraer nombre del término
  go$Termino_GO <- trimws(
    sub("^GO:[0-9]+\\s+", "", go$Pathway)
  )

  # Extraer genes representativos
  genes_split <- strsplit(trimws(go$Genes), "[,;[:space:]]+")

  go$Genes_representativos <- sapply(
    genes_split,
    function(x) {
      x <- unique(x[x != ""])
      paste(head(x, 8), collapse = "; ")
    }
  )

  # Filtrar términos significativos
  go_filtered <- go[
    !is.na(go$`Enrichment FDR`) &
      go$`Enrichment FDR` < FDR_cutoff &
      !is.na(go$nGenes) &
      go$nGenes >= min_genes,
  ]

  # Ordenar términos para el reporte
  go_filtered <- go_filtered[
    order(
      go_filtered$`Enrichment FDR`,
      -go_filtered$nGenes,
      -go_filtered$`Fold Enrichment`
    ),
  ]

  # Conservar términos principales
  go_filtered <- head(go_filtered, top_n)

  # Construir tabla final
  if (nrow(go_filtered) > 0) {
    go_report <- go_filtered[, c(
      "GO_ID",
      "Termino_GO",
      "Enrichment FDR",
      "nGenes",
      "Pathway Genes",
      "Fold Enrichment",
      "Genes_representativos"
    )]

    colnames(go_report) <- c(
      "GO_ID",
      "Termino_GO",
      "FDR",
      "Genes_en_lista",
      "Genes_en_termino",
      "Fold_enrichment",
      "Genes_representativos"
    )

    # Redondear valores para presentación
    go_report$FDR <- signif(go_report$FDR, 3)
    go_report$Fold_enrichment <- round(go_report$Fold_enrichment, 1)
  } else {
    # Crear tabla vacía si no hay términos significativos
    go_report <- data.frame(
      GO_ID = character(),
      Termino_GO = character(),
      FDR = numeric(),
      Genes_en_lista = numeric(),
      Genes_en_termino = numeric(),
      Fold_enrichment = numeric(),
      Genes_representativos = character()
    )
  }

  # Agregar información del contraste
  go_report <- cbind(
    contrast = contrast,
    direction = direction,
    Categoria_GO = "Biological Process",
    go_report
  )

  # Guardar tabla filtrada
  write.csv(
    go_report,
    output_file,
    row.names = FALSE
  )

  cat("\nTabla generada:", contrast, direction, "\n")
  cat("Archivo:", output_file, "\n")
  cat("Términos seleccionados:", nrow(go_report), "\n")

  return(go_report)
}

# Lista para guardar resultados
go_reports <- list()

# Procesar todos los contrastes y direcciones
for (contrast_name in lactate_contrasts) {
  for (direction in directions) {
    # Definir archivos de entrada y salida
    input_file <- file.path(
      raw_dir,
      paste0(contrast_name, "_", direction, "_BP.csv")
    )

    output_file <- file.path(
      filtered_dir,
      paste0(contrast_name, "_", direction, "_BP_terms_for_report.csv")
    )

    # Filtrar tabla
    go_reports[[paste0(contrast_name, "_", direction)]] <- filter_go_bp_table(
      input_file = input_file,
      output_file = output_file,
      contrast = contrast_name,
      direction = direction
    )
  }
}

# Combinar todas las tablas
go_combined <- do.call(rbind, go_reports)

# Guardar tabla combinada
write.csv(
  go_combined,
  out_combined,
  row.names = FALSE
)

cat("\nFiltrado de GO Biological Process terminado correctamente.\n")
cat("\nTabla combinada:\n")
cat(out_combined, "\n")
