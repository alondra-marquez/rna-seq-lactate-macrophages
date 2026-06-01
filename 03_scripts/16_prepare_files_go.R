# 16_prepare_files_go.R
#
# Author: Alondra Márquez
# Date: 2026-05-31
#
# Description:
#   Prepara listas de genes para anotación funcional con GO Biological Process.
#   Usa los resultados completos de DESeq2 de los contrastes asociados con
#   ácido láctico.
#
#   Para cada contraste se generan dos listas:
#     - Genes UP
#     - Genes DOWN
#
#   También se genera un archivo background con los genes evaluados en DESeq2.
#
#   No genera archivos para STRING ni GSEA.
#
# Input:
#   - 04_results/11_deseq2_dge/all_results/
#
# Output:
#   - 04_results/12_go_bp_inputs/gene_lists/
#   - 04_results/12_go_bp_inputs/go_bp_input_summary.csv
# Use:
#   Rscript 03_scripts/16_prepare_files_go.R

# Rutas de entrada y salida
results_dir <- "04_results/11_deseq2_dge/all_results"
out_dir <- "04_results/12_go_bp_inputs"
out_gene_lists <- file.path(out_dir, "gene_lists")

# Crear carpetas de salida
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_gene_lists, recursive = TRUE, showWarnings = FALSE)

# Umbrales usados en el análisis diferencial
FDR_cutoff <- 0.01
LFC_cutoff <- 0.5

# Contrastes principales del proyecto
lactate_contrasts <- c(
  "Lac_8h_vs_M1_8h",
  "Lac_16h_vs_M1_16h",
  "Lac_24h_vs_M1_24h"
)

# Función para guardar listas de genes
write_gene_list <- function(x, file) {
  x <- unique(x)
  x <- x[!is.na(x) & x != ""]
  writeLines(x, con = file)
}

# Lista para guardar el resumen
summary_list <- list()

# Objeto para guardar el background
background_genes <- NULL

# Procesar cada contraste
for (contrast_name in lactate_contrasts) {
  message("contraste: ", contrast_name)

  # Archivo de resultados de DESeq2
  input_file <- file.path(
    results_dir,
    paste0(contrast_name, "_all.csv")
  )

  # Verificar que exista el archivo
  if (!file.exists(input_file)) {
    stop(paste("No se encontró el archivo:", input_file))
  }

  # Cargar tabla de resultados
  res <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Revisar columnas necesarias
  expected_cols <- c(
    "gene_id",
    "log2FoldChange",
    "padj"
  )

  missing_cols <- setdiff(expected_cols, colnames(res))

  if (length(missing_cols) > 0) {
    stop(paste(
      "Faltan columnas en la tabla:",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Limpiar identificadores Ensembl
  res$gene_id <- sub("\\..*$", "", res$gene_id)

  # Eliminar genes sin identificador válido
  res <- res[!is.na(res$gene_id) & res$gene_id != "", ]

  # Guardar background a partir del primer contraste
  if (is.null(background_genes)) {
    background_genes <- unique(res$gene_id)
  }

  # Recalcular clasificación UP/DOWN con los umbrales definidos
  res$regulation_check <- "NO"

  res$regulation_check[
    !is.na(res$padj) &
      res$padj < FDR_cutoff &
      res$log2FoldChange > LFC_cutoff
  ] <- "UP"

  res$regulation_check[
    !is.na(res$padj) &
      res$padj < FDR_cutoff &
      res$log2FoldChange < -LFC_cutoff
  ] <- "DOWN"

  # Separar genes por dirección
  up <- res[res$regulation_check == "UP", ]
  down <- res[res$regulation_check == "DOWN", ]

  # Definir nombres de salida
  up_file <- file.path(
    out_gene_lists,
    paste0(contrast_name, "_UP_ensembl.txt")
  )

  down_file <- file.path(
    out_gene_lists,
    paste0(contrast_name, "_DOWN_ensembl.txt")
  )

  # Guardar listas para GO
  write_gene_list(up$gene_id, up_file)
  write_gene_list(down$gene_id, down_file)

  # Guardar resumen del contraste
  summary_list[[paste0(contrast_name, "_UP")]] <- data.frame(
    contrast = contrast_name,
    direction = "UP",
    FDR = FDR_cutoff,
    LFC = LFC_cutoff,
    genes_in_list = length(unique(up$gene_id)),
    output_file = up_file
  )

  summary_list[[paste0(contrast_name, "_DOWN")]] <- data.frame(
    contrast = contrast_name,
    direction = "DOWN",
    FDR = FDR_cutoff,
    LFC = LFC_cutoff,
    genes_in_list = length(unique(down$gene_id)),
    output_file = down_file
  )
}

# Guardar background para GO
background_file <- file.path(
  out_gene_lists,
  "background_filterByExpr_ensembl.txt"
)

write_gene_list(background_genes, background_file)

# Crear tabla resumen
summary_df <- do.call(rbind, summary_list)

summary_df$background_genes <- length(unique(background_genes))
summary_df$background_file <- background_file

# Guardar resumen para el reporte
write.csv(
  summary_df,
  file.path(out_dir, "go_bp_input_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

message("Listas para GO Biological Process generadas correctamente.")
