# 11_deseq2_dge.R
#
# Author: Alondra Márquez
# Date: 2026-05-31
#
# Description:
#   Realiza el análisis de expresión diferencial con DESeq2 a partir de la
#   matriz de conteos génicos generada con featureCounts.
#
#   El análisis usa una variable combinada de condición experimental, que
#   integra el tiempo de estimulación y la presencia o ausencia de ácido
#   láctico.
#
#   El modelo usado es:
#     design = ~ 0 + condition
#
#   Este diseño permite comparar directamente condiciones específicas del
#   experimento. Se evalúan dos grupos de contrastes:
#     1. Activación M1 respecto a BMDM 0 h.
#     2. Efecto del ácido láctico respecto al control M1 del mismo tiempo.
#
#   El script realiza los siguientes pasos:
#     - Carga la metadata y la matriz de conteos.
#     - Verifica que las muestras coincidan entre ambos archivos.
#     - Filtra genes con baja expresión usando filterByExpr.
#     - Ajusta el modelo de DESeq2 con la variable condition.
#     - Extrae resultados para todos los contrastes definidos.
#     - Clasifica genes como UP, DOWN o NO según FDR y log2FC.
#     - Guarda tablas completas y tablas de genes diferencialmente expresados.
#     - Guarda una matriz VST para visualizaciones posteriores.
#
# Input:
#   - 01_metadata/metadata_dge.tsv
#   - 04_results/09_featurecounts/gene_counts_strand_2.csv
#
# Output:
#   - 04_results/11_deseq2_dge/filter_summary.csv
#   - 04_results/11_deseq2_dge/deseq2_coefficients.txt
#   - 04_results/11_deseq2_dge/dge_summary_by_contrast.csv
#   - 04_results/11_deseq2_dge/vst_matrix.rds
#   - 04_results/11_deseq2_dge/metadata_ordered.tsv
#   - 04_results/11_deseq2_dge/all_results/
#   - 04_results/11_deseq2_dge/significant_results/
#
# Requirements:
#   - DESeq2
#   - edgeR
#   - limma
#   - dplyr
#   - tibble
#   - readr
#
# Use:
#   Rscript 03_scripts/14_deseq2_dge.R

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tibble)
  library(readr)
})

# Rutas de entrada
metadata_file <- "01_metadata/metadata_dge.tsv"
counts_file <- "04_results/09_featurecounts/gene_counts_strand_2.csv"

# Rutas de salida
out_dir <- "04_results/11_deseq2_dge"

# Crear carpetas de salida
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(
  file.path(out_dir, "all_results"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(out_dir, "significant_results"),
  recursive = TRUE,
  showWarnings = FALSE
)

# Archivo de correspondencia gene_id - gene_name
gene_map_file <- "02_data/reference/mm39-gencode-M36-gene_id-gene_name.txt"

# Cargar tabla de anotación de genes
gene_map <- read.table(
  gene_map_file,
  header = FALSE,
  col.names = c("gene_id_clean", "gene_name"),
  stringsAsFactors = FALSE
)

# Eliminar versiones si existen en el archivo de anotación
gene_map$gene_id_clean <- sub("\\..*$", "", gene_map$gene_id_clean)

# Umbrales para clasificar genes diferencialmente expresados
fdr_cutoff <- 0.01
lfc_cutoff <- 0.5

# Cargar metadata
metadata <- read.delim(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Definir orden de condiciones
metadata$condition <- factor(
  metadata$condition,
  levels = c(
    "BMDM_0h",
    "M1_4h",
    "M1_8h",
    "M1_16h",
    "M1_24h",
    "M1_Lac_8h",
    "M1_Lac_16h",
    "M1_Lac_24h"
  )
)

# Definir variables usadas en visualizaciones
metadata$time <- factor(
  metadata$time,
  levels = c("0h", "4h", "8h", "16h", "24h")
)

metadata$lactate <- factor(metadata$lactate, levels = c("no_LA", "LA"))

# Usar sample_id como nombre de fila
rownames(metadata) <- metadata$sample_id

# Cargar matriz de conteos
counts_raw <- read.csv(counts_file, check.names = FALSE)

# Separar gene_id y conteos
gene_id <- counts_raw[[1]]
counts <- counts_raw[, -1, drop = FALSE]

rownames(counts) <- gene_id

# Convertir conteos a matriz numérica entera
counts[] <- lapply(counts, as.numeric)
counts <- as.matrix(counts)
counts <- round(counts)

# Limpiar nombres de columnas si conservan el sufijo del BAM
colnames(counts) <- sub("_Aligned.sortedByCoord.out.bam$", "", colnames(counts))

# Verificar que no haya valores problemáticos
if (any(is.na(counts))) {
  stop("La matriz de conteos contiene NA.")
}

if (any(counts < 0)) {
  stop("La matriz de conteos contiene valores negativos.")
}

# Verificar correspondencia entre conteos y metadata
missing_in_metadata <- setdiff(colnames(counts), rownames(metadata))
missing_in_counts <- setdiff(rownames(metadata), colnames(counts))

if (length(missing_in_metadata) > 0) {
  stop(
    paste(
      "Muestras en conteos pero no en metadata:",
      paste(missing_in_metadata, collapse = ", ")
    )
  )
}

# Reordenar metadata según la matriz de conteos
metadata <- metadata[colnames(counts), , drop = FALSE]

# Verificar que el orden sea idéntico
if (!all(colnames(counts) == rownames(metadata))) {
  stop("El orden de muestras no coincide entre conteos y metadata.")
}

# Guardar metadata ordenada para visualizaciones posteriores
write_tsv(metadata, file.path(out_dir, "metadata_ordered.tsv"))

# Crear matriz de diseño para filterByExpr y contrastes
design_matrix <- model.matrix(~ 0 + condition, data = metadata)

# Filtrar genes con baja expresión
keep <- filterByExpr(counts, design = design_matrix)
counts_filtered <- counts[keep, ]

# Guardar resumen del filtrado
filter_summary <- data.frame(
  genes_before_filter = nrow(counts),
  genes_after_filter = sum(keep),
  genes_removed = nrow(counts) - sum(keep)
)

write.csv(
  filter_summary,
  file = file.path(out_dir, "filter_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Crear objeto DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ 0 + condition
)

# Ajustar modelo
dds <- DESeq(dds)

# Guardar nombres de coeficientes del modelo
writeLines(
  resultsNames(dds),
  con = file.path(out_dir, "deseq2_coefficients.txt")
)

# Definir contrastes
contrast_matrix <- makeContrasts(
  M1_4h_vs_BMDM_0h = conditionM1_4h - conditionBMDM_0h,
  M1_8h_vs_BMDM_0h = conditionM1_8h - conditionBMDM_0h,
  M1_16h_vs_BMDM_0h = conditionM1_16h - conditionBMDM_0h,
  M1_24h_vs_BMDM_0h = conditionM1_24h - conditionBMDM_0h,
  Lac_8h_vs_M1_8h = conditionM1_Lac_8h - conditionM1_8h,
  Lac_16h_vs_M1_16h = conditionM1_Lac_16h - conditionM1_16h,
  Lac_24h_vs_M1_24h = conditionM1_Lac_24h - conditionM1_24h,
  levels = design_matrix
)

# Verificar que los contrastes coincidan con los coeficientes de DESeq2
if (!all(rownames(contrast_matrix) == resultsNames(dds))) {
  stop(
    "Los nombres de la matriz de contrastes no coinciden con resultsNames(dds)."
  )
}

# Clasificar genes según FDR y log2FC
classify_genes <- function(res_df, fdr_cutoff, lfc_cutoff) {
  res_df %>%
    mutate(
      regulation = case_when(
        !is.na(padj) & padj < fdr_cutoff & log2FoldChange > lfc_cutoff ~ "UP",
        !is.na(padj) &
          padj < fdr_cutoff &
          log2FoldChange < -lfc_cutoff ~ "DOWN",
        TRUE ~ "NO"
      )
    )
}

# Ejecutar contrastes
summary_list <- list()

for (contrast_name in colnames(contrast_matrix)) {
  # Extraer resultados del contraste
  res <- results(
    dds,
    contrast = contrast_matrix[, contrast_name],
    alpha = fdr_cutoff
  )

  # Convertir resultados a tabla
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(
      gene_id_clean = sub("\\..*$", "", gene_id)
    ) %>%
    left_join(
      gene_map,
      by = "gene_id_clean"
    ) %>%
    classify_genes(fdr_cutoff, lfc_cutoff) %>%
    select(
      gene_id,
      gene_name,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      regulation,
      everything()
    ) %>%
    arrange(padj)

  # Guardar tabla completa
  write_csv(
    res_df,
    file.path(out_dir, "all_results", paste0(contrast_name, "_all.csv"))
  )

  # Guardar genes diferencialmente expresados
  sig_df <- res_df %>%
    filter(regulation %in% c("UP", "DOWN"))

  write_csv(
    sig_df,
    file.path(out_dir, "significant_results", paste0(contrast_name, "_DEG.csv"))
  )

  # Guardar resumen del contraste
  summary_list[[contrast_name]] <- data.frame(
    contrast = contrast_name,
    tested_genes = sum(!is.na(res_df$padj)),
    up_genes = sum(res_df$regulation == "UP"),
    down_genes = sum(res_df$regulation == "DOWN"),
    total_deg = sum(res_df$regulation %in% c("UP", "DOWN"))
  )
}

# Guardar tabla resumen de todos los contrastes
dge_summary <- bind_rows(summary_list)

write_csv(dge_summary, file.path(out_dir, "dge_summary_by_contrast.csv"))

# Aplicar VST para visualizaciones posteriores
vsd <- vst(dds, blind = FALSE)
vst_matrix <- assay(vsd)

# Guardar matriz VST
saveRDS(vst_matrix, file.path(out_dir, "vst_matrix.rds"))
