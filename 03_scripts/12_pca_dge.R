# 10_pca_dge.R
#
# Author: Alondra Márquez
# Date: 2026-05-31
#
# Description:
#   Verifica la correspondencia entre metadata y matriz de conteos,
#   filtra genes con baja expresión, aplica VST y genera un PCA
#   para revisar la estructura global de las muestras.
#
# Use:
#   Rscript 03_scripts/13_pca_dge.R

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

# Rutas de entrada
metadata_file <- "01_metadata/metadata_dge.tsv"
counts_file <- "04_results/09_featurecounts/gene_counts_strand_2.csv"

# Rutas de salida
fig_dir <- "05_figures/10_dge_qc"
stats_dir <- "04_results/10_dge_qc"

# Crear carpetas de salida
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

# Cargar metadata
metadata <- read.delim(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Definir el orden de las condiciones
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

# Definir sample_id como nombre de fila
rownames(metadata) <- metadata$sample_id

# Cargar matriz de conteos
counts_raw <- read.csv(
  counts_file,
  check.names = FALSE
)

# Separar identificadores de genes y conteos
gene_id <- counts_raw[[1]]
counts <- counts_raw[, -1, drop = FALSE]
rownames(counts) <- gene_id

# Convertir conteos a matriz numérica entera
counts[] <- lapply(counts, as.numeric)
counts <- as.matrix(counts)
counts <- round(counts)

# Limpiar nombres de columnas si conservan el sufijo del BAM
colnames(counts) <- sub(
  "_Aligned.sortedByCoord.out.bam$",
  "",
  colnames(counts)
)

# Verificar que no haya valores que puedan ser problematicos
if (any(is.na(counts))) {
  stop("La matriz de conteos contiene NA después de la conversión.")
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
      "Estas muestras están en la matriz de conteos pero no en la metadata:",
      paste(missing_in_metadata, collapse = ", ")
    )
  )
}

if (length(missing_in_counts) > 0) {
  stop(
    paste(
      "Estas muestras están en la metadata pero no en la matriz de conteos:",
      paste(missing_in_counts, collapse = ", ")
    )
  )
}

# Reordenar metadata según la matriz de conteos
metadata <- metadata[colnames(counts), , drop = FALSE]

# Confirmar que el orden sea idéntico
if (!all(colnames(counts) == rownames(metadata))) {
  stop("El orden de las muestras no coincide entre conteos y metadata.")
}

# Guardar tabla de distribución de condiciones
condition_table <- as.data.frame(table(metadata$condition))
colnames(condition_table) <- c("condition", "n_samples")

write.csv(
  condition_table,
  file = file.path(stats_dir, "condition_table.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Crear matriz de diseño
design_matrix <- model.matrix(~ 0 + condition, data = metadata)

# Filtrar genes con baja expresión
keep <- filterByExpr(counts, design = design_matrix)

cat("Genes antes de filtrar:", nrow(counts), "\n")
cat("Genes después de filtrar:", sum(keep), "\n")

counts_filtered <- counts[keep, ]

# Guardar resumen del filtrado
filter_summary <- data.frame(
  genes_before_filter = nrow(counts),
  genes_after_filter = sum(keep),
  genes_removed = nrow(counts) - sum(keep)
)

write.csv(
  filter_summary,
  file = file.path(stats_dir, "filter_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Crear objeto DESeq2 para VST y PCA
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ 0 + condition
)

# Aplicar VST
vsd <- vst(dds)

# Guardar matriz VST para visualizaciones posteriores
vst_matrix <- assay(vsd)

write.csv(
  vst_matrix,
  file = file.path(stats_dir, "vst_matrix.csv"),
  quote = FALSE
)

# Calcular PCA
pca <- prcomp(t(vst_matrix))

percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)

metadata_plot <- metadata[rownames(pca$x), , drop = FALSE]

# Evitar duplicar sample_id porque ya se agregará manualmente al PCA
metadata_plot$sample_id <- NULL

pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  metadata_plot,
  check.names = FALSE
)

# Guardar coordenadas del PCA
write.csv(
  pca_df,
  file = file.path(stats_dir, "pca_coordinates.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Guardar resumen del PCA
pca_summary <- data.frame(
  PC1_percent = percent_var[1],
  PC2_percent = percent_var[2]
)

write.csv(
  pca_summary,
  file = file.path(stats_dir, "pca_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Generar PCA
p <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = condition)
) +
  geom_point(size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = sample_id),
    size = 2.4,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  labs(
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)"),
    color = "Condición"
  ) +
  theme_classic(base_size = 10, base_family = "Liberation Serif") +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Guardar figura
ggsave(
  filename = file.path(fig_dir, "PCA_DGE_global.png"),
  plot = p,
  width = 7,
  height = 5,
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, "PCA_DGE_global.pdf"),
  plot = p,
  width = 7,
  height = 5,
  bg = "white",
  device = cairo_pdf
)
