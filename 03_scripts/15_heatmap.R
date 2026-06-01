# 15_heatmap.R
#
# Author: Alondra Márquez
# Date: 2026-05-31
#
# Description:
#   Genera un heatmap enfocado en los genes más significativos asociados al
#   ácido láctico. Usa la matriz VST generada en el análisis de DESeq2 y las
#   tablas de resultados de los contrastes de lactato.
#
#   Los genes se seleccionan a partir de los contrastes:
#     - Lac_8h_vs_M1_8h
#     - Lac_16h_vs_M1_16h
#     - Lac_24h_vs_M1_24h
#
#   Para cada contraste se toman los 10 genes diferencialmente expresados más
#   significativos. Después se calcula un Z-score por gen para visualizar
#   cambios relativos de expresión entre muestras.
#
#   Para hacer la figura más clara, el heatmap solo incluye las condiciones
#   directamente usadas en los contrastes de lactato:
#     - M1_8h y M1_Lac_8h
#     - M1_16h y M1_Lac_16h
#     - M1_24h y M1_Lac_24h
#
# Input:
#   - 04_results/11_deseq2_dge/vst_matrix.rds
#   - 04_results/11_deseq2_dge/metadata_ordered.tsv
#   - 04_results/11_deseq2_dge/all_results/
#
# Output:
#   - 05_figures/11_deseq2_dge/heatmap_top_lactate_genes.pdf
#   - 04_results/11_deseq2_dge/top_lactate_genes_heatmap.csv
#
# Requirements:
#   - dplyr
#   - readr
#   - ComplexHeatmap
#   - circlize
#   - grid
#
# Use:
#   Rscript 03_scripts/15_heatmap.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# Rutas de entrada
out_dir <- "04_results/11_deseq2_dge"
results_dir <- file.path(out_dir, "all_results")
vst_file <- file.path(out_dir, "vst_matrix.rds")
metadata_file <- file.path(out_dir, "metadata_ordered.tsv")

# Ruta de salida
fig_dir <- "05_figures/11_deseq2_dge"

# Crear carpeta de salida
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Contrastes principales del proyecto
lactate_contrasts <- c(
  "Lac_8h_vs_M1_8h",
  "Lac_16h_vs_M1_16h",
  "Lac_24h_vs_M1_24h"
)

# Condiciones que se mostrarán en el heatmap
lactate_conditions <- c(
  "M1_8h",
  "M1_Lac_8h",
  "M1_16h",
  "M1_Lac_16h",
  "M1_24h",
  "M1_Lac_24h"
)

# Paleta de colores para condición
condition_cols <- c(
  "BMDM_0h" = "#D9D9D9",
  "M1_4h" = "#BDBDBD",
  "M1_8h" = "#80B1D3",
  "M1_Lac_8h" = "#B3CDE3",
  "M1_16h" = "#66C2A5",
  "M1_Lac_16h" = "#B2DF8A",
  "M1_24h" = "#FC8D62",
  "M1_Lac_24h" = "#F4A6C8"
)

# Paleta de colores para tiempo
time_cols <- c(
  "0h" = "#D9D9D9",
  "4h" = "#BDBDBD",
  "8h" = "#80B1D3",
  "16h" = "#66C2A5",
  "24h" = "#FC8D62"
)

# Paleta de colores para tratamiento con ácido láctico
lactate_cols <- c(
  "no_LA" = "#BDBDBD",
  "LA" = "#C51B7D"
)

# Cargar matriz VST
vst_matrix <- readRDS(vst_file)

# Cargar metadata ordenada
metadata <- read_tsv(metadata_file, show_col_types = FALSE)

# Convertir metadata a data.frame para usar rownames
metadata <- as.data.frame(metadata)

# Usar sample_id como nombre de fila
rownames(metadata) <- metadata$sample_id

# Reordenar metadata según la matriz VST
metadata <- metadata[colnames(vst_matrix), , drop = FALSE]

# Definir variables para anotaciones
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

metadata$time <- factor(
  metadata$time,
  levels = c("0h", "4h", "8h", "16h", "24h")
)

metadata$lactate <- factor(
  metadata$lactate,
  levels = c("no_LA", "LA")
)

# Conservar solo muestras relevantes para lactato
keep_samples <- metadata$condition %in% lactate_conditions
metadata <- metadata[keep_samples, , drop = FALSE]
vst_matrix <- vst_matrix[, rownames(metadata), drop = FALSE]

# Eliminar niveles no usados después del filtrado
metadata$condition <- factor(
  as.character(metadata$condition),
  levels = lactate_conditions
)

metadata$time <- droplevels(metadata$time)
metadata$lactate <- droplevels(metadata$lactate)

# Conservar solo colores de categorías presentes
condition_cols_use <- condition_cols[levels(metadata$condition)]
time_cols_use <- time_cols[levels(metadata$time)]
lactate_cols_use <- lactate_cols[levels(metadata$lactate)]

# Leer resultados de contrastes de lactato
top_lactate_table <- bind_rows(
  lapply(lactate_contrasts, function(contrast_name) {
    file <- file.path(
      results_dir,
      paste0(contrast_name, "_all.csv")
    )

    if (!file.exists(file)) {
      stop(paste0("No se encontró el archivo: ", file))
    }

    read_csv(file, show_col_types = FALSE) %>%
      filter(regulation %in% c("UP", "DOWN")) %>%
      arrange(padj) %>%
      slice_head(n = 10) %>%
      mutate(contrast = contrast_name)
  })
) %>%
  distinct(gene_id, .keep_all = TRUE)

# Extraer genes seleccionados
top_lactate_genes <- top_lactate_table$gene_id

# Conservar solo genes presentes en la matriz VST
top_lactate_genes <- top_lactate_genes[
  top_lactate_genes %in% rownames(vst_matrix)
]

# Detener el script si no hay suficientes genes
if (length(top_lactate_genes) < 2) {
  stop("No hay suficientes genes para generar el heatmap de lactato.")
}

# Extraer matriz para el heatmap
heatmap_matrix <- vst_matrix[top_lactate_genes, , drop = FALSE]

# Calcular Z-score por gen
heatmap_matrix_z <- t(scale(t(heatmap_matrix)))

# Eliminar genes con valores no finitos
finite_rows <- apply(
  heatmap_matrix_z,
  1,
  function(x) all(is.finite(x))
)

heatmap_matrix_z <- heatmap_matrix_z[finite_rows, , drop = FALSE]

# Reordenar tabla de genes según el heatmap
top_lactate_table <- top_lactate_table[
  match(rownames(heatmap_matrix_z), top_lactate_table$gene_id),
]

# Usar gene_name como etiqueta si está disponible
row_labels <- top_lactate_table$gene_name

# Usar gene_id cuando no haya gene_name
row_labels[is.na(row_labels) | row_labels == ""] <- top_lactate_table$gene_id[
  is.na(row_labels) | row_labels == ""
]

# Crear anotación superior
ha <- HeatmapAnnotation(
  condition = metadata$condition,
  time = metadata$time,
  lactate = metadata$lactate,
  col = list(
    condition = condition_cols_use,
    time = time_cols_use,
    lactate = lactate_cols_use
  ),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    condition = list(title = "Condición"),
    time = list(title = "Tiempo"),
    lactate = list(title = "Lactato")
  )
)

# Crear heatmap
ht <- Heatmap(
  heatmap_matrix_z,
  name = "Z",
  col = colorRamp2(
    c(-2, 0, 2),
    c("#3B4CC0", "white", "#D73027")
  ),
  top_annotation = ha,
  column_split = metadata$time,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  row_labels = row_labels,
  row_names_gp = gpar(fontsize = 6),
  show_column_names = FALSE,
  heatmap_legend_param = list(
    title = "Z",
    legend_height = unit(2.2, "cm")
  )
)

# Guardar heatmap en PDF
pdf(
  file = file.path(fig_dir, "heatmap_top_lactate_genes.pdf"),
  width = 8,
  height = 5.5
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

# Guardar genes incluidos en el heatmap
write_csv(
  top_lactate_table,
  file.path(out_dir, "top_lactate_genes_heatmap.csv")
)
