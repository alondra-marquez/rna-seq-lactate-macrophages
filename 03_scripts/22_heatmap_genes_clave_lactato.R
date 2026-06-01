#!/usr/bin/env Rscript

# 22_heatmap_genes_clave_pheatmap.R
#
# Objetivo:
#   Generar un heatmap clásico de genes clave asociados al ácido láctico
#   usando valores VST promedio por condición.
#
#   El escalado por fila (row scaling) se hace dentro de pheatmap,
#   por lo que no se necesita calcular Z-score manualmente.
#
# Entradas:
#   - 04_results/10_dge_qc/vst_matrix.csv
#   - 04_results/11_deseq2_dge/metadata_ordered.tsv
#   - 02_data/reference/gencode.vM36.basic.annotation.gtf
#
# Salidas:
#   - Tabla de expresión promedio por condición
#   - Heatmap PNG y PDF

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(pheatmap)
  library(RColorBrewer)
})


vst_file <- "04_results/10_dge_qc/vst_matrix.csv"
metadata_file <- "04_results/11_deseq2_dge/metadata_ordered.tsv"
gtf_file <- "02_data/reference/gencode.vM36.basic.annotation.gtf"

out_dir <- "04_results/12_candidate_gene_heatmap"
fig_dir <- "05_figures/12_candidate_gene_heatmap"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# Genes clave

genes_clave <- tibble::tribble(
  ~gene_name , ~grupo                     ,
  "Nos2"     , "Inflamación M1"           ,
  "Il6"      , "Inflamación M1"           ,
  "Il1b"     , "Inflamación M1"           ,
  "Cxcl10"   , "Inflamación M1"           ,
  "Arg1"     , "Homeostático / M2-like"   ,
  "Vegfa"    , "Homeostático / M2-like"   ,
  "Hcar2"    , "Metabolitos / regulación" ,
  "Dusp1"    , "Metabolitos / regulación" ,
  "Hspa1b"   , "Estrés celular"           ,
  "Hspa1l"   , "Estrés celular"
)

orden_genes <- genes_clave$gene_name

metadata <- read_tsv(metadata_file, show_col_types = FALSE) %>%
  mutate(
    condition = factor(
      condition,
      levels = c(
        "BMDM_0h",
        "M1_4h",
        "M1_8h",
        "M1_Lac_8h",
        "M1_16h",
        "M1_Lac_16h",
        "M1_24h",
        "M1_Lac_24h"
      )
    ),
    time = factor(time, levels = c("0h", "4h", "8h", "16h", "24h")),
    lactate = factor(lactate, levels = c("no_LA", "LA"))
  )

conditions_to_plot <- c(
  "M1_8h",
  "M1_Lac_8h",
  "M1_16h",
  "M1_Lac_16h",
  "M1_24h",
  "M1_Lac_24h"
)

metadata_plot <- metadata %>%
  filter(condition %in% conditions_to_plot)


#leer GTF y mapear IDs a gene_name

gtf_cols <- c(
  "seqname",
  "source",
  "feature",
  "start",
  "end",
  "score",
  "strand",
  "frame",
  "attribute"
)

gtf <- read_tsv(
  gtf_file,
  comment = "#",
  col_names = gtf_cols,
  show_col_types = FALSE
)

gene_map <- gtf %>%
  filter(feature == "gene") %>%
  transmute(
    gene_id = str_match(attribute, 'gene_id "([^"]+)"')[, 2],
    gene_name = str_match(attribute, 'gene_name "([^"]+)"')[, 2]
  ) %>%
  mutate(
    gene_id_clean = str_remove(gene_id, "\\.[0-9]+$")
  ) %>%
  distinct(gene_id_clean, .keep_all = TRUE)


# Leer matriz VST

vst <- read.csv(
  vst_file,
  row.names = 1,
  check.names = FALSE
)

vst_tbl <- vst %>%
  rownames_to_column("gene_id") %>%
  mutate(
    gene_id_clean = str_remove(gene_id, "\\.[0-9]+$")
  ) %>%
  left_join(
    gene_map %>% select(gene_id_clean, gene_name),
    by = "gene_id_clean"
  )

# Conservar genes clave

vst_genes_clave <- vst_tbl %>%
  filter(gene_name %in% genes_clave$gene_name)

genes_no_encontrados <- genes_clave %>%
  filter(!gene_name %in% vst_genes_clave$gene_name)

write_csv(
  genes_no_encontrados,
  file.path(out_dir, "genes_clave_no_encontrados.csv")
)

# Si hay duplicados por símbolo, conservar el de mayor expresión media
vst_genes_clave_unique <- vst_genes_clave %>%
  mutate(
    mean_vst_global = rowMeans(across(all_of(metadata$sample_id)), na.rm = TRUE)
  ) %>%
  group_by(gene_name) %>%
  arrange(desc(mean_vst_global), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

# Calcular VST promedio por condición

vst_long <- vst_genes_clave_unique %>%
  select(gene_name, all_of(metadata_plot$sample_id)) %>%
  pivot_longer(
    cols = all_of(metadata_plot$sample_id),
    names_to = "sample_id",
    values_to = "vst"
  ) %>%
  left_join(metadata_plot, by = "sample_id")

vst_means <- vst_long %>%
  group_by(gene_name, condition) %>%
  summarise(
    mean_vst = mean(vst, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  vst_means,
  file.path(out_dir, "genes_clave_mean_vst_by_condition_long.csv")
)


# Convertir a matriz para heatmap

heatmap_df <- vst_means %>%
  mutate(
    gene_name = factor(gene_name, levels = orden_genes),
    condition = factor(condition, levels = conditions_to_plot)
  ) %>%
  arrange(gene_name, condition) %>%
  pivot_wider(
    names_from = condition,
    values_from = mean_vst
  ) %>%
  as.data.frame()

rownames(heatmap_df) <- heatmap_df$gene_name
heatmap_df$gene_name <- NULL

heatmap_mat <- as.matrix(heatmap_df)

# Reordenar filas exactamente como en orden_genes
heatmap_mat <- heatmap_mat[orden_genes, conditions_to_plot]

write.csv(
  heatmap_mat,
  file.path(out_dir, "genes_clave_mean_vst_matrix.csv"),
  quote = FALSE
)


#  Anotaciones
annotation_col <- data.frame(
  Tiempo = c("8h", "8h", "16h", "16h", "24h", "24h"),
  Lactato = c("no_LA", "LA", "no_LA", "LA", "no_LA", "LA")
)
rownames(annotation_col) <- conditions_to_plot

annotation_row <- genes_clave %>%
  mutate(
    gene_name = factor(gene_name, levels = orden_genes)
  ) %>%
  arrange(gene_name) %>%
  as.data.frame()

rownames(annotation_row) <- annotation_row$gene_name
annotation_row$gene_name <- NULL

# Colores de anotación
ann_colors <- list(
  Tiempo = c(
    "8h" = "#9ecae1",
    "16h" = "#66c2a4",
    "24h" = "#fc8d62"
  ),
  Lactato = c(
    "no_LA" = "#bdbdbd",
    "LA" = "#c51b8a"
  ),
  grupo = c(
    "Inflamación M1" = "#f4a3a8",
    "Homeostático / M2-like" = "#92c5de",
    "Metabolitos / regulación" = "#b2df8a",
    "Estrés celular" = "#fdbf6f"
  )
)

# Heatmap

# Paleta clásica azul-blanco-rojo
heat_colors <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)

png(
  filename = file.path(fig_dir, "heatmap_genes_clave_lactato_pheatmap.png"),
  width = 2200,
  height = 1500,
  res = 220
)

pheatmap(
  heatmap_mat,
  scale = "row", # z-score por gen
  color = heat_colors,
  cluster_rows = FALSE, # mantener orden narrativo
  cluster_cols = FALSE, # mantener orden temporal
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  border_color = NA,
  fontsize = 11,
  fontsize_row = 12,
  fontsize_col = 11,
  angle_col = 45,
  cellwidth = 35,
  cellheight = 28,
  main = "Genes clave asociados al ácido láctico"
)

dev.off()

pdf(
  file = file.path(fig_dir, "heatmap_genes_clave_lactato_pheatmap.pdf"),
  width = 10,
  height = 7
)

pheatmap(
  heatmap_mat,
  scale = "row",
  color = heat_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  border_color = NA,
  fontsize = 11,
  fontsize_row = 12,
  fontsize_col = 11,
  angle_col = 45,
  cellwidth = 35,
  cellheight = 28,
  main = "Genes clave asociados al ácido láctico"
)

dev.off()
