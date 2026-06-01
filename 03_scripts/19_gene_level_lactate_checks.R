#!/usr/bin/env Rscript

# 19_gene_level_lactate_checks.R
#
# Objetivo:
#   Revisar genes específicos asociados con activación M1, respuesta a lactato
#   y genes homeostáticos/M2-like reportados en el contexto biológico del dataset.
#
#   Este script NO redefine los resultados de DESeq2. Solo agrega una revisión
#   biológica dirigida para interpretar genes representativos y contrastarlos
#   con el modelo del artículo original.
#
# Entradas:
#   - Resultados completos de DESeq2 para contrastes de lactato.
#   - Matriz VST.
#   - Metadata ordenada.
#   - GTF GENCODE vM36.
#
# Salidas:
#   - Top genes UP/DOWN por contraste.
#   - Tabla de genes candidatos.
#   - Heatmap de genes candidatos por condición.
#   - Figura de log2FC de genes candidatos por contraste.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(forcats)
})


base_dir <- "."

deseq_dir <- file.path(base_dir, "04_results/11_deseq2_dge/all_results")
vst_file <- file.path(base_dir, "04_results/10_dge_qc/vst_matrix.csv")
metadata_file <- file.path(
  base_dir,
  "04_results/11_deseq2_dge/metadata_ordered.tsv"
)
gtf_file <- file.path(
  base_dir,
  "02_data/reference/gencode.vM36.basic.annotation.gtf"
)

out_dir <- file.path(base_dir, "04_results/11_gene_level_checks")
fig_dir <- file.path(base_dir, "05_figures/11_gene_level_checks")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Parámetros generales

fdr_cutoff <- 0.01
lfc_cutoff <- 0.5
top_n <- 10

lactate_contrasts <- c(
  "Lac_8h_vs_M1_8h",
  "Lac_16h_vs_M1_16h",
  "Lac_24h_vs_M1_24h"
)

# Genes dirigidos para revisión biológica.
candidate_genes <- tibble::tribble(
  ~gene_name , ~category                   , ~biological_reason                                                                        ,
  "Arg1"     , "Homeostático / M2-like"    , "Gen destacado en el modelo de lactato/lactilación; metabolismo de arginina y reparación" ,
  "Vegfa"    , "Homeostático / M2-like"    , "Gen M2-like reportado como sensible a lactato en el contexto del artículo"               ,
  "Nos2"     , "M1 inflamatorio"           , "Marcador de activación M1 y producción de óxido nítrico"                                 ,
  "Tnf"      , "M1 inflamatorio"           , "Citocina proinflamatoria asociada a activación M1"                                       ,
  "Il6"      , "M1 inflamatorio"           , "Citocina proinflamatoria asociada a activación M1"                                       ,
  "Il1b"     , "M1 inflamatorio"           , "Citocina inflamatoria asociada a respuesta innata"                                       ,
  "Cxcl10"   , "M1 inflamatorio"           , "Quimiocina inducida en activación inflamatoria"                                          ,
  "Hcar2"    , "Respuesta a metabolitos"   , "Receptor asociado a señales metabólicas; apareció en genes representativos"              ,
  "Dusp1"    , "Regulación / estrés"       , "Regulador de señalización MAPK; apareció entre genes destacados"                         ,
  "Hspa1b"   , "Estrés celular"            , "Proteína de choque térmico; posible respuesta a estrés celular"                          ,
  "Hspa1l"   , "Estrés celular"            , "Proteína de choque térmico; posible respuesta a estrés celular"                          ,
  "Cd9"      , "Estado celular / membrana" , "Gen representativo observado en el heatmap del análisis"
)

write_csv(candidate_genes, file.path(out_dir, "candidate_gene_set.csv"))

#  Crear mapeo gene_id_clean -> gene_name desde GTF

message("Leyendo GTF para completar nombres de genes...")

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
    gene_name_gtf = str_match(attribute, 'gene_name "([^"]+)"')[, 2],
    gene_type = str_match(attribute, 'gene_type "([^"]+)"')[, 2]
  ) %>%
  mutate(
    gene_id_clean = str_remove(gene_id, "\\.[0-9]+$")
  ) %>%
  distinct(gene_id_clean, .keep_all = TRUE)

write_csv(gene_map, file.path(out_dir, "gencode_vM36_gene_map.csv"))


# Leer resultados de DESeq2 y completar gene_name

read_deseq_result <- function(contrast_name) {
  file <- file.path(deseq_dir, paste0(contrast_name, "_all.csv"))

  if (!file.exists(file)) {
    stop("No se encontró el archivo: ", file)
  }

  read_csv(file, show_col_types = FALSE) %>%
    mutate(
      contrast = contrast_name,
      gene_id_clean = ifelse(
        is.na(gene_id_clean),
        str_remove(gene_id, "\\.[0-9]+$"),
        gene_id_clean
      )
    ) %>%
    left_join(
      gene_map %>% select(gene_id_clean, gene_name_gtf, gene_type),
      by = "gene_id_clean"
    ) %>%
    mutate(
      gene_name_final = case_when(
        !is.na(gene_name) & gene_name != "" & gene_name != "NA" ~ gene_name,
        !is.na(gene_name_gtf) & gene_name_gtf != "" ~ gene_name_gtf,
        TRUE ~ gene_id_clean
      ),
      regulation_strict = case_when(
        !is.na(padj) & padj < fdr_cutoff & log2FoldChange > lfc_cutoff ~ "UP",
        !is.na(padj) &
          padj < fdr_cutoff &
          log2FoldChange < -lfc_cutoff ~ "DOWN",
        TRUE ~ "NS"
      )
    )
}

all_lactate_results <- bind_rows(
  lapply(lactate_contrasts, read_deseq_result)
)

write_csv(
  all_lactate_results,
  file.path(out_dir, "all_lactate_contrasts_annotated.csv")
)


# top genes UP y DOWN por contraste

top_up_down <- all_lactate_results %>%
  filter(
    regulation_strict %in% c("UP", "DOWN"),
    !is.na(padj)
  ) %>%
  group_by(contrast, regulation_strict) %>%
  arrange(padj, desc(abs(log2FoldChange)), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup() %>%
  select(
    contrast,
    regulation = regulation_strict,
    gene_name = gene_name_final,
    gene_id,
    gene_id_clean,
    gene_type,
    baseMean,
    log2FoldChange,
    padj
  )

write_csv(
  top_up_down,
  file.path(out_dir, "top10_UP_DOWN_lactate_contrasts.csv")
)


#  Genes candidatos en contrastes de lactato

candidate_results <- all_lactate_results %>%
  filter(gene_name_final %in% candidate_genes$gene_name) %>%
  group_by(contrast, gene_name_final) %>%
  arrange(padj, desc(abs(log2FoldChange)), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  left_join(candidate_genes, by = c("gene_name_final" = "gene_name")) %>%
  mutate(
    significant = case_when(
      !is.na(padj) &
        padj < fdr_cutoff &
        abs(log2FoldChange) > lfc_cutoff ~ "DEG",
      !is.na(padj) ~ "No_DEG",
      TRUE ~ "Not_tested"
    )
  ) %>%
  select(
    contrast,
    gene_name = gene_name_final,
    category,
    biological_reason,
    gene_id,
    gene_id_clean,
    gene_type,
    baseMean,
    log2FoldChange,
    padj,
    regulation = regulation_strict,
    significant
  )

# Asegurar que aparezcan candidatos ausentes.
found_candidates <- unique(candidate_results$gene_name)
missing_candidates <- candidate_genes %>%
  filter(!gene_name %in% found_candidates)

write_csv(
  candidate_results,
  file.path(out_dir, "candidate_genes_lactate_contrasts.csv")
)

write_csv(
  missing_candidates,
  file.path(out_dir, "candidate_genes_not_found.csv")
)

# una fila por gen y columnas por contraste.
candidate_wide <- candidate_results %>%
  select(
    gene_name,
    category,
    contrast,
    log2FoldChange,
    padj,
    regulation,
    significant
  ) %>%
  pivot_wider(
    names_from = contrast,
    values_from = c(log2FoldChange, padj, regulation, significant)
  )

write_csv(
  candidate_wide,
  file.path(out_dir, "candidate_genes_lactate_contrasts_wide.csv")
)


#Heatmap dirigido usando matriz VST

message("Generando heatmap dirigido de genes candidatos...")

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
    gene_map %>% select(gene_id_clean, gene_name_gtf),
    by = "gene_id_clean"
  ) %>%
  mutate(
    gene_name = gene_name_gtf
  )

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

candidate_vst <- vst_tbl %>%
  filter(gene_name %in% candidate_genes$gene_name)

# Si hay duplicados por símbolo, conservar el transcrito/gen con mayor expresión media.
candidate_vst_unique <- candidate_vst %>%
  mutate(
    mean_vst = rowMeans(across(all_of(metadata$sample_id)), na.rm = TRUE)
  ) %>%
  group_by(gene_name) %>%
  arrange(desc(mean_vst), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

vst_long <- candidate_vst_unique %>%
  select(gene_id, gene_id_clean, gene_name, all_of(metadata$sample_id)) %>%
  pivot_longer(
    cols = all_of(metadata$sample_id),
    names_to = "sample_id",
    values_to = "vst"
  ) %>%
  left_join(metadata, by = "sample_id")

candidate_condition_means <- vst_long %>%
  group_by(gene_name, condition) %>%
  summarise(
    mean_vst = mean(vst, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(gene_name) %>%
  mutate(
    z = ifelse(
      sd(mean_vst, na.rm = TRUE) == 0,
      0,
      (mean_vst - mean(mean_vst, na.rm = TRUE)) / sd(mean_vst, na.rm = TRUE)
    )
  ) %>%
  ungroup() %>%
  left_join(candidate_genes, by = "gene_name")

write_csv(
  candidate_condition_means,
  file.path(out_dir, "candidate_gene_mean_vst_by_condition.csv")
)

p_heatmap <- ggplot(
  candidate_condition_means,
  aes(x = condition, y = fct_rev(factor(gene_name)), fill = z)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    name = "Z-score"
  ) +
  labs(
    title = "Expresión relativa de genes candidatos",
    subtitle = "Promedio VST por condición, estandarizado por gen",
    x = "Condición",
    y = "Gen"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  file.path(fig_dir, "candidate_gene_heatmap_mean_vst_by_condition.png"),
  p_heatmap,
  width = 9,
  height = 5.5,
  dpi = 300
)
