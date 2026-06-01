#!/usr/bin/env Rscript

# 23_prepare_string_lists_24h.R
#
# Objetivo:
#   Preparar listas de genes UP y DOWN del contraste Lac_24h_vs_M1_24h
#   para análisis en STRING.
#
#   Se generan listas top 100 y top 200 ordenadas por:
#   1) padj más bajo
#   2) mayor |log2FoldChange|
#   3) mayor baseMean
#
# Entrada:
#   - 04_results/11_deseq2_dge/all_results/Lac_24h_vs_M1_24h_all.csv
#   - 02_data/reference/gencode.vM36.basic.annotation.gtf
#
# Salida:
#   - Listas .txt para STRING
#   - Tablas .csv con información completa

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})


deseq_file <- "04_results/11_deseq2_dge/all_results/Lac_24h_vs_M1_24h_all.csv"
gtf_file <- "02_data/reference/gencode.vM36.basic.annotation.gtf"

out_dir <- "04_results/13_string_lists"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#Parámetros

fdr_cutoff <- 0.01
lfc_cutoff <- 0.5
top_n_values <- c(100, 200)

contrast_name <- "Lac_24h_vs_M1_24h"

# Leer GTF para completar gene_name y gene_type

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

#Leer resultados DESeq2 y anotar

res <- read_csv(deseq_file, show_col_types = FALSE) %>%
  mutate(
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
    gene_symbol = case_when(
      !is.na(gene_name) & gene_name != "" & gene_name != "NA" ~ gene_name,
      !is.na(gene_name_gtf) & gene_name_gtf != "" ~ gene_name_gtf,
      TRUE ~ NA_character_
    ),
    regulation_strict = case_when(
      !is.na(padj) & padj < fdr_cutoff & log2FoldChange > lfc_cutoff ~ "UP",
      !is.na(padj) & padj < fdr_cutoff & log2FoldChange < -lfc_cutoff ~ "DOWN",
      TRUE ~ "NS"
    ),
    abs_log2FC = abs(log2FoldChange),
    rank_score = -log10(padj) * abs_log2FC
  )

# Guardar resultado anotado completo
write_csv(
  res,
  file.path(out_dir, paste0(contrast_name, "_annotated_for_string.csv"))
)

# Filtrar genes adecuados para STRING

# Para STRING tiene sentido priorizar genes protein_coding, porque la red
# representa proteínas o asociaciones funcionales entre proteínas.

res_string <- res %>%
  filter(
    regulation_strict %in% c("UP", "DOWN"),
    !is.na(gene_symbol),
    gene_symbol != "",
    gene_type == "protein_coding"
  ) %>%
  distinct(gene_symbol, .keep_all = TRUE)


# Generar listas UP/DOWN top 100 y top 200

make_string_list <- function(data, direction, n_top) {
  top_genes <- data %>%
    filter(regulation_strict == direction) %>%
    arrange(padj, desc(abs_log2FC), desc(baseMean)) %>%
    slice_head(n = n_top)

  # Tabla completa para revisar
  csv_file <- file.path(
    out_dir,
    paste0(contrast_name, "_", direction, "_top", n_top, "_table.csv")
  )

  write_csv(
    top_genes %>%
      select(
        gene_symbol,
        gene_id,
        gene_id_clean,
        gene_type,
        baseMean,
        log2FoldChange,
        padj,
        regulation = regulation_strict,
        abs_log2FC,
        rank_score
      ),
    csv_file
  )

  # Lista simple para STRING: un gen por línea
  txt_file <- file.path(
    out_dir,
    paste0(contrast_name, "_", direction, "_top", n_top, "_STRING.txt")
  )

  write_lines(top_genes$gene_symbol, txt_file)

  message("Generado: ", txt_file)
}

for (n_top in top_n_values) {
  make_string_list(res_string, "UP", n_top)
  make_string_list(res_string, "DOWN", n_top)
}

summary_table <- res_string %>%
  count(regulation_strict, name = "n_protein_coding_genes")

write_csv(
  summary_table,
  file.path(out_dir, paste0(contrast_name, "_STRING_summary.csv"))
)

message("Listas para STRING generadas en: ", out_dir)
