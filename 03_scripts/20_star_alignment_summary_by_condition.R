#!/usr/bin/env Rscript

# 20_star_alignment_summary_by_condition.R
#
# Genera una tabla resumida de métricas de alineamiento STAR por condición.
# Útil para presentación y para revisar si alguna condición tuvo peor alineamiento.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# =========================
# Rutas
# =========================

metadata_file <- "04_results/11_deseq2_dge/metadata_ordered.tsv"
star_file <- "04_results/07_star_pe_stats/star_pe_stats.tsv"

out_dir <- "04_results/07_star_pe_stats"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# Cargar metadata
# =========================

metadata <- read_tsv(metadata_file, show_col_types = FALSE) %>%
  mutate(
    sample_id = as.character(sample_id),
    condition = factor(
      condition,
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
  )

# =========================
# Cargar métricas STAR
# =========================
# Se lee sin encabezado porque en tu archivo el encabezado aparece fusionado
# en la parte de multi_reads / multi_percent.

star_raw <- read_tsv(
  star_file,
  col_names = FALSE,
  show_col_types = FALSE
)

star_data <- star_raw[-1, ]

colnames(star_data) <- c(
  "sample_id",
  "input_reads",
  "unique_reads",
  "unique_percent",
  "multi_reads",
  "multi_percent"
)

star_stats <- star_data %>%
  mutate(
    sample_id = as.character(sample_id),
    input_reads = as.numeric(input_reads),
    unique_reads = as.numeric(unique_reads),
    unique_percent = as.numeric(unique_percent),
    multi_reads = as.numeric(multi_reads),
    multi_percent = as.numeric(multi_percent)
  )

# =========================
# Integrar STAR + metadata
# =========================

star_annotated <- star_stats %>%
  left_join(metadata, by = "sample_id") %>%
  select(
    sample_id,
    condition,
    time,
    treatment,
    lactate,
    input_reads,
    unique_reads,
    unique_percent,
    multi_reads,
    multi_percent
  )

write_csv(
  star_annotated,
  file.path(out_dir, "star_alignment_metrics_with_metadata.csv")
)

# =========================
# Resumen por condición
# =========================

star_summary_by_condition <- star_annotated %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean_input_reads = round(mean(input_reads, na.rm = TRUE)),
    mean_unique_percent = round(mean(unique_percent, na.rm = TRUE), 2),
    sd_unique_percent = round(sd(unique_percent, na.rm = TRUE), 2),
    mean_multi_percent = round(mean(multi_percent, na.rm = TRUE), 2),
    sd_multi_percent = round(sd(multi_percent, na.rm = TRUE), 2),
    min_unique_percent = round(min(unique_percent, na.rm = TRUE), 2),
    max_unique_percent = round(max(unique_percent, na.rm = TRUE), 2),
    .groups = "drop"
  )

write_csv(
  star_summary_by_condition,
  file.path(out_dir, "star_alignment_summary_by_condition.csv")
)

# =========================
# Tabla simplificada para diapositiva
# =========================

star_summary_for_slide <- star_summary_by_condition %>%
  mutate(
    unique_alignment = paste0(
      mean_unique_percent,
      " ± ",
      sd_unique_percent,
      "%"
    ),
    multimapping = paste0(mean_multi_percent, " ± ", sd_multi_percent, "%"),
    unique_range = paste0(min_unique_percent, "–", max_unique_percent, "%")
  ) %>%
  select(
    Condición = condition,
    n,
    `Alineamiento único` = unique_alignment,
    `Multimapeo` = multimapping,
    `Rango alineamiento único` = unique_range
  )

write_csv(
  star_summary_for_slide,
  file.path(out_dir, "star_alignment_summary_for_slide.csv")
)

message("Tablas generadas en: ", out_dir)
message(
  "Archivo principal: star_alignment_summary_for_slide.csv"
)
