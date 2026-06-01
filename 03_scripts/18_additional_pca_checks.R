#!/usr/bin/env Rscript

# 17_additional_pca_checks.R
#
# Objetivo:
#   Generar PCAs adicionales y revisiones técnicas para evaluar si la
#   estructura global de las muestras se asocia principalmente con variables
#   biológicas o con métricas técnicas disponibles.
#
# Entradas:
#   - 04_results/10_dge_qc/vst_matrix.csv
#   - 04_results/11_deseq2_dge/metadata_ordered.tsv
#   - 01_metadata/GSE115354_sra_table.csv
#   - 04_results/07_star_pe_stats/star_pe_stats.tsv
#   - 04_results/09_featurecounts/featurecounts_strand_2_stats.csv
#
# Salidas:
#   - Tablas integradas de PCA + metadata + métricas técnicas
#   - PCAs coloreadas por condición, tiempo, lactato y métricas técnicas
#   - Boxplots de métricas técnicas por condición

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

vst_file <- file.path(base_dir, "04_results/10_dge_qc/vst_matrix.csv")
metadata_file <- file.path(
  base_dir,
  "04_results/11_deseq2_dge/metadata_ordered.tsv"
)
sra_file <- file.path(base_dir, "01_metadata/GSE115354_sra_table.csv")
star_file <- file.path(
  base_dir,
  "04_results/07_star_pe_stats/star_pe_stats.tsv"
)
fc_file <- file.path(
  base_dir,
  "04_results/09_featurecounts/featurecounts_strand_2_stats.csv"
)

out_dir <- file.path(base_dir, "04_results/10_dge_qc_additional_pca")
fig_dir <- file.path(base_dir, "05_figures/10_dge_qc_additional_pca")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Cargar matriz VST

vst <- read.csv(
  vst_file,
  row.names = 1,
  check.names = FALSE
)

# La matriz viene como genes x muestras; para PCA necesitamos muestras x genes.
vst_t <- t(as.matrix(vst))


# Cargar metadata biológica
metadata <- read_tsv(metadata_file, show_col_types = FALSE) %>%
  mutate(
    sample_id = as.character(sample_id),
    time = factor(time, levels = c("0h", "4h", "8h", "16h", "24h")),
    lactate = factor(lactate, levels = c("no_LA", "LA")),
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
    ),
    srr_number = as.numeric(str_remove(sample_id, "SRR"))
  )

# Verificar correspondencia entre matriz y metadata.
samples_vst <- rownames(vst_t)

if (!all(samples_vst %in% metadata$sample_id)) {
  missing_meta <- setdiff(samples_vst, metadata$sample_id)
  stop(
    "Hay muestras en VST sin metadata: ",
    paste(missing_meta, collapse = ", ")
  )
}

metadata <- metadata %>%
  filter(sample_id %in% samples_vst) %>%
  arrange(match(sample_id, samples_vst))

if (!identical(metadata$sample_id, samples_vst)) {
  stop("El orden de metadata y matriz VST no coincide.")
}

# PCA

pca <- prcomp(vst_t, center = TRUE, scale. = FALSE)

pca_var <- (pca$sdev^2) / sum(pca$sdev^2)
pca_summary <- tibble(
  PC = paste0("PC", seq_along(pca_var)),
  variance_explained = pca_var,
  percent_explained = round(pca_var * 100, 2)
)

pca_coords <- as.data.frame(pca$x) %>%
  rownames_to_column("sample_id") %>%
  as_tibble()


# 5. Cargar métricas STAR

star_raw <- read_tsv(
  star_file,
  col_names = FALSE,
  show_col_types = FALSE
)

# Quitar la fila de encabezado original.
star_data <- star_raw[-1, ]

if (ncol(star_data) == 6) {
  colnames(star_data) <- c(
    "sample_id",
    "input_reads",
    "unique_reads",
    "unique_percent",
    "multi_reads",
    "multi_percent"
  )
} else {
  stop(
    "El archivo STAR no tiene 6 columnas después de quitar el encabezado. Revisar formato."
  )
}

star_stats <- star_data %>%
  mutate(
    sample_id = as.character(sample_id),
    input_reads = as.numeric(input_reads),
    unique_reads = as.numeric(unique_reads),
    unique_percent = as.numeric(unique_percent),
    multi_reads = as.numeric(multi_reads),
    multi_percent = as.numeric(multi_percent)
  )


# Cargar estadísticas de featureCounts

fc_stats_wide <- read_csv(fc_file, show_col_types = FALSE)

fc_stats_long <- fc_stats_wide %>%
  pivot_longer(
    cols = -Status,
    names_to = "bam_file",
    values_to = "count"
  ) %>%
  mutate(
    sample_id = str_extract(bam_file, "SRR[0-9]+"),
    count = as.numeric(count)
  )

fc_summary <- fc_stats_long %>%
  group_by(sample_id) %>%
  summarise(
    assigned_reads = count[Status == "Assigned"][1],
    total_featurecounts = sum(count, na.rm = TRUE),
    assigned_percent_fc = 100 * assigned_reads / total_featurecounts,
    unassigned_no_features = ifelse(
      any(Status == "Unassigned_NoFeatures"),
      count[Status == "Unassigned_NoFeatures"][1],
      NA_real_
    ),
    .groups = "drop"
  )


# Revisar metadata SRA

sra <- read_csv(sra_file, show_col_types = FALSE)

sra_rnaseq <- sra %>%
  filter(`Assay Type` == "RNA-Seq") %>%
  transmute(
    sample_id = Run,
    assay_type = `Assay Type`,
    instrument = Instrument,
    platform = Platform,
    library_layout = LibraryLayout,
    library_selection = LibrarySelection,
    library_source = LibrarySource,
    center_name = `Center Name`,
    bioproject = BioProject,
    sra_study = `SRA Study`,
    biosample = BioSample,
    avg_spot_len = AvgSpotLen,
    bases = Bases,
    bytes = Bytes,
    release_date = ReleaseDate,
    create_date = create_date,
    sra_time = `time_after_lps/ifng_treatment`,
    sra_treatment = treatment
  )

# Resumen para saber qué variables SRA son constantes o variables.
sra_variable_summary <- sra_rnaseq %>%
  summarise(
    across(
      .cols = c(
        instrument,
        platform,
        library_layout,
        library_selection,
        library_source,
        center_name,
        bioproject,
        sra_study
      ),
      .fns = ~ n_distinct(.x, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "sra_variable",
    values_to = "n_distinct_values"
  )

write_csv(
  sra_variable_summary,
  file.path(out_dir, "sra_variable_summary.csv")
)


# Integrar PCA + metadata + métricas técnicas

pca_annot <- pca_coords %>%
  left_join(metadata, by = "sample_id") %>%
  left_join(star_stats, by = "sample_id") %>%
  left_join(fc_summary, by = "sample_id") %>%
  left_join(sra_rnaseq, by = "sample_id") %>%
  mutate(
    log10_input_reads = log10(input_reads),
    log10_assigned_reads = log10(assigned_reads)
  )

write_csv(
  pca_annot,
  file.path(out_dir, "pca_metadata_qc_integrated.csv")
)

write_csv(
  pca_summary,
  file.path(out_dir, "pca_summary_additional.csv")
)


# Funciones de graficación

plot_pca_discrete <- function(
  data,
  color_var,
  title,
  pcx = "PC1",
  pcy = "PC2"
) {
  x_lab <- paste0(
    pcx,
    " (",
    pca_summary$percent_explained[pca_summary$PC == pcx],
    "%)"
  )
  y_lab <- paste0(
    pcy,
    " (",
    pca_summary$percent_explained[pca_summary$PC == pcy],
    "%)"
  )

  ggplot(
    data,
    aes(x = .data[[pcx]], y = .data[[pcy]], color = .data[[color_var]])
  ) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text(
      aes(label = sample_id),
      size = 2.4,
      vjust = -0.7,
      show.legend = FALSE
    ) +
    labs(
      title = title,
      x = x_lab,
      y = y_lab,
      color = color_var
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

plot_pca_continuous <- function(
  data,
  color_var,
  title,
  pcx = "PC1",
  pcy = "PC2"
) {
  x_lab <- paste0(
    pcx,
    " (",
    pca_summary$percent_explained[pca_summary$PC == pcx],
    "%)"
  )
  y_lab <- paste0(
    pcy,
    " (",
    pca_summary$percent_explained[pca_summary$PC == pcy],
    "%)"
  )

  ggplot(
    data,
    aes(x = .data[[pcx]], y = .data[[pcy]], color = .data[[color_var]])
  ) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text(
      aes(label = sample_id),
      size = 2.4,
      vjust = -0.7,
      show.legend = FALSE
    ) +
    labs(
      title = title,
      x = x_lab,
      y = y_lab,
      color = color_var
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

save_plot <- function(plot, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

# PCAs por variables biológicas
p_condition <- plot_pca_discrete(
  pca_annot,
  "condition",
  "PCA coloreada por condición experimental"
)

p_time <- plot_pca_discrete(
  pca_annot,
  "time",
  "PCA coloreada por tiempo"
)

p_lactate <- plot_pca_discrete(
  pca_annot,
  "lactate",
  "PCA coloreada por presencia de ácido láctico"
)

p_treatment <- plot_pca_discrete(
  pca_annot,
  "treatment",
  "PCA coloreada por tratamiento"
)

p_pc34_condition <- plot_pca_discrete(
  pca_annot,
  "condition",
  "PCA PC3 vs PC4 coloreada por condición",
  pcx = "PC3",
  pcy = "PC4"
)

save_plot(p_condition, "PCA_by_condition.png")
save_plot(p_time, "PCA_by_time.png")
save_plot(p_lactate, "PCA_by_lactate.png")
save_plot(p_treatment, "PCA_by_treatment.png")
save_plot(p_pc34_condition, "PCA_PC3_PC4_by_condition.png")

# PCAs por métricas técnicas

p_srr_order <- plot_pca_continuous(
  pca_annot,
  "srr_number",
  "PCA coloreada por orden de accesión SRR"
)

p_input_reads <- plot_pca_continuous(
  pca_annot,
  "log10_input_reads",
  "PCA coloreada por número de lecturas de entrada"
)

p_unique <- plot_pca_continuous(
  pca_annot,
  "unique_percent",
  "PCA coloreada por porcentaje de alineamiento único"
)

p_multi <- plot_pca_continuous(
  pca_annot,
  "multi_percent",
  "PCA coloreada por porcentaje de multimapeo"
)

p_assigned <- plot_pca_continuous(
  pca_annot,
  "log10_assigned_reads",
  "PCA coloreada por lecturas asignadas"
)

p_assigned_pct <- plot_pca_continuous(
  pca_annot,
  "assigned_percent_fc",
  "PCA coloreada por porcentaje de asignación featureCounts"
)

save_plot(p_srr_order, "PCA_by_SRR_order.png")
save_plot(p_input_reads, "PCA_by_input_reads.png")
save_plot(p_unique, "PCA_by_unique_percent.png")
save_plot(p_multi, "PCA_by_multi_percent.png")
save_plot(p_assigned, "PCA_by_assigned_reads.png")
save_plot(p_assigned_pct, "PCA_by_assigned_percent_featureCounts.png")

# Boxplots de métricas técnicas por condición

technical_long <- pca_annot %>%
  select(
    sample_id,
    condition,
    input_reads,
    unique_percent,
    multi_percent,
    assigned_reads,
    assigned_percent_fc
  ) %>%
  pivot_longer(
    cols = c(
      input_reads,
      unique_percent,
      multi_percent,
      assigned_reads,
      assigned_percent_fc
    ),
    names_to = "metric",
    values_to = "value"
  )

p_box <- ggplot(technical_long, aes(x = condition, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  labs(
    title = "Métricas técnicas por condición",
    x = "Condición",
    y = "Valor"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

save_plot(p_box, "technical_metrics_by_condition.png", width = 10, height = 8)

# Resumen técnico por condición

technical_summary_by_condition <- pca_annot %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean_input_reads = mean(input_reads, na.rm = TRUE),
    sd_input_reads = sd(input_reads, na.rm = TRUE),
    mean_unique_percent = mean(unique_percent, na.rm = TRUE),
    sd_unique_percent = sd(unique_percent, na.rm = TRUE),
    mean_multi_percent = mean(multi_percent, na.rm = TRUE),
    sd_multi_percent = sd(multi_percent, na.rm = TRUE),
    mean_assigned_reads = mean(assigned_reads, na.rm = TRUE),
    sd_assigned_reads = sd(assigned_reads, na.rm = TRUE),
    mean_assigned_percent_fc = mean(assigned_percent_fc, na.rm = TRUE),
    sd_assigned_percent_fc = sd(assigned_percent_fc, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  technical_summary_by_condition,
  file.path(out_dir, "technical_summary_by_condition.csv")
)

message("Análisis adicional de PCA terminado.")
message("Tablas en: ", out_dir)
message("Figuras en: ", fig_dir)
