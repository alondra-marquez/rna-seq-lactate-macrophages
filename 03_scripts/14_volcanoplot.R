# 14_volcanoplot.R
#
# Author: Alondra Márquez
# Date: 2026-05-31
#
# Description:
#   Genera volcano plots a partir de las tablas completas producidas por
#   DESeq2. El script crea un volcano plot por contraste y un collage con los
#   tres contrastes principales asociados al ácido láctico.
#
# Input:
#   - 04_results/11_deseq2_dge/all_results/
#
# Output:
#   - 05_figures/11_deseq2_dge/volcano_*.png
#   - 05_figures/11_deseq2_dge/volcano_lactate_collage.png
#
# Requirements:
#   - dplyr
#   - readr
#   - ggplot2
#   - patchwork
#
# Use:
#   Rscript 03_scripts/14_volcanoplot.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

# Rutas de entrada y salida
results_dir <- "04_results/11_deseq2_dge/all_results"
fig_dir <- "05_figures/11_deseq2_dge"

# Crear carpeta de salida
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Umbrales usados en DESeq2
fdr_cutoff <- 0.01
lfc_cutoff <- 0.5

# Colores para volcano plots
vpcolors <- c(
  "NO" = "gray65",
  "DOWN" = "#25BFC1",
  "UP" = "#F08CB4"
)

# Orden de contrastes
contrast_order <- c(
  "M1_4h_vs_BMDM_0h",
  "M1_8h_vs_BMDM_0h",
  "M1_16h_vs_BMDM_0h",
  "M1_24h_vs_BMDM_0h",
  "Lac_8h_vs_M1_8h",
  "Lac_16h_vs_M1_16h",
  "Lac_24h_vs_M1_24h"
)

# Contrastes principales del proyecto
lactate_contrasts <- c(
  "Lac_8h_vs_M1_8h",
  "Lac_16h_vs_M1_16h",
  "Lac_24h_vs_M1_24h"
)

# Detectar archivos de resultados
result_files <- list.files(
  results_dir,
  pattern = "_all\\.csv$",
  full.names = TRUE
)

if (length(result_files) == 0) {
  stop("No se encontraron archivos de resultados completos de DESeq2.")
}

# Generar volcano plot
make_volcano <- function(res_df, contrast_name) {
  # Preparar datos para graficar
  volcano_df <- res_df %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange)
    ) %>%
    mutate(
      padj_plot = pmax(padj, 1e-300),
      neglog10padj = -log10(padj_plot),
      regulation = factor(regulation, levels = c("NO", "DOWN", "UP"))
    )

  # Construir volcano plot
  ggplot(
    volcano_df,
    aes(
      x = log2FoldChange,
      y = neglog10padj,
      color = regulation
    )
  ) +
    geom_point(size = 0.7, alpha = 0.75) +
    scale_color_manual(values = vpcolors, drop = FALSE) +
    geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      color = "gray30",
      linetype = "dashed",
      linewidth = 0.35
    ) +
    geom_hline(
      yintercept = -log10(fdr_cutoff),
      color = "gray30",
      linetype = "dashed",
      linewidth = 0.35
    ) +
    coord_cartesian(
      xlim = c(-10, 10),
      ylim = c(0, 30)
    ) +
    labs(
      title = contrast_name,
      x = "log2 Fold Change",
      y = "-log10 FDR",
      color = "Clasificación"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "gray95", color = "gray70"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      legend.position = "right"
    )
}

# Lista para guardar plots
volcano_plots <- list()

# Generar volcano plot para cada contraste
for (file in result_files) {
  contrast_name <- basename(file)
  contrast_name <- sub("_all\\.csv$", "", contrast_name)

  # Cargar resultados
  res_df <- read_csv(file, show_col_types = FALSE)

  # Verificar columnas necesarias
  if (!all(c("log2FoldChange", "padj", "regulation") %in% colnames(res_df))) {
    stop(paste0("Faltan columnas necesarias en: ", file))
  }

  # Generar volcano plot
  p <- make_volcano(res_df, contrast_name)
  volcano_plots[[contrast_name]] <- p

  # Guardar volcano individual en PNG
  ggsave(
    filename = file.path(fig_dir, paste0("volcano_", contrast_name, ".png")),
    plot = p,
    width = 5.8,
    height = 4.5,
    dpi = 400,
    bg = "white"
  )
}
# Seleccionar volcano plots de ácido láctico
available_lactate <- lactate_contrasts[
  lactate_contrasts %in% names(volcano_plots)
]

if (length(available_lactate) == 0) {
  stop("No se encontraron volcano plots para contrastes de lactato.")
}

# Crear collage de volcano plots de ácido láctico
volcano_lactate_plots <- lapply(
  available_lactate,
  function(x) {
    volcano_plots[[x]]
  }
)

names(volcano_lactate_plots) <- available_lactate

# Unir plots y conservar una sola leyenda compartida
volcano_lactate_collage <- wrap_plots(
  volcano_lactate_plots,
  ncol = 3
) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

# Guardar collage en PNG
ggsave(
  filename = file.path(fig_dir, "volcano_lactate_collage.png"),
  plot = volcano_lactate_collage,
  width = 12,
  height = 4,
  dpi = 400,
  bg = "white"
)
