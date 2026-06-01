#!/usr/bin/env Rscript

# 21_featurecounts_summary_by_condition.R
#
# Genera una tabla resumida de asignación featureCounts por condición.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

metadata_file <- "04_results/11_deseq2_dge/metadata_ordered.tsv"
fc_file <- "04_results/09_featurecounts/featurecounts_strand_2_stats.csv"

out_dir <- "04_results/09_featurecounts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

fc_wide <- read_csv(fc_file, show_col_types = FALSE)

fc_long <- fc_wide %>%
  pivot_longer(
    cols = -Status,
    names_to = "bam_file",
    values_to = "count"
  ) %>%
  mutate(
    sample_id = str_extract(bam_file, "SRR[0-9]+"),
    count = as.numeric(count)
  )

fc_by_sample <- fc_long %>%
  group_by(sample_id) %>%
  summarise(
    assigned = count[Status == "Assigned"][1],
    unassigned_no_features = ifelse(
      any(Status == "Unassigned_NoFeatures"),
      count[Status == "Unassigned_NoFeatures"][1],
      NA_real_
    ),
    total = sum(count, na.rm = TRUE),
    assigned_percent = 100 * assigned / total,
    unassigned_no_features_percent = 100 * unassigned_no_features / total,
    .groups = "drop"
  ) %>%
  left_join(metadata, by = "sample_id")

write_csv(
  fc_by_sample,
  file.path(out_dir, "featurecounts_assignment_by_sample.csv")
)

fc_summary_by_condition <- fc_by_sample %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean_assigned = mean(assigned, na.rm = TRUE),
    sd_assigned = sd(assigned, na.rm = TRUE),
    mean_assigned_percent = mean(assigned_percent, na.rm = TRUE),
    sd_assigned_percent = sd(assigned_percent, na.rm = TRUE),
    mean_unassigned_no_features_percent = mean(
      unassigned_no_features_percent,
      na.rm = TRUE
    ),
    sd_unassigned_no_features_percent = sd(
      unassigned_no_features_percent,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    `Fragmentos asignados` = paste0(
      round(mean_assigned / 1e6, 2),
      " ± ",
      round(sd_assigned / 1e6, 2),
      " M"
    ),
    `% asignación` = paste0(
      round(mean_assigned_percent, 2),
      " ± ",
      round(sd_assigned_percent, 2),
      "%"
    ),
    `% no asignado a genes` = paste0(
      round(mean_unassigned_no_features_percent, 2),
      " ± ",
      round(sd_unassigned_no_features_percent, 2),
      "%"
    )
  ) %>%
  select(
    Condición = condition,
    n,
    `Fragmentos asignados`,
    `% asignación`,
    `% no asignado a genes`
  )

write_csv(
  fc_summary_by_condition,
  file.path(out_dir, "featurecounts_summary_for_slide.csv")
)
