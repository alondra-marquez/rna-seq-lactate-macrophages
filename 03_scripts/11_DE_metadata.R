# 10_prepare_metadata_dge.R
#
# Author: Alondra Márquez
# Date: 2026-05-31
#
# Description:
#   Prepara la metadata para el análisis de expresión diferencial.
#   Usa la metadata filtrada y genera una columna condition para DESeq2.
#
# Use:
#   Rscript 03_scripts/10_prepare_metadata_dge.R

# Cargar metadata
metadata <- read.delim(
  "01_metadata/samples_metadata.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Crear ID de muestra compatible con la matriz de conteos
metadata$sample_id <- metadata$run_accession

# Identificar si la muestra tiene ácido láctico
metadata$lactate <- ifelse(
  grepl("lactic acid", metadata$treatment, ignore.case = TRUE),
  "LA",
  "no_LA"
)
# Crear condición combinando tiempo y tratamiento
metadata$condition <- ifelse(
  metadata$lactate == "LA",
  paste0("M1_Lac_", metadata$time),
  paste0("M1_", metadata$time)
)

# Renombrar la condición 0 h para usarla como referencia inicial
metadata$condition[
  metadata$time == "0h" & metadata$lactate == "no_LA"
] <- "BMDM_0h"

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

metadata$time <- factor(
  metadata$time,
  levels = c("0h", "4h", "8h", "16h", "24h")
)

metadata$lactate <- factor(
  metadata$lactate,
  levels = c("no_LA", "LA")
)

# Seleccionar columnas útiles para DGE
metadata_dge <- metadata[, c(
  "sample_id",
  "geo_accession",
  "time",
  "treatment",
  "lactate",
  "condition"
)]

# Revisar distribución final
print(table(metadata_dge$condition))

# Guardar metadata final
write.table(
  metadata_dge,
  file = "01_metadata/metadata_dge.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
