# 09_featurecounts_star_pe.R
#
# Author: Alondra Márquez
# Date: 2026-05-30
#
# Description:
#   Ejecuta la cuantificación génica con featureCounts a partir de archivos BAM
#   paired-end generados con STAR. El script usa una anotación GTF y permite
#   definir el valor de strandSpecific para evaluar la orientación de la librería.
#
# Use:
#   Rscript 03_scripts/08_featurecounts_star_pe.R \
#     04_results/06_star_pe \
#     02_data/reference/gencode.vM36.basic.annotation.gtf \
#     04_results/08_featurecounts_star_pe \
#     2 \
#     6

library(Rsubread)
# Argumentos de entrada
args <- commandArgs(trailingOnly = TRUE)

input_dir <- args[1]
gtf_file <- args[2]
out_dir <- args[3]
strand_value <- as.integer(args[4])
threads <- as.integer(args[5])

# Crear carpeta de salida si no existe
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Buscar archivos BAM generados por STAR
bam_files <- list.files(
  input_dir,
  pattern = "Aligned.sortedByCoord.out.bam$",
  recursive = TRUE,
  full.names = TRUE
)

# Ordenar archivos para mantener el mismo orden
bam_files <- sort(bam_files)

# Extraer nombres de muestra a partir de los BAM
sample_names <- sub("_Aligned.sortedByCoord.out.bam$", "", basename(bam_files))

# Ejecutar featureCounts
fc <- featureCounts(
  files = bam_files,
  nthreads = threads,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  largestOverlap = TRUE,
  isPairedEnd = TRUE,
  requireBothEndsMapped = TRUE,
  strandSpecific = strand_value
)

# Asignar nombres de muestra a la matriz de conteos
colnames(fc$counts) <- sample_names

# Guardar matriz de conteos génicos
write.csv(
  fc$counts,
  file = file.path(out_dir, paste0("gene_counts_strand_", strand_value, ".csv"))
)

# Guardar estadisticas de asignacion
write.csv(
  fc$stat,
  file = file.path(
    out_dir,
    paste0("featurecounts_strand_", strand_value, "_stats.csv")
  ),
  row.names = FALSE
)
