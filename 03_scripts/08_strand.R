# 08_test_featurecounts_strand.R
#
# Author: Alondra Márquez
# Date: 2026-05-30
#
# Description:
#   Prueba tres configuraciones de strandedness en featureCounts
#   para determinar cuál asigna mejor las lecturas alineadas con STAR.
#
# Use:
#   Rscript 03_scripts/08_test_featurecounts_strand.R \
#     04_results/06_star_pe \
#     02_data/reference/gencode.vM36.basic.annotation.gtf \
#     04_results/08_featurecounts_strand_test \
#     6

library(Rsubread)

# Argumentos de entrada
args <- commandArgs(trailingOnly = TRUE)

input_dir <- args[1]
gtf_file <- args[2]
out_dir <- args[3]
threads <- as.integer(args[4])

# Crear carpeta de salida
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Buscar archivos BAM de STAR
bam_files <- list.files(
  input_dir,
  pattern = "Aligned.sortedByCoord.out.bam$",
  recursive = TRUE,
  full.names = TRUE
)

bam_files <- sort(bam_files)

# Probar las tres configuraciones de strandedness
for (strand in c(0, 1, 2)) {
  fc <- featureCounts(
    files = bam_files,
    annot.ext = gtf_file,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    isPairedEnd = TRUE,
    requireBothEndsMapped = TRUE,
    strandSpecific = strand,
    nthreads = threads
  )

  # Guardar estadisticas de asignacion
  write.csv(
    fc$stat,
    file = file.path(
      out_dir,
      paste0("featurecounts_strand_", strand, "_stats.csv")
    ),
    row.names = FALSE
  )
}
