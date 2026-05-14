# ===========================================================================
# US-003: Download BRCA Somatic Mutation data from TCGA via TCGAbiolinks
# 项目：BRCA多组学数据挖掘
# 日期：2026-05-14
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE, timeout = 1200)

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(maftools)
})

cat("\n========== US-003: Download BRCA Somatic Mutation Data ==========\n\n")

OUT_DIR <- "D:/Users/Desktop/R_Work/data/public_data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Download MAF via GDCquery_Maf ----
cat("Step 1: Downloading TCGA-BRCA MAF (Mutect2 pipeline)...\n")

# Try Mutect2 first, fall back to varscan or muse
maf_file <- file.path(OUT_DIR, "brca_mutations.maf")
rds_file <- file.path(OUT_DIR, "brca_mutations.rds")

maf_download_success <- FALSE

tryCatch({
  maf_query <- GDCquery_Maf(
    tumor     = "BRCA",
    pipelines = "mutect2",
    save.csv  = TRUE,
    directory = OUT_DIR
  )
  maf_download_success <- TRUE
  cat("  Successfully downloaded Mutect2 MAF\n")
}, error = function(e) {
  cat(sprintf("  Mutect2 failed: %s\n", e$message))
  cat("  Trying varscan2 pipeline...\n")
  tryCatch({
    maf_query <<- GDCquery_Maf(
      tumor     = "BRCA",
      pipelines = "varscan2",
      save.csv  = TRUE,
      directory = OUT_DIR
    )
    maf_download_success <<- TRUE
    cat("  Successfully downloaded varscan2 MAF\n")
  }, error = function(e2) {
    cat(sprintf("  varscan2 also failed: %s\n", e2$message))
    cat("  Trying muse pipeline...\n")
    tryCatch({
      maf_query <<- GDCquery_Maf(
        tumor     = "BRCA",
        pipelines = "muse",
        save.csv  = TRUE,
        directory = OUT_DIR
      )
      maf_download_success <<- TRUE
      cat("  Successfully downloaded muse MAF\n")
    }, error = function(e3) {
      cat(sprintf("  All pipelines failed: %s\n", e3$message))
    })
  })
})

# ---- 2. If direct MAF download failed, try GDCquery approach ----
if (!maf_download_success) {
  cat("\nStep 1b: Trying alternative GDCquery approach for mutations...\n")

  tryCatch({
    mut_query <- GDCquery(
      project       = "TCGA-BRCA",
      data.category = "Simple Nucleotide Variation",
      data.type     = "Masked Somatic Mutation",
      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
      access        = "open"
    )

    if (nrow(mut_query$results[[1]]) > 0) {
      cat(sprintf("  Found %d mutation files\n", nrow(mut_query$results[[1]])))

      GDCdownload(mut_query, method = "api", files.per.chunk = 20)

      maf_data <- GDCprepare(mut_query)
      saveRDS(maf_data, rds_file)

      cat(sprintf("  Mutation data saved: %d variants\n", nrow(maf_data)))
      maf_download_success <- TRUE
    }
  }, error = function(e) {
    cat(sprintf("  Alternative approach also failed: %s\n", e$message))
  })
}

# ---- 3. Process MAF with maftools ----
if (maf_download_success && exists("maf_query")) {
  cat("\nStep 2: Processing MAF with maftools...\n")

  maf <- read.maf(
    maf     = maf_query,
    verbose = TRUE
  )

  cat(sprintf("  MAF summary: %d samples, %d genes, %d variants\n",
              length(unique(maf@data$Tumor_Sample_Barcode)),
              length(unique(maf@data$Hugo_Symbol)),
              nrow(maf@data)))

  # Save as RDS for easy loading
  saveRDS(maf, rds_file)

  # ---- 4. Mutation summary statistics ----
  cat("\nStep 3: Computing mutation summary statistics...\n")

  # Top mutated genes
  gene_summary <- maf@data %>%
    dplyr::count(Hugo_Symbol, sort = TRUE) %>%
    dplyr::rename(N_Mutations = n)

  cat(sprintf("  Top 10 mutated genes:\n"))
  for (i in 1:min(10, nrow(gene_summary))) {
    cat(sprintf("    %s: %d mutations\n", gene_summary$Hugo_Symbol[i],
                gene_summary$N_Mutations[i]))
  }

  # Variant type distribution
  var_summary <- maf@data %>%
    dplyr::count(Variant_Classification, sort = TRUE)
  cat(sprintf("\n  Variant types:\n"))
  for (i in 1:min(8, nrow(var_summary))) {
    cat(sprintf("    %s: %d (%.1f%%)\n", var_summary$Variant_Classification[i],
                var_summary$n[i], 100 * var_summary$n[i] / sum(var_summary$n)))
  }

  # Variant type distribution table
  write.csv(var_summary, file.path(OUT_DIR, "brca_variant_types.csv"), row.names = FALSE)
  write.csv(gene_summary, file.path(OUT_DIR, "brca_mutated_genes.csv"), row.names = FALSE)

  # ---- 5. Download log ----
  log_lines <- c(
    paste("Download date:", Sys.time()),
    paste("Data type: Somatic Mutation (MAF format)"),
    paste("Pipeline: Mutect2 / Varscan2 / Muse"),
    paste("Samples:", length(unique(maf@data$Tumor_Sample_Barcode))),
    paste("Genes mutated:", length(unique(maf@data$Hugo_Symbol))),
    paste("Total variants:", nrow(maf@data)),
    paste("Top mutated gene:", gene_summary$Hugo_Symbol[1],
          sprintf("(%d mutations)", gene_summary$N_Mutations[1]))
  )
  writeLines(log_lines, file.path(OUT_DIR, "brca_mutation_download_log.txt"))

} else {
  cat("\nWARNING: MAF data could not be downloaded or processed.\n")
  cat("Will proceed with available data in integration step.\n")

  log_lines <- c(
    paste("Download date:", Sys.time()),
    "Status: FAILED - No MAF data could be retrieved",
    "Note: Mutation analysis will be skipped or use alternative source"
  )
  writeLines(log_lines, file.path(OUT_DIR, "brca_mutation_download_log.txt"))
}

cat("\n========== US-003 Complete: Mutation Download ==========\n")
