#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(oligo)
})

cat("== 03b_rma_from_CEL started ==\n")
cat("Time:", format(Sys.time()), "\n")

proj_dir <- getwd()
cel_dir <- file.path(proj_dir, "data", "CEL")
log_dir <- file.path(proj_dir, "logs")
out_dir <- file.path(proj_dir, "data")

dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(log_dir, "03b_rma_from_CEL.log")
zz <- file(log_file, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")

cat("Project:", proj_dir, "\n")
cat("CEL dir:", cel_dir, "\n\n")

cel_files <- list.celfiles(cel_dir, full.names = TRUE)
cat("Found CEL files:", length(cel_files), "\n")
if (length(cel_files) == 0) stop("No CEL files found under ", cel_dir)

cat("\nReading CEL files...\n")
raw <- read.celfiles(cel_files)

cat("Running RMA normalization (oligo)...\n")
eset <- rma(raw, target = "core")

saveRDS(eset, file.path(out_dir, "GSE84571_eset_rma.rds"))
write.csv(exprs(eset), file.path(out_dir, "GSE84571_expr_rma_log2.csv"))

cat("\nSaved RMA outputs\n")
cat("== DONE ==\n")

sink(type = "output"); sink(type = "message"); close(zz)
