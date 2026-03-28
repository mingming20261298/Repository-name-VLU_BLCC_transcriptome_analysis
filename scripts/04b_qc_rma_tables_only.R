suppressPackageStartupMessages({
  library(data.table)
  library(stats)
})

proj <- Sys.getenv("PROJ", unset="~/JTM_GSE84571")
in_probe <- file.path(proj, "data", "GSE84571_expr_rma_log2.csv")
in_gene  <- file.path(proj, "data", "GSE84571_expr_rma_gene_log2.csv")
out_dir  <- file.path(proj, "results")
logfile  <- file.path(proj, "logs", "04b_qc_rma_tables_only.log")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(proj, "logs"), showWarnings=FALSE, recursive=TRUE)

sink(logfile, split=TRUE)
cat("== 04b_qc_rma_tables_only started ==\n")
cat("Time:", format(Sys.time()), "\n")
cat("Project:", proj, "\n\n")

if (!file.exists(in_probe)) stop("Missing: data/GSE84571_expr_rma_log2.csv (run 03b first).")
if (!file.exists(in_gene))  stop("Missing: data/GSE84571_expr_rma_gene_log2.csv (run 03b first).")

probe <- fread(in_probe, data.table=FALSE, check.names=FALSE)
gene  <- fread(in_gene,  data.table=FALSE, check.names=FALSE)

# first col is ID
probe_id_col <- colnames(probe)[1]
gene_id_col  <- colnames(gene)[1]

probe_mat <- as.matrix(probe[, -1, drop=FALSE])
rownames(probe_mat) <- probe[[probe_id_col]]
gene_mat  <- as.matrix(gene[, -1, drop=FALSE])
rownames(gene_mat)  <- gene[[gene_id_col]]

samples <- colnames(probe_mat)
cat("Probe matrix:", nrow(probe_mat), "x", ncol(probe_mat), "\n")
cat("Gene  matrix:", nrow(gene_mat),  "x", ncol(gene_mat),  "\n")

# ---------- QC summary per sample ----------
qc_summary <- data.frame(
  sample = samples,
  n_probes = nrow(probe_mat),
  mean = colMeans(probe_mat, na.rm=TRUE),
  sd   = apply(probe_mat, 2, sd, na.rm=TRUE),
  median = apply(probe_mat, 2, median, na.rm=TRUE),
  iqr    = apply(probe_mat, 2, IQR, na.rm=TRUE),
  min    = apply(probe_mat, 2, min, na.rm=TRUE),
  max    = apply(probe_mat, 2, max, na.rm=TRUE),
  na_frac = colMeans(is.na(probe_mat)),
  stringsAsFactors = FALSE
)

# ---------- Sample-sample correlation (gene-level recommended) ----------
# use gene-level for robustness and speed
cor_mat <- cor(gene_mat, use="pairwise.complete.obs", method="pearson")
cor_df <- as.data.frame(as.table(cor_mat))
colnames(cor_df) <- c("sample1","sample2","correlation")

# Average correlation per sample
avg_cor <- sapply(samples, function(s){
  mean(cor_mat[s, setdiff(samples, s)], na.rm=TRUE)
})
qc_summary$avg_sample_correlation_gene <- as.numeric(avg_cor[qc_summary$sample])

# ---------- PCA coordinates (gene-level) ----------
# center/scale genes across samples
x <- t(scale(t(gene_mat)))  # z-score per gene
x[is.na(x)] <- 0
pca <- prcomp(t(x), center=TRUE, scale.=FALSE)
pca_coords <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PC3 = if (ncol(pca$x) >= 3) pca$x[,3] else NA_real_,
  stringsAsFactors = FALSE
)

# ---------- Simple outlier flags (tables only) ----------
# criteria: low avg correlation OR extreme PC1/PC2 (|z|>3)
z_pc1 <- as.numeric(scale(pca_coords$PC1))
z_pc2 <- as.numeric(scale(pca_coords$PC2))
pc_out <- (abs(z_pc1) > 3) | (abs(z_pc2) > 3)

cor_out <- qc_summary$avg_sample_correlation_gene < (median(qc_summary$avg_sample_correlation_gene) - 0.03)

outliers <- data.frame(
  sample = qc_summary$sample,
  avg_cor_gene = qc_summary$avg_sample_correlation_gene,
  flag_low_avg_cor = as.logical(cor_out),
  flag_extreme_pca = as.logical(pc_out[match(qc_summary$sample, pca_coords$sample)]),
  flagged = as.logical(cor_out) | as.logical(pc_out[match(qc_summary$sample, pca_coords$sample)]),
  stringsAsFactors = FALSE
)

# ---------- Write outputs ----------
fwrite(qc_summary, file.path(out_dir, "QC_RMA_sample_summary.csv"))
fwrite(cor_df,     file.path(out_dir, "QC_RMA_sample_correlation_long.csv"))
fwrite(pca_coords, file.path(out_dir, "QC_RMA_PCA_coordinates.csv"))
fwrite(outliers,   file.path(out_dir, "QC_RMA_outlier_flags.csv"))

cat("\nSaved QC tables to:", out_dir, "\n")
cat(" - QC_RMA_sample_summary.csv\n - QC_RMA_sample_correlation_long.csv\n - QC_RMA_PCA_coordinates.csv\n - QC_RMA_outlier_flags.csv\n")

cat("\nQuick checks:\n")
cat("Avg correlation range:", paste(range(outliers$avg_cor_gene, na.rm=TRUE), collapse=" ~ "), "\n")
cat("Flagged samples:", paste(outliers$sample[outliers$flagged], collapse=", "), "\n")

cat("\n== DONE 04b_qc_rma_tables_only ==\n")
sink()
