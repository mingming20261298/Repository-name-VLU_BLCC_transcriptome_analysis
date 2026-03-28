suppressPackageStartupMessages({
  library(readr)
})

# 1) 自动找 MAS5 log2 表达矩阵（优先用 csv）
csv1 <- "data/GSE84571_expr_mas5_log2.csv"
csv2 <- "data/GSE84571_expr_mas5_log2.tsv"
rds1 <- "data/GSE84571_eset_mas5.rds"

expr <- NULL

if (file.exists(csv1)) {
  message("Reading: ", csv1)
  dat <- readr::read_csv(csv1, show_col_types = FALSE)
  # 第一列一般是 probe_id
  rn <- dat[[1]]
  mat <- as.matrix(dat[,-1])
  rownames(mat) <- rn
  expr <- mat
} else if (file.exists(csv2)) {
  message("Reading: ", csv2)
  dat <- readr::read_tsv(csv2, show_col_types = FALSE)
  rn <- dat[[1]]
  mat <- as.matrix(dat[,-1])
  rownames(mat) <- rn
  expr <- mat
} else if (file.exists(rds1)) {
  message("Reading: ", rds1)
  eset <- readRDS(rds1)
  # ExpressionSet / matrix 兼容
  if ("exprs" %in% methods::slotNames(eset)) {
    expr <- Biobase::exprs(eset)
  } else if (is.matrix(eset)) {
    expr <- eset
  } else {
    stop("Unknown object in RDS: ", rds1)
  }
} else {
  stop("Cannot find MAS5 expression file: ", csv1, " or ", csv2, " or ", rds1)
}

# 2) 基本 QC 数值表（不画图）
expr <- apply(expr, 2, as.numeric)  # 确保是数值
colnames(expr) <- colnames(expr)

qc <- data.frame(
  sample = colnames(expr),
  n_probes = nrow(expr),
  n_na = colSums(is.na(expr)),
  mean = colMeans(expr, na.rm=TRUE),
  sd = apply(expr, 2, sd, na.rm=TRUE),
  median = apply(expr, 2, median, na.rm=TRUE),
  IQR = apply(expr, 2, IQR, na.rm=TRUE),
  min = apply(expr, 2, min, na.rm=TRUE),
  max = apply(expr, 2, max, na.rm=TRUE),
  stringsAsFactors = FALSE
)

dir.create("results", showWarnings = FALSE)
write.csv(qc, "results/QC_MAS5_metrics.csv", row.names = FALSE)
message("Saved: results/QC_MAS5_metrics.csv")

# 3) 额外给一个样本相关性矩阵（同样不画图）
cor_mat <- cor(expr, use="pairwise.complete.obs", method="pearson")
write.csv(cor_mat, "results/QC_MAS5_sample_correlation.csv")
message("Saved: results/QC_MAS5_sample_correlation.csv")

message("QC (no-plot) finished.")
