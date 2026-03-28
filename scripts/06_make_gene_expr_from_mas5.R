suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

expr_file <- "data/GSE84571_expr_mas5_log2.csv"
annot_file <- "data/GPL570_platform_table.tsv"
out_file <- "data/GSE84571_expr_gene_log2.csv"

if (!file.exists(expr_file)) stop("❌ not found: ", expr_file)
if (!file.exists(annot_file)) stop("❌ not found: ", annot_file)

message("📄 Reading expression matrix: ", expr_file)
expr <- read_csv(expr_file, show_col_types = FALSE)
# 期望第一列是 probe_id 或者类似列名
# 兼容：如果第一列不是 probe_id，则取第一列作为 probe_id
probe_col <- names(expr)[1]
names(expr)[1] <- "probe_id"

message("📄 Reading GPL570 annotation table: ", annot_file)
gpl <- read_tsv(annot_file, show_col_types = FALSE)

# 尝试找到 symbol 列
pick <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}
col_symbol <- pick(gpl, c("Gene Symbol", "GENE_SYMBOL", "Symbol", "SYMBOL"))
if (is.na(col_symbol)) stop("❌ Cannot find Gene Symbol column in GPL table")

annot <- gpl %>%
  transmute(probe_id = ID, SYMBOL = .data[[col_symbol]]) %>%
  filter(!is.na(SYMBOL), SYMBOL != "")

# 合并注释到表达矩阵
expr2 <- expr %>%
  inner_join(annot, by = "probe_id")

# 计算每个探针的平均表达（用于选代表探针）
sample_cols <- setdiff(names(expr2), c("probe_id","SYMBOL"))
expr2$AveExpr_probe <- rowMeans(as.matrix(expr2[, sample_cols]), na.rm = TRUE)

# 每个 SYMBOL 选 AveExpr 最大的 probe
best <- expr2 %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = AveExpr_probe, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(SYMBOL, all_of(sample_cols))

write_csv(best, out_file)
message("✅ Saved gene-level expression matrix: ", out_file)
message("ℹ️ Genes: ", nrow(best), "  Samples: ", length(sample_cols))
