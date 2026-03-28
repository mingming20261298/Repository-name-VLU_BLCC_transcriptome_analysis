suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(fgsea)
})

# ---- 输入 / 输出
diff_file <- "results/Table_interaction_all_annotated.csv"
out_file  <- "results/GSEA_Hallmark_interaction.csv"

stopifnot(file.exists(diff_file))

# ---- 读差异结果（用 logFC 排序）
tab <- read_csv(diff_file, show_col_types = FALSE) %>%
  filter(!is.na(SYMBOL), SYMBOL != "")

# gene-level ranking（同一基因取 |logFC| 最大的探针）
rank_df <- tab %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
  ungroup()

ranks <- rank_df$logFC
names(ranks) <- rank_df$SYMBOL
ranks <- sort(ranks, decreasing = TRUE)

# ---- Hallmark gene sets（用 msigdbr 已下载缓存）
ms <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
hallmark <- split(ms$gene_symbol, ms$gs_name)

# ---- fgsea
set.seed(1)
fg <- fgsea(
  pathways = hallmark,
  stats    = ranks,
  minSize  = 15,
  maxSize  = 500,
  nperm    = 10000
)

fg_out <- fg %>%
  arrange(padj, desc(abs(NES)))

write_csv(fg_out, out_file)
message("✅ Saved: ", out_file)
message("ℹ️ Significant (FDR<0.10): ", sum(fg_out$padj < 0.10, na.rm=TRUE))
