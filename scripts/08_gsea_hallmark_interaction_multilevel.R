suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(fgsea)
  library(msigdbr)
})

diff_file <- "results/Table_interaction_all_annotated.csv"
out_csv   <- "results/GSEA_Hallmark_interaction_multilevel.csv"
out_rds   <- "results/GSEA_Hallmark_interaction_multilevel.rds"

stopifnot(file.exists(diff_file))

tab <- read_csv(diff_file, show_col_types = FALSE) %>%
  filter(!is.na(SYMBOL), SYMBOL != "")

rank_df <- tab %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
  ungroup()

ranks <- rank_df$logFC
names(ranks) <- rank_df$SYMBOL
ranks <- sort(ranks, decreasing = TRUE)

ms <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(ms$gene_symbol, ms$gs_name)

set.seed(1)
fg <- fgseaMultilevel(
  pathways = hallmark,
  stats    = ranks,
  minSize  = 15,
  maxSize  = 500
)

saveRDS(fg, out_rds)

fg_out <- fg %>% arrange(padj, desc(abs(NES)))
write_csv(fg_out, out_csv)

message("✅ Saved CSV: ", out_csv)
message("✅ Saved RDS: ", out_rds)
message("ℹ️ Significant (FDR<0.10): ", sum(fg_out$padj < 0.10, na.rm=TRUE))
