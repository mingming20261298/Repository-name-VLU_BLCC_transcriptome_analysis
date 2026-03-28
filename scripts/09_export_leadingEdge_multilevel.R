suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

in_rds   <- "results/GSEA_Hallmark_interaction_multilevel.rds"
out_long <- "results/GSEA_Hallmark_multilevel_leadingEdge_long.csv"
out_wide <- "results/GSEA_Hallmark_multilevel_leadingEdge_byPathway.csv"

stopifnot(file.exists(in_rds))

fg <- readRDS(in_rds) |> as_tibble()

# 只取显著通路
sig <- fg %>% filter(padj < 0.10) %>% arrange(padj)

# ---- 手动展开 leadingEdge（不用 tidyr::unnest）
le_long <- do.call(
  rbind,
  lapply(seq_len(nrow(sig)), function(i){
    data.frame(
      pathway = sig$pathway[i],
      pval    = sig$pval[i],
      padj    = sig$padj[i],
      NES     = sig$NES[i],
      size    = sig$size[i],
      gene    = sig$leadingEdge[[i]],
      stringsAsFactors = FALSE
    )
  })
)

write_csv(le_long, out_long)

# ---- 每条通路汇总
le_wide <- le_long %>%
  group_by(pathway) %>%
  summarise(
    padj = min(padj),
    NES  = NES[which.min(padj)],
    n_leadingEdge = n_distinct(gene),
    leadingEdge_genes = paste(unique(gene), collapse=";"),
    .groups = "drop"
  ) %>%
  arrange(padj)

write_csv(le_wide, out_wide)

message("✅ Saved: ", out_long)
message("✅ Saved: ", out_wide)
message("ℹ️ Pathways: ", n_distinct(le_long$pathway))
message("ℹ️ Leading-edge gene entries: ", nrow(le_long))
