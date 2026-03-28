suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

gsea_file <- "results/GSEA_Hallmark_interaction.csv"
out_long  <- "results/GSEA_Hallmark_leadingEdge_long.csv"
out_wide  <- "results/GSEA_Hallmark_leadingEdge_byPathway.csv"

stopifnot(file.exists(gsea_file))

g <- read_csv(gsea_file, show_col_types = FALSE)

# fgsea 的 leadingEdge 可能是 list 列，也可能为空字符串（你这里是空）
# 如果为空，我们先提示你换 fgseaMultilevel 以输出 leadingEdge
if (!"leadingEdge" %in% names(g)) stop("❌ leadingEdge column not found.")
if (all(is.na(g$leadingEdge)) || all(g$leadingEdge == "")) {
  stop("❌ leadingEdge is empty in your GSEA output. Next step: rerun fgsea with fgseaMultilevel (no nperm) to get leadingEdge.")
}

sig <- g %>% filter(padj < 0.10) %>% arrange(padj)

le_long <- sig %>%
  select(pathway, padj, NES, leadingEdge) %>%
  unnest(leadingEdge) %>%
  rename(gene = leadingEdge)

write_csv(le_long, out_long)

le_wide <- le_long %>%
  group_by(pathway) %>%
  summarise(leadingEdge_genes = paste(unique(gene), collapse=";"), .groups="drop")

write_csv(le_wide, out_wide)

message("✅ Saved: ", out_long)
message("✅ Saved: ", out_wide)
message("ℹ️ Pathways (FDR<0.10): ", n_distinct(le_long$pathway))
message("ℹ️ Total leading-edge gene entries: ", nrow(le_long))
