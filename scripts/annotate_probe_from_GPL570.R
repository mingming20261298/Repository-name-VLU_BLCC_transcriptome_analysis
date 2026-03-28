suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

diff_file <- "results/Table_interaction_all.csv"
gpl_tsv   <- "data/GPL570_platform_table.tsv"
out_file  <- "results/Table_interaction_all_annotated.csv"

if (!file.exists(diff_file)) stop("❌ not found: ", diff_file)
if (!file.exists(gpl_tsv))   stop("❌ not found: ", gpl_tsv)

diff <- read_csv(diff_file, show_col_types = FALSE)
if (!"probe_id" %in% names(diff)) stop("❌ probe_id column not found in results/Table_interaction_all.csv")

message("📄 Reading local GPL570 platform table: ", gpl_tsv)
gpl_table <- read_tsv(gpl_tsv, show_col_types = FALSE)

# 兼容不同列名
pick <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}
col_symbol   <- pick(gpl_table, c("Gene Symbol", "GENE_SYMBOL", "Symbol", "SYMBOL"))
col_genename <- pick(gpl_table, c("Gene Title", "GENE_TITLE", "Gene title", "GENENAME"))
col_entrez   <- pick(gpl_table, c("ENTREZ_GENE_ID", "ENTREZID", "Entrez Gene", "ENTREZ_GENE"))

annot <- tibble(
  probe_id = gpl_table$ID,
  SYMBOL   = if (!is.na(col_symbol))   gpl_table[[col_symbol]]   else NA_character_,
  GENENAME = if (!is.na(col_genename)) gpl_table[[col_genename]] else NA_character_,
  ENTREZID = if (!is.na(col_entrez))   gpl_table[[col_entrez]]   else NA_character_
)

res_annot <- diff %>% left_join(annot, by = "probe_id")
write_csv(res_annot, out_file)

message("✅ Saved: ", out_file)
message("ℹ️ Probes with SYMBOL: ", sum(!is.na(res_annot$SYMBOL) & res_annot$SYMBOL != ""))
