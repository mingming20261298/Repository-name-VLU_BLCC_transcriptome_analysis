suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
})

cat("== MCP-counter DiD started ==\n")

## expr: genes x samples (log2)
expr <- read.csv("data/GSE84571_expr_gene_log2.csv",
                 row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(expr)

meta <- read.csv("data/sample_sheet.csv", stringsAsFactors = FALSE)
stopifnot(all(c("sample_id","arm","time","patient") %in% colnames(meta)))

# match your encoding: arm is BLCC/Control; time is 0/1
meta$arm     <- factor(meta$arm, levels=c("Control","BLCC"))
meta$time    <- factor(meta$time, levels=c(0,1))
meta$patient <- factor(meta$patient)

expr_mat <- expr_mat[, meta$sample_id]

# ensure gene symbols as rownames (MCP-counter uses HUGO symbols)
genes <- rownames(expr_mat)
dup <- duplicated(genes)
if (any(dup)) {
  expr_mat <- expr_mat[!dup, , drop=FALSE]
}

# Install MCPcounter if missing
if (!requireNamespace("MCPcounter", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  BiocManager::install("MCPcounter", update=FALSE, ask=FALSE)
}

library(MCPcounter)

# MCPcounter expects HUGO symbols; our rownames are gene symbols already
cat("Running MCP-counter...\n")
mcp <- MCPcounter.estimate(expr_mat, featuresType="HUGO_symbols")  # celltypes x samples

mcp_df <- as.data.frame(t(mcp))  # samples x celltypes
mcp_df$sample_id <- rownames(mcp_df)
write.csv(mcp_df, "results/MCPcounter_scores_all_samples.csv", row.names=FALSE)

cell_types <- rownames(mcp)

out <- lapply(cell_types, function(cell) {
  df <- meta %>% dplyr::select(sample_id, arm, time, patient)
  df$score <- as.numeric(mcp[cell, df$sample_id])

  fit <- lm(score ~ arm * time + patient, data=df)
  cf <- summary(fit)$coefficients

  term <- "armBLCC:time1"
  if (!(term %in% rownames(cf))) return(data.frame(cell_type=cell, logFC=NA_real_, t=NA_real_, P.Value=NA_real_))
  data.frame(cell_type=cell,
             logFC=unname(cf[term,"Estimate"]),
             t=unname(cf[term,"t value"]),
             P.Value=unname(cf[term,"Pr(>|t|)"]))
})

res <- bind_rows(out)
res$FDR <- p.adjust(res$P.Value, method="BH")

write.csv(res, "results/MCPcounter_DiD_interaction_all.csv", row.names=FALSE)
write.csv(res[!is.na(res$FDR) & res$FDR < 0.10, ],
          "results/MCPcounter_DiD_interaction_sig_FDR0.10.csv", row.names=FALSE)

cat("== DONE ==\n")
cat("Saved:\n")
cat("  results/MCPcounter_scores_all_samples.csv\n")
cat("  results/MCPcounter_DiD_interaction_all.csv\n")
cat("  results/MCPcounter_DiD_interaction_sig_FDR0.10.csv\n")
