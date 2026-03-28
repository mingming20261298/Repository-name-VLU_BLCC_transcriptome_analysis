suppressPackageStartupMessages({
  library(xCell)
  library(dplyr)
})

cat("== xCell DiD analysis started ==\n")

## 1) load gene-level expression (log2), genes x samples
expr <- read.csv("data/GSE84571_expr_gene_log2.csv",
                 row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(expr)

## 2) load sample sheet (must contain: sample_id, arm, time, patient)
meta <- read.csv("data/sample_sheet.csv", stringsAsFactors = FALSE)
stopifnot(all(c("sample_id","arm","time","patient") %in% colnames(meta)))

meta$arm     <- factor(meta$arm, levels = c("Control","BLCC"))
meta$time <- factor(meta$time, levels = c(0,1))
meta$patient <- factor(meta$patient)

## reorder columns to match meta
expr_mat <- expr_mat[, meta$sample_id]

## 3) run xCell
cat("Running xCell...\n")
xcell_scores <- xCellAnalysis(expr_mat)   # celltypes x samples

## write per-sample scores: samples x celltypes
xcell_df <- as.data.frame(t(xcell_scores))
xcell_df$sample_id <- rownames(xcell_df)
write.csv(xcell_df, "results/xCell_scores_all_samples.csv", row.names = FALSE)

## 4) DiD model for each cell type using lm()
cell_types <- rownames(xcell_scores)

out <- lapply(cell_types, function(cell) {
  df <- meta %>% dplyr::select(sample_id, arm, time, patient)

  ## score per sample
  df$score <- as.numeric(xcell_scores[cell, df$sample_id])

  fit <- lm(score ~ arm * time + patient, data = df)
  cf  <- summary(fit)$coefficients

  # handle missing interaction term gracefully
  term <- "armBLCC:time1"
  if (!(term %in% rownames(cf))) {
    term <- "armBLCC:timeweek1"  # fallback for other codings
  }
  if (!(term %in% rownames(cf))) {
    return(data.frame(cell_type=cell, logFC=NA_real_, t=NA_real_, P.Value=NA_real_, stringsAsFactors=FALSE))
  }


  if (!("armBLCC:timeweek1" %in% rownames(cf))) {
  }

  data.frame(
    cell_type = cell,
    logFC = unname(cf[term,"Estimate"]),
    t = unname(cf[term,"t value"]),
    P.Value = unname(cf[term,"Pr(>|t|)"]),
    stringsAsFactors = FALSE
  )
})

res <- bind_rows(out)
res$FDR <- p.adjust(res$P.Value, method = "BH")

write.csv(res, "results/xCell_DiD_interaction_all.csv", row.names = FALSE)
write.csv(res[res$FDR < 0.10, ],
          "results/xCell_DiD_interaction_sig_FDR0.10.csv",
          row.names = FALSE)

cat("== DONE ==\n")
cat("Saved:\n")
cat("  results/xCell_scores_all_samples.csv\n")
cat("  results/xCell_DiD_interaction_all.csv\n")
cat("  results/xCell_DiD_interaction_sig_FDR0.10.csv\n")
