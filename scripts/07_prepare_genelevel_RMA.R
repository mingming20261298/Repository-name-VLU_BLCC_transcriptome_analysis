suppressPackageStartupMessages({
  library(Biobase)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(matrixStats)
})

dir.create("derived_rma", showWarnings=FALSE, recursive=TRUE)

# 1) load RMA eset (probes x samples)
eset <- readRDS("data/GSE84571_eset_rma.rds")
expr_probe <- exprs(eset)
stopifnot(ncol(expr_probe)==18)

# 2) sample sheet from colnames
samp <- colnames(expr_probe)
arm <- ifelse(grepl("_C[0-9]+_", samp), "Control",
       ifelse(grepl("_T[0-9]+_", samp), "BLCC", NA))
time <- ifelse(grepl("Week0", samp, ignore.case=TRUE), "Week0",
        ifelse(grepl("Week1", samp, ignore.case=TRUE), "Week1", NA))
mm <- regexpr("_[CT][0-9]+_", samp)
patient <- substr(samp, mm + 1, mm + attr(mm,"match.length") - 2)

ss <- data.frame(sample=samp, arm=arm, time=time, patient=patient, stringsAsFactors=FALSE)
stopifnot(all(!is.na(ss$arm)), all(!is.na(ss$time)), all(!is.na(ss$patient)))
write.csv(ss, "derived_rma/sample_sheet_rma.csv", row.names=FALSE)

# 3) probe -> SYMBOL
probes <- rownames(expr_probe)
sym <- mapIds(hgu133plus2.db, keys=probes, keytype="PROBEID", column="SYMBOL")
sym <- as.character(sym)

keep <- !is.na(sym) & sym != ""
expr_probe2 <- expr_probe[keep, , drop=FALSE]
sym2 <- sym[keep]

# 4) aggregate probes per gene (median across probes)
# build index by symbol
ord <- order(sym2)
expr_probe2 <- expr_probe2[ord,]
sym2 <- sym2[ord]

# split-apply median (fast)
u <- unique(sym2)
expr_gene <- matrix(NA_real_, nrow=length(u), ncol=ncol(expr_probe2),
                    dimnames=list(u, colnames(expr_probe2)))

# group boundaries
starts <- c(1, which(sym2[-1] != sym2[-length(sym2)]) + 1)
ends   <- c(starts[-1]-1, length(sym2))

for(i in seq_along(u)){
  idx <- starts[i]:ends[i]
  if (length(idx)==1) {
    expr_gene[i,] <- expr_probe2[idx,]
  } else {
    expr_gene[i,] <- colMedians(expr_probe2[idx,,drop=FALSE], na.rm=TRUE)
  }
}

write.csv(expr_gene, "derived_rma/expr_gene_rma_log2.csv")
saveRDS(expr_gene, "derived_rma/expr_gene_rma_log2.rds")

cat("DONE gene-level RMA: genes=", nrow(expr_gene), " samples=", ncol(expr_gene), "\n", sep="")
