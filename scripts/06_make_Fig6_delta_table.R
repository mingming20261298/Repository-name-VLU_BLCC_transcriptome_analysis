#!/usr/bin/env Rscript

# ===============================
# Input files
# ===============================
EXPR_FILE <- "data/GSE84571_expr_mas5_log2.csv"
PHENO_FILE <- "data/pheno_fig6.csv"
SIG_GENE_FILE <- "results/signature_genes.txt"
OUT_FILE <- "results/Fig6_delta_table.csv"

# ===============================
# Load phenotype
# ===============================
pheno <- read.csv(PHENO_FILE, stringsAsFactors = FALSE, check.names = FALSE)

pheno$time <- tolower(pheno$time)
pheno$time[grepl("0", pheno$time)] <- "week0"
pheno$time[grepl("1", pheno$time)] <- "week1"

# ===============================
# Load expression matrix
# ===============================
expr_raw <- read.csv(EXPR_FILE, check.names = FALSE)
colnames(expr_raw)[1] <- "feature"
rownames(expr_raw) <- expr_raw$feature
expr_raw$feature <- NULL
expr <- as.matrix(expr_raw)

# ===============================
# Signature score (probe-level)
# ===============================
sig <- readLines(SIG_GENE_FILE)
sig <- sig[sig %in% rownames(expr)]

if(length(sig) < 5){
  stop("Too few signature probes matched expression matrix.")
}

sig_expr <- expr[sig, , drop = FALSE]

# z-score per probe
sig_z <- t(scale(t(sig_expr)))
sig_score <- colMeans(sig_z, na.rm = TRUE)

# ===============================
# Sample-level table
# ===============================
sample_tbl <- data.frame(
  sample_id = pheno$sample_id,
  patient_id = pheno$patient_id,
  arm = pheno$arm,
  time = pheno$time,
  signature = sig_score[pheno$sample_id],
  stringsAsFactors = FALSE
)

# ===============================
# Patient-level delta (week1 - week0)
# ===============================
patients <- unique(sample_tbl$patient_id)
out <- data.frame()

for(p in patients){
  tmp <- sample_tbl[sample_tbl$patient_id == p, ]
  if(!all(c("week0","week1") %in% tmp$time)) next

  for(a in unique(tmp$arm)){
    tmp2 <- tmp[tmp$arm == a, ]
    if(!all(c("week0","week1") %in% tmp2$time)) next

    v0 <- tmp2$signature[tmp2$time == "week0"][1]
    v1 <- tmp2$signature[tmp2$time == "week1"][1]

    out <- rbind(out, data.frame(
      patient_id = p,
      arm = a,
      delta_signature = v1 - v0,
      stringsAsFactors = FALSE
    ))
  }
}

# ===============================
# High / Low grouping (median split)
# ===============================
thr <- median(out$delta_signature, na.rm = TRUE)
out$group <- ifelse(out$delta_signature >= thr,
                    "signature_high", "signature_low")

# ===============================
# Output
# ===============================
write.csv(out, OUT_FILE, row.names = FALSE)
cat("Saved:", OUT_FILE, "\n")
print(out)
