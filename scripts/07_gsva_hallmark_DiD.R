suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(limma)
  library(GSVA)
})

expr_gene <- "data/GSE84571_expr_gene_log2.csv"
pheno_file <- "data/sample_sheet.csv"   # 你项目里通常已有；如果名字不同我再按你的实际文件改
out_file <- "results/GSVA_Hallmark_interaction.csv"

if (!file.exists(expr_gene)) stop("❌ not found: ", expr_gene)
if (!file.exists(pheno_file)) stop("❌ not found: ", pheno_file, " (please tell me your sample sheet path)")

# ---- 1) 读 gene expression
dat <- read_csv(expr_gene, show_col_types = FALSE)
genes <- dat$SYMBOL
mat <- as.matrix(dat[, -1])
rownames(mat) <- genes

# ---- 2) 读样本信息（必须包含：sample_id, arm, time, patient）
pheno <- read_csv(pheno_file, show_col_types = FALSE)

need <- c("sample_id","arm","time","patient")
miss <- setdiff(need, names(pheno))
if (length(miss) > 0) stop("❌ sample_sheet missing columns: ", paste(miss, collapse=", "))

# 对齐顺序
mat <- mat[, pheno$sample_id, drop=FALSE]

# ---- 3) Hallmark gene sets
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  stop("❌ msigdbr not installed. Install it or tell me to use an offline gene set file.")
}
ms <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
gs <- split(ms$gene_symbol, ms$gs_name)

# ---- 4) GSVA
message("🧪 Running GSVA (Hallmark)...")
gsva_es <- GSVA::gsva(mat, gs, method = "gsva", kcdf = "Gaussian", verbose = FALSE)

# ---- 5) DiD model (arm × time interaction) with patient blocking
# arm: Control/Treat, time: Week0/Week1 (以你的 sample_sheet 为准)
pheno$arm  <- factor(pheno$arm)
pheno$time <- factor(pheno$time)
pheno$patient <- factor(pheno$patient)

design <- model.matrix(~ arm * time, data = pheno)

corfit <- duplicateCorrelation(gsva_es, design, block = pheno$patient)
fit <- lmFit(gsva_es, design, block = pheno$patient, correlation = corfit$consensus)
fit <- eBayes(fit)

# 交互项系数名自动寻找：包含 arm 与 time 的冒号项
coef_name <- grep("arm.*:time|time.*:arm", colnames(design), value = TRUE)
if (length(coef_name) != 1) {
  stop("❌ Cannot uniquely identify interaction coefficient in design: ", paste(colnames(design), collapse=" | "))
}

tab <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("Pathway")

write_csv(tab, out_file)
message("✅ Saved: ", out_file)
message("ℹ️ Significant (FDR<0.10): ", sum(tab$adj.P.Val < 0.10, na.rm=TRUE))
