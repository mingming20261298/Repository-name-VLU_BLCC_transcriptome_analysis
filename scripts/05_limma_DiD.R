#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
})

# -----------------------------
# Paths
# -----------------------------
proj_dir <- normalizePath("~/.")  # will be overridden below safely
proj_dir <- normalizePath("~/JTM_GSE84571")
data_dir <- file.path(proj_dir, "data")
res_dir  <- file.path(proj_dir, "results")
fig_dir  <- file.path(proj_dir, "figs")
log_dir  <- file.path(proj_dir, "logs")

dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(log_dir, "05_limma_DiD.log")
zz <- file(log_file, open="wt")
sink(zz); sink(zz, type="message")

cat("== 05_limma_DiD started ==\n")
cat("Time:", format(Sys.time()), "\n\n")

# -----------------------------
# Load ExpressionSet (MAS5)
# -----------------------------
expr_path <- file.path(data_dir, "GSE84571_expr_mas5_log2.csv")
if (!file.exists(expr_path)) stop("Missing file: ", expr_path)
expr <- read.csv(expr_path, row.names=1, check.names=FALSE)
expr <- as.matrix(expr)
eset_path <- file.path(data_dir, "GSE84571_eset_mas5.rds")

cat("Expression matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")

# -----------------------------
# Load sample sheet
# -----------------------------
ss_path <- file.path(data_dir, "sample_sheet.csv")
if (!file.exists(ss_path)) stop("Missing file: ", ss_path)

ss <- read.csv(ss_path, stringsAsFactors = FALSE, check.names = FALSE)
cat("Sample sheet rows:", nrow(ss), "cols:", ncol(ss), "\n")
cat("Sample sheet columns:", paste(colnames(ss), collapse=", "), "\n\n")

# Heuristic: find a column that matches sample names in ExpressionSet
eset_samp <- colnames(expr)

# candidate columns likely to contain sample identifiers
cand_cols <- intersect(colnames(ss), c("sample", "Sample", "sample_id", "SampleID", "geo_accession",
                                      "GEO", "gsm", "GSM", "filename", "file", "cel", "CEL",
                                      "title", "name"))
if (length(cand_cols) == 0) cand_cols <- colnames(ss)

best_col <- NA
best_hit <- -Inf
for (cc in cand_cols) {
  v <- as.character(ss[[cc]])
  # normalize: remove .CEL/.cel/.gz and path
  v2 <- basename(v)
  v2 <- sub("\\.gz$", "", v2, ignore.case=TRUE)
  v2 <- sub("\\.CEL$", "", v2, ignore.case=TRUE)

  hit <- sum(v2 %in% eset_samp) + sum(v %in% eset_samp)
  if (hit > best_hit) { best_hit <- hit; best_col <- cc }
}

if (is.na(best_col) || best_hit <= 0) {
  cat("WARNING: Could not directly match sample sheet to ExpressionSet by any column.\n")
  cat("ExpressionSet sample names example:\n")
  cat(paste(head(eset_samp, 5), collapse="\n"), "\n\n")
  cat("Sample sheet head:\n")
  print(head(ss))
  stop("Please ensure sample_sheet.csv contains a column matching ExpressionSet colnames (CEL basenames or GSM IDs).")
}

cat("Matched sample identifier column:", best_col, " (hits=", best_hit, ")\n\n")

sid <- as.character(ss[[best_col]])
sid <- basename(sid)
sid <- sub("\\.gz$", "", sid, ignore.case=TRUE)

# Keep only rows that exist in ExpressionSet
keep <- sid %in% eset_samp
if (sum(keep) == 0) stop("No sample sheet entries matched ExpressionSet sample names.")
ss2 <- ss[keep, , drop=FALSE]
sid2 <- sid[keep]

# Reorder to ExpressionSet column order
ord <- match(eset_samp, sid2)
if (any(is.na(ord))) {
  missing <- eset_samp[is.na(ord)]
  cat("WARNING: some ExpressionSet samples not found in sample sheet:\n")
  cat(paste(missing, collapse="\n"), "\n")
  stop("Sample sheet incomplete: missing some ExpressionSet samples.")
}
ss2 <- ss2[ord, , drop=FALSE]
stopifnot(identical(sid2[ord], eset_samp))

# -----------------------------
# Build design variables: treatment, time, patient
# -----------------------------
# Try to use existing columns if present; otherwise infer from sample names / description.
get_col <- function(df, keys) {
  k <- intersect(colnames(df), keys)
  if (length(k) > 0) return(df[[k[1]]])
  return(NULL)
}

treat <- get_col(ss2, c("arm","Arm","group","Group","treatment","Treatment","condition","Condition"))
time  <- get_col(ss2, c("time","Time","week","Week","timepoint","Timepoint"))
patient <- get_col(ss2, c("patient","Patient","subject","Subject","id","ID","pair","Pair"))

# Fallback inference from sample name (ExpressionSet colnames)
sname <- eset_samp

infer_treat <- function(x) {
  x0 <- tolower(x)
  if (grepl("blcc|apligraf|treat", x0)) return("BLCC")
  if (grepl("control|std|standard", x0)) return("Control")
  return(NA_character_)
}
infer_time <- function(x) {
  x0 <- tolower(x)
  if (grepl("week1|w1|1wk|day7|d7", x0)) return("Week1")
  if (grepl("week0|w0|baseline|day0|d0", x0)) return("Week0")
  return(NA_character_)
}

if (is.null(treat)) treat <- vapply(sname, infer_treat, character(1))
if (is.null(time))  time  <- vapply(sname, infer_time,  character(1))

# patient: if not provided, try to parse a trailing C1/C2... or patient-like token
if (is.null(patient)) {
  # common pattern: *_C1_* etc
  m <- regexpr("C[0-9]+", sname, ignore.case=TRUE)
  patient <- ifelse(m > 0, toupper(regmatches(sname, m)), NA_character_)
}
# If still NA, assign sequential (NOT ideal, but avoids crash; user should fix sheet)
if (any(is.na(patient))) {
  cat("WARNING: patient ID missing for some samples; please add patient column in sample_sheet.csv for proper pairing.\n")
  # create a dummy patient id per sample (will effectively remove pairing benefit)
  patient[is.na(patient)] <- paste0("P", seq_len(sum(is.na(patient))))
}

# Normalize labels
treat <- ifelse(tolower(treat) %in% c("blcc","apligraf","treatment","treat"), "BLCC", treat)
treat <- ifelse(tolower(treat) %in% c("control","std","standard"), "Control", treat)
time  <- ifelse(grepl("1", as.character(time)), "Week1", time)
time  <- ifelse(grepl("0", as.character(time)), "Week0", time)

treat <- factor(treat, levels=c("Control","BLCC"))
time  <- factor(time, levels=c("Week0","Week1"))
patient <- factor(patient)

# sanity check
cat("treatment table:\n"); print(table(treat, useNA="ifany")); cat("\n")
cat("time table:\n"); print(table(time, useNA="ifany")); cat("\n")
cat("patient table (n unique):", nlevels(patient), "\n\n")

if (any(is.na(treat)) || any(is.na(time))) {
  cat("ERROR: Could not infer treatment/time for all samples.\n")
  cat("Please ensure sample_sheet.csv has columns like arm/group and week/time.\n")
  cat("Head of mapping table:\n")
  print(data.frame(sample=sname, treat=treat, time=time, patient=patient)[1:10, ])
  stop("Missing treatment/time labels.")
}

# -----------------------------
# limma DiD model with pairing (duplicateCorrelation)
# -----------------------------
design <- model.matrix(~ treat * time)   # includes treat, time, interaction
colnames(design) <- make.names(colnames(design))
cat("Design matrix columns:", paste(colnames(design), collapse=", "), "\n\n")

# estimate correlation within patient
corfit <- duplicateCorrelation(expr, design, block=patient)
cat("Estimated within-patient correlation:", corfit$consensus, "\n\n")

fit <- lmFit(expr, design, block=patient, correlation=corfit$consensus)
fit <- eBayes(fit)

# coefficients:
# (Intercept) : Control Week0
# treatBLCC   : BLCC vs Control at Week0
# timeWeek1   : Week1 vs Week0 in Control
# treatBLCC:timeWeek1 : DiD interaction
coef_name <- grep("treatBLCC.*timeWeek1", colnames(design), value=TRUE)
if (length(coef_name) != 1) {
  cat("ERROR: cannot find interaction coefficient.\n")
  cat("Available coefficients:\n"); print(colnames(design))
  stop("Interaction term not found.")
}

cat("Interaction coefficient:", coef_name, "\n\n")

tab_int <- topTable(fit, coef=coef_name, number=Inf, sort.by="P")
tab_int$probe_id <- rownames(tab_int)

# Save full and filtered results
out_full <- file.path(res_dir, "Table_interaction_all.csv")
write.csv(tab_int, out_full, row.names=FALSE)
cat("Saved:", out_full, "\n")

# choose a main threshold
fdr_cut <- 0.05
lfc_cut <- log2(1.2)  # modest threshold for microarray; adjust as needed
tab_sig <- subset(tab_int, adj.P.Val < fdr_cut & abs(logFC) >= lfc_cut)
out_sig <- file.path(res_dir, "Table_interaction_sig_FDR0.05_LFC1.2.csv")
write.csv(tab_sig, out_sig, row.names=FALSE)
cat("Saved:", out_sig, " (n=", nrow(tab_sig), ")\n\n")

# -----------------------------
# Figures (base R: no ggplot2 dependency)
# -----------------------------
# Volcano
vol_png <- file.path(fig_dir, "Fig_DiD_volcano_interaction.png")
png(vol_png, width=1600, height=1200, res=200)
with(tab_int, {
  x <- logFC
  y <- -log10(P.Value)
  plot(x, y, pch=16, cex=0.6, main="Interaction (treatment×time) Volcano",
       xlab="log2FC (DiD: ΔBLCC − ΔControl)", ylab="-log10(P)")
  abline(v=c(-lfc_cut, lfc_cut), lty=2)
  abline(h=-log10(0.05), lty=2)
  sig <- (adj.P.Val < fdr_cut & abs(logFC) >= lfc_cut)
  points(x[sig], y[sig], pch=16, cex=0.6)
})
dev.off()
cat("Saved:", vol_png, "\n")

# Heatmap of top 50 probes by interaction
topN <- 50
top_probes <- head(tab_int$probe_id, topN)
mat <- expr[top_probes, , drop=FALSE]
# z-score by probe
mat_z <- t(scale(t(mat)))
hm_png <- file.path(fig_dir, "Fig_DiD_heatmap_top50.png")
png(hm_png, width=1800, height=1400, res=200)
par(mar=c(6,6,4,2))
heatmap(mat_z, Rowv=NA, Colv=NA, scale="none",
        col=colorRampPalette(c("navy","white","firebrick3"))(101),
        margins=c(8,10),
        main="Top 50 probes (interaction) z-score heatmap")
dev.off()
cat("Saved:", hm_png, "\n")

# Paired line plot for a representative top gene (top1 probe)
top1 <- tab_int$probe_id[1]
y <- expr[top1, ]
# compute within-patient delta for each arm for this probe
df <- data.frame(sample=sname, patient=patient, treat=treat, time=time, y=y, stringsAsFactors=FALSE)
# order by arm, then patient, then time
df <- df[order(df$treat, df$patient, df$time), ]

pair_png <- file.path(fig_dir, "Fig_DiD_paired_lines_top1probe.png")
png(pair_png, width=1600, height=1200, res=200)
par(mar=c(6,6,4,2))
plot(NA, xlim=c(0.8,2.2), ylim=range(df$y, na.rm=TRUE),
     xaxt="n", xlab="", ylab="Expression (MAS5)", main=paste0("Paired change: ", top1))
axis(1, at=c(1,2), labels=c("Week0","Week1"))

# draw per patient lines separately for each arm with slight x offsets
xpos <- ifelse(df$time=="Week0", 1, 2)
offset <- ifelse(df$treat=="Control", -0.06, 0.06)
x <- xpos + offset

# connect within patient for each arm
for (arm in levels(treat)) {
  sub <- df[df$treat==arm, ]
  pats <- unique(sub$patient)
  for (p in pats) {
    s2 <- sub[sub$patient==p, ]
    if (nrow(s2)==2) lines(x[x%in%(xpos+offset)][df$patient==p & df$treat==arm], s2$y)
  }
}
points(x, df$y, pch=16)
legend("topleft", legend=c("Control","BLCC"), bty="n")
dev.off()
cat("Saved:", pair_png, "\n\n")

cat("== 05_limma_DiD finished ==\n")
cat("Time:", format(Sys.time()), "\n")

sink(type="message"); sink()
close(zz)
