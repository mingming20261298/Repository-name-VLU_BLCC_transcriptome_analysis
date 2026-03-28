# ============================================================
# Fig6 (full) - colored, reference-style syntax (base R)
# Inputs:
#   data/GSE84571_expr_mas5_log2.csv   (probe_id + sample columns *.CEL)
#   data/pheno_fig6.csv                (sample_id with .CEL, patient_id, arm, time week0/week1)
#   results/Fig6_delta_table.csv       (patient-level delta_signature + group)
#   results/signature_genes.txt        (probe list)
# Output:
#   figs/Fig6_full.png / figs/Fig6_full.pdf
# ============================================================

# ---------- packages (only annotation) ----------
need <- c("AnnotationDbi","hgu133plus2.db")
ok <- sapply(need, requireNamespace, quietly=TRUE)
if(!all(ok)){
  stop("Missing packages: ", paste(need[!ok], collapse=", "),
       "\nInstall in R env (conda): conda install -c bioconda bioconductor-annotationdbi bioconductor-hgu133plus2.db")
}
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(hgu133plus2.db)
})

# ---------- load data ----------
expr_raw <- read.csv("data/GSE84571_expr_mas5_log2.csv", check.names=FALSE)
colnames(expr_raw)[1] <- "probe_id"
rownames(expr_raw) <- expr_raw$probe_id
expr_raw$probe_id <- NULL
expr <- as.matrix(expr_raw)

ph <- read.csv("data/pheno_fig6.csv", stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("sample_id","patient_id","arm","time") %in% colnames(ph)))

# ensure sample_id exist in expr cols
miss <- setdiff(ph$sample_id, colnames(expr))
if(length(miss) > 0){
  stop("Some pheno sample_id not found in expression columns. Example: ", paste(head(miss,5), collapse=", "))
}

# patient-level delta signature + group
dtab <- read.csv("results/Fig6_delta_table.csv", check.names=FALSE)
dtab$group <- as.character(dtab$group)

# use same group mapping as Fig6_delta_table (signature_high/low)
dtab$group <- factor(dtab$group, levels=c("signature_low","signature_high"))

# colors (match reference: pink vs green)
COL_LOW  <- "#2ECC71"  # green
COL_HIGH <- "#FF4FA3"  # pink

# ---------- helpers ----------
z_by_row <- function(m){
  # z-score per row (probe), across samples
  t(scale(t(m)))
}

violin_poly <- function(y, xcenter, fill, width=0.40){
  y <- y[is.finite(y)]
  if(length(y) < 2) return()
  d <- density(y, na.rm=TRUE)
  d$y <- d$y / max(d$y) * width
  polygon(c(xcenter - d$y, rev(xcenter + d$y)),
          c(d$x, rev(d$x)),
          border=NA, col=fill)
}

panel_violin <- function(df, ycol, main, ylab=NULL, show_p=TRUE){
  df <- df[is.finite(df[[ycol]]) & !is.na(df$group), ]
  df$group <- factor(df$group, levels=c("signature_low","signature_high"))

  pval <- NA_real_
  if(show_p){
    pval <- tryCatch(wilcox.test(df[[ycol]] ~ df$group)$p.value, error=function(e) NA_real_)
  }
  ptxt <- ifelse(is.na(pval), "Wilcoxon p=NA", paste0("Wilcoxon p=", format(pval, digits=3)))

  ylim <- range(df[[ycol]], na.rm=TRUE)
  plot(1, type="n", xlim=c(0.4,2.6), ylim=ylim,
       xaxt="n", xlab="", ylab=ifelse(is.null(ylab), ycol, ylab),
       main=main, cex.main=1.0, cex.lab=0.95)
  box(bty="l")
  axis(1, at=c(1,2), labels=c("risk low","risk high"), las=1, tick=FALSE, cex.axis=0.9)

  y_low  <- df[[ycol]][df$group=="signature_low"]
  y_high <- df[[ycol]][df$group=="signature_high"]

  violin_poly(y_low,  1, fill=COL_LOW)
  violin_poly(y_high, 2, fill=COL_HIGH)

  bp <- boxplot(df[[ycol]] ~ df$group, plot=FALSE)
  bxp(bp, add=TRUE, at=c(1,2), boxwex=0.22,
      outline=FALSE, axes=FALSE,
      border="black", lwd=1.0,
      col="white")

  set.seed(1)
  xj <- ifelse(df$group=="signature_low", 1, 2) + runif(nrow(df), -0.08, 0.08)
  points(xj, df[[ycol]], pch=16, cex=0.7, col="black")

  if(show_p) mtext(ptxt, side=3, line=0.1, adj=0.02, cex=0.85)
}

# map gene symbols -> probes (hgu133plus2)
sym2probes <- function(symbols){
  symbols <- unique(symbols[nzchar(symbols)])
  m <- AnnotationDbi::select(hgu133plus2.db, keys=symbols, keytype="SYMBOL", columns=c("PROBEID","SYMBOL"))
  m <- m[!is.na(m$PROBEID), ]
  unique(m$PROBEID)
}

# compute sample-level module score: mean z(probe) across probe set
module_score <- function(probes, expr_z){
  probes <- intersect(probes, rownames(expr_z))
  if(length(probes) < 3) return(rep(NA_real_, ncol(expr_z)))
  colMeans(expr_z[probes, , drop=FALSE], na.rm=TRUE)
}

# compute patient-level delta (week1-week0) for a sample-level vector
delta_by_patient <- function(sample_vec, ph){
  # sample_vec named by sample_id
  out <- data.frame()
  pts <- unique(ph$patient_id)
  for(p in pts){
    tmp <- ph[ph$patient_id==p, ]
    if(!all(c("week0","week1") %in% tmp$time)) next
    for(a in unique(tmp$arm)){
      tmp2 <- tmp[tmp$arm==a, ]
      if(!all(c("week0","week1") %in% tmp2$time)) next
      s0 <- tmp2$sample_id[tmp2$time=="week0"][1]
      s1 <- tmp2$sample_id[tmp2$time=="week1"][1]
      v0 <- sample_vec[s0]; v1 <- sample_vec[s1]
      out <- rbind(out, data.frame(patient_id=p, arm=a, delta=v1-v0, stringsAsFactors=FALSE))
    }
  }
  out
}

# ---------- prepare z-scored expression ----------
expr_z <- z_by_row(expr)

# ---------- define 5 wound-repair modules (gene symbols) ----------
# (这些是“论文可解释”的模块，不依赖GSVA文件)
GS_REPAIR <- c("KRT14","KRT16","KRT6A","KRT6B","ITGA6","ITGB4","EPCAM","TP63","JUN","FOS")
GS_ECM    <- c("COL1A1","COL1A2","COL3A1","FN1","TNC","VCAN","MMP9","TIMP1","SPP1","LGALS3")
GS_ANGIO  <- c("KDR","FLT1","PECAM1","VWF","ENG","TEK","ANGPT1","ANGPT2","VEGFA","ESM1")
GS_INFLAM <- c("IL1B","CXCL8","CXCL2","CCL2","TNF","PTGS2","S100A8","S100A9","NFKBIA","IRF1")
GS_KERAT  <- c("IVL","LOR","FLG","DSG1","DSC1","KRT1","KRT10","SPRR1A","SPRR2D","TGM1")

P_REPAIR <- sym2probes(GS_REPAIR)
P_ECM    <- sym2probes(GS_ECM)
P_ANGIO  <- sym2probes(GS_ANGIO)
P_INFLAM <- sym2probes(GS_INFLAM)
P_KERAT  <- sym2probes(GS_KERAT)

# sample-level module scores
sc_repair <- module_score(P_REPAIR, expr_z); names(sc_repair) <- colnames(expr_z)
sc_ecm    <- module_score(P_ECM,    expr_z); names(sc_ecm)    <- colnames(expr_z)
sc_angio  <- module_score(P_ANGIO,  expr_z); names(sc_angio)  <- colnames(expr_z)
sc_inflam <- module_score(P_INFLAM, expr_z); names(sc_inflam) <- colnames(expr_z)
sc_kerat  <- module_score(P_KERAT,  expr_z); names(sc_kerat)  <- colnames(expr_z)

# patient-level deltas for modules
d_repair <- delta_by_patient(sc_repair, ph); names(d_repair)[3] <- "delta_repair"
d_ecm    <- delta_by_patient(sc_ecm,    ph); names(d_ecm)[3]    <- "delta_ecm"
d_angio  <- delta_by_patient(sc_angio,  ph); names(d_angio)[3]  <- "delta_angio"
d_inflam <- delta_by_patient(sc_inflam, ph); names(d_inflam)[3] <- "delta_inflam"
d_kerat  <- delta_by_patient(sc_kerat,  ph); names(d_kerat)[3]  <- "delta_kerat"

# merge into one patient table + bring group from dtab
m1 <- merge(dtab, d_repair, by=c("patient_id","arm"), all.x=TRUE)
m1 <- merge(m1, d_ecm,    by=c("patient_id","arm"), all.x=TRUE)
m1 <- merge(m1, d_angio,  by=c("patient_id","arm"), all.x=TRUE)
m1 <- merge(m1, d_inflam, by=c("patient_id","arm"), all.x=TRUE)
m1 <- merge(m1, d_kerat,  by=c("patient_id","arm"), all.x=TRUE)

# composite "repair program" (for panel G)
m1$delta_repair_comp <- rowMeans(m1[,c("delta_repair","delta_ecm","delta_angio","delta_kerat")], na.rm=TRUE)

# ---------- Panel I genes (delta expression by patient, symbol->probe, mean if multi-probe) ----------
GENE_I <- c("COL1A1","COL3A1","FN1","MMP9","TIMP1","KRT14","KRT16","IVL","IL1B","S100A9")
P_I <- AnnotationDbi::select(hgu133plus2.db, keys=GENE_I, keytype="SYMBOL", columns=c("PROBEID","SYMBOL"))
P_I <- P_I[!is.na(P_I$PROBEID), ]

# sample-level gene expression (z) aggregated by SYMBOL
gene_sample_scores <- list()
for(sym in unique(P_I$SYMBOL)){
  probes <- unique(P_I$PROBEID[P_I$SYMBOL==sym])
  probes <- intersect(probes, rownames(expr_z))
  if(length(probes) < 1) next
  gene_sample_scores[[sym]] <- colMeans(expr_z[probes, , drop=FALSE], na.rm=TRUE)
}

# patient-level delta for each gene
gene_delta <- data.frame(patient_id=m1$patient_id, arm=m1$arm, stringsAsFactors=FALSE)
for(sym in names(gene_sample_scores)){
  dd <- delta_by_patient(gene_sample_scores[[sym]], ph)
  names(dd)[3] <- paste0("d_", sym)
  gene_delta <- merge(gene_delta, dd, by=c("patient_id","arm"), all.x=TRUE)
}
gene_delta <- merge(gene_delta, m1[,c("patient_id","arm","group")], by=c("patient_id","arm"), all.x=TRUE)

# ---------- Panel H correlation heatmap (delta_signature + module deltas) ----------
Hvars <- c("delta_signature","delta_repair","delta_ecm","delta_angio","delta_inflam","delta_kerat")
Hmat <- as.matrix(m1[,Hvars])
# spearman correlation (pairwise)
C <- suppressWarnings(cor(Hmat, use="pairwise.complete.obs", method="spearman"))

panel_heatmap <- function(C, main="Correlation (Spearman)"){
  n <- ncol(C)
  # draw heatmap
  image(1:n, 1:n, t(C[n:1, ]), axes=FALSE, xlab="", ylab="", main=main)
  box(bty="l")
  axis(1, at=1:n, labels=colnames(C), las=2, cex.axis=0.6)
  axis(2, at=1:n, labels=rev(rownames(C)), las=2, cex.axis=0.6)
  # write values
  for(i in 1:n){
    for(j in 1:n){
      val <- C[i,j]
      text(j, n - i + 1, labels=sprintf("%.2f", val), cex=0.55, col="black")
    }
  }
}

# ---------- Panel I: grouped boxplots for genes ----------
panel_genes_box <- function(gene_delta, genes){
  # prepare matrices by group
  gene_delta$group <- factor(gene_delta$group, levels=c("signature_low","signature_high"))
  genes <- genes[genes %in% colnames(gene_delta)]
  if(length(genes)==0){
    plot.new(); text(0.5,0.5,"Panel I: no gene deltas available"); return()
  }

  # build data list in alternating order: low/high for each gene
  dat <- list()
  cols <- c()
  for(g in genes){
    ylow  <- gene_delta[[g]][gene_delta$group=="signature_low"]
    yhigh <- gene_delta[[g]][gene_delta$group=="signature_high"]
    dat[[length(dat)+1]] <- ylow
    dat[[length(dat)+1]] <- yhigh
    cols <- c(cols, adjustcolor(COL_LOW,0.65), adjustcolor(COL_HIGH,0.65))
  }
  at <- seq(1, 2*length(genes), by=1)
  bp <- boxplot(dat, plot=FALSE)

  plot(1, type="n",
       xlim=c(0.2, 2*length(genes)+0.8),
       ylim=range(unlist(dat), na.rm=TRUE),
       xaxt="n", xlab="", ylab="Δ expression (z, week1-week0)",
       main="Key markers (Δ)")
  box(bty="l")

  # draw boxes
  bxp(bp, add=TRUE, at=at, boxwex=0.6, outline=FALSE,
      border="black", lwd=0.8, col=cols, axes=FALSE)

  # jitter points
  set.seed(1)
  for(k in seq_along(dat)){
    yy <- dat[[k]]
    xx <- rep(at[k], length(yy)) + runif(length(yy), -0.12, 0.12)
    points(xx, yy, pch=16, cex=0.45, col=adjustcolor("black",0.45))
  }

  # x labels centered per gene
  gene_centers <- seq(1.5, 2*length(genes)-0.5, by=2)
  axis(1, at=gene_centers, labels=sub("^d_","",genes), las=2, tick=FALSE, cex.axis=0.75)

  # legend
  legend("topright", legend=c("risk low","risk high"),
         fill=c(adjustcolor(COL_LOW,0.65), adjustcolor(COL_HIGH,0.65)),
         border="black", bty="n", cex=0.8)
}

# ---------- layout to mimic reference ----------
# Layout matrix:
# Row1: A B C
# Row2: D E (blank)
# Row3: H (big left) F G (right)
# Row4: I across
lay <- matrix(c(
  1, 2, 3,
  4, 5, 0,
  6, 6, 7,
  6, 6, 8,
  9, 9, 9
), nrow=5, byrow=TRUE)

# device function
draw_fig6 <- function(file, pdf=FALSE){
  if(pdf) pdf(file, width=8.2, height=10.0) else png(file, width=2200, height=2700, res=300, type="cairo-png")
  par(oma=c(1.5,1.2,1.2,0.5))
  layout(lay, widths=c(1,1,1), heights=c(1,1,1.05,1.05,1.35))

  # A–E (module deltas)
  par(mar=c(4.2,4.2,2.0,0.6))
  panel_violin(m1, "delta_repair", "A  Re-epithelialization", ylab="Δ score")
  panel_violin(m1, "delta_ecm",    "B  ECM remodeling",      ylab="Δ score")
  panel_violin(m1, "delta_angio",  "C  Angiogenesis",        ylab="Δ score")
  panel_violin(m1, "delta_inflam", "D  Inflammation program",ylab="Δ score")
  panel_violin(m1, "delta_kerat",  "E  Keratin/barrier",     ylab="Δ score")

  # blank cell
  plot.new()

  # H (correlation heatmap) - large
  par(mar=c(6.5,5.2,2.0,0.6))
  panel_heatmap(C, main="H  Module correlation (Δ)")

  # F / G (distribution panels)
  par(mar=c(4.2,4.2,2.0,0.6))
  panel_violin(m1, "delta_signature", "F  Δ Signature score", ylab="Δ score")
  panel_violin(m1, "delta_repair_comp","G  Δ Repair composite", ylab="Δ score")

  # I (genes)
  par(mar=c(7.2,4.2,2.0,0.6))
  gene_cols <- paste0("d_", GENE_I)
  panel_genes_box(gene_delta, gene_cols)

  mtext("Fig. 6  Translational signature stratification (BLCC VLU)", outer=TRUE, side=3, line=0.2, cex=1.1)
  dev.off()
}

draw_fig6("figs/Fig6_full.png", pdf=FALSE)
draw_fig6("figs/Fig6_full.pdf", pdf=TRUE)

cat("Saved:\n  figs/Fig6_full.png\n  figs/Fig6_full.pdf\n")
