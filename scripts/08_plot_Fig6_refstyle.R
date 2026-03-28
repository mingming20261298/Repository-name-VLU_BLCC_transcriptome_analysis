# Fig6 ref-style (pink/green violin + white box) -> PDF then convert to PNG

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(hgu133plus2.db)
})

# -------- inputs --------
expr_raw <- read.csv("data/GSE84571_expr_mas5_log2.csv", check.names=FALSE)
colnames(expr_raw)[1] <- "probe_id"
rownames(expr_raw) <- expr_raw$probe_id
expr_raw$probe_id <- NULL
expr <- as.matrix(expr_raw)

ph <- read.csv("data/pheno_fig6.csv", stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("sample_id","patient_id","arm","time") %in% colnames(ph)))
stopifnot(all(ph$sample_id %in% colnames(expr)))

dtab <- read.csv("results/Fig6_delta_table.csv", check.names=FALSE)
dtab$group <- factor(dtab$group, levels=c("signature_low","signature_high"))

COL_LOW  <- "#2ECC71"  # green
COL_HIGH <- "#FF4FA3"  # pink

# -------- helpers --------
z_by_row <- function(m) t(scale(t(m)))
expr_z <- z_by_row(expr)

sym2probes <- function(symbols){
  symbols <- unique(symbols[nzchar(symbols)])
  m <- AnnotationDbi::select(hgu133plus2.db, keys=symbols, keytype="SYMBOL", columns=c("PROBEID","SYMBOL"))
  m <- m[!is.na(m$PROBEID), ]
  unique(m$PROBEID)
}

module_score <- function(probes, expr_z){
  probes <- intersect(probes, rownames(expr_z))
  if(length(probes) < 3) return(rep(NA_real_, ncol(expr_z)))
  colMeans(expr_z[probes, , drop=FALSE], na.rm=TRUE)
}

delta_by_patient <- function(sample_vec, ph){
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
      out <- rbind(out, data.frame(patient_id=p, arm=a, delta=sample_vec[s1]-sample_vec[s0], stringsAsFactors=FALSE))
    }
  }
  out
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

panel_violin <- function(df, ycol, main){
  df <- df[is.finite(df[[ycol]]) & !is.na(df$group), ]
  df$group <- factor(df$group, levels=c("signature_low","signature_high"))

  pval <- tryCatch(wilcox.test(df[[ycol]] ~ df$group)$p.value, error=function(e) NA_real_)
  ptxt <- ifelse(is.na(pval), "Wilcoxon p=NA", paste0("Wilcoxon p=", format(pval, digits=3)))

  plot(1, type="n", xlim=c(0.4,2.6),
       ylim=range(df[[ycol]], na.rm=TRUE),
       xaxt="n", xlab="", ylab="Δ score",
       main=main, cex.main=1.0, cex.lab=0.95)
  box(bty="l")
  axis(1, at=c(1,2), labels=c("risk low","risk high"), tick=FALSE, cex.axis=0.85)

  y_low  <- df[[ycol]][df$group=="signature_low"]
  y_high <- df[[ycol]][df$group=="signature_high"]

  # solid-like fill (avoid alpha bugs) but still soft
  violin_poly(y_low,  1, fill=COL_LOW)
  violin_poly(y_high, 2, fill=COL_HIGH)

  bp <- boxplot(df[[ycol]] ~ df$group, plot=FALSE)
  bxp(bp, add=TRUE, at=c(1,2), boxwex=0.22,
      outline=FALSE, axes=FALSE,
      border="black", lwd=1.0,
      col="white")

  set.seed(1)
  xj <- ifelse(df$group=="signature_low", 1, 2) + runif(nrow(df), -0.08, 0.08)
  points(xj, df[[ycol]], pch=16, cex=0.65, col="black")

  mtext(ptxt, side=3, line=0.05, adj=0.02, cex=0.8)
}

panel_heatmap <- function(C, main="H"){
  n <- ncol(C)
  image(1:n, 1:n, t(C[n:1, ]), axes=FALSE, xlab="", ylab="", main=main)
  box(bty="l")
  axis(1, at=1:n, labels=colnames(C), las=2, cex.axis=0.55)
  axis(2, at=1:n, labels=rev(rownames(C)), las=2, cex.axis=0.55)
  for(i in 1:n) for(j in 1:n){
    text(j, n-i+1, sprintf("%.2f", C[i,j]), cex=0.5)
  }
}

# -------- modules (可解释、对标“语法”) --------
GS_REPAIR <- c("KRT14","KRT16","KRT6A","ITGA6","ITGB4","EPCAM","TP63","JUN","FOS","DSG1")
GS_ECM    <- c("COL1A1","COL1A2","COL3A1","FN1","TNC","VCAN","MMP9","TIMP1","SPP1","LGALS3")
GS_ANGIO  <- c("KDR","FLT1","PECAM1","VWF","ENG","TEK","ANGPT1","ANGPT2","VEGFA","ESM1")
GS_INFLAM <- c("IL1B","CXCL8","CCL2","TNF","PTGS2","S100A8","S100A9","NFKBIA","IRF1","CXCL2")
GS_KERAT  <- c("IVL","LOR","FLG","KRT1","KRT10","SPRR1A","SPRR2D","TGM1","CLDN1","DSC1")

P_REPAIR <- sym2probes(GS_REPAIR)
P_ECM    <- sym2probes(GS_ECM)
P_ANGIO  <- sym2probes(GS_ANGIO)
P_INFLAM <- sym2probes(GS_INFLAM)
P_KERAT  <- sym2probes(GS_KERAT)

sc_repair <- module_score(P_REPAIR, expr_z); names(sc_repair) <- colnames(expr_z)
sc_ecm    <- module_score(P_ECM,    expr_z); names(sc_ecm)    <- colnames(expr_z)
sc_angio  <- module_score(P_ANGIO,  expr_z); names(sc_angio)  <- colnames(expr_z)
sc_inflam <- module_score(P_INFLAM, expr_z); names(sc_inflam) <- colnames(expr_z)
sc_kerat  <- module_score(P_KERAT,  expr_z); names(sc_kerat)  <- colnames(expr_z)

d_repair <- delta_by_patient(sc_repair, ph); names(d_repair)[3] <- "A_repair"
d_ecm    <- delta_by_patient(sc_ecm,    ph); names(d_ecm)[3]    <- "B_ecm"
d_angio  <- delta_by_patient(sc_angio,  ph); names(d_angio)[3]  <- "C_angio"
d_inflam <- delta_by_patient(sc_inflam, ph); names(d_inflam)[3] <- "D_inflam"
d_kerat  <- delta_by_patient(sc_kerat,  ph); names(d_kerat)[3]  <- "E_kerat"

m1 <- merge(dtab, d_repair, by=c("patient_id","arm"), all.x=TRUE)
m1 <- merge(m1, d_ecm,    by=c("patient_id","arm"), all.x=TRUE)
m1 <- merge(m1, d_angio,  by=c("patient_id","arm"), all.x=TRUE)
m1 <- merge(m1, d_inflam, by=c("patient_id","arm"), all.x=TRUE)
m1 <- merge(m1, d_kerat,  by=c("patient_id","arm"), all.x=TRUE)

m1$G_comp <- rowMeans(m1[,c("A_repair","B_ecm","C_angio","E_kerat")], na.rm=TRUE)

# correlation heatmap (H)
Hvars <- c("delta_signature","A_repair","B_ecm","C_angio","D_inflam","E_kerat")
Cmat <- suppressWarnings(cor(as.matrix(m1[,Hvars]), use="pairwise.complete.obs", method="spearman"))

# Panel I genes
GENE_I <- c("COL1A1","COL3A1","FN1","MMP9","TIMP1","KRT14","KRT16","IVL","IL1B","S100A9")
P_I <- AnnotationDbi::select(hgu133plus2.db, keys=GENE_I, keytype="SYMBOL", columns=c("PROBEID","SYMBOL"))
P_I <- P_I[!is.na(P_I$PROBEID), ]

gene_sample_scores <- list()
for(sym in unique(P_I$SYMBOL)){
  probes <- intersect(unique(P_I$PROBEID[P_I$SYMBOL==sym]), rownames(expr_z))
  if(length(probes) < 1) next
  gene_sample_scores[[sym]] <- colMeans(expr_z[probes, , drop=FALSE], na.rm=TRUE)
}
gene_delta <- data.frame(patient_id=m1$patient_id, arm=m1$arm, group=m1$group, stringsAsFactors=FALSE)
for(sym in names(gene_sample_scores)){
  dd <- delta_by_patient(gene_sample_scores[[sym]], ph)
  names(dd)[3] <- paste0("d_", sym)
  gene_delta <- merge(gene_delta, dd, by=c("patient_id","arm"), all.x=TRUE)
}

panel_I <- function(gene_delta){
  gene_delta$group <- factor(gene_delta$group, levels=c("signature_low","signature_high"))
  genes <- paste0("d_", GENE_I)
  genes <- genes[genes %in% colnames(gene_delta)]
  if(length(genes)==0){ plot.new(); text(0.5,0.5,"I: no genes"); return() }

  plot(1, type="n",
       xlim=c(0.2, 2*length(genes)+0.8),
       ylim=range(as.matrix(gene_delta[,genes]), na.rm=TRUE),
       xaxt="n", xlab="", ylab="Δ expression (z)",
       main="I  Key genes (Δ)")
  box(bty="l")

  at <- 1:(2*length(genes))
  dat <- list(); cols <- c()
  for(g in genes){
    dat[[length(dat)+1]] <- gene_delta[[g]][gene_delta$group=="signature_low"]
    dat[[length(dat)+1]] <- gene_delta[[g]][gene_delta$group=="signature_high"]
    cols <- c(cols, COL_LOW, COL_HIGH)
  }
  bp <- boxplot(dat, plot=FALSE)
  bxp(bp, add=TRUE, at=at, boxwex=0.55, outline=FALSE, axes=FALSE,
      border="black", col=cols)

  set.seed(1)
  for(k in seq_along(dat)){
    yy <- dat[[k]]
    xx <- rep(at[k], length(yy)) + runif(length(yy), -0.12, 0.12)
    points(xx, yy, pch=16, cex=0.40, col="black")
  }
  centers <- seq(1.5, 2*length(genes)-0.5, by=2)
  axis(1, at=centers, labels=sub("^d_","",genes), las=2, tick=FALSE, cex.axis=0.75)
  legend("topright", legend=c("risk low","risk high"), fill=c(COL_LOW, COL_HIGH), bty="n", cex=0.8)
}

# -------- layout mimic reference --------
lay <- matrix(c(
  1, 2, 3,
  4, 5, 0,
  6, 6, 7,
  6, 6, 8,
  9, 9, 9
), nrow=5, byrow=TRUE)

out_pdf <- "figs/Fig6_refstyle.pdf"
out_png <- "figs/Fig6_refstyle.png"

# Use cairo_pdf for stable vector output in headless
cairo_pdf(out_pdf, width=8.2, height=10.0)
par(oma=c(1.2,1.0,1.2,0.5))
layout(lay)

par(mar=c(4.0,4.2,2.0,0.5))
panel_violin(m1, "A_repair", "A  Re-epithelialization")
panel_violin(m1, "B_ecm",    "B  ECM remodeling")
panel_violin(m1, "C_angio",  "C  Angiogenesis")
panel_violin(m1, "D_inflam", "D  Inflammation")
panel_violin(m1, "E_kerat",  "E  Keratin/barrier")

plot.new()

par(mar=c(6.0,5.0,2.0,0.5))
panel_heatmap(Cmat, main="H  Module correlation")

par(mar=c(4.0,4.2,2.0,0.5))
panel_violin(m1, "delta_signature", "F  Δ Signature score")
panel_violin(m1, "G_comp",          "G  Δ Repair composite")

par(mar=c(7.0,4.2,2.0,0.5))
panel_I(gene_delta)

mtext("Fig. 6  Signature stratification (BLCC VLU) — reference-style", outer=TRUE, side=3, line=0.2, cex=1.05)
dev.off()

cat("Saved PDF:", out_pdf, "\n")
