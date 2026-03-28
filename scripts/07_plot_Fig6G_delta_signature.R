df <- read.csv("results/Fig6_delta_table.csv", check.names = FALSE)
df <- df[is.finite(df$delta_signature) & !is.na(df$group), ]
df$group <- factor(df$group, levels=c("signature_low","signature_high"))

# p-value
pval <- tryCatch(wilcox.test(delta_signature ~ group, data=df)$p.value, error=function(e) NA_real_)
ptxt <- ifelse(is.na(pval), "Wilcoxon p=NA", paste0("Wilcoxon p=", format(pval, digits=3)))

# violin helper
violin_one <- function(y, xcenter, width=0.35){
  y <- y[is.finite(y)]
  if(length(y) < 2) return()
  d <- density(y, na.rm=TRUE)
  d$y <- d$y / max(d$y) * width
  polygon(c(xcenter - d$y, rev(xcenter + d$y)),
          c(d$x, rev(d$x)),
          border=NA, col=rgb(0,0,0,0.12))
}

draw_one <- function(device="png"){
  if(device=="png"){
    png("figs/Fig6G_delta_signature.png", width=1600, height=1200, res=200)
  } else {
    pdf("figs/Fig6G_delta_signature.pdf", width=6.8, height=5.2)
  }

  par(mar=c(6,5,3,1))

  ylim <- range(df$delta_signature, na.rm=TRUE)
  plot(1, type="n",
       xlim=c(0.5,2.5), ylim=ylim,
       xaxt="n", xlab="",
       ylab=expression(Delta*" Signature score (week1 - week0)"),
       cex.lab=1.2)

  axis(1, at=c(1,2), labels=c("signature_low","signature_high"), las=2)

  # violins
  violin_one(df$delta_signature[df$group=="signature_low"], 1)
  violin_one(df$delta_signature[df$group=="signature_high"], 2)

  # boxplot overlay via bxp (most robust)
  bp <- boxplot(delta_signature ~ group, data=df, plot=FALSE)
  bxp(bp, add=TRUE, at=c(1,2), boxwex=0.25,
      outline=FALSE, border="black", col=rgb(0,0,0,0.20), axes=FALSE)

  # jitter points
  set.seed(1)
  xj <- ifelse(df$group=="signature_low", 1, 2) + runif(nrow(df), -0.08, 0.08)
  points(xj, df$delta_signature, pch=16, cex=1.0, col=rgb(0,0,0,0.55))

  mtext(ptxt, side=3, line=0.5, adj=0.02, cex=1.0)

  dev.off()
}

draw_one("png")
draw_one("pdf")

cat("Saved figs/Fig6G_delta_signature.png and .pdf\n")
