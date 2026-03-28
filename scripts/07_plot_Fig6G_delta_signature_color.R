df <- read.csv("results/Fig6_delta_table.csv", check.names = FALSE)
df <- df[is.finite(df$delta_signature) & !is.na(df$group), ]
df$group <- factor(df$group, levels=c("signature_low","signature_high"))

# colors (match "pink vs green" syntax)
COL_LOW  <- "#2ECC71"  # green
COL_HIGH <- "#FF4FA3"  # pink

# p-value
pval <- tryCatch(wilcox.test(delta_signature ~ group, data=df)$p.value, error=function(e) NA_real_)
ptxt <- ifelse(is.na(pval), "Wilcoxon p=NA", paste0("Wilcoxon p=", format(pval, digits=3)))

# violin polygon
violin_poly <- function(y, xcenter, fill, width=0.40){
  y <- y[is.finite(y)]
  if(length(y) < 2) return()
  d <- density(y, na.rm=TRUE)
  d$y <- d$y / max(d$y) * width
  polygon(c(xcenter - d$y, rev(xcenter + d$y)),
          c(d$x, rev(d$x)),
          border=NA, col=fill)
}

draw_panel <- function(fn, w=6.6, h=4.6, is_pdf=FALSE){
  if(is_pdf) pdf(fn, width=w, height=h) else png(fn, width=1600, height=1100, res=220)

  par(mar=c(5,4.8,2.2,1), xaxs="i", yaxs="i")
  # white background, clean axes
  plot(1, type="n", xlim=c(0.4,2.6),
       ylim=range(df$delta_signature, na.rm=TRUE),
       xaxt="n", xlab="",
       ylab=expression(Delta*" Signature score (week1 - week0)"),
       cex.lab=1.15)
  box(bty="l")

  axis(1, at=c(1,2), labels=c("risk low","risk high"), las=1, tick=FALSE)

  y_low  <- df$delta_signature[df$group=="signature_low"]
  y_high <- df$delta_signature[df$group=="signature_high"]

  # translucent fills like paper
  violin_poly(y_low,  1, fill=adjustcolor(COL_LOW,  alpha.f=0.35))
  violin_poly(y_high, 2, fill=adjustcolor(COL_HIGH, alpha.f=0.35))

  # white box overlay via bxp (robust)
  bp <- boxplot(delta_signature ~ group, data=df, plot=FALSE)
  bxp(bp, add=TRUE, at=c(1,2), boxwex=0.22,
      outline=FALSE, axes=FALSE,
      border="black", lwd=1.2,
      col="white")

  # jitter points (darker)
  set.seed(1)
  xj <- ifelse(df$group=="signature_low", 1, 2) + runif(nrow(df), -0.08, 0.08)
  points(xj, df$delta_signature, pch=16, cex=0.9, col=adjustcolor("black", alpha.f=0.55))

  # p-value on top-left (like paper)
  mtext(ptxt, side=3, line=0.2, adj=0.02, cex=1.0)

  dev.off()
}

draw_panel("figs/Fig6G_delta_signature_color.png", is_pdf=FALSE)
draw_panel("figs/Fig6G_delta_signature_color.pdf", is_pdf=TRUE)

cat("Saved:\n  figs/Fig6G_delta_signature_color.png\n  figs/Fig6G_delta_signature_color.pdf\n")
