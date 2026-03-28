suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

setwd("~/JTM_GSE84571")
dir.create("figs", showWarnings = FALSE, recursive = TRUE)

# Install ggrepel if missing (for non-overlapping labels)
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel", repos="https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(ggrepel))

# load expression
expr <- fread("data/GSE84571_expr_mas5_log2.csv")
mat <- as.matrix(expr[, -1, with=FALSE])
rownames(mat) <- expr$probe_id

# load sample sheet
ss <- fread("data/meta/sample_sheet.csv")
ss[, short_label := paste0(patient, "_", ifelse(time=="Week0","W0","W1"))]
ss[, cel_base := basename(cel_file)]

# match columns
col_base <- basename(colnames(mat))
idx <- match(ss$cel_base, col_base)
if (any(is.na(idx))) stop("Cannot match matrix columns with sample_sheet cel_file basenames.")
mat <- mat[, idx, drop=FALSE]
colnames(mat) <- ss$short_label

group <- paste(ss$arm, ss$time, sep="_")
group_levels <- c("Control_Week0","Control_Week1","BLCC_Week0","BLCC_Week1")
group <- factor(group, levels = group_levels)

cols <- c("Control_Week0"="#1f77b4",
          "Control_Week1"="#6baed6",
          "BLCC_Week0"="#d62728",
          "BLCC_Week1"="#ff9896")

# --- Boxplot (base is fine)
plot_box <- function() {
  par(mar=c(7,4,3,1))
  boxplot(mat, outline=FALSE, las=2,
          main="MAS5 log2 expression (colored by group)",
          col=cols[as.character(group)],
          border="grey30")
  legend("topleft", legend=names(cols), fill=cols, cex=0.8, bty="n")
}

# --- PCA with ggrepel (non-overlapping labels)
pca <- prcomp(t(mat), scale. = FALSE)
dfp <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  label = colnames(mat),
  group = group
)

p_pca <- ggplot(dfp, aes(PC1, PC2, color=group, label=label)) +
  geom_point(size=3) +
  ggrepel::geom_text_repel(
    size=3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.3
  ) +
  scale_color_manual(values=cols) +
  labs(title="PCA (MAS5 log2)", x="PC1", y="PC2") +
  theme_bw() +
  theme(legend.position="left")

# Save: PDF + SVG + PNG
# 1) PDF
pdf("figs/Fig_QC_MAS5_boxplot_PCA_color.pdf", width=12, height=5)
par(mfrow=c(1,2))
plot_box()
plot.new();  # placeholder for ggplot in base device
dev.off()

# Better: export combined as ggplot for PCA, base for boxplot -> use two files:
# (keeps it simple and clean for publication)

# Boxplot files
pdf("figs/Fig_QC_MAS5_boxplot_color.pdf", width=12, height=5)
plot_box()
dev.off()

svg("figs/Fig_QC_MAS5_boxplot_color.svg", width=12, height=5)
plot_box()
dev.off()

png("figs/Fig_QC_MAS5_boxplot_color.png", width=12, height=5, units="in", res=300)
plot_box()
dev.off()

# PCA files (ggplot)
ggsave("figs/Fig_QC_MAS5_PCA_color_repel.pdf", plot=p_pca, width=12, height=5)
ggsave("figs/Fig_QC_MAS5_PCA_color_repel.svg", plot=p_pca, width=12, height=5)
ggsave("figs/Fig_QC_MAS5_PCA_color_repel.png", plot=p_pca, width=12, height=5, dpi=300)

cat("Saved:\n")
cat(" - figs/Fig_QC_MAS5_boxplot_color.(pdf/svg/png)\n")
cat(" - figs/Fig_QC_MAS5_PCA_color_repel.(pdf/svg/png)\n")
