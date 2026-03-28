source("~/JTM_GSE84571/scripts/figs/_common.R")
EXPR_FILE <- "~/JTM_GSE84571/data/GSE84571_expr_rma_log2.csv"
META_FILE <- "~/JTM_GSE84571/data/sample_sheet.csv"
stopifnot(file.exists(EXPR_FILE), file.exists(META_FILE))
expr <- read_csv(EXPR_FILE, show_col_types = FALSE)
meta <- read_csv(META_FILE, show_col_types = FALSE)
id_col <- colnames(expr)[1]
mat <- as.matrix(expr[,-1]); rownames(mat) <- expr[[id_col]]
sample_col <- pick_col(meta, c("sample_id","SampleID","gsm","GSM","sample","Sample"))
arm_col    <- pick_col(meta, c("arm","group","treatment","Arm"))
time_col   <- pick_col(meta, c("time","Time","week","Week","timepoint"))
stopifnot(!is.na(sample_col), !is.na(arm_col), !is.na(time_col))
keep <- intersect(colnames(mat), meta[[sample_col]])
mat <- mat[, keep, drop=FALSE]
meta2 <- meta %>% mutate(sample=.data[[sample_col]], arm=.data[[arm_col]], time=.data[[time_col]]) %>%
  filter(sample %in% keep) %>% distinct(sample, .keep_all=TRUE)
set.seed(1); probe_n <- min(5000, nrow(mat)); idx <- sample(seq_len(nrow(mat)), probe_n)
df_den <- data.frame(value=as.vector(mat[idx,]), sample=rep(colnames(mat), each=probe_n)) %>%
  left_join(meta2 %>% select(sample, arm, time), by="sample")
p_den <- ggplot(df_den, aes(x=value, group=sample)) + geom_density(alpha=0.15) +
  facet_grid(arm ~ time) + theme_bw() + labs(x="log2 expression (RMA)", y="Density", title="QC: Density (RMA)")
save_plot(p_den, "Fig2A_QC_density_RMA", w=10, h=6)
p_box <- ggplot(df_den, aes(x=sample, y=value)) + geom_boxplot(outlier.size=0.2) + theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="", y="log2 expression", title="QC: Boxplot (RMA)") +
  facet_grid(arm ~ time, scales="free_x", space="free_x")
save_plot(p_box, "Fig2B_QC_boxplot_RMA", w=12, h=6)
pc <- prcomp(t(mat), scale.=TRUE)
df_pca <- data.frame(pc$x[,1:2], sample=rownames(pc$x)) %>% left_join(meta2 %>% select(sample, arm, time), by="sample")
p_pca <- ggplot(df_pca, aes(x=PC1, y=PC2)) + geom_point(aes(shape=time), size=3) + facet_wrap(~arm) +
  theme_bw() + labs(title="QC: PCA (RMA)", x="PC1", y="PC2")
save_plot(p_pca, "Fig2C_QC_PCA", w=8, h=5)
if (requireNamespace("uwot", quietly=TRUE)) {
  um <- uwot::umap(t(mat), n_neighbors=10, min_dist=0.3, metric="cosine")
  df_um <- data.frame(UMAP1=um[,1], UMAP2=um[,2], sample=colnames(mat)) %>%
    left_join(meta2 %>% select(sample, arm, time), by="sample")
  p_um <- ggplot(df_um, aes(x=UMAP1, y=UMAP2)) + geom_point(aes(shape=time), size=3) + facet_wrap(~arm) +
    theme_bw() + labs(title="QC: UMAP (RMA)", x="UMAP1", y="UMAP2")
  save_plot(p_um, "Fig2D_QC_UMAP", w=8, h=5)
} else {
  message("uwot not installed; skipped UMAP.")
}
