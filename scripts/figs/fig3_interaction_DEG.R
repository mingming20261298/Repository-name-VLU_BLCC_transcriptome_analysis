source("~/JTM_GSE84571/scripts/figs/_common.R")
RES_FILE  <- "~/JTM_GSE84571/results/Table_interaction_all.csv"
META_FILE <- "~/JTM_GSE84571/data/sample_sheet.csv"
EXPR_FILE <- "~/JTM_GSE84571/data/GSE84571_expr_rma_log2.csv"
stopifnot(file.exists(RES_FILE), file.exists(META_FILE), file.exists(EXPR_FILE))
res <- read_csv(RES_FILE, show_col_types = FALSE)
logfc_col <- pick_col(res, c("logFC","log2FC","log_fc","LFC"))
padj_col  <- pick_col(res, c("adj.P.Val","FDR","padj","qval","q_value"))
stopifnot(!is.na(logfc_col), !is.na(padj_col))
res2 <- res %>% mutate(logFC=to_num(.data[[logfc_col]]), FDR=to_num(.data[[padj_col]]),
                       negLog10FDR=-log10(pmax(FDR,1e-300)))
p_vol <- ggplot(res2, aes(x=logFC, y=negLog10FDR)) + geom_point(size=1.2, alpha=0.7) +
  theme_bw() + labs(title="Interaction (arm×time): Volcano", x="log2FC (interaction)", y="-log10(FDR)")
save_plot(p_vol, "Fig3A_interaction_volcano", w=7, h=5)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
if (!requireNamespace("ComplexHeatmap", quietly=TRUE)) BiocManager::install("ComplexHeatmap", ask=FALSE, update=FALSE)
library(ComplexHeatmap)
expr <- read_csv(EXPR_FILE, show_col_types = FALSE)
meta <- read_csv(META_FILE, show_col_types = FALSE)
id_col <- colnames(expr)[1]
mat <- as.matrix(expr[,-1]); rownames(mat) <- expr[[id_col]]
sample_col <- pick_col(meta, c("sample_id","SampleID","gsm","GSM","sample","Sample"))
arm_col    <- pick_col(meta, c("arm","group","treatment","Arm"))
time_col   <- pick_col(meta, c("time","Time","week","Week","timepoint"))
patient_col<- pick_col(meta, c("patient","Patient","subject","Subject","id","ID"))
stopifnot(!is.na(sample_col), !is.na(arm_col), !is.na(time_col), !is.na(patient_col))
meta2 <- meta %>% mutate(sample=.data[[sample_col]], arm=.data[[arm_col]], time=.data[[time_col]], patient=.data[[patient_col]]) %>%
  filter(sample %in% colnames(mat)) %>% distinct(sample, .keep_all=TRUE)
top_ids <- res2 %>% arrange(FDR) %>% slice_head(n=50) %>%
  { if ("probe_id" %in% colnames(.)) .$probe_id else .[[1]] } %>% as.character()
top_ids <- intersect(top_ids, rownames(mat))
hm <- mat[top_ids, meta2$sample, drop=FALSE]
hm_z <- t(scale(t(hm)))
ha <- HeatmapAnnotation(Arm=meta2$arm, Time=meta2$time, Patient=meta2$patient)
pdf("~/JTM_GSE84571/figs/Fig3B_interaction_top50_heatmap.pdf", width=9, height=7)
Heatmap(hm_z, name="Z", top_annotation=ha, show_row_names=TRUE, show_column_names=FALSE,
        column_split=paste(meta2$arm, meta2$time, sep="_"))
dev.off()
message("Saved: Fig3B_interaction_top50_heatmap.pdf")
