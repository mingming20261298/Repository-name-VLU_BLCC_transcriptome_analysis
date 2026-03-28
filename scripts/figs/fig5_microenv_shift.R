source("~/JTM_GSE84571/scripts/figs/_common.R")
SCORE_FILE <- "~/JTM_GSE84571/results/MCPcounter_scores.csv"
META_FILE  <- "~/JTM_GSE84571/data/sample_sheet.csv"
stopifnot(file.exists(SCORE_FILE), file.exists(META_FILE))
sc <- read_csv(SCORE_FILE, show_col_types = FALSE)
meta <- read_csv(META_FILE, show_col_types = FALSE)
sample_col <- pick_col(meta, c("sample_id","SampleID","gsm","GSM","sample","Sample"))
arm_col    <- pick_col(meta, c("arm","group","treatment","Arm"))
time_col   <- pick_col(meta, c("time","Time","week","Week","timepoint"))
patient_col<- pick_col(meta, c("patient","Patient","subject","Subject","id","ID"))
stopifnot(!is.na(sample_col), !is.na(arm_col), !is.na(time_col), !is.na(patient_col))
meta2 <- meta %>% mutate(sample=.data[[sample_col]], arm=.data[[arm_col]], time=.data[[time_col]], patient=.data[[patient_col]]) %>%
  distinct(sample, .keep_all=TRUE)
if (all(c("sample","celltype","score") %in% colnames(sc))) {
  long <- sc %>% left_join(meta2, by="sample") %>% filter(!is.na(arm), !is.na(time))
} else {
  cell_col <- colnames(sc)[1]
  long <- sc %>% pivot_longer(-all_of(cell_col), names_to="sample", values_to="score") %>%
    rename(celltype=all_of(cell_col)) %>% left_join(meta2, by="sample") %>% filter(!is.na(arm), !is.na(time))
}
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
if (!requireNamespace("ComplexHeatmap", quietly=TRUE)) BiocManager::install("ComplexHeatmap", ask=FALSE, update=FALSE)
library(ComplexHeatmap)
mat <- long %>% select(celltype, sample, score) %>% pivot_wider(names_from=sample, values_from=score) %>% as.data.frame()
rownames(mat) <- mat$celltype; mat$celltype <- NULL; mat <- as.matrix(mat)
keep <- intersect(colnames(mat), meta2$sample)
mat <- mat[, keep, drop=FALSE]; meta3 <- meta2 %>% filter(sample %in% keep)
ha <- HeatmapAnnotation(Arm=meta3$arm, Time=meta3$time)
pdf("~/JTM_GSE84571/figs/Fig5A_microenv_heatmap.pdf", width=10, height=7)
Heatmap(mat, name="Score", top_annotation=ha, show_column_names=FALSE)
dev.off()
message("Saved: Fig5A_microenv_heatmap.pdf")
