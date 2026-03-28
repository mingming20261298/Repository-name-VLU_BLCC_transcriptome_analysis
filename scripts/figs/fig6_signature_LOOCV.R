source("~/JTM_GSE84571/scripts/figs/_common.R")
DEG_FILE  <- "~/JTM_GSE84571/results/Table_interaction_all.csv"
EXPR_FILE <- "~/JTM_GSE84571/data/GSE84571_expr_rma_log2.csv"
META_FILE <- "~/JTM_GSE84571/data/sample_sheet.csv"
stopifnot(file.exists(DEG_FILE), file.exists(EXPR_FILE), file.exists(META_FILE))
deg <- read_csv(DEG_FILE, show_col_types = FALSE)
expr <- read_csv(EXPR_FILE, show_col_types = FALSE)
meta <- read_csv(META_FILE, show_col_types = FALSE)
logfc_col <- pick_col(deg, c("logFC","log2FC","log_fc","LFC"))
padj_col  <- pick_col(deg, c("adj.P.Val","FDR","padj","qval"))
id_deg    <- pick_col(deg, c("probe_id","ProbeID","ID","id"))
stopifnot(!is.na(logfc_col), !is.na(padj_col))
id_expr <- colnames(expr)[1]
mat <- as.matrix(expr[,-1]); rownames(mat) <- expr[[id_expr]]
sample_col <- pick_col(meta, c("sample_id","SampleID","gsm","GSM","sample","Sample"))
arm_col    <- pick_col(meta, c("arm","group","treatment","Arm"))
time_col   <- pick_col(meta, c("time","Time","week","Week","timepoint"))
patient_col<- pick_col(meta, c("patient","Patient","subject","Subject","id","ID"))
stopifnot(!is.na(sample_col), !is.na(arm_col), !is.na(time_col), !is.na(patient_col))
meta2 <- meta %>% mutate(sample=.data[[sample_col]], arm=.data[[arm_col]], time=.data[[time_col]], patient=.data[[patient_col]]) %>%
  filter(sample %in% colnames(mat)) %>% distinct(sample, .keep_all=TRUE)
deg2 <- deg %>% mutate(id=if(!is.na(id_deg)) as.character(.data[[id_deg]]) else as.character(.[[1]]),
                       logFC=to_num(.data[[logfc_col]]), FDR=to_num(.data[[padj_col]])) %>%
  filter(!is.na(logFC), !is.na(FDR)) %>% arrange(FDR)
sig_ids <- deg2 %>% filter(FDR < 0.10) %>% slice_head(n=50) %>% pull(id)
sig_ids <- intersect(sig_ids, rownames(mat))
stopifnot(length(sig_ids) >= 5)
w <- deg2 %>% filter(id %in% sig_ids) %>% select(id, logFC)
w_vec <- w$logFC; names(w_vec) <- w$id
sub <- mat[sig_ids, meta2$sample, drop=FALSE]
sub_z <- t(scale(t(sub)))
score <- colSums(sub_z * w_vec[sig_ids], na.rm=TRUE)
df_score <- data.frame(sample=names(score), score=as.numeric(score)) %>%
  left_join(meta2 %>% select(sample, arm, time, patient), by="sample")
dlt <- df_score %>% group_by(arm, patient) %>%
  summarise(w0=mean(score[as.character(time) %in% c("week0","0","Week0","W0")], na.rm=TRUE),
            w1=mean(score[as.character(time) %in% c("week1","1","Week1","W1")], na.rm=TRUE),
            delta=w1-w0, .groups="drop")
p_delta <- ggplot(dlt, aes(x=arm, y=delta, group=patient)) +
  geom_hline(yintercept=0, linetype=2) + geom_line(alpha=0.6) + geom_point(size=2) +
  theme_bw() + labs(title="Signature: patient-level Δscore (week1-week0)", x="", y="Δ score")
save_plot(p_delta, "Fig6_signature_patient_delta", w=6.5, h=4)
