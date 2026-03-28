suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2)
})

res_f  <- "results/Table_interaction_all.csv"
meta_f <- "data/sample_sheet.csv"
f_rma  <- "data/GSE84571_expr_rma_log2.csv"
f_mas5 <- "data/GSE84571_expr_mas5_log2.csv"
expr_f <- if (file.exists(f_rma)) f_rma else if (file.exists(f_mas5)) f_mas5 else NA
stopifnot(file.exists(res_f), file.exists(meta_f), !is.na(expr_f))
dir.create("figs", showWarnings=FALSE, recursive=TRUE)

pick <- function(df, cands){ hit <- cands[cands %in% colnames(df)]; if(length(hit)==0) NA_character_ else hit[1] }
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# ---- interaction ----
tab <- read_csv(res_f, show_col_types=FALSE)
logfc_col <- pick(tab, c("logFC","log2FC","log2FoldChange"))
fdr_col   <- pick(tab, c("adj.P.Val","FDR","padj"))
probe_col <- pick(tab, c("probe_id","ProbeID","PROBEID","ID","probe","...1"))
tab2 <- tab %>%
  transmute(probe=as.character(.data[[probe_col]]),
            logFC=as_num(.data[[logfc_col]]),
            FDR=as_num(.data[[fdr_col]])) %>%
  filter(!is.na(logFC), !is.na(FDR))

# ---- expression ----
expr <- read_csv(expr_f, show_col_types=FALSE)
id_col <- colnames(expr)[1]
expr[-1] <- lapply(expr[-1], as_num)
mat <- as.matrix(expr[,-1])
rownames(mat) <- as.character(expr[[id_col]])

# ---- metadata ----
meta <- read_csv(meta_f, show_col_types=FALSE)
sample_col  <- pick(meta, c("sample","Sample","sample_id","SampleID","gsm","GSM"))
arm_col     <- pick(meta, c("arm","Arm","group","Group","treatment","Treatment"))
time_col    <- pick(meta, c("time","Time","week","Week","timepoint","Timepoint"))
patient_col <- pick(meta, c("patient","Patient","subject","Subject","id","ID"))
meta2 <- meta %>%
  transmute(sample=as.character(.data[[sample_col]]),
            arm=as.character(.data[[arm_col]]),
            time=as.character(.data[[time_col]]),
            patient=as.character(.data[[patient_col]])) %>%
  filter(sample %in% colnames(mat)) %>%
  mutate(time2=tolower(time),
         time2=str_replace_all(time2,"week","w"),
         time2=ifelse(time2%in%c("0","w0","week0"),"w0",time2),
         time2=ifelse(time2%in%c("1","w1","week1"),"w1",time2))

# ---- probe-based signature: top 40 probes by |logFC| among FDR<0.10 ----
FDR_TH <- 0.10
topk <- 40
sig_probes <- tab2 %>%
  filter(FDR < FDR_TH, probe %in% rownames(mat)) %>%
  arrange(desc(abs(logFC))) %>%
  slice_head(n=topk)

if(nrow(sig_probes) < 10) stop("Too few probes for signature; loosen FDR_TH or check probe IDs match.")

up_p   <- sig_probes %>% filter(logFC > 0) %>% pull(probe) %>% unique()
down_p <- sig_probes %>% filter(logFC < 0) %>% pull(probe) %>% unique()
if(length(up_p) < 3 || length(down_p) < 3) stop("Too few up/down probes; loosen thresholds.")

# z per probe across samples
long <- tibble(probe=rep(sig_probes$probe, times=length(meta2$sample)),
               sample=rep(meta2$sample, each=nrow(sig_probes)),
               expr=as_num(as.vector(mat[sig_probes$probe, meta2$sample, drop=FALSE])))
long <- long %>%
  group_by(probe) %>%
  mutate(z=(expr-mean(expr,na.rm=TRUE))/sd(expr,na.rm=TRUE)) %>%
  ungroup()

score_df <- long %>%
  mutate(set=case_when(probe %in% up_p ~ "up",
                       probe %in% down_p ~ "down",
                       TRUE ~ NA_character_)) %>%
  filter(!is.na(set)) %>%
  group_by(sample, set) %>%
  summarize(m=mean(z,na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=set, values_from=m) %>%
  mutate(signature=up-down) %>%
  select(sample, signature)

delta_df <- meta2 %>%
  left_join(score_df, by="sample") %>%
  filter(time2 %in% c("w0","w1")) %>%
  select(patient, arm, time2, signature) %>%
  pivot_wider(names_from=time2, values_from=signature) %>%
  filter(!is.na(w0), !is.na(w1)) %>%
  mutate(delta=w1-w0)

if(nrow(delta_df)==0) stop("delta_df is empty: check pairing (patient has both w0 & w1).")

write_csv(delta_df, "figs/Fig3C_patientDelta_table.csv")

# arm color auto-map（不写死 BLCC/Control）
arms <- unique(delta_df$arm)
cols <- setNames(c("#E64B35","#4DBBD5","#7E6148","#3C5488")[seq_along(arms)], arms)

delta_df2 <- delta_df %>%
  group_by(arm) %>%
  arrange(desc(delta), .by_group=TRUE) %>%
  mutate(order=row_number()) %>%
  ungroup() %>%
  mutate(x=interaction(arm, order, drop=TRUE))

ymax <- max(abs(delta_df2$delta), na.rm=TRUE)

pC <- ggplot(delta_df2, aes(x=x, y=delta, color=arm)) +
  geom_hline(yintercept=0, linewidth=0.6) +
  geom_segment(aes(xend=x, y=0, yend=delta), linewidth=0.9) +
  geom_point(size=3.2) +
  scale_color_manual(values=cols) +
  coord_cartesian(ylim=c(-ymax, ymax)) +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.title=element_blank()) +
  labs(x="Patients (sorted within arm)",
       y=expression(Delta~"signature (week1 - week0)"))

ggsave("figs/Fig3C_only.png", pC, width=8.5, height=3.6, dpi=300)
ggsave("figs/Fig3C_only.pdf", pC, width=8.5, height=3.6, device=cairo_pdf)
message("Saved figs/Fig3C_only.(png/pdf) and figs/Fig3C_patientDelta_table.csv")
