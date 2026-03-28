#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggrepel)
  library(pheatmap); library(RColorBrewer); library(patchwork)
})

# --------------------
# paths
# --------------------
res_f  <- "results/Table_interaction_all.csv"
meta_f <- "data/sample_sheet.csv"
f_rma  <- "data/GSE84571_expr_rma_log2.csv"
f_mas5 <- "data/GSE84571_expr_mas5_log2.csv"
expr_f <- if (file.exists(f_rma)) f_rma else if (file.exists(f_mas5)) f_mas5 else NA

stopifnot(file.exists(res_f), file.exists(meta_f))
if (is.na(expr_f)) stop("Missing expression matrix: ", f_rma, " or ", f_mas5)

dir.create("figs", showWarnings=FALSE, recursive=TRUE)

# --------------------
# helpers
# --------------------
pick <- function(df, cands){
  hit <- cands[cands %in% colnames(df)]
  if(length(hit)==0) NA_character_ else hit[1]
}
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# ============================================================
# 1) Interaction results
# ============================================================
tab <- read_csv(res_f, show_col_types=FALSE, progress=FALSE)

logfc_col <- pick(tab, c("logFC","log2FC","log2FoldChange"))
fdr_col   <- pick(tab, c("adj.P.Val","FDR","padj"))
p_col     <- pick(tab, c("P.Value","pvalue","pval"))
probe_col <- pick(tab, c("probe_id","ProbeID","PROBEID","ID","probe","...1","X"))

if(any(is.na(c(logfc_col,fdr_col,probe_col)))) {
  stop("Result table missing needed columns. Found: ", paste(colnames(tab), collapse=", "))
}

tab2 <- tab %>%
  dplyr::transmute(
    probe = as.character(.data[[probe_col]]),
    logFC = as_num(.data[[logfc_col]]),
    FDR   = as_num(.data[[fdr_col]]),
    PVAL  = if(!is.na(p_col)) as_num(.data[[p_col]]) else NA_real_
  ) %>%
  dplyr::filter(!is.na(logFC), !is.na(FDR)) %>%
  dplyr::mutate(neglogFDR = -log10(pmax(FDR, 1e-300)))

# ============================================================
# 2) Probe -> Gene symbol (optional; if fails, volcano labels fall back to probe)
# ============================================================
gene_map <- NULL
have_anno <- requireNamespace("AnnotationDbi", quietly=TRUE) &&
             requireNamespace("hgu133plus2.db", quietly=TRUE)

if (have_anno) {
  suppressPackageStartupMessages({
    library(AnnotationDbi); library(hgu133plus2.db)
  })
  raw_map <- AnnotationDbi::select(
    hgu133plus2.db,
    keys = unique(tab2$probe),
    columns = c("SYMBOL"),
    keytype = "PROBEID"
  ) %>% as.data.frame()

  key_col <- if ("PROBEID" %in% colnames(raw_map)) "PROBEID" else colnames(raw_map)[1]
  gene_map <- raw_map %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      probe = as.character(.data[[key_col]]),
      gene  = as.character(SYMBOL)
    ) %>%
    dplyr::filter(!is.na(gene), gene!="") %>%
    dplyr::distinct(probe, .keep_all=TRUE)
} else {
  gene_map <- tibble::tibble(probe=unique(tab2$probe), gene=NA_character_)
}

tab2 <- tab2 %>%
  dplyr::left_join(gene_map, by="probe") %>%
  dplyr::mutate(label_gene = dplyr::if_else(is.na(gene) | gene=="", probe, gene))

# ============================================================
# Panel A: Volcano (Up/Down/NS)
# ============================================================
FDR_TH <- 0.10
LFC_TH <- 0.50
TOP_LAB <- 10

tab2 <- tab2 %>%
  dplyr::mutate(
    sig = (FDR < FDR_TH) & (abs(logFC) >= LFC_TH),
    direction = dplyr::case_when(
      sig & logFC > 0 ~ "Up",
      sig & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    ),
    rank_score = neglogFDR * abs(logFC)
  )

lab_df <- tab2 %>%
  dplyr::filter(sig) %>%
  dplyr::arrange(dplyr::desc(rank_score)) %>%
  dplyr::distinct(label_gene, .keep_all=TRUE) %>%
  dplyr::slice_head(n=TOP_LAB)

pA <- ggplot(tab2, aes(logFC, neglogFDR)) +
  geom_point(aes(color=direction), size=1.4, alpha=0.9) +
  geom_vline(xintercept=c(-LFC_TH, LFC_TH), linetype="dashed", linewidth=0.35) +
  geom_hline(yintercept=-log10(FDR_TH), linetype="dashed", linewidth=0.35) +
  ggrepel::geom_text_repel(
    data=lab_df,
    aes(label=label_gene),
    size=3.2,
    max.overlaps=Inf,
    box.padding=0.25,
    point.padding=0.20,
    min.segment.length=0,
    segment.size=0.25,
    seed=1
  ) +
  scale_color_manual(values=c(Up="#E64B35", Down="#4DBBD5", NS="grey75")) +
  theme_classic(base_size=12) +
  theme(legend.position="top", legend.title=element_blank()) +
  labs(x="Interaction (ΔBLCC − ΔControl, log2FC)", y="-log10(FDR)")

# ============================================================
# 3) Expression + metadata
# ============================================================
expr <- read_csv(expr_f, show_col_types=FALSE, progress=FALSE)
id_col <- colnames(expr)[1]
expr[-1] <- lapply(expr[-1], as_num)
mat <- as.matrix(expr[,-1])
rownames(mat) <- as.character(expr[[id_col]])

meta <- read_csv(meta_f, show_col_types=FALSE, progress=FALSE)
sample_col  <- pick(meta, c("sample","Sample","sample_id","SampleID","gsm","GSM"))
arm_col     <- pick(meta, c("arm","Arm","group","Group","treatment","Treatment"))
time_col    <- pick(meta, c("time","Time","week","Week","timepoint","Timepoint"))
patient_col <- pick(meta, c("patient","Patient","subject","Subject","id","ID"))

if(any(is.na(c(sample_col,arm_col,time_col,patient_col)))) {
  stop("sample_sheet.csv must contain sample/arm/time/patient. Found: ", paste(colnames(meta), collapse=", "))
}

meta2 <- meta %>%
  dplyr::transmute(
    sample  = as.character(.data[[sample_col]]),
    arm     = as.character(.data[[arm_col]]),
    time    = as.character(.data[[time_col]]),
    patient = as.character(.data[[patient_col]])
  ) %>%
  dplyr::filter(sample %in% colnames(mat)) %>%
  dplyr::mutate(
    time2 = tolower(time),
    time2 = str_replace_all(time2, "week", "w"),
    time2 = ifelse(time2 %in% c("0","w0","week0"), "w0", time2),
    time2 = ifelse(time2 %in% c("1","w1","week1"), "w1", time2),
    grp   = paste0(arm, "_", time2)
  )

# clinical ordering: Control_w0 -> Control_w1 -> BLCC_w0 -> BLCC_w1 (按字母排序可能不符合，强制用出现顺序)
arm_levels <- unique(meta2$arm)
time_levels <- c("w0","w1")
meta2$arm <- factor(meta2$arm, levels=arm_levels)
meta2$time2 <- factor(meta2$time2, levels=time_levels)

meta_ord <- meta2 %>% dplyr::arrange(arm, time2, patient, sample)
mat_ord  <- mat[, meta_ord$sample, drop=FALSE]

# ============================================================
# Panel B: Heatmap (Top 40 probes; rownames show gene if available)
# ============================================================
topN <- 40
top_probes <- tab2 %>%
  dplyr::arrange(FDR, dplyr::desc(abs(logFC))) %>%
  dplyr::filter(probe %in% rownames(mat_ord)) %>%
  dplyr::slice_head(n=topN) %>%
  dplyr::pull(probe)

if(length(top_probes) < 10) stop("Too few probes overlap with expression matrix rownames.")

hm <- mat_ord[top_probes, , drop=FALSE]
hm_z <- t(scale(t(hm))); hm_z[is.na(hm_z)] <- 0

# heatmap row labels: gene symbol if present else probe
row_labels <- tab2 %>%
  dplyr::filter(probe %in% top_probes) %>%
  dplyr::select(probe, label_gene) %>%
  dplyr::distinct(probe, .keep_all=TRUE)

lab_vec <- setNames(row_labels$label_gene, row_labels$probe)
rownames(hm_z) <- lab_vec[rownames(hm_z)]
rownames(hm_z)[is.na(rownames(hm_z))] <- top_probes[is.na(rownames(hm_z))]

ann <- meta_ord %>% dplyr::transmute(arm=as.character(arm), time=as.character(time2)) %>% as.data.frame()
rownames(ann) <- meta_ord$sample

arm_cols <- setNames(colorRampPalette(brewer.pal(8,"Set2"))(length(unique(ann$arm))), unique(ann$arm))
time_cols <- setNames(colorRampPalette(brewer.pal(8,"Pastel1"))(length(unique(ann$time))), unique(ann$time))

grp_blocks <- meta_ord %>% dplyr::mutate(block=paste0(as.character(arm), "_", as.character(time2)))
gaps <- cumsum(table(grp_blocks$block))

ph <- pheatmap(
  hm_z,
  annotation_col = ann,
  annotation_colors = list(arm=arm_cols, time=time_cols),
  show_colnames = FALSE,
  show_rownames = TRUE,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  border_color = NA,
  gaps_col = gaps,
  legend = FALSE
)
pB <- patchwork::wrap_elements(ph[[4]])

# ============================================================
# Panel C: patient-level Δ signature (ALWAYS works; probe-based signature)
# signature = mean(z(up probes)) - mean(z(down probes))
# ============================================================
# choose signature probes from interaction (FDR<0.10, top |logFC|)
sig_topk <- 40
sig_tbl <- tab2 %>%
  dplyr::filter(FDR < FDR_TH, probe %in% rownames(mat)) %>%
  dplyr::arrange(dplyr::desc(abs(logFC))) %>%
  dplyr::slice_head(n=sig_topk)

if(nrow(sig_tbl) < 10) stop("Too few significant probes for signature; loosen FDR_TH or check IDs.")

up_p   <- sig_tbl %>% dplyr::filter(logFC > 0) %>% dplyr::pull(probe) %>% unique()
down_p <- sig_tbl %>% dplyr::filter(logFC < 0) %>% dplyr::pull(probe) %>% unique()
if(length(up_p) < 3 || length(down_p) < 3) stop("Too few up/down probes; loosen thresholds.")

# z-score per probe across samples (only samples used)
use_samples <- meta2$sample
m_sig <- mat[c(sig_tbl$probe), use_samples, drop=FALSE]
# ensure numeric
m_sig <- apply(m_sig, 2, as_num)
rownames(m_sig) <- sig_tbl$probe
m_sig_z <- t(scale(t(m_sig)))
m_sig_z[is.na(m_sig_z)] <- 0

score_df <- tibble::tibble(
  sample = use_samples,
  up = colMeans(m_sig_z[up_p, , drop=FALSE]),
  down = colMeans(m_sig_z[down_p, , drop=FALSE])
) %>%
  dplyr::mutate(signature = up - down) %>%
  dplyr::select(sample, signature)

delta_df <- meta2 %>%
  dplyr::left_join(score_df, by="sample") %>%
  dplyr::filter(time2 %in% c("w0","w1")) %>%
  dplyr::select(patient, arm, time2, signature) %>%
  tidyr::pivot_wider(names_from=time2, values_from=signature) %>%
  dplyr::filter(!is.na(w0), !is.na(w1)) %>%
  dplyr::mutate(delta = w1 - w0)

if(nrow(delta_df)==0) stop("delta_df is empty (pairing failed): check patient has both w0 and w1.")

write_csv(delta_df, "figs/Fig3C_patientDelta_table.csv")

# pretty lollipop waterfall (single panel, no facet)
delta_df2 <- delta_df %>%
  dplyr::group_by(arm) %>%
  dplyr::arrange(dplyr::desc(delta), .by_group=TRUE) %>%
  dplyr::mutate(order=row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(x = interaction(arm, order, drop=TRUE))

ymax <- max(abs(delta_df2$delta), na.rm=TRUE)
ymax <- ifelse(is.finite(ymax) && ymax>0, ymax, 1)

# arm colors (auto)
arms <- unique(as.character(delta_df2$arm))
base_cols <- c("#E64B35", "#4DBBD5", "#7E6148", "#3C5488", "#00A087", "#F39B7F")
arm_col_map <- setNames(base_cols[seq_along(arms)], arms)

pC <- ggplot(delta_df2, aes(x=x, y=delta, color=arm)) +
  geom_hline(yintercept=0, linewidth=0.6, color="black") +
  geom_segment(aes(xend=x, y=0, yend=delta), linewidth=0.95) +
  geom_point(size=3.3) +
  scale_color_manual(values=arm_col_map) +
  coord_cartesian(ylim=c(-ymax, ymax)) +
  theme_classic(base_size=12) +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.position="top",
    legend.title=element_blank()
  ) +
  labs(x="Patients (sorted within arm)", y=expression(Delta~"signature (week1 - week0)"))

# ============================================================
# Assemble Fig3 (A|B on top; C emphasized below)
# ============================================================
fig <- (pA | pB) / pC +
  plot_layout(heights=c(1, 1.15), widths=c(1, 1.15))

ggsave("figs/Fig3_ONE.pdf", fig, width=13, height=9, device=cairo_pdf)
ggsave("figs/Fig3_ONE.png", fig, width=13, height=9, dpi=300)

message("✅ Saved: figs/Fig3_ONE.(png/pdf) and figs/Fig3C_patientDelta_table.csv")
