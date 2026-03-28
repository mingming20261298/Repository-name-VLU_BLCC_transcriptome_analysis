suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggrepel)
  library(pheatmap); library(RColorBrewer); library(patchwork)
})

# ---- avoid select/rename masking if AnnotationDbi loaded elsewhere ----
`%>%` <- dplyr::`%>%`

# ---- paths ----
res_f  <- "~/JTM_GSE84571/results/Table_interaction_all.csv"
meta_f <- "~/JTM_GSE84571/data/sample_sheet.csv"
f_rma  <- "~/JTM_GSE84571/data/GSE84571_expr_rma_log2.csv"
f_mas5 <- "~/JTM_GSE84571/data/GSE84571_expr_mas5_log2.csv"
expr_f <- if (file.exists(f_rma)) f_rma else if (file.exists(f_mas5)) f_mas5 else NA

stopifnot(file.exists(res_f), file.exists(meta_f))
if(is.na(expr_f)) stop("Missing expression matrix.")

# ---- helpers ----
pick <- function(df, cands){
  hit <- cands[cands %in% colnames(df)]
  if(length(hit)==0) NA_character_ else hit[1]
}
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# =========================
# 1) Interaction table
# =========================
tab <- read_csv(res_f, show_col_types=FALSE, progress=FALSE)
logfc_col <- pick(tab, c("logFC","log2FC","log2FoldChange"))
fdr_col   <- pick(tab, c("adj.P.Val","FDR","padj"))
p_col     <- pick(tab, c("P.Value","pvalue","pval"))
probe_col <- pick(tab, c("probe_id","ProbeID","ID","probe","PROBEID","...1","X"))
if(any(is.na(c(logfc_col,fdr_col,probe_col)))) stop("Missing key columns in interaction table.")

tab2 <- tab %>%
  dplyr::transmute(
    probe = as.character(.data[[probe_col]]),
    logFC = as_num(.data[[logfc_col]]),
    FDR   = as_num(.data[[fdr_col]]),
    PVAL  = if(!is.na(p_col)) as_num(.data[[p_col]]) else NA_real_
  ) %>%
  dplyr::filter(!is.na(logFC), !is.na(FDR)) %>%
  dplyr::mutate(neglogFDR = -log10(pmax(FDR, 1e-300)))

# =========================
# 2) Probe -> Gene SYMBOL (robust, no masking)
# =========================
if(!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager", repos="https://cloud.r-project.org")
}
if(!requireNamespace("AnnotationDbi", quietly=TRUE)) BiocManager::install("AnnotationDbi", update=FALSE, ask=FALSE)
if(!requireNamespace("hgu133plus2.db", quietly=TRUE)) BiocManager::install("hgu133plus2.db", update=FALSE, ask=FALSE)

suppressPackageStartupMessages({
  library(AnnotationDbi); library(hgu133plus2.db)
})

raw_map <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = unique(tab2$probe),
  columns = c("SYMBOL"),
  keytype = "PROBEID"
) %>% as.data.frame()

key_col <- "PROBEID"
if(!(key_col %in% colnames(raw_map))) key_col <- colnames(raw_map)[1]

map <- raw_map %>%
  as_tibble() %>%
  dplyr::transmute(
    probe  = as.character(.data[[key_col]]),
    gene   = as.character(SYMBOL)
  ) %>%
  dplyr::filter(!is.na(gene), gene!="") %>%
  dplyr::distinct(probe, .keep_all=TRUE)

tab2 <- tab2 %>%
  dplyr::left_join(map, by="probe")

# =========================
# Panel A: Volcano (pretty)
# =========================
FDR_TH <- 0.10
LFC_TH <- 0.50
TOP_N  <- 10

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
  dplyr::filter(sig, !is.na(gene)) %>%
  dplyr::arrange(dplyr::desc(rank_score)) %>%
  dplyr::distinct(gene, .keep_all=TRUE) %>%
  dplyr::slice_head(n=TOP_N)

pA <- ggplot(tab2, aes(logFC, neglogFDR)) +
  geom_point(aes(color=direction), size=1.35, alpha=0.85) +
  geom_vline(xintercept=c(-LFC_TH, LFC_TH), linetype="dashed", linewidth=0.35) +
  geom_hline(yintercept=-log10(FDR_TH), linetype="dashed", linewidth=0.35) +
  ggrepel::geom_text_repel(
    data=lab_df,
    aes(label=gene),
    size=3.2,
    box.padding=0.25,
    point.padding=0.2,
    min.segment.length=0,
    segment.size=0.25,
    seed=1,
    max.overlaps=Inf
  ) +
  scale_color_manual(values=c(Up="#D94841", Down="#2B8CBE", NS="grey75")) +
  theme_classic(base_size=12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(6,6,6,6)
  ) +
  labs(x="Interaction (ΔBLCC − ΔControl, log2FC)", y="-log10(FDR)")

# =========================
# 3) Expression + metadata
# =========================
expr <- read_csv(expr_f, show_col_types=FALSE, progress=FALSE)
id_col <- colnames(expr)[1]
expr[-1] <- lapply(expr[-1], as_num)
mat <- as.matrix(expr[,-1])
rownames(mat) <- as.character(expr[[id_col]])
mat <- mat[, colSums(!is.na(mat))>0, drop=FALSE]

meta <- read_csv(meta_f, show_col_types=FALSE, progress=FALSE)
sample_col  <- pick(meta, c("sample","Sample","sample_id","SampleID","gsm","GSM"))
arm_col     <- pick(meta, c("arm","Arm","group","Group","treatment","Treatment"))
time_col    <- pick(meta, c("time","Time","week","Week","timepoint","Timepoint"))
patient_col <- pick(meta, c("patient","Patient","subject","Subject","id","ID"))
if(any(is.na(c(sample_col,arm_col,time_col,patient_col)))) stop("sample_sheet missing required cols")

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
    grp = paste0(arm, "_", time2)
  )

meta_ord <- meta2 %>% dplyr::arrange(arm, time2, patient, sample)
mat_ord  <- mat[, meta_ord$sample, drop=FALSE]

# =========================
# Panel B: Heatmap (clean)
# =========================
topN <- 40
top_probes <- tab2 %>%
  dplyr::arrange(FDR, dplyr::desc(abs(logFC))) %>%
  dplyr::filter(probe %in% rownames(mat_ord)) %>%
  dplyr::slice_head(n=topN) %>%
  dplyr::pull(probe)

hm <- mat_ord[top_probes, , drop=FALSE]
hm_z <- t(scale(t(hm))); hm_z[is.na(hm_z)] <- 0

ann <- meta_ord %>%
  dplyr::transmute(arm=arm, time=time2) %>%
  as.data.frame()
rownames(ann) <- meta_ord$sample

arm_levels  <- unique(meta_ord$arm)
time_levels <- unique(meta_ord$time2)
arm_cols <- setNames(colorRampPalette(brewer.pal(8,"Set2"))(length(arm_levels)), arm_levels)
time_cols <- setNames(colorRampPalette(brewer.pal(8,"Pastel1"))(length(time_levels)), time_levels)

grp_blocks <- meta_ord %>% dplyr::mutate(block=paste0(arm,"_",time2))
gaps <- cumsum(table(grp_blocks$block))

ph <- pheatmap(
  hm_z,
  annotation_col = ann,
  annotation_colors = list(arm=arm_cols, time=time_cols),
  show_colnames = FALSE,
  show_rownames = FALSE,
  cluster_cols  = FALSE,
  cluster_rows  = TRUE,
  border_color  = NA,
  gaps_col      = gaps,
  legend        = FALSE,
  annotation_legend = TRUE
)
pB <- patchwork::wrap_elements(ph[[4]])

# =========================
# Panel C: Δ waterfall (facet by arm, cleaner)
# signature = mean(z(up genes)) - mean(z(down genes))
# =========================
sig_gene <- tab2 %>%
  dplyr::filter(FDR < FDR_TH, abs(logFC) >= LFC_TH, !is.na(gene)) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(logFC = mean(logFC, na.rm=TRUE), .groups="drop")

up_genes   <- sig_gene %>% dplyr::filter(logFC > 0) %>% dplyr::pull(gene) %>% unique()
down_genes <- sig_gene %>% dplyr::filter(logFC < 0) %>% dplyr::pull(gene) %>% unique()
if(length(up_genes) < 2 || length(down_genes) < 2) stop("Too few up/down genes; loosen thresholds.")

probe2gene <- map %>% dplyr::select(probe, gene)

long_g <- as.data.frame(mat) %>%
  tibble::rownames_to_column("probe") %>%
  dplyr::inner_join(probe2gene, by="probe") %>%
  dplyr::filter(gene %in% unique(c(up_genes, down_genes))) %>%
  tidyr::pivot_longer(-c(probe,gene), names_to="sample", values_to="expr") %>%
  dplyr::mutate(expr = as_num(expr)) %>%
  dplyr::group_by(gene, sample) %>%
  dplyr::summarize(expr_g = median(expr, na.rm=TRUE), .groups="drop") %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(z = (expr_g - mean(expr_g, na.rm=TRUE))/sd(expr_g, na.rm=TRUE)) %>%
  dplyr::ungroup()

score_df <- long_g %>%
  dplyr::mutate(set = dplyr::case_when(
    gene %in% up_genes ~ "up",
    gene %in% down_genes ~ "down",
    TRUE ~ NA_character_
  )) %>%
  dplyr::filter(!is.na(set)) %>%
  dplyr::group_by(sample, set) %>%
  dplyr::summarize(m = mean(z, na.rm=TRUE), .groups="drop") %>%
  tidyr::pivot_wider(names_from=set, values_from=m) %>%
  dplyr::mutate(signature = up - down) %>%
  dplyr::select(sample, signature)

delta_df <- meta2 %>%
  dplyr::left_join(score_df, by="sample") %>%
  dplyr::filter(!is.na(signature), time2 %in% c("w0","w1")) %>%
  dplyr::select(patient, arm, time2, signature) %>%
  tidyr::pivot_wider(names_from=time2, values_from=signature) %>%
  dplyr::filter(!is.na(w0), !is.na(w1)) %>%
  dplyr::mutate(delta = w1 - w0)

write_csv(delta_df, "~/JTM_GSE84571/figs/Fig3C_patientDelta_table.csv")

delta_df2 <- delta_df %>%
  dplyr::group_by(arm) %>%
  dplyr::arrange(dplyr::desc(delta), .by_group=TRUE) %>%
  dplyr::mutate(idx = row_number()) %>%
  dplyr::ungroup()

# symmetric y for nicer look
ymax <- max(abs(delta_df2$delta), na.rm=TRUE)
ymax <- ifelse(is.finite(ymax), ymax, 1)

pC <- ggplot(delta_df2, aes(x=idx, y=delta)) +
  geom_col(width=0.85) +
  geom_hline(yintercept=0, linewidth=0.5) +
  facet_wrap(~arm, nrow=1, scales="free_x") +
  coord_cartesian(ylim=c(-ymax, ymax)) +
  theme_classic(base_size=12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face="bold"),
    plot.margin = margin(6,6,6,6)
  ) +
  labs(x="Patients (sorted)", y="Δ signature (week1 − week0)")

# =========================
# Assemble
# =========================
fig <- (pA | pB) / pC +
  plot_layout(heights=c(1, 0.9), widths=c(1, 1.15))

ggsave("~/JTM_GSE84571/figs/Fig3_interaction_pretty.pdf", fig, width=13, height=9, device=cairo_pdf)
ggsave("~/JTM_GSE84571/figs/Fig3_interaction_pretty.png", fig, width=13, height=9, dpi=300)

message("✅ Saved:",
        "\n - ~/JTM_GSE84571/figs/Fig3_interaction_pretty.pdf",
        "\n - ~/JTM_GSE84571/figs/Fig3_interaction_pretty.png",
        "\n - ~/JTM_GSE84571/figs/Fig3C_patientDelta_table.csv")
