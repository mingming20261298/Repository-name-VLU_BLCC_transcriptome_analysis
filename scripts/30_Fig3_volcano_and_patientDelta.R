suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(ggplot2); library(ggrepel)
})

# -------------------------
# I/O
# -------------------------
res_f  <- "results/Table_interaction_all.csv"
meta_f <- "data/sample_sheet.csv"
f_rma  <- "data/GSE84571_expr_rma_log2.csv"
f_mas5 <- "data/GSE84571_expr_mas5_log2.csv"
expr_f <- if (file.exists(f_rma)) f_rma else if (file.exists(f_mas5)) f_mas5 else NA

out_dir <- "results/fig3"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

stopifnot(file.exists(res_f), file.exists(meta_f))
if (is.na(expr_f)) stop("No expression matrix found: ", f_rma, " or ", f_mas5)

# -------------------------
# helper: pick column name
# -------------------------
pick_col <- function(df, candidates){
  hit <- candidates[candidates %in% colnames(df)]
  if(length(hit)==0) return(NA_character_)
  hit[1]
}

# -------------------------
# 1) Load interaction results
# -------------------------
x <- read_csv(res_f, show_col_types = FALSE)

# Ensure probe id column exists
if (!("probe_id" %in% colnames(x))) {
  # sometimes exported as first unnamed column
  if (any(grepl("^\\.{3}1$", colnames(x))) || any(grepl("^X$", colnames(x)))) {
    c1 <- if ("...1" %in% colnames(x)) "...1" else "X"
    x <- x %>% rename(probe_id = all_of(c1))
  } else {
    stop("Cannot find probe_id column in ", res_f)
  }
}

# Numeric coercion (safe)
num_cols <- intersect(c("logFC","AveExpr","t","P.Value","adj.P.Val","B"), colnames(x))
for(cc in num_cols) x[[cc]] <- suppressWarnings(as.numeric(x[[cc]]))

# -------------------------
# 2) Probe -> Gene SYMBOL mapping (GPL570 / hgu133plus2)
# -------------------------
symbol_map <- NULL
suppressWarnings({
  ok <- requireNamespace("AnnotationDbi", quietly=TRUE) &&
        requireNamespace("hgu133plus2.db", quietly=TRUE)
})

if (ok) {
  suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(hgu133plus2.db)
  })
  symbol_map <- AnnotationDbi::select(
    hgu133plus2.db,
    keys     = unique(as.character(x$probe_id)),
    columns  = c("SYMBOL"),
    keytype  = "PROBEID"
  ) %>%
    as_tibble() %>%
    rename(probe_id = PROBEID) %>%
    mutate(
      probe_id = as.character(probe_id),
      SYMBOL   = as.character(SYMBOL)
    ) %>%
    filter(!is.na(SYMBOL), SYMBOL != "")
} else {
  message("⚠️ Missing AnnotationDbi/hgu133plus2.db. Volcano will fall back to probe_id labels.")
  symbol_map <- tibble(probe_id = as.character(unique(x$probe_id)),
                       SYMBOL   = NA_character_)
}

x2 <- x %>%
  mutate(probe_id = as.character(probe_id)) %>%
  left_join(symbol_map, by="probe_id") %>%
  mutate(label_gene = ifelse(is.na(SYMBOL) | SYMBOL=="", probe_id, SYMBOL))

# -------------------------
# 3) Fig3A: Volcano (Up/Down coloring + label top gene symbols)
#    thresholds: FDR < 0.10, |logFC| >= 0.5  (你可改)
# -------------------------
FDR_TH <- 0.10
LFC_TH <- 0.50
TOP_N  <- 12   # 仅标注 top N（按 |logFC|，在显著基因里）

x2 <- x2 %>%
  mutate(
    neglog10p = -log10(pmax(P.Value, 1e-300)),
    sig = ifelse(!is.na(`adj.P.Val`) & `adj.P.Val` < FDR_TH & abs(logFC) >= LFC_TH, TRUE, FALSE),
    direction = case_when(
      sig & logFC > 0 ~ "Up",
      sig & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

lab_df <- x2 %>%
  filter(sig) %>%
  arrange(desc(abs(logFC))) %>%
  distinct(label_gene, .keep_all=TRUE) %>%
  slice_head(n = TOP_N)

p_vol <- ggplot(x2, aes(x=logFC, y=neglog10p)) +
  geom_point(aes(color=direction), size=1.2, alpha=0.85) +
  geom_vline(xintercept=c(-LFC_TH, LFC_TH), linetype="dashed", linewidth=0.4) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth=0.4) +
  ggrepel::geom_text_repel(
    data=lab_df,
    aes(label=label_gene),
    size=3.2,
    max.overlaps=Inf,
    box.padding=0.35,
    point.padding=0.25,
    min.segment.length=0
  ) +
  scale_color_manual(values=c(Up="#D94841", Down="#2B8CBE", NS="grey70")) +
  labs(x="log2 fold-change (DiD interaction)", y="-log10(P-value)") +
  theme_classic(base_size=12) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    plot.title = element_blank()
  )

ggsave(file.path(out_dir, "Fig3A_volcano.png"), p_vol, width=6.4, height=4.8, dpi=300)
ggsave(file.path(out_dir, "Fig3A_volcano.pdf"), p_vol, width=6.4, height=4.8)

# -------------------------
# 4) Fig3C: patient-level Δ waterfall (week1 - week0)
#    我这里用“交互项显著基因”构建一个简单 signature：
#      score = mean(z(up genes)) - mean(z(down genes))
# -------------------------
meta <- read_csv(meta_f, show_col_types = FALSE)

sample_col  <- pick_col(meta, c("sample","Sample","gsm","GSM","geo_accession","accession"))
arm_col     <- pick_col(meta, c("arm","Arm","group","Group","treatment","Treatment"))
time_col    <- pick_col(meta, c("time","Time","week","Week","timepoint","Timepoint"))
patient_col <- pick_col(meta, c("patient","Patient","subject","Subject","id","ID"))

if(any(is.na(c(sample_col,arm_col,time_col,patient_col)))) {
  stop("sample_sheet.csv must contain columns for sample/arm/time/patient. Found: ",
       paste(colnames(meta), collapse=", "))
}

meta2 <- meta %>%
  transmute(
    sample  = as.character(.data[[sample_col]]),
    arm     = as.character(.data[[arm_col]]),
    time    = as.character(.data[[time_col]]),
    patient = as.character(.data[[patient_col]])
  ) %>%
  mutate(
    time2 = tolower(time),
    time2 = str_replace_all(time2, "week", "w"),
    time2 = ifelse(time2 %in% c("0","w0","week0"), "w0", time2),
    time2 = ifelse(time2 %in% c("1","w1","week1"), "w1", time2)
  )

# expression: probes x samples, first col = probe id
mat <- read_csv(expr_f, show_col_types = FALSE)
pid_col <- colnames(mat)[1]
mat <- mat %>% rename(probe_id = all_of(pid_col))
mat$probe_id <- as.character(mat$probe_id)

# keep samples present
keep_samps <- intersect(meta2$sample, colnames(mat))
if(length(keep_samps) < 6) stop("Too few matched samples between sample_sheet and expression matrix.")
meta2 <- meta2 %>% filter(sample %in% keep_samps)
mat2 <- mat %>% select(probe_id, all_of(keep_samps))

# signature genes from interaction table
sig_tbl <- x2 %>%
  filter(sig) %>%
  mutate(gene = ifelse(is.na(SYMBOL) | SYMBOL=="", NA_character_, SYMBOL)) %>%
  filter(!is.na(gene)) %>%
  group_by(gene) %>%
  summarize(logFC = mean(logFC, na.rm=TRUE), .groups="drop")

# if mapping missing, fallback to top probes directly
if(nrow(sig_tbl) < 6) {
  message("⚠️ Too few mapped significant genes; fallback to probe-based signature (top 40 probes by |logFC|).")
  top_probes <- x2 %>%
    filter(sig) %>%
    arrange(desc(abs(logFC))) %>%
    slice_head(n=40) %>%
    pull(probe_id) %>%
    unique()
  sig_up   <- top_probes[1:floor(length(top_probes)/2)]
  sig_down <- setdiff(top_probes, sig_up)

  long <- mat2 %>%
    filter(probe_id %in% c(sig_up, sig_down)) %>%
    pivot_longer(-probe_id, names_to="sample", values_to="expr") %>%
    mutate(expr = as.numeric(expr))
  # z-score across samples within probe
  long <- long %>%
    group_by(probe_id) %>%
    mutate(z = (expr - mean(expr, na.rm=TRUE))/sd(expr, na.rm=TRUE)) %>%
    ungroup()

  score_df <- long %>%
    mutate(set = ifelse(probe_id %in% sig_up, "up", "down")) %>%
    group_by(sample, set) %>%
    summarize(m = mean(z, na.rm=TRUE), .groups="drop") %>%
    pivot_wider(names_from=set, values_from=m) %>%
    mutate(signature = up - down) %>%
    select(sample, signature)

} else {
  # gene-level collapse: map probes->gene then median
  probe2gene <- symbol_map %>%
    filter(!is.na(SYMBOL), SYMBOL!="") %>%
    distinct(probe_id, SYMBOL)

  long <- mat2 %>%
    inner_join(probe2gene, by="probe_id") %>%
    filter(SYMBOL %in% sig_tbl$gene) %>%
    pivot_longer(-c(probe_id, SYMBOL), names_to="sample", values_to="expr") %>%
    mutate(expr = as.numeric(expr)) %>%
    group_by(SYMBOL, sample) %>%
    summarize(expr_g = median(expr, na.rm=TRUE), .groups="drop")

  # z-score per gene across samples
  long <- long %>%
    group_by(SYMBOL) %>%
    mutate(z = (expr_g - mean(expr_g, na.rm=TRUE))/sd(expr_g, na.rm=TRUE)) %>%
    ungroup()

  up_genes   <- sig_tbl %>% filter(logFC > 0) %>% pull(gene) %>% unique()
  down_genes <- sig_tbl %>% filter(logFC < 0) %>% pull(gene) %>% unique()

  # require at least a few genes each side
  if(length(up_genes) < 2 || length(down_genes) < 2) {
    message("⚠️ Up/Down genes too few; using all sig genes as 'up' set (mean z).")
    score_df <- long %>%
      group_by(sample) %>%
      summarize(signature = mean(z, na.rm=TRUE), .groups="drop")
  } else {
    score_df <- long %>%
      mutate(set = case_when(
        SYMBOL %in% up_genes ~ "up",
        SYMBOL %in% down_genes ~ "down",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(set)) %>%
      group_by(sample, set) %>%
      summarize(m = mean(z, na.rm=TRUE), .groups="drop") %>%
      pivot_wider(names_from=set, values_from=m) %>%
      mutate(signature = up - down) %>%
      select(sample, signature)
  }
}

# patient-level delta
delta_df <- meta2 %>%
  left_join(score_df, by="sample") %>%
  filter(!is.na(signature), time2 %in% c("w0","w1")) %>%
  select(patient, arm, time2, signature) %>%
  pivot_wider(names_from=time2, values_from=signature) %>%
  filter(!is.na(w0), !is.na(w1)) %>%
  mutate(delta = w1 - w0) %>%
  arrange(desc(delta))

write_csv(delta_df, file.path(out_dir, "Fig3C_patientDelta_table.csv"))

# waterfall plot, sorted within arm (更像临床试验)
delta_df2 <- delta_df %>%
  group_by(arm) %>%
  arrange(desc(delta), .by_group=TRUE) %>%
  mutate(rank_in_arm = row_number()) %>%
  ungroup() %>%
  mutate(x = paste0(arm, "_", rank_in_arm))

p_wf <- ggplot(delta_df2, aes(x=reorder(x, delta), y=delta, fill=arm)) +
  geom_col(width=0.85) +
  geom_hline(yintercept=0, linewidth=0.4) +
  labs(x="Patients (sorted within arm)", y="Δ signature score (week1 - week0)") +
  theme_classic(base_size=12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    plot.title = element_blank()
  )

ggsave(file.path(out_dir, "Fig3C_patientDelta_waterfall.png"), p_wf, width=7.2, height=4.8, dpi=300)
ggsave(file.path(out_dir, "Fig3C_patientDelta_waterfall.pdf"), p_wf, width=7.2, height=4.8)

message("✅ Saved:",
        "\n - ", file.path(out_dir, "Fig3A_volcano.png"),
        "\n - ", file.path(out_dir, "Fig3A_volcano.pdf"),
        "\n - ", file.path(out_dir, "Fig3C_patientDelta_waterfall.png"),
        "\n - ", file.path(out_dir, "Fig3C_patientDelta_waterfall.pdf"),
        "\n - ", file.path(out_dir, "Fig3C_patientDelta_table.csv")
)
