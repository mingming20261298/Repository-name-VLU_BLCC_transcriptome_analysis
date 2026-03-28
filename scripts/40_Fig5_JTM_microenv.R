#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggrepel); library(patchwork)
  library(grid)
})

dir.create("figs", showWarnings=FALSE, recursive=TRUE)
dir.create("figs/Fig5_tables", showWarnings=FALSE, recursive=TRUE)

meta_f   <- "data/sample_sheet.csv"
gsva_f   <- "results/GSVA_Hallmark_interaction.csv"
xcell_f  <- "results/xCell_scores_all_samples.csv"
mcp_f    <- "results/MCPcounter_scores_all_samples.csv"

stopifnot(file.exists(meta_f))
stopifnot(file.exists(gsva_f))
stopifnot(file.exists(xcell_f))
stopifnot(file.exists(mcp_f))

# -------------------------
# helpers
# -------------------------
pick <- function(df, cands){ hit <- cands[cands %in% colnames(df)]; if(length(hit)==0) NA_character_ else hit[1] }
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
to_time2 <- function(x){
  x2 <- tolower(as.character(x))
  x2 <- str_replace_all(x2, "week", "w")
  x2 <- ifelse(x2 %in% c("0","w0","week0"), "w0", x2)
  x2 <- ifelse(x2 %in% c("1","w1","week1"), "w1", x2)
  x2
}
short_nm <- function(x, n=34){
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^REACTOME_", "", x)
  x <- gsub("_", " ", x)
  x <- str_squish(x)
  str_trunc(x, width=n, side="right", ellipsis="…")
}

# 样本名标准化：解决 GSM / “GSM123” / GSM-123 / GSM.123 / 空格等
norm_id <- function(x){
  x <- as.character(x)
  x <- str_replace_all(x, '"', "")
  x <- str_replace_all(x, "'", "")
  x <- str_trim(x)
  x <- tolower(x)
  x <- str_replace_all(x, "[\\.\\-\\s]+", "_")
  x <- str_replace_all(x, "__+", "_")
  x
}

# JTM colors
COL <- list(
  ctrl = "#3182BD",
  trt  = "#CB181D",
  cat = c(
    "Repair/Epithelial"    = "#D73027",
    "Inflammation/Immune"  = "#4575B4",
    "Stromal/ECM/Vascular" = "#1A9850",
    "Other"                = "#7F7F7F"
  )
)

cat_prog <- function(p){
  p2 <- tolower(p)
  if (str_detect(p2, "keratin|epiderm|epitheli|barrier|cornif|wound|migration|angiogen|vegf|tgf|egf|fgf|pdgf|mitotic|cell cycle")) return("Repair/Epithelial")
  if (str_detect(p2, "extracellular matrix|ecm|collagen|integrin|matrix|focal adhesion|wnt|beta-catenin|notch|vascular|endothelial")) return("Stromal/ECM/Vascular")
  if (str_detect(p2, "inflamm|interferon|tnf|nf-kb|il-|chemok|immune|complement|allograft|t cell|b cell|neutroph|macroph")) return("Inflammation/Immune")
  "Other"
}

# reader
read_score_matrix <- function(f, feat_name){
  x <- read_csv(f, show_col_types=FALSE, progress=FALSE)
  fc <- colnames(x)[1]
  x <- x %>% rename(!!feat_name := all_of(fc))
  x[-1] <- lapply(x[-1], as_num)
  # normalize sample colnames (except feature col)
  cn <- colnames(x)
  cn2 <- cn
  cn2[-1] <- norm_id(cn[-1])
  colnames(x) <- cn2
  x
}

# -------------------------
# meta
# -------------------------
meta <- read_csv(meta_f, show_col_types=FALSE, progress=FALSE)
sample_col  <- pick(meta, c("sample","Sample","sample_id","SampleID","gsm","GSM"))
arm_col     <- pick(meta, c("arm","Arm","group","Group","treatment","Treatment"))
time_col    <- pick(meta, c("time","Time","week","Week","timepoint","Timepoint"))
patient_col <- pick(meta, c("patient","Patient","subject","Subject","id","ID"))
if(any(is.na(c(sample_col,arm_col,time_col,patient_col)))) {
  stop("sample_sheet.csv must contain sample/arm/time/patient columns. Found:\n",
       paste(colnames(meta), collapse=", "))
}

meta2 <- meta %>%
  transmute(
    sample_raw  = as.character(.data[[sample_col]]),
    sample      = norm_id(.data[[sample_col]]),
    arm         = as.character(.data[[arm_col]]),
    time        = as.character(.data[[time_col]]),
    patient     = as.character(.data[[patient_col]])
  ) %>%
  mutate(time2 = factor(to_time2(time), levels=c("w0","w1"))) %>%
  filter(time2 %in% c("w0","w1"))

# guess BLCC vs Control
arm_names <- unique(meta2$arm)
arm_lower <- tolower(arm_names)
arm_blcc <- if(any(str_detect(arm_lower, "blcc|apligraf|treat"))) arm_names[which(str_detect(arm_lower, "blcc|apligraf|treat"))[1]] else arm_names[1]
arm_ctrl <- if(any(str_detect(arm_lower, "control|ctrl|standard|soc|compression|placebo"))) arm_names[which(str_detect(arm_lower, "control|ctrl|standard|soc|compression|placebo"))[1]] else setdiff(arm_names, arm_blcc)[1]

meta2 <- meta2 %>%
  mutate(arm = factor(arm, levels=c(arm_ctrl, arm_blcc)))

# -------------------------
# load matrices (colnames already normalized)
# -------------------------
gsva  <- read_score_matrix(gsva_f,  "program")
xcell <- read_score_matrix(xcell_f, "cell")
mcp   <- read_score_matrix(mcp_f,   "cell")

# -------------------------
# choose sample sets with maximum overlap (do not force all 4)
# -------------------------
meta_s <- meta2$sample
gsva_s <- colnames(gsva)[-1]
xcel_s <- colnames(xcell)[-1]
mcp_s  <- colnames(mcp)[-1]

ov_gsva <- intersect(meta_s, gsva_s)
ov_xcel <- intersect(meta_s, xcel_s)
ov_mcp  <- intersect(meta_s, mcp_s)

cat("Matched samples:\n")
cat(" meta∩gsva =", length(ov_gsva), "\n")
cat(" meta∩xcell=", length(ov_xcel), "\n")
cat(" meta∩mcp  =", length(ov_mcp), "\n")

if(length(ov_gsva) < 6) stop("❌ meta 与 GSVA 样本匹配 < 6，请先检查 GSVA_Hallmark_interaction.csv 的列名是不是 GSM。")

use_xcell <- length(ov_xcel) >= 6
use_mcp   <- length(ov_mcp)  >= 6

# for A/D (program-level): only need meta + gsva
keep_prog <- ov_gsva
meta_prog <- meta2 %>% filter(sample %in% keep_prog)

# for cell panels: prefer MCP for composition; use whichever matches
keep_cell <- NULL
cell_source <- NULL
if(use_mcp) { keep_cell <- ov_mcp; cell_source <- "MCP" }
if(!use_mcp && use_xcell) { keep_cell <- ov_xcel; cell_source <- "xCell" }

# -------------------------
# build long tables
# -------------------------
to_long <- function(df, feat){
  df %>% pivot_longer(cols=-all_of(feat), names_to="sample", values_to="score")
}

# program matrix pathway x sample -> long
gsva_long <- to_long(gsva %>% select(program, all_of(keep_prog)), "program") %>%
  left_join(meta_prog %>% select(sample, arm, time2, patient), by="sample")

make_delta <- function(long_df, feat){
  long_df %>%
    select(arm, time2, patient, !!sym(feat), score) %>%
    pivot_wider(names_from=time2, values_from=score) %>%
    filter(!is.na(w0), !is.na(w1)) %>%
    mutate(delta = w1 - w0)
}
did_from_delta <- function(df, feat){
  df %>%
    group_by(.data[[feat]]) %>%
    summarize(
      effect = mean(delta[arm==arm_blcc], na.rm=TRUE) - mean(delta[arm==arm_ctrl], na.rm=TRUE),
      p = tryCatch(t.test(delta[arm==arm_blcc], delta[arm==arm_ctrl])$p.value, error=function(e) NA_real_),
      .groups="drop"
    ) %>%
    mutate(FDR=p.adjust(p,"BH"), neglogFDR=-log10(pmax(FDR,1e-300)))
}

gsva_delta <- make_delta(gsva_long, "program")
gsva_did   <- did_from_delta(gsva_delta, "program") %>%
  mutate(category=vapply(program, cat_prog, character(1)),
         label=short_nm(program, 36))

write_csv(gsva_did, "figs/Fig5_tables/GSVA_DiD.csv")

# healing shift ordering (week1 sample ordering)
top_sig <- gsva_did %>% arrange(FDR, desc(abs(effect))) %>% slice_head(n=40) %>% pull(program)
hm <- gsva %>% filter(program %in% top_sig) %>% column_to_rownames("program") %>% as.matrix()
hm <- hm[, keep_prog, drop=FALSE]
hm_z <- t(scale(t(hm))); hm_z[is.na(hm_z)] <- 0
prog_cat <- tibble(program=rownames(hm_z), category=vapply(rownames(hm_z), cat_prog, character(1)))

shift_df <- as_tibble(hm_z, rownames="program") %>%
  pivot_longer(-program, names_to="sample", values_to="z") %>%
  left_join(prog_cat, by="program") %>%
  mutate(set = case_when(
    category %in% c("Repair/Epithelial","Stromal/ECM/Vascular") ~ "repair",
    category=="Inflammation/Immune" ~ "inflam",
    TRUE ~ NA_character_
  )) %>% filter(!is.na(set)) %>%
  group_by(sample, set) %>% summarize(m=mean(z, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=set, values_from=m) %>%
  mutate(healing_shift = repair - inflam)

meta_prog <- meta_prog %>% left_join(shift_df %>% select(sample, healing_shift), by="sample")
write_csv(meta_prog, "figs/Fig5_tables/meta_prog_with_shift.csv")

# ============================================================
# Panel A (program delta jitter + mean±CI) — 18 programs
# ============================================================
selA <- gsva_did %>%
  arrange(FDR, desc(abs(effect))) %>%
  group_by(category) %>% slice_head(n=8) %>% ungroup() %>%
  slice_head(n=18) %>% pull(program)

A_df <- gsva_delta %>%
  filter(program %in% selA) %>%
  left_join(gsva_did %>% select(program, FDR, effect, category, label), by="program") %>%
  mutate(label = factor(label, levels=rev(unique(label[match(selA, program)]))),
         arm = factor(as.character(arm), levels=c(arm_ctrl, arm_blcc)))

A_sum <- A_df %>%
  group_by(label, arm) %>%
  summarize(
    m=mean(delta, na.rm=TRUE),
    se=sd(delta, na.rm=TRUE)/sqrt(sum(!is.na(delta))),
    lo=m-1.96*se, hi=m+1.96*se,
    .groups="drop"
  )

stars <- function(q) ifelse(q<0.001,"***",ifelse(q<0.01,"**",ifelse(q<0.05,"*","")))
A_star <- gsva_did %>%
  filter(program %in% selA) %>%
  transmute(label=factor(short_nm(program,36), levels=levels(A_df$label)),
            star=stars(FDR))

xmaxA <- max(A_df$delta, na.rm=TRUE)

pA <- ggplot(A_df, aes(x=delta, y=label)) +
  geom_vline(xintercept=0, linewidth=0.3, color="grey60") +
  geom_point(aes(color=arm), size=1.8,
             position=position_jitter(height=0.12, width=0),
             alpha=0.9) +
  geom_errorbarh(data=A_sum, aes(xmin=lo, xmax=hi, y=label, color=arm),
                 height=0.18, linewidth=0.7) +
  geom_point(data=A_sum, aes(x=m, y=label, color=arm), size=2.5) +
  geom_text(data=A_star, aes(x=xmaxA, y=label, label=star),
            inherit.aes=FALSE, hjust=-0.2, size=4) +
  scale_color_manual(values=setNames(c(COL$ctrl, COL$trt), c(arm_ctrl, arm_blcc))) +
  theme_classic(base_size=12) +
  theme(legend.position="top", legend.title=element_blank(),
        axis.title.y=element_blank(),
        plot.margin=margin(5,8,5,5)) +
  labs(x="Δ program score (week1 − week0)", y=NULL)

write_csv(A_df, "figs/Fig5_tables/PanelA_program_delta.csv")

# ============================================================
# Panel B/C (cell) — auto downgrade if not enough samples
# ============================================================
make_placeholder <- function(txt){
  ggplot() + theme_void() +
    annotate("text", x=0, y=0, label=txt, size=4.2, hjust=0.5, vjust=0.5) +
    xlim(-1,1) + ylim(-1,1)
}

pB <- NULL; pC <- NULL

if(is.null(keep_cell) || length(keep_cell) < 6) {
  pB <- make_placeholder("Panel B: cell scores sample IDs do not match meta\n(need >=6 matched samples)")
  pC <- make_placeholder("Panel C: composition needs matched cell scores\n(fix sample IDs then rerun)")
} else {
  if(cell_source=="MCP"){
    cell_mat <- mcp %>% select(cell, all_of(keep_cell))
    cell_long <- to_long(cell_mat, "cell") %>%
      left_join(meta2 %>% filter(sample %in% keep_cell) %>% select(sample, arm, time2, patient), by="sample") %>%
      mutate(method="MCPcounter")
  } else {
    cell_mat <- xcell %>% select(cell, all_of(keep_cell))
    cell_long <- to_long(cell_mat, "cell") %>%
      left_join(meta2 %>% filter(sample %in% keep_cell) %>% select(sample, arm, time2, patient), by="sample") %>%
      mutate(method="xCell")
  }

  cell_delta <- make_delta(cell_long, "cell")
  cell_did <- did_from_delta(cell_delta, "cell") %>%
    mutate(label=short_nm(cell, 30),
           method=cell_source)

  write_csv(cell_did, "figs/Fig5_tables/Cell_DiD_used.csv")

  # Panel B: top 12
  selB <- cell_did %>% arrange(FDR, desc(abs(effect))) %>% slice_head(n=12) %>% pull(cell)
  B_df <- cell_delta %>%
    filter(cell %in% selB) %>%
    left_join(cell_did %>% select(cell, label, FDR), by="cell") %>%
    mutate(
      arm = factor(as.character(arm), levels=c(arm_ctrl, arm_blcc)),
      feat = factor(label, levels=rev(unique(label[match(selB, cell)])))
    )

  pB <- ggplot(B_df, aes(x=delta, y=feat, color=arm)) +
    geom_vline(xintercept=0, linewidth=0.3, color="grey60") +
    geom_boxplot(aes(group=interaction(arm, feat)),
                 width=0.55, alpha=0.10, outlier.shape=NA,
                 linewidth=0.45, color="black") +
    geom_point(size=1.6, alpha=0.9,
               position=position_jitter(height=0.12, width=0)) +
    scale_color_manual(values=setNames(c(COL$ctrl, COL$trt), c(arm_ctrl, arm_blcc))) +
    theme_classic(base_size=12) +
    theme(legend.position="top", legend.title=element_blank(),
          axis.title.y=element_blank(),
          plot.margin=margin(5,8,5,5)) +
    labs(x=paste0("Δ cell score (week1 − week0) [", cell_source, "]"), y=NULL)

  # Panel C: week1 composition stacked (top 8 mean)
  w1_mean <- cell_long %>%
    filter(time2=="w1") %>%
    group_by(cell) %>% summarize(m=mean(score, na.rm=TRUE), .groups="drop") %>%
    arrange(desc(m))
  comp_cells <- head(w1_mean$cell, 8)

  C_df <- cell_long %>%
    filter(time2=="w1", cell %in% comp_cells) %>%
    left_join(meta2 %>% select(sample, healing_shift) %>% distinct(), by="sample") %>%
    mutate(score_pos=pmax(score,0))

  C_rel <- C_df %>%
    group_by(sample, arm, healing_shift) %>%
    mutate(total=sum(score_pos, na.rm=TRUE),
           rel=ifelse(total>0, score_pos/total, 0)) %>%
    ungroup() %>%
    mutate(cell_short=factor(short_nm(cell,20), levels=short_nm(comp_cells,20)))

  ord_samples <- meta2 %>% filter(time2=="w1", sample %in% keep_cell) %>%
    arrange(desc(healing_shift), arm, patient, sample) %>% pull(sample)
  C_rel$sample <- factor(C_rel$sample, levels=ord_samples)

  cell_levels <- levels(C_rel$cell_short)
  pal_cells <- setNames(c(
    "#1b9e77","#d95f02","#7570b3","#e7298a",
    "#66a61e","#e6ab02","#a6761d","#666666"
  )[seq_along(cell_levels)], cell_levels)

  pC <- ggplot(C_rel, aes(x=sample, y=rel, fill=cell_short)) +
    geom_col(width=0.95) +
    facet_grid(. ~ arm, scales="free_x", space="free_x") +
    theme_classic(base_size=12) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.position="right",
      legend.title=element_blank(),
      strip.background=element_blank(),
      strip.text=element_text(face="bold"),
      plot.margin=margin(5,8,5,5)
    ) +
    scale_fill_manual(values=pal_cells) +
    labs(x="Week1 samples (ordered by healing shift)", y="Relative composition")

  write_csv(C_rel, "figs/Fig5_tables/PanelC_composition_week1.csv")
}

# ============================================================
# Panel D: summary effect (program + cell if exists)
# ============================================================
D_prog <- gsva_did %>% arrange(FDR, desc(abs(effect))) %>% slice_head(n=12) %>%
  transmute(feature=short_nm(program,36), effect=effect, category=category)

D_all <- D_prog
if(exists("cell_did") && !is.null(cell_did) && nrow(cell_did) > 0){
  D_cell <- cell_did %>% arrange(FDR, desc(abs(effect))) %>% slice_head(n=10) %>%
    transmute(feature=paste0(short_nm(cell,30), " [", cell_source, "]"),
              effect=effect,
              category=ifelse(effect>=0,"Stromal/ECM/Vascular","Inflammation/Immune"))
  D_all <- bind_rows(D_prog, D_cell)
}

D_all <- D_all %>%
  mutate(cat2=ifelse(category %in% names(COL$cat), category, "Other"),
         feature=factor(feature, levels=rev(unique(feature))))

pD <- ggplot(D_all, aes(x=effect, y=feature, color=cat2)) +
  geom_vline(xintercept=0, linewidth=0.3, color="grey60") +
  geom_segment(aes(x=0, xend=effect, yend=feature), linewidth=0.95) +
  geom_point(size=2.7) +
  scale_color_manual(values=COL$cat) +
  theme_classic(base_size=12) +
  theme(legend.position="top", legend.title=element_blank(),
        axis.title.y=element_blank(),
        plot.margin=margin(5,8,5,5)) +
  labs(x="Effect size (ΔBLCC − ΔControl)", y=NULL)

write_csv(D_all, "figs/Fig5_tables/PanelD_effect_summary.csv")

# ============================================================
# Assemble storyboard with row labels
# ============================================================
row1_lab <- wrap_elements(grid::textGrob("Program-level remodeling", rot=90,
                                        gp=grid::gpar(fontface="bold", fontsize=12)))
row2_lab <- wrap_elements(grid::textGrob("Cellular microenvironment shift", rot=90,
                                        gp=grid::gpar(fontface="bold", fontsize=12)))

row1 <- (row1_lab | (pA | pB)) + plot_layout(widths=c(0.06, 0.94))
row2 <- (row2_lab | (pC | pD)) + plot_layout(widths=c(0.06, 0.94))

fig5 <- row1 / row2 +
  plot_layout(heights=c(1,1)) +
  plot_annotation(tag_levels="A") &
  theme(plot.tag = element_text(face="bold", size=14),
        plot.tag.position = c(0.01, 0.98))

ggsave("figs/Fig5_JTM.pdf", fig5, width=14, height=10.5, device="pdf")
ggsave("figs/Fig5_JTM.png", fig5, width=14, height=10.5, dpi=300)

cat("✅ Saved: figs/Fig5_JTM.(png/pdf)\n")
cat("arm_ctrl=", arm_ctrl, " arm_blcc=", arm_blcc, "\n")
cat("Cell panels source:", ifelse(is.null(cell_source), "NONE", cell_source), "\n")
