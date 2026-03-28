#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggrepel); library(patchwork)
})

meta_f <- "data/sample_sheet.csv"
cand_score <- c(
  "results_gsva/GSVA_Reactome_scores.csv",
  "results_gsva/GSVA_Hallmark_scores.csv",
  "results_gsva/GSVA_scores.csv",
  "results_gsva/GSVA_Reactome_matrix.csv",
  "results_gsva/GSVA_Hallmark_matrix.csv",
  "results_gsva/GSVA_matrix.csv"
)
cand_did <- c(
  "results_gsva/GSVA_Reactome_DiD.csv",
  "results_gsva/GSVA_Hallmark_DiD.csv",
  "results/GSVA_Reactome_DiD.csv",
  "results/GSVA_Hallmark_DiD.csv"
)

stopifnot(file.exists(meta_f))
score_f <- cand_score[file.exists(cand_score)][1]
did_f   <- cand_did[file.exists(cand_did)][1]
if (is.na(score_f)) stop("❌ Missing GSVA per-sample matrix in results_gsva/.")

dir.create("figs", showWarnings=FALSE, recursive=TRUE)
dir.create("figs/Fig4_tables", showWarnings=FALSE, recursive=TRUE)

# ---------- palettes ----------
col_pathway <- c(
  "Repair"        = "#D73027",
  "Inflammation"  = "#4575B4",
  "ECM/Wnt"       = "#1A9850",
  "Other"         = "#7F7F7F"
)

pick <- function(df, cands){ hit <- cands[cands %in% colnames(df)]; if(length(hit)==0) NA_character_ else hit[1] }
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
to_time2 <- function(x){
  x2 <- tolower(as.character(x))
  x2 <- str_replace_all(x2, "week", "w")
  x2 <- ifelse(x2 %in% c("0","w0","week0"), "w0", x2)
  x2 <- ifelse(x2 %in% c("1","w1","week1"), "w1", x2)
  x2
}
short_path <- function(x, n=44){
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^REACTOME_", "", x)
  x <- gsub("_", " ", x)
  x <- str_squish(x)
  str_trunc(x, width=n, side="right", ellipsis="…")
}
cat_path <- function(p){
  p2 <- tolower(p)
  if (str_detect(p2, "keratin|epiderm|cornif|skin|epitheli|barrier|differentiation|wound")) return("Repair")
  if (str_detect(p2, "vegf|angiogen|growth factor|fgf|tgf|egf|pdgf|migration")) return("Repair")
  if (str_detect(p2, "wnt|beta-catenin|notch|hedgehog")) return("ECM/Wnt")
  if (str_detect(p2, "extracellular matrix|ecm|collagen|integrin|matrix|focal adhesion")) return("ECM/Wnt")
  if (str_detect(p2, "inflamm|interferon|il-|tnf|nf-kb|chemok|immune|t cell|b cell|neutroph|macroph|complement")) return("Inflammation")
  "Other"
}

# ---------- meta ----------
meta <- read_csv(meta_f, show_col_types=FALSE, progress=FALSE)
sample_col  <- pick(meta, c("sample","Sample","sample_id","SampleID","gsm","GSM"))
arm_col     <- pick(meta, c("arm","Arm","group","Group","treatment","Treatment"))
time_col    <- pick(meta, c("time","Time","week","Week","timepoint","Timepoint"))
patient_col <- pick(meta, c("patient","Patient","subject","Subject","id","ID"))
if(any(is.na(c(sample_col,arm_col,time_col,patient_col)))) {
  stop("❌ sample_sheet.csv must contain sample/arm/time/patient columns.")
}
meta2 <- meta %>%
  transmute(sample=as.character(.data[[sample_col]]),
            arm=as.character(.data[[arm_col]]),
            time=as.character(.data[[time_col]]),
            patient=as.character(.data[[patient_col]])) %>%
  mutate(time2=to_time2(time)) %>%
  filter(time2 %in% c("w0","w1"))

arm_names <- unique(meta2$arm)
arm_lower <- tolower(arm_names)
arm_blcc <- if(any(str_detect(arm_lower, "blcc|apligraf|treat"))) arm_names[which(str_detect(arm_lower, "blcc|apligraf|treat"))[1]] else arm_names[1]
arm_ctrl <- if(any(str_detect(arm_lower, "control|ctrl|standard|soc|compression|placebo"))) arm_names[which(str_detect(arm_lower, "control|ctrl|standard|soc|compression|placebo"))[1]] else setdiff(arm_names, arm_blcc)[1]

# ---------- GSVA matrix ----------
gs <- read_csv(score_f, show_col_types=FALSE, progress=FALSE)
path_col <- colnames(gs)[1]
gs <- gs %>% rename(pathway = all_of(path_col))
gs[-1] <- lapply(gs[-1], as_num)

keep_samps <- intersect(meta2$sample, colnames(gs))
if(length(keep_samps) < 6) stop("❌ Too few matched samples between GSVA and sample_sheet.")
meta2 <- meta2 %>% filter(sample %in% keep_samps)

mat <- as.matrix(gs[, keep_samps, drop=FALSE]); rownames(mat) <- gs$pathway

meta2 <- meta2 %>%
  mutate(time2=factor(time2, levels=c("w0","w1")),
         arm=factor(arm, levels=unique(arm)))
meta_ord <- meta2 %>% arrange(arm, time2, patient, sample)
mat_ord <- mat[, meta_ord$sample, drop=FALSE]

# ---------- long & delta ----------
wide <- tibble(sample=colnames(mat_ord)) %>%
  left_join(meta_ord, by="sample") %>%
  cbind(t(mat_ord)) %>%
  as_tibble()

long_scores <- wide %>%
  pivot_longer(cols = -c(sample, arm, time, time2, patient),
               names_to="pathway", values_to="score")

delta_long <- long_scores %>%
  select(arm, time2, patient, pathway, score) %>%
  pivot_wider(names_from=time2, values_from=score) %>%
  filter(!is.na(w0), !is.na(w1)) %>%
  mutate(delta = w1 - w0)

# ---------- DiD stats ----------
did_tbl <- NULL
if(!is.na(did_f)) {
  did_raw <- read_csv(did_f, show_col_types=FALSE, progress=FALSE)
  path_c <- pick(did_raw, c("pathway","Pathway","set","term","NAME"))
  eff_c  <- pick(did_raw, c("effect","DiD","estimate","coef","delta","diff"))
  fdr_c  <- pick(did_raw, c("FDR","adj.P.Val","padj","qvalue"))
  if(!any(is.na(c(path_c, eff_c, fdr_c)))) {
    did_tbl <- did_raw %>%
      transmute(pathway=as.character(.data[[path_c]]),
                effect=as_num(.data[[eff_c]]),
                FDR=as_num(.data[[fdr_c]])) %>%
      filter(!is.na(effect), !is.na(FDR))
  }
}
if(is.null(did_tbl)) {
  did_tbl <- delta_long %>%
    group_by(pathway) %>%
    summarize(
      effect = mean(delta[arm==arm_blcc], na.rm=TRUE) - mean(delta[arm==arm_ctrl], na.rm=TRUE),
      p = tryCatch(t.test(delta[arm==arm_blcc], delta[arm==arm_ctrl])$p.value, error=function(e) NA_real_),
      .groups="drop"
    ) %>%
    mutate(FDR=p.adjust(p, method="BH")) %>%
    select(-p)
}

did_tbl <- did_tbl %>%
  filter(pathway %in% rownames(mat_ord)) %>%
  mutate(
    category = vapply(pathway, cat_path, character(1)),
    neglogFDR = -log10(pmax(FDR, 1e-300)),
    pathway_short = short_path(pathway, 44)
  )

write_csv(did_tbl, "figs/Fig4_tables/DiD_pathway_table.csv")

# =========================
# Panel A: Forest (Top 28)
# =========================
TOP_A <- 28
A_df <- did_tbl %>%
  arrange(FDR, desc(abs(effect))) %>%
  slice_head(n=TOP_A) %>%
  mutate(pathway_short = factor(pathway_short, levels=rev(pathway_short)))

pA <- ggplot(A_df, aes(y=pathway_short, x=effect, color=category)) +
  geom_vline(xintercept=0, linewidth=0.35, color="grey55") +
  geom_segment(aes(x=0, xend=effect, yend=pathway_short), linewidth=0.95) +
  geom_point(size=2.7) +
  scale_color_manual(values=col_pathway) +
  theme_classic(base_size=12) +
  theme(legend.position="top", legend.title=element_blank(),
        axis.title.y=element_blank()) +
  labs(x="DiD effect size (ΔBLCC − ΔControl)")

# =========================
# Panel B: Volcano (label Top 8)
# =========================
labB <- did_tbl %>%
  filter(FDR < 0.05) %>%
  arrange(desc(abs(effect) * neglogFDR)) %>%
  slice_head(n=8)

pB <- ggplot(did_tbl, aes(effect, neglogFDR, color=category)) +
  geom_vline(xintercept=0, linewidth=0.3, color="grey55") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth=0.3, color="grey55") +
  geom_point(size=2.0, alpha=0.9) +
  ggrepel::geom_text_repel(
    data=labB, aes(label=pathway_short),
    size=3.0, box.padding=0.25, point.padding=0.15,
    segment.size=0.25, max.overlaps=Inf
  ) +
  scale_color_manual(values=col_pathway) +
  theme_classic(base_size=12) +
  theme(legend.position="top", legend.title=element_blank()) +
  labs(x="DiD effect size (ΔBLCC − ΔControl)", y="-log10(FDR)")

# =========================
# Panel C: Patient Δ waterfall (2 pathways × 2 arms)
# =========================
pick_one <- function(cat, keys){
  x <- did_tbl %>% filter(category==cat) %>% arrange(FDR, desc(abs(effect)))
  for(k in keys){
    hit <- x %>% filter(str_detect(toupper(pathway), k))
    if(nrow(hit)>0) return(hit$pathway[1])
  }
  if(nrow(x)>0) return(x$pathway[1])
  NA_character_
}
pw_repair <- pick_one("Repair", c("KERATIN","EPIDERM","ANGIO","VEGF","WOUND","MIGRAT"))
pw_inflam <- pick_one("Inflammation", c("INFLAM","TNF","INTERFERON","COMPLEMENT","IL","NF"))
selC <- unique(na.omit(c(pw_repair, pw_inflam)))
if(length(selC) < 1) selC <- did_tbl %>% arrange(FDR, desc(abs(effect))) %>% slice_head(n=1) %>% pull(pathway)

C_df <- delta_long %>%
  filter(pathway %in% selC) %>%
  mutate(
    pathway_short = short_path(pathway, 36),
    arm = factor(as.character(arm), levels=c(arm_ctrl, arm_blcc))
  ) %>%
  group_by(pathway_short, arm) %>%
  arrange(desc(delta), .by_group=TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

pC <- ggplot(C_df, aes(x=rank, y=delta, color=arm)) +
  geom_hline(yintercept=0, linewidth=0.35, color="grey55") +
  geom_segment(aes(xend=rank, y=0, yend=delta), linewidth=0.95) +
  geom_point(size=2.6) +
  facet_grid(pathway_short ~ arm, scales="free_x", space="free_x") +
  theme_classic(base_size=12) +
  theme(
    legend.position="none",
    strip.background=element_blank(),
    strip.text=element_text(face="bold", size=10),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  ) +
  scale_color_manual(values=setNames(c("#3182BD", "#CB181D"), c(arm_ctrl, arm_blcc))) +
  labs(x="Patients (sorted within arm)", y=expression(Delta~"GSVA (week1 - week0)"))

write_csv(C_df, "figs/Fig4_tables/PanelC_delta_waterfall_table.csv")
write_lines(selC, "figs/Fig4_tables/PanelC_pathways_used.txt")

# =========================
# Panel D: LOPO robustness（改：不同颜色，按通路类别上色）
# =========================
patients <- sort(unique(meta_ord$patient))
topD <- did_tbl %>% arrange(FDR, desc(abs(effect))) %>% slice_head(n=12) %>% pull(pathway)

lopo_prop <- lapply(topD, function(path){
  dfp <- delta_long %>% filter(pathway==path)
  if(nrow(dfp)==0) return(NA_real_)
  full_eff <- did_tbl$effect[match(path, did_tbl$pathway)]
  v <- sapply(patients, function(p){
    df2 <- dfp %>% filter(patient != p)
    if(length(unique(df2$arm))<2) return(NA_real_)
    mean(df2$delta[df2$arm==arm_blcc], na.rm=TRUE) - mean(df2$delta[df2$arm==arm_ctrl], na.rm=TRUE)
  })
  v <- v[is.finite(v)]
  if(length(v)==0 || !is.finite(full_eff) || full_eff==0) return(NA_real_)
  mean(sign(v) == sign(full_eff))
})

D_df <- did_tbl %>%
  filter(pathway %in% topD) %>%
  transmute(
    pathway,
    pathway_short = short_path(pathway, 44),
    category,
    prop_consistent = unlist(lopo_prop)[match(pathway, topD)]
  ) %>%
  arrange(desc(prop_consistent), pathway_short) %>%
  mutate(pathway_short = factor(pathway_short, levels=rev(pathway_short)))

write_csv(D_df, "figs/Fig4_tables/LOPO_robustness_table.csv")

pD <- ggplot(D_df, aes(x=prop_consistent, y=pathway_short, color=category)) +
  geom_vline(xintercept=0.8, linetype="dashed", linewidth=0.3, color="grey55") +
  geom_segment(aes(x=0, xend=prop_consistent, yend=pathway_short),
               linewidth=0.95, alpha=0.95) +
  geom_point(size=2.6) +
  scale_color_manual(values=col_pathway) +
  scale_x_continuous(limits=c(0,1)) +
  theme_classic(base_size=12) +
  theme(
    legend.position="top",
    legend.title=element_blank(),
    axis.title.y=element_blank()
  ) +
  labs(x="Proportion consistent (LOPO)")

# =========================
# Assemble
# =========================
left  <- pA / pD + plot_layout(heights=c(1.15, 0.85))
right <- pB / pC + plot_layout(heights=c(1.0, 1.0))
fig <- (left | right) + plot_layout(widths=c(1.05, 1.15))

ggsave("figs/Fig4_paperStyle_colorD.pdf", fig, width=14, height=10.5, device="pdf")
ggsave("figs/Fig4_paperStyle_colorD.png", fig, width=14, height=10.5, dpi=300)

cat("✅ Saved: figs/Fig4_paperStyle_colorD.(png/pdf)\n")
cat("Using GSVA matrix:", score_f, "\n")
cat("Using DiD file:", did_f, "\n")
cat("arm_ctrl=", arm_ctrl, " arm_blcc=", arm_blcc, "\n")
cat("Panel C pathways:\n", paste(selC, collapse="\n"), "\n")
