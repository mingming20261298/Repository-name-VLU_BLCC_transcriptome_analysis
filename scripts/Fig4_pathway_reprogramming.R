#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(pheatmap); library(RColorBrewer)
  library(patchwork); library(ggrepel)
})
dir.create("figs", showWarnings=FALSE, recursive=TRUE)
dir.create("figs/Fig4_tables", showWarnings=FALSE, recursive=TRUE)

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

cat("Using score_f:", score_f, "\n")
cat("Using did_f:", did_f, "\n")

if (is.na(score_f)) {
  stop("❌ 没找到 GSVA per-sample scores 矩阵（Panel A/C 必须要它）。\n",
       "请确认以下任一文件存在（rows=pathway, cols=samples, 第一列=pathway）：\n",
       paste(cand_score, collapse="\n"))
}

pick <- function(df, cands){ hit <- cands[cands %in% colnames(df)]; if(length(hit)==0) NA_character_ else hit[1] }
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
to_time2 <- function(x){
  x2 <- tolower(as.character(x))
  x2 <- str_replace_all(x2, "week", "w")
  x2 <- ifelse(x2 %in% c("0","w0","week0"), "w0", x2)
  x2 <- ifelse(x2 %in% c("1","w1","week1"), "w1", x2)
  x2
}
cat_path <- function(p){
  p2 <- tolower(p)
  if (str_detect(p2, "keratin|epiderm|cornif|skin|epitheli|barrier")) return("Repair/Keratin")
  if (str_detect(p2, "vegf|angiogen|growth factor|fgf|tgf|egf|pdgf")) return("Repair/Growth")
  if (str_detect(p2, "wnt|beta-catenin|notch|hedgehog|extracellular matrix|ecm|collagen|integrin|matrix")) return("ECM/Wnt")
  if (str_detect(p2, "inflamm|interferon|il-|tnf|nf-kb|chemok|immune|t cell|b cell|neutroph|macroph")) return("Inflammation/Immune")
  "Other"
}

# ---- meta
meta <- read_csv(meta_f, show_col_types=FALSE, progress=FALSE)
sample_col  <- pick(meta, c("sample","Sample","sample_id","SampleID","gsm","GSM"))
arm_col     <- pick(meta, c("arm","Arm","group","Group","treatment","Treatment"))
time_col    <- pick(meta, c("time","Time","week","Week","timepoint","Timepoint"))
patient_col <- pick(meta, c("patient","Patient","subject","Subject","id","ID"))
if(any(is.na(c(sample_col,arm_col,time_col,patient_col)))) {
  stop("❌ sample_sheet.csv 缺 sample/arm/time/patient 列。现有列：\n", paste(colnames(meta), collapse=", "))
}
meta2 <- meta %>%
  transmute(sample=as.character(.data[[sample_col]]),
            arm=as.character(.data[[arm_col]]),
            time=as.character(.data[[time_col]]),
            patient=as.character(.data[[patient_col]])) %>%
  mutate(time2=to_time2(time)) %>%
  filter(time2 %in% c("w0","w1"))

# ---- GSVA score matrix (pathway x sample)
gs <- read_csv(score_f, show_col_types=FALSE, progress=FALSE)
path_col <- colnames(gs)[1]
gs <- gs %>% rename(pathway = all_of(path_col))
gs[-1] <- lapply(gs[-1], as_num)

keep_samps <- intersect(meta2$sample, colnames(gs))
if(length(keep_samps) < 6) stop("❌ GSVA matrix 与 sample_sheet 匹配样本太少。")
meta2 <- meta2 %>% filter(sample %in% keep_samps)
mat <- as.matrix(gs[, keep_samps, drop=FALSE]); rownames(mat) <- gs$pathway

meta2 <- meta2 %>% mutate(time2=factor(time2, levels=c("w0","w1")))
arm_levels <- unique(meta2$arm)
meta2 <- meta2 %>% mutate(arm=factor(arm, levels=arm_levels))
meta_ord <- meta2 %>% arrange(arm, time2, patient, sample)
mat_ord <- mat[, meta_ord$sample, drop=FALSE]

# ---- patient delta (week1-week0)
wide <- as_tibble(t(mat_ord), rownames="sample") %>% left_join(meta_ord, by="sample")
delta_long <- wide %>%
  pivot_longer(cols = -c(sample, arm, time2, patient, time), names_to="pathway", values_to="score") %>%
  select(arm, time2, patient, pathway, score) %>%
  pivot_wider(names_from=time2, values_from=score) %>%
  filter(!is.na(w0), !is.na(w1)) %>%
  mutate(delta = w1 - w0)

# determine BLCC vs Control names
arm_names <- as.character(unique(meta_ord$arm))
arm_lower <- tolower(arm_names)
arm_blcc <- if(any(str_detect(arm_lower,"blcc|apligraf|treat"))) arm_names[which(str_detect(arm_lower,"blcc|apligraf|treat"))[1]] else arm_names[1]
arm_ctrl <- if(any(str_detect(arm_lower,"control|ctrl|standard|soc|compression|placebo"))) arm_names[which(str_detect(arm_lower,"control|ctrl|standard|soc|compression|placebo"))[1]] else setdiff(arm_names, arm_blcc)[1]
cat("arm_ctrl=", arm_ctrl, " arm_blcc=", arm_blcc, "\n")

# ---- DiD table
did_tbl <- NULL
if(!is.na(did_f)) {
  did_raw <- read_csv(did_f, show_col_types=FALSE, progress=FALSE)
  path_c <- pick(did_raw, c("pathway","Pathway","set","gene_set","gs_name","term","NAME"))
  eff_c  <- pick(did_raw, c("effect","DiD","Did","estimate","coef","delta","diff"))
  fdr_c  <- pick(did_raw, c("FDR","adj.P.Val","padj","qvalue","q.val"))
  if(any(is.na(c(path_c, eff_c, fdr_c)))) {
    did_tbl <- NULL
  } else {
    did_tbl <- did_raw %>% transmute(pathway=as.character(.data[[path_c]]),
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
    ) %>% mutate(FDR=p.adjust(p, method="BH")) %>% select(-p)
}

did_tbl <- did_tbl %>%
  mutate(category=vapply(pathway, cat_path, character(1)),
         neglogFDR=-log10(pmax(FDR, 1e-300))) %>%
  filter(pathway %in% rownames(mat_ord))

write_csv(did_tbl, "figs/Fig4_tables/DiD_pathway_table.csv")

# ---- Panel A heatmap top 40
topM <- 40
top_paths <- did_tbl %>% arrange(FDR, desc(abs(effect))) %>% slice_head(n=topM) %>% pull(pathway)
hm <- mat_ord[top_paths, , drop=FALSE]
hm_z <- t(scale(t(hm))); hm_z[is.na(hm_z)] <- 0

ann <- meta_ord %>% transmute(arm=as.character(arm), time=as.character(time2)) %>% as.data.frame()
rownames(ann) <- meta_ord$sample
arm_cols <- setNames(colorRampPalette(brewer.pal(8,"Set2"))(length(unique(ann$arm))), unique(ann$arm))
time_cols <- setNames(colorRampPalette(brewer.pal(8,"Pastel1"))(length(unique(ann$time))), unique(ann$time))
grp_blocks <- meta_ord %>% mutate(block=paste0(as.character(arm), "_", as.character(time2)))
gaps <- cumsum(table(grp_blocks$block))

ph <- pheatmap(hm_z, annotation_col=ann,
               annotation_colors=list(arm=arm_cols, time=time_cols),
               show_colnames=FALSE, show_rownames=TRUE,
               cluster_cols=FALSE, cluster_rows=TRUE,
               border_color=NA, gaps_col=gaps, legend=FALSE)
pA <- patchwork::wrap_elements(ph[[4]])

# ---- Panel B pathway-volcano
did_tbl <- did_tbl %>%
  mutate(color_group = case_when(
    category %in% c("Repair/Keratin","Repair/Growth") ~ "Repair",
    category %in% c("Inflammation/Immune") ~ "Inflammation",
    category %in% c("ECM/Wnt") ~ "ECM/Wnt",
    TRUE ~ "Other"
  ))
col_map <- c(Repair="#E64B35", Inflammation="#4DBBD5", `ECM/Wnt`="#00A087", Other="grey70")
labB <- did_tbl %>% filter(FDR < 0.05) %>% arrange(desc(abs(effect)*neglogFDR)) %>% slice_head(n=10)

pB <- ggplot(did_tbl, aes(effect, neglogFDR)) +
  geom_point(aes(color=color_group), size=2.0, alpha=0.9) +
  geom_vline(xintercept=0, linewidth=0.4) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth=0.35) +
  ggrepel::geom_text_repel(data=labB, aes(label=pathway),
                           size=3, box.padding=0.25, point.padding=0.2, segment.size=0.25,
                           max.overlaps=Inf) +
  scale_color_manual(values=col_map) +
  theme_classic(base_size=12) +
  theme(legend.position="top", legend.title=element_blank()) +
  labs(x="DiD effect size (ΔBLCC − ΔControl)", y="-log10(FDR)")

# ---- Panel C patient-level paired lines (5–8 pathways)
prefs <- c("WNT","KERATIN","VEGF","ANGIO","INFLAM","INTERFERON","TNF","ECM","EXTRACELLULAR MATRIX","EPIDERM","TGF","FGF")
pick_core <- function(paths, prefs){
  out <- character(0)
  for(k in prefs){
    hit <- paths[str_detect(toupper(paths), k)]
    if(length(hit)>0) out <- c(out, hit[1])
  }
  unique(out)
}
core <- pick_core(did_tbl$pathway, prefs)
core <- core[core %in% delta_long$pathway]
if(length(core) < 5){
  core <- did_tbl %>% arrange(FDR, desc(abs(effect))) %>% slice_head(n=8) %>% pull(pathway)
} else core <- core[1:min(8, length(core))]

delta_core <- wide %>%
  pivot_longer(cols = -c(sample, arm, time2, patient, time), names_to="pathway", values_to="score") %>%
  filter(pathway %in% core) %>%
  mutate(time2=factor(time2, levels=c("w0","w1")),
         arm=factor(as.character(arm), levels=c(arm_ctrl, arm_blcc)))

pC <- ggplot(delta_core, aes(time2, score, group=patient, color=arm)) +
  geom_line(linewidth=0.7, alpha=0.85) +
  geom_point(size=1.9, alpha=0.95) +
  geom_boxplot(aes(group=interaction(arm,time2)),
               width=0.45, alpha=0.12, outlier.shape=NA, linewidth=0.4) +
  facet_wrap(~pathway, ncol=4, scales="free_y") +
  scale_color_manual(values=setNames(c("#4DBBD5","#E64B35"), c(arm_ctrl, arm_blcc))) +
  theme_classic(base_size=12) +
  theme(strip.background=element_blank(),
        strip.text=element_text(face="bold", size=9),
        legend.position="top", legend.title=element_blank()) +
  labs(x="", y="GSVA score")

# ---- Panel D LOPO robustness
patients <- sort(unique(meta_ord$patient))
top_robust <- did_tbl %>% arrange(FDR, desc(abs(effect))) %>% slice_head(n=15) %>% pull(pathway)

lopo_prop <- lapply(top_robust, function(path){
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

rob_tbl <- tibble(pathway=top_robust, prop_consistent=unlist(lopo_prop)) %>%
  arrange(desc(prop_consistent))
write_csv(rob_tbl, "figs/Fig4_tables/LOPO_robustness_table.csv")

pD <- ggplot(rob_tbl, aes(x=reorder(pathway, prop_consistent), y=prop_consistent)) +
  geom_col(width=0.75, fill="grey35") +
  coord_flip() +
  theme_classic(base_size=12) +
  theme(axis.title.y = element_blank()) +
  labs(y="Proportion consistent (LOPO)")

# ---- Assemble & Save (no cairo dependency)
fig4 <- pA / (pB | pC) / pD + plot_layout(heights=c(1.05, 1.25, 0.85))

ggsave("figs/Fig4_pathway_reprogramming.pdf", fig4, width=14, height=12, device="pdf")
ggsave("figs/Fig4_pathway_reprogramming.png", fig4, width=14, height=12, dpi=300)

cat("✅ Saved figs/Fig4_pathway_reprogramming.(png/pdf)\n")
cat("Core pathways used (Panel C):\n"); cat(paste(core, collapse="\n"), "\n")
