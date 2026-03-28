# ============================================================
# Fig 5 Pretty Layout (A-F) | Robust + auto-install + logging
# ============================================================

# ---- auto install missing packages ----
need <- c("readr","dplyr","tidyr","ggplot2","patchwork","scales")
to_install <- need[!vapply(need, requireNamespace, logical(1), quietly=TRUE)]
if(length(to_install) > 0){
  message("Installing missing packages: ", paste(to_install, collapse=", "))
  install.packages(to_install, repos="https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(scales)
})

# ---- paths check ----
stopifnot(file.exists("data/sample_sheet.csv"))
stopifnot(file.exists("results/GSVA_Hallmark_scores_wide.csv"))
stopifnot(file.exists("results/xCell_scores_wide.csv"))
stopifnot(file.exists("results/Fig5D_GSVA_DiD_LOPO.csv"))
dir.create("figs", showWarnings=FALSE)

cat("▶ Loaded packages + checked input files OK\n")

# -------------------------
# global theme
# -------------------------
theme_jtm <- theme_classic(base_size=12) +
  theme(
    axis.text.y = element_text(size=9),
    axis.text.x = element_text(size=9),
    legend.title = element_blank(),
    plot.title = element_text(face="bold", size=12),
    plot.margin = margin(6, 10, 6, 6)
  )

# group colors (fixed, consistent)
col_group <- c("Control"="#2E86C1","BLCC"="#C0392B")

# -------------------------
# meta
# -------------------------
meta <- read_csv("data/sample_sheet.csv", show_col_types=FALSE)
sample_col <- names(meta)[grep("sample|gsm", names(meta), ignore.case=TRUE)][1]
arm_col    <- names(meta)[grep("arm|group|treat", names(meta), ignore.case=TRUE)][1]
if(is.na(sample_col) || is.na(arm_col)){
  stop("❌ sample_sheet.csv must contain columns like sample/GSM and arm/group/treat.\n",
       "Found columns: ", paste(names(meta), collapse=", "))
}
meta2 <- meta %>% transmute(sample=.data[[sample_col]], group=.data[[arm_col]])
cat("▶ Meta columns: sample=", sample_col, " group=", arm_col, "\n")

# normalize group names (optional but helps legend consistency)
meta2 <- meta2 %>%
  mutate(group = as.character(group)) %>%
  mutate(group = ifelse(tolower(group) %in% c("control","ctrl"), "Control", group)) %>%
  mutate(group = ifelse(grepl("blcc|apligraf", tolower(group)), "BLCC", group))

# -------------------------
# helper: clean box + light points, NO legend per panel
# -------------------------
plot_box_clean <- function(df, features, xlab, title){
  df2 <- df %>% filter(feature %in% features) %>%
    mutate(feature=factor(feature, levels=rev(features)))

  ggplot(df2, aes(x=score, y=feature, fill=group)) +
    geom_boxplot(width=0.55, outlier.shape=NA, alpha=0.85,
                 position=position_dodge(width=0.68)) +
    geom_point(aes(color=group),
               position=position_jitterdodge(jitter.width=0.06, dodge.width=0.68),
               alpha=0.30, size=1.05, show.legend=FALSE) +
    scale_fill_manual(values=col_group, limits=c("BLCC","Control")) +
    scale_color_manual(values=col_group, limits=c("BLCC","Control")) +
    labs(x=xlab, y=NULL, title=title) +
    theme_jtm +
    theme(legend.position="none")
}

# -------------------------
# GSVA long
# -------------------------
cat("▶ Loading GSVA wide → long...\n")
gsva <- read_csv("results/GSVA_Hallmark_scores_wide.csv", show_col_types=FALSE) %>%
  pivot_longer(-feature, names_to="sample", values_to="score") %>%
  inner_join(meta2, by="sample")

# choose features (you can edit)
repair_prog <- c(
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_TGF_BETA_SIGNALING"
)

inflam_prog <- c(
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_UV_RESPONSE_DN"
)

pA <- plot_box_clean(gsva, repair_prog, "GSVA score", "Repair-associated programs")
pB <- plot_box_clean(gsva, inflam_prog, "GSVA score", "Inflammatory / stress programs")

ggsave("figs/Fig5A_pretty.pdf", pA, width=5.1, height=3.6)
ggsave("figs/Fig5B_pretty.pdf", pB, width=5.1, height=3.6)
cat("✔ Saved Fig5A_pretty / Fig5B_pretty\n")

# -------------------------
# xCell long
# -------------------------
cat("▶ Loading xCell wide → long...\n")
xcell <- read_csv("results/xCell_scores_wide.csv", show_col_types=FALSE) %>%
  pivot_longer(-feature, names_to="sample", values_to="score") %>%
  inner_join(meta2, by="sample")

immune_cells  <- c("Neutrophils","Macrophages","CD8+ T-cells","Tregs")
stromal_cells <- c("Fibroblasts","Endothelial cells","Keratinocytes")

pC <- plot_box_clean(xcell, immune_cells,  "xCell score", "Immune cell infiltration")
pD <- plot_box_clean(xcell, stromal_cells, "xCell score", "Stromal / epithelial cells")

ggsave("figs/Fig5C_pretty.pdf", pC, width=5.1, height=3.2)
ggsave("figs/Fig5D_pretty.pdf", pD, width=5.1, height=3.2)
cat("✔ Saved Fig5C_pretty / Fig5D_pretty\n")

# -------------------------
# Panel E: composition (group-average, HORIZONTAL, wide)
# -------------------------
cat("▶ Building Panel E (group-average composition)...\n")
comp_cells <- c("Neutrophils","Macrophages","CD8+ T-cells","Fibroblasts","Endothelial cells")

dfE <- xcell %>%
  filter(feature %in% comp_cells) %>%
  group_by(group, sample) %>%
  mutate(prop = score / sum(score)) %>%
  ungroup() %>%
  group_by(group, feature) %>%
  summarise(prop = mean(prop), .groups="drop") %>%
  mutate(group = factor(group, levels=c("Control","BLCC")))

pE <- ggplot(dfE, aes(x=prop, y=group, fill=feature)) +
  geom_bar(stat="identity", width=0.55) +
  scale_x_continuous(labels=percent_format(accuracy=1)) +
  labs(x="Average relative composition", y=NULL, title="Microenvironment composition") +
  theme_jtm +
  theme(legend.position="right")

ggsave("figs/Fig5E_pretty.pdf", pE, width=7.2, height=2.4)
cat("✔ Saved Fig5E_pretty\n")

# -------------------------
# Panel F: Effect + LOPO (margin-safe)
# -------------------------
cat("▶ Building Panel F (DiD + LOPO)...\n")
dfF <- read_csv("results/Fig5D_GSVA_DiD_LOPO.csv", show_col_types=FALSE) %>%
  arrange(desc(abs(DiD))) %>% slice(1:12) %>%
  mutate(feature=factor(feature, levels=rev(feature)),
         label=paste0(consistent,"/",N))

pF <- ggplot(dfF, aes(x=DiD, y=feature, fill=DiD>0)) +
  geom_col(width=0.7) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50") +
  geom_text(aes(label=label),
            hjust=ifelse(dfF$DiD>0,-0.15,1.15), size=3) +
  scale_fill_manual(values=c("TRUE"="#C0392B","FALSE"="#2E86C1"), guide="none") +
  labs(x="DiD effect", y=NULL, title="Effect summary & LOPO stability") +
  theme_jtm +
  theme(plot.margin = margin(6, 20, 6, 6)) +
  coord_cartesian(xlim=c(min(dfF$DiD)*1.35, max(dfF$DiD)*1.35))

ggsave("figs/Fig5F_pretty.pdf", pF, width=5.1, height=3.8)
cat("✔ Saved Fig5F_pretty\n")

# -------------------------
# One shared legend (ONLY)
# -------------------------
cat("▶ Creating shared legend...\n")
leg <- ggplot(meta2, aes(x=group, fill=group)) +
  geom_bar() +
  scale_fill_manual(values=col_group, limits=c("BLCC","Control")) +
  theme_void() +
  theme(legend.position="bottom", legend.title=element_blank())

# -------------------------
# Pretty layout (key!)
# Top: 2x2
# Bottom: E wide, F narrow
# Legend: one at bottom
# Patchwork auto tag A-F
# -------------------------
cat("▶ Combining layout...\n")

top <- (pA | pB) / (pC | pD) + plot_layout(heights=c(1,1))
bottom <- (pE | pF) + plot_layout(widths=c(1.35, 1))

fig5 <- (top / bottom / wrap_elements(leg)) +
  plot_annotation(tag_levels="A") +
  plot_layout(heights=c(1.9, 1.05, 0.15))

ggsave("figs/Fig5_pretty.pdf", fig5, width=11.5, height=12)
ggsave("figs/Fig5_pretty.png", fig5, width=11.5, height=12, dpi=300)

cat("✅ DONE. Output:\n")
cat("  figs/Fig5_pretty.pdf\n")
cat("  figs/Fig5_pretty.png\n")
