# ============================
# Fig5 FINAL — FIXED LAYOUT
# Single output only
# ============================

need <- c("readr","dplyr","tidyr","ggplot2","patchwork","scales")
to_install <- need[!vapply(need, requireNamespace, logical(1), quietly=TRUE)]
if(length(to_install)>0) install.packages(to_install, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(scales)
})

# ---- style ----
theme_jtm <- theme_classic(base_size=12) +
  theme(
    axis.text.y = element_text(size=9),
    axis.text.x = element_text(size=9),
    plot.title  = element_text(face="bold", size=12),
    legend.title = element_blank(),
    plot.margin = margin(6, 10, 6, 6)
  )

col_group <- c(Control="#2E86C1", BLCC="#C0392B")

# ---- meta ----
meta <- read_csv("data/sample_sheet.csv", show_col_types=FALSE)
sample_col <- names(meta)[grep("sample|gsm", names(meta), ignore.case=TRUE)][1]
arm_col    <- names(meta)[grep("arm|group|treat", names(meta), ignore.case=TRUE)][1]

meta2 <- meta %>%
  transmute(sample=.data[[sample_col]],
            group=.data[[arm_col]]) %>%
  mutate(group=ifelse(grepl("blcc|apligraf",tolower(group)),"BLCC","Control"))

# ---- helper ----
plot_box <- function(df, feats, title, xlab){
  df %>% filter(feature %in% feats) %>%
    mutate(feature=factor(feature, levels=rev(feats))) %>%
    ggplot(aes(score, feature, fill=group)) +
    geom_boxplot(width=0.55, outlier.shape=NA,
                 position=position_dodge(0.65)) +
    geom_point(position=position_jitterdodge(0.05,0.65),
               size=1, alpha=0.3, show.legend=FALSE) +
    scale_fill_manual(values=col_group) +
    labs(title=title, x=xlab, y=NULL) +
    theme_jtm +
    theme(legend.position="none")
}

# ---- GSVA ----
gsva <- read_csv("results/GSVA_Hallmark_scores_wide.csv", show_col_types=FALSE) %>%
  pivot_longer(-feature, names_to="sample", values_to="score") %>%
  inner_join(meta2, by="sample")

pA <- plot_box(gsva,
  c("HALLMARK_ANGIOGENESIS",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_TGF_BETA_SIGNALING"),
  "Repair-associated programs","GSVA score")

pB <- plot_box(gsva,
  c("HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MYC_TARGETS_V2",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_IL2_STAT5_SIGNALING"),
  "Inflammatory / stress programs","GSVA score")

# ---- xCell ----
xcell <- read_csv("results/xCell_scores_wide.csv", show_col_types=FALSE) %>%
  pivot_longer(-feature, names_to="sample", values_to="score") %>%
  inner_join(meta2, by="sample")

pC <- plot_box(xcell,
  c("Neutrophils","Macrophages","CD8+ T-cells","Tregs"),
  "Immune cell infiltration","xCell score")

pD <- plot_box(xcell,
  c("Fibroblasts","Endothelial cells","Keratinocytes"),
  "Stromal / epithelial cells","xCell score")

# ---- E: composition (group mean only) ----
dfE <- xcell %>%
  filter(feature %in% c("Neutrophils","Macrophages","CD8+ T-cells",
                        "Fibroblasts","Endothelial cells")) %>%
  group_by(group, sample) %>%
  mutate(p = score/sum(score)) %>%
  ungroup() %>%
  group_by(group, feature) %>%
  summarise(p=mean(p), .groups="drop")

pE <- ggplot(dfE, aes(p, group, fill=feature)) +
  geom_col(width=0.6) +
  scale_x_continuous(labels=percent) +
  labs(title="Microenvironment composition",
       x="Relative proportion", y=NULL) +
  theme_jtm

# ---- F: DiD + LOPO + LEGEND ----
dfF <- read_csv("results/Fig5D_GSVA_DiD_LOPO.csv", show_col_types=FALSE) %>%
  arrange(desc(abs(DiD))) %>% slice(1:10) %>%
  mutate(feature=factor(feature, levels=rev(feature)))

pF <- ggplot(dfF, aes(DiD, feature, fill=DiD>0)) +
  geom_col(width=0.7) +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_fill_manual(values=c("#2E86C1","#C0392B"),
                    labels=c("Control ↓","BLCC ↑")) +
  labs(title="Effect summary (DiD + LOPO)",
       x="Difference-in-differences", y=NULL) +
  theme_jtm +
  theme(legend.position="bottom")

# ---- FINAL LAYOUT: fixed 2×3 ----
fig5 <- (pA | pB | pC) /
        (pD | pE | pF) +
        plot_annotation(tag_levels="A")

ggsave("figs/Fig5_FINAL.pdf", fig5, width=14, height=9)
ggsave("figs/Fig5_FINAL.png", fig5, width=14, height=9, dpi=300)

cat("✅ Fig5 saved:\n  figs/Fig5_FINAL.pdf\n  figs/Fig5_FINAL.png\n")
