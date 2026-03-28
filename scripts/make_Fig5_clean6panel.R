suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

theme_jtm <- theme_classic(base_size=12) +
  theme(
    axis.text.y = element_text(size=9),
    axis.text.x = element_text(size=9),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face="bold")
  )
col_group <- c("Control"="#2E86C1","BLCC"="#C0392B")

# ---- meta ----
meta <- read_csv("data/sample_sheet.csv", show_col_types=FALSE)
sample_col <- names(meta)[grep("sample|gsm", names(meta), ignore.case=TRUE)][1]
arm_col    <- names(meta)[grep("arm|group|treat", names(meta), ignore.case=TRUE)][1]
meta2 <- meta %>% transmute(sample=.data[[sample_col]], group=.data[[arm_col]])

# ---- helper: mean ± 95%CI ----
summ_ci <- function(df){
  df %>%
    group_by(feature, group) %>%
    summarise(
      mean = mean(score, na.rm=TRUE),
      se   = sd(score, na.rm=TRUE)/sqrt(sum(is.finite(score))),
      .groups="drop"
    ) %>%
    mutate(ci=1.96*se, xmin=mean-ci, xmax=mean+ci)
}

plot_lollipop <- function(df_long, features, xlab, title){
  df <- df_long %>% filter(feature %in% features) %>%
    mutate(feature=factor(feature, levels=rev(features)))

  ss <- summ_ci(df)

  ggplot() +
    geom_point(data=df, aes(x=score, y=feature, color=group),
               position=position_jitterdodge(jitter.width=0.03, dodge.width=0.55),
               alpha=0.30, size=1.1) +
    geom_errorbarh(data=ss, aes(y=feature, xmin=xmin, xmax=xmax, color=group),
                   position=position_dodge(width=0.55), height=0.18, linewidth=0.7) +
    geom_point(data=ss, aes(x=mean, y=feature, color=group),
               position=position_dodge(width=0.55), size=2.4) +
    scale_color_manual(values=col_group) +
    labs(x=xlab, y=NULL, title=title) +
    theme_jtm
}

# ---- GSVA long ----
gsva <- read_csv("results/GSVA_Hallmark_scores_wide.csv", show_col_types=FALSE) %>%
  pivot_longer(-feature, names_to="sample", values_to="score") %>%
  inner_join(meta2, by="sample")

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

pA <- plot_lollipop(gsva, repair_prog, "GSVA score", "A  Repair-associated programs")
pB <- plot_lollipop(gsva, inflam_prog, "GSVA score", "B  Inflammatory / stress programs")

ggsave("figs/Fig5A_clean.pdf", pA, width=5.2, height=3.6)
ggsave("figs/Fig5B_clean.pdf", pB, width=5.2, height=3.6)

# ---- xCell long ----
xcell <- read_csv("results/xCell_scores_wide.csv", show_col_types=FALSE) %>%
  pivot_longer(-feature, names_to="sample", values_to="score") %>%
  inner_join(meta2, by="sample")

immune_cells  <- c("Neutrophils","Macrophages","CD8+ T-cells","Tregs")
stromal_cells <- c("Fibroblasts","Endothelial cells","Keratinocytes")

pC <- plot_lollipop(xcell, immune_cells,  "xCell score", "C  Immune cell infiltration")
pD <- plot_lollipop(xcell, stromal_cells, "xCell score", "D  Stromal / epithelial cells")

ggsave("figs/Fig5C_clean.pdf", pC, width=5.2, height=3.2)
ggsave("figs/Fig5D_clean.pdf", pD, width=5.2, height=3.2)

# ---- Panel E: group-average composition (clean, no sample stacking) ----
comp_cells <- c("Neutrophils","Macrophages","CD8+ T-cells","Fibroblasts","Endothelial cells")
dfE <- xcell %>%
  filter(feature %in% comp_cells) %>%
  group_by(group, sample) %>%
  mutate(prop = score / sum(score)) %>%
  ungroup() %>%
  group_by(group, feature) %>%
  summarise(prop=mean(prop), .groups="drop") %>%
  mutate(group=factor(group, levels=c("Control","BLCC")))

pE <- ggplot(dfE, aes(x=prop, y=group, fill=feature)) +
  geom_bar(stat="identity", width=0.6) +
  scale_x_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(x="Average relative composition", y=NULL, title="E  Microenvironment composition") +
  theme_jtm +
  theme(legend.position="right")

ggsave("figs/Fig5E_clean.pdf", pE, width=7.2, height=2.3)

# ---- Panel F: DiD + LOPO ----
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
  labs(x="Difference-in-differences effect", y=NULL, title="F  Effect & LOPO stability") +
  theme_jtm

ggsave("figs/Fig5F_clean.pdf", pF, width=5.2, height=3.6)

# ---- Combine (clean layout) ----
fig5 <- (pA | pB) / (pC | pD) / (pE | pF) + plot_layout(heights=c(1,1,0.85))
ggsave("figs/Fig5_clean6panel.pdf", fig5, width=10, height=12)
ggsave("figs/Fig5_clean6panel.png", fig5, width=10, height=12, dpi=300)

cat("✅ DONE. Outputs in figs/: Fig5_clean6panel.pdf/png + panels A-F pdf\n")
