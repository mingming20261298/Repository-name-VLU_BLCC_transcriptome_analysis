# ============================================================
# Fig5 (JTM-like): A/B/C mean±CI + faint points; D/E diverging COLOR bars
# Output: figs/Fig5_JTM_like.pdf/png (ONE FIG ONLY)
# Uses patient-level delta (w1-w0) if patient+time exist in meta.
# ============================================================

need <- c("readr","dplyr","tidyr","stringr","ggplot2","patchwork","scales")
to_install <- need[!vapply(need, requireNamespace, logical(1), quietly=TRUE)]
if(length(to_install)>0) install.packages(to_install, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(patchwork); library(scales)
})

dir.create("figs", showWarnings=FALSE)

# -----------------------------
# helpers
# -----------------------------
pick_col <- function(df, patterns){
  nms <- names(df)
  for(p in patterns){
    hit <- grep(p, nms, ignore.case=TRUE, value=TRUE)
    if(length(hit)>0) return(hit[1])
  }
  NA_character_
}

norm_time <- function(x){
  x <- tolower(as.character(x))
  x <- str_replace_all(x, "week", "w")
  x <- str_replace_all(x, "[^a-z0-9]+", "")
  x <- ifelse(x %in% c("0","w0","t0","baseline","pre"), "w0", x)
  x <- ifelse(x %in% c("1","w1","t1","post"), "w1", x)
  x
}

mean_ci_df <- function(x){
  x <- x[is.finite(x)]
  n <- length(x)
  if(n < 2) return(data.frame(y=mean(x), ymin=NA_real_, ymax=NA_real_))
  m <- mean(x)
  se <- sd(x)/sqrt(n)
  ci <- qt(0.975, df=n-1) * se
  data.frame(y=m, ymin=m-ci, ymax=m+ci)
}

# -----------------------------
# style
# -----------------------------
col_group <- c(Control="#2E86C1", BLCC="#C0392B")
col_dir   <- c(NEG="#2E86C1", POS="#C0392B")   # diverging bars: blue / red

theme_jtm <- theme_classic(base_size=12) +
  theme(
    plot.title = element_text(face="bold", size=12),
    axis.title = element_text(size=11),
    axis.text.y = element_text(size=9),
    axis.text.x = element_text(size=9),
    legend.title = element_blank(),
    legend.position = "right",
    plot.margin = margin(6, 10, 6, 6)
  )

# -----------------------------
# load meta
# -----------------------------
meta <- read_csv("data/sample_sheet.csv", show_col_types=FALSE)
sample_col <- pick_col(meta, c("^sample$","sample","gsm"))
arm_col    <- pick_col(meta, c("^arm$","arm","group","treat"))
time_col   <- pick_col(meta, c("^time$","time","week","timepoint"))
pt_col     <- pick_col(meta, c("^patient$","patient","subject","id"))

if(any(is.na(c(sample_col, arm_col)))){
  stop("❌ sample_sheet.csv needs sample/GSM and arm/group/treat columns. Found: ",
       paste(names(meta), collapse=", "))
}

meta2 <- meta %>%
  transmute(
    sample = as.character(.data[[sample_col]]),
    arm_raw = as.character(.data[[arm_col]]),
    time_raw = if(!is.na(time_col)) as.character(.data[[time_col]]) else NA_character_,
    patient_raw = if(!is.na(pt_col)) as.character(.data[[pt_col]]) else NA_character_
  ) %>%
  mutate(
    group = ifelse(grepl("blcc|apligraf", tolower(arm_raw)), "BLCC", "Control"),
    time  = ifelse(is.na(time_raw), NA_character_, norm_time(time_raw)),
    patient = ifelse(is.na(patient_raw), NA_character_, as.character(patient_raw))
  )

have_delta <- all(!is.na(meta2$time)) && all(meta2$time %in% c("w0","w1")) && all(!is.na(meta2$patient))

cat("▶ patient-level delta mode:", have_delta, "\n")

# -----------------------------
# load GSVA & xCell
# -----------------------------
stopifnot(file.exists("results/GSVA_Hallmark_scores_wide.csv"))
stopifnot(file.exists("results/xCell_scores_wide.csv"))

gsva_w <- read_csv("results/GSVA_Hallmark_scores_wide.csv", show_col_types=FALSE)
xcell_w <- read_csv("results/xCell_scores_wide.csv", show_col_types=FALSE)

gsva_l <- gsva_w %>% pivot_longer(-feature, names_to="sample", values_to="score") %>%
  inner_join(meta2 %>% select(sample, group, time, patient), by="sample")
xcell_l <- xcell_w %>% pivot_longer(-feature, names_to="sample", values_to="score") %>%
  inner_join(meta2 %>% select(sample, group, time, patient), by="sample")

make_delta <- function(df_long){
  if(have_delta){
    df_long %>%
      filter(time %in% c("w0","w1")) %>%
      select(group, patient, feature, time, score) %>%
      pivot_wider(names_from=time, values_from=score) %>%
      mutate(delta = w1 - w0) %>%
      filter(is.finite(delta)) %>%
      transmute(group, id=patient, feature, value=delta)
  } else {
    df_long %>%
      transmute(group, id=sample, feature, value=score) %>%
      filter(is.finite(value))
  }
}

gsva_d <- make_delta(gsva_l)
xcell_d <- make_delta(xcell_l)

# -----------------------------
# choose features for A/B/C
# -----------------------------
choose_top <- function(df, k=6){
  tmp <- df %>% group_by(feature, group) %>% summarise(m=mean(value), .groups="drop") %>%
    pivot_wider(names_from=group, values_from=m) %>%
    mutate(diff = BLCC - Control) %>% arrange(desc(abs(diff)))
  tmp$feature[1:min(k, nrow(tmp))]
}

A_prog <- c(
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_KERATINIZATION",
  "HALLMARK_WOUND_HEALING"
)
A_prog <- A_prog[A_prog %in% unique(gsva_d$feature)]
if(length(A_prog) < 4) A_prog <- choose_top(gsva_d, k=6)

B_prog <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_COMPLEMENT"
)
B_prog <- B_prog[B_prog %in% unique(gsva_d$feature)]
if(length(B_prog) < 4) B_prog <- choose_top(gsva_d, k=6)

C_cells <- c("Neutrophils","Macrophages","CD8+ T-cells","Tregs","T cells","Fibroblasts","Endothelial cells")
C_cells <- C_cells[C_cells %in% unique(xcell_d$feature)]
if(length(C_cells) < 4) C_cells <- choose_top(xcell_d, k=8)

# -----------------------------
# A/B/C: mean±CI + faint points
# -----------------------------
plot_mean_ci <- function(df, features, title, xlab){
  d <- df %>% filter(feature %in% features) %>%
    mutate(feature = factor(feature, levels=rev(features)))

  xr <- range(d$value, finite=TRUE)
  pad <- 0.08 * diff(xr); if(!is.finite(pad) || pad==0) pad <- 0.2

  ggplot(d, aes(x=value, y=feature)) +
    geom_vline(xintercept=0, linetype="dashed", color="grey60") +
    geom_point(color="black", alpha=0.22, size=1.1,
               position=position_jitter(height=0.12, width=0)) +
    stat_summary(aes(color=group),
                 fun.data = mean_ci_df,
                 geom="pointrange",
                 position=position_dodge(width=0.55),
                 size=0.55, linewidth=0.55) +
    stat_summary(aes(color=group),
                 fun=mean, geom="point",
                 position=position_dodge(width=0.55),
                 size=2.2) +
    scale_color_manual(values=col_group) +
    coord_cartesian(xlim=c(xr[1]-pad, xr[2]+pad)) +
    labs(title=title, x=xlab, y=NULL) +
    theme_jtm
}

xlab_delta <- if(have_delta) "Δ score (week1 − week0)" else "Score"

pA <- plot_mean_ci(gsva_d, A_prog, "A  Repair-associated programs", xlab_delta)
pB <- plot_mean_ci(gsva_d, B_prog, "B  Inflammatory programs", xlab_delta)
pC <- plot_mean_ci(xcell_d, C_cells, "C  Immune / stromal cell signals", xlab_delta)

# -----------------------------
# D/E: diverging COLOR bars
# -----------------------------
stopifnot(file.exists("results/Fig5D_GSVA_DiD_LOPO.csv"))
did_g <- read_csv("results/Fig5D_GSVA_DiD_LOPO.csv", show_col_types=FALSE)

eff_col <- pick_col(did_g, c("^did$","did","effect","estimate","logfc"))
feat_col <- pick_col(did_g, c("^feature$","feature","pathway"))

did_g2 <- did_g %>%
  transmute(feature = as.character(.data[[feat_col]]),
            effect  = as.numeric(.data[[eff_col]]),
            consistent = if("consistent" %in% names(did_g)) as.integer(consistent) else NA_integer_,
            N = if("N" %in% names(did_g)) as.integer(N) else NA_integer_) %>%
  filter(is.finite(effect)) %>%
  mutate(dir = ifelse(effect >= 0, "POS", "NEG"),
         label = ifelse(!is.na(consistent) & !is.na(N), paste0(consistent,"/",N), NA_character_)) %>%
  arrange(desc(abs(effect))) %>%
  slice(1:18) %>%
  mutate(feature=factor(feature, levels=rev(feature)))

pD <- ggplot(did_g2, aes(x=effect, y=feature, fill=dir)) +
  geom_col(width=0.72) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50") +
  geom_text(aes(label=label),
            hjust=ifelse(did_g2$effect>0,-0.15,1.15),
            size=3, na.rm=TRUE) +
  scale_fill_manual(values=col_dir, guide="none") +
  labs(title="D  Program-level DiD summary (top effects)",
       x="DiD effect (BLCC − Control)", y=NULL) +
  theme_jtm +
  theme(plot.margin = margin(6, 18, 6, 6)) +
  coord_cartesian(xlim=c(min(did_g2$effect)*1.35, max(did_g2$effect)*1.35))

cell_eff <- xcell_d %>%
  group_by(feature, group) %>%
  summarise(m=mean(value), .groups="drop") %>%
  pivot_wider(names_from=group, values_from=m) %>%
  mutate(effect = BLCC - Control) %>%
  filter(is.finite(effect)) %>%
  mutate(dir = ifelse(effect >= 0, "POS", "NEG")) %>%
  arrange(desc(abs(effect))) %>%
  slice(1:18) %>%
  mutate(feature=factor(feature, levels=rev(feature)))

pE <- ggplot(cell_eff, aes(x=effect, y=feature, fill=dir)) +
  geom_col(width=0.72) +
  geom_vline(xintercept=0, linetype="dashed", color="grey50") +
  scale_fill_manual(values=col_dir, guide="none") +
  labs(title="E  Cell-signal shift summary (top effects)",
       x=if(have_delta) "Δ(week1−week0): BLCC − Control" else "BLCC − Control", y=NULL) +
  theme_jtm +
  theme(plot.margin = margin(6, 12, 6, 6)) +
  coord_cartesian(xlim=c(min(cell_eff$effect)*1.35, max(cell_eff$effect)*1.35))

# -----------------------------
# Layout: 3 stacked + 2 bottom (like reference)
# -----------------------------
fig5 <- (pA / pB / pC / (pD | pE)) +
  plot_layout(heights=c(1.05, 1.05, 1.05, 1.15), widths=c(1,1))

ggsave("figs/Fig5_JTM_like.pdf", fig5, width=12.5, height=12)
ggsave("figs/Fig5_JTM_like.png", fig5, width=12.5, height=12, dpi=350)

cat("✅ Saved:\n  figs/Fig5_JTM_like.pdf\n  figs/Fig5_JTM_like.png\n")
