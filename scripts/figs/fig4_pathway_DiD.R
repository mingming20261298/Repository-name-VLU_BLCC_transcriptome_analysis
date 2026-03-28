source("~/JTM_GSE84571/scripts/figs/_common.R")
INFILE <- "~/JTM_GSE84571/results/GSVA_Reactome_DiD.csv"
stopifnot(file.exists(INFILE))
gs <- read_csv(INFILE, show_col_types = FALSE)
path_col <- pick_col(gs, c("pathway","Pathway","set","term","NAME"))
eff_col  <- pick_col(gs, c("estimate","logFC","coef","beta","Effect"))
fdr_col  <- pick_col(gs, c("FDR","adj.P.Val","padj","qval"))
stopifnot(!is.na(path_col), !is.na(eff_col), !is.na(fdr_col))
gs2 <- gs %>% mutate(Pathway=.data[[path_col]], Effect=to_num(.data[[eff_col]]), FDR=to_num(.data[[fdr_col]])) %>% arrange(FDR)
topk <- 10
up <- gs2 %>% filter(Effect>0) %>% slice_head(n=topk)
dn <- gs2 %>% filter(Effect<0) %>% slice_head(n=topk)
plot_df <- bind_rows(up, dn) %>% mutate(Pathway=str_trunc(Pathway, 60)) %>%
  arrange(Effect) %>% mutate(Pathway=factor(Pathway, levels=Pathway))
p <- ggplot(plot_df, aes(x=Effect, y=Pathway)) +
  geom_vline(xintercept=0, linetype=2) +
  geom_point(aes(size=-log10(pmax(FDR,1e-300))), alpha=0.9) +
  theme_bw() + labs(title="Pathway shift (GSVA DiD: arm×time)", x="Effect", y="") +
  guides(size=guide_legend(title="-log10(FDR)"))
save_plot(p, "Fig4_pathway_DiD_top", w=10, h=6)
