suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr); library(stringr); library(scales)
})
dir.create("~/JTM_GSE84571/figs", showWarnings=FALSE, recursive=TRUE)
dir.create("~/JTM_GSE84571/results", showWarnings=FALSE, recursive=TRUE)
save_plot <- function(p, stem, w=8, h=5, dpi=300){
  ggsave(file.path("~/JTM_GSE84571/figs", paste0(stem,".pdf")), p, width=w, height=h, device=cairo_pdf)
  ggsave(file.path("~/JTM_GSE84571/figs", paste0(stem,".png")), p, width=w, height=h, dpi=dpi)
  message("Saved: ", stem)
}
pick_col <- function(df, candidates){ hit <- candidates[candidates %in% colnames(df)]; if(length(hit)==0) NA_character_ else hit[1] }
to_num <- function(x) suppressWarnings(as.numeric(x))
