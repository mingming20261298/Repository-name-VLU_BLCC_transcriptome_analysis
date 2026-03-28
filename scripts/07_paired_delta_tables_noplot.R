suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
})
dir.create("results/paired", showWarnings = FALSE, recursive = TRUE)

tab_int <- readr::read_csv("results/Table_interaction_all.csv", show_col_types = FALSE)

expr_candidates <- c(
  "data/GSE84571_expr_rma_log2.csv",
  "data/GSE84571_expr_mas5_log2.csv",
  "data/GSE84571_expr_log2.csv",
  "data/expr_log2.csv"
)
expr_path <- expr_candidates[file.exists(expr_candidates)][1]
if(is.na(expr_path)) stop("Cannot find expression matrix in data/.")
expr <- readr::read_csv(expr_path, show_col_types = FALSE)
expr_mat <- as.matrix(expr[,-1]); rownames(expr_mat) <- expr[[1]]

ss_candidates <- c("data/sample_sheet.csv","data/GSE84571_sample_sheet.csv","data/targets.csv","data/phenodata.csv")
ss_path <- ss_candidates[file.exists(ss_candidates)][1]
if(is.na(ss_path)) stop("Cannot find sample sheet in data/. Put as data/sample_sheet.csv.")
pheno <- readr::read_csv(ss_path, show_col_types = FALSE)
names(pheno) <- tolower(names(pheno))

if(!("sample_id" %in% names(pheno))){
  alt <- intersect(names(pheno), c("sample","sampleid","gsm","geo_accession","array","cel"))
  if(length(alt)==0) stop("sample_id missing.")
  pheno <- pheno %>% rename(sample_id = all_of(alt[1]))
}
if(!("arm" %in% names(pheno))){
  alt <- intersect(names(pheno), c("group","treatment","condition","study_arm"))
  if(length(alt)==0) stop("arm missing.")
  pheno <- pheno %>% rename(arm = all_of(alt[1]))
}
if(!("time" %in% names(pheno))){
  alt <- intersect(names(pheno), c("week","timepoint","tp"))
  if(length(alt)==0) stop("time missing.")
  pheno <- pheno %>% rename(time = all_of(alt[1]))
}
if(!("patient" %in% names(pheno))){
  alt <- intersect(names(pheno), c("subject","patientid","id","pair","individual"))
  if(length(alt)==0) stop("patient missing.")
  pheno <- pheno %>% rename(patient = all_of(alt[1]))
}

pheno <- pheno %>%
  mutate(
    sample_id = as.character(sample_id),
    patient = as.character(patient),
    arm = ifelse(grepl("blcc|apligraf", arm, ignore.case=TRUE), "BLCC",
                 ifelse(grepl("control|std|standard", arm, ignore.case=TRUE), "Control", as.character(arm))),
    time = ifelse(grepl("0", time), "week0",
                  ifelse(grepl("1", time), "week1", as.character(time)))
  )

common <- intersect(colnames(expr_mat), pheno$sample_id)
if(length(common) < 10) stop("Too few matching samples; fix sample_id.")
expr_mat <- expr_mat[, common, drop=FALSE]
pheno <- pheno %>% filter(sample_id %in% common) %>% arrange(match(sample_id, common))
stopifnot(all(pheno$sample_id == colnames(expr_mat)))

if(!("adj.p.val" %in% tolower(names(tab_int)))) stop("Need adj.P.Val in interaction table.")
if(!("probe_id" %in% names(tab_int))) stop("Need probe_id in interaction table.")

# normalize adj.P.Val column name
if(!("adj.P.Val" %in% names(tab_int))){
  nm <- names(tab_int)[tolower(names(tab_int))=="adj.p.val"][1]
  tab_int <- tab_int %>% rename(adj.P.Val = all_of(nm))
}

sig <- tab_int %>% filter(adj.P.Val < 0.10) %>% arrange(adj.P.Val)
if(nrow(sig) == 0) sig <- tab_int %>% arrange(adj.P.Val)
sig23 <- sig %>% slice_head(n=23)
readr::write_csv(sig23, "results/paired/Top23_interaction_genes.csv")

probes <- intersect(sig23$probe_id, rownames(expr_mat))
if(length(probes) < 5) stop("Too few selected probes found in expression matrix. Check probe_id mapping.")

df <- as.data.frame(t(expr_mat[probes, , drop=FALSE]))
df$sample_id <- rownames(df)

long <- df %>%
  left_join(pheno, by="sample_id") %>%
  pivot_longer(cols = all_of(probes), names_to="probe_id", values_to="expr") %>%
  mutate(time = factor(time, levels=c("week0","week1")),
         arm  = factor(arm, levels=c("Control","BLCC")))

delta <- long %>%
  select(patient, arm, time, probe_id, expr) %>%
  pivot_wider(names_from=time, values_from=expr) %>%
  mutate(delta = week1 - week0) %>%
  arrange(probe_id, arm, patient)

readr::write_csv(delta, "results/paired/Top23_patient_level_week0_week1_delta.csv")

delta_sum <- delta %>%
  group_by(probe_id, arm) %>%
  summarise(n=n(), mean_delta=mean(delta, na.rm=TRUE), sd_delta=sd(delta, na.rm=TRUE), .groups="drop")

readr::write_csv(delta_sum, "results/paired/Top23_delta_summary_by_arm.csv")
message("DONE. CSV saved to results/paired/")
