suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
})

ip <- function(pkgs){
  for(p in pkgs){
    if(!requireNamespace(p, quietly=TRUE)){
      install.packages(p, repos="https://cloud.r-project.org")
    }
  }
}
ip(c("dplyr","tidyr","readr","tibble","stringr"))

if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
bp <- function(pkgs){
  for(p in pkgs){
    if(!requireNamespace(p, quietly=TRUE)){
      BiocManager::install(p, ask=FALSE, update=FALSE)
    }
  }
}
bp(c("fgsea","msigdbr","GSVA","limma","hgu133plus2.db","AnnotationDbi","org.Hs.eg.db",
     "progeny","dorothea","viper"))

suppressPackageStartupMessages({
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(GSVA)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(org.Hs.eg.db)
  library(progeny)
  library(dorothea)
  library(viper)
  library(stringr)
  library(tibble)
})

dir.create("results/mechanism", showWarnings = FALSE, recursive = TRUE)

tab_int_path <- "results/Table_interaction_all.csv"
stopifnot(file.exists(tab_int_path))
tab_int <- readr::read_csv(tab_int_path, show_col_types = FALSE)

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
  if(length(alt)==0) stop("sample_id not found.")
  pheno <- pheno %>% rename(sample_id = all_of(alt[1]))
}
if(!("arm" %in% names(pheno))){
  alt <- intersect(names(pheno), c("group","treatment","condition","study_arm"))
  if(length(alt)==0) stop("arm not found.")
  pheno <- pheno %>% rename(arm = all_of(alt[1]))
}
if(!("time" %in% names(pheno))){
  alt <- intersect(names(pheno), c("week","timepoint","tp"))
  if(length(alt)==0) stop("time not found.")
  pheno <- pheno %>% rename(time = all_of(alt[1]))
}
if(!("patient" %in% names(pheno))){
  alt <- intersect(names(pheno), c("subject","patientid","id","pair","individual"))
  if(length(alt)==0) stop("patient not found.")
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
if(length(common) < 10){
  message("Expr colnames head: ", paste(head(colnames(expr_mat)), collapse=", "))
  message("Sample sheet head: ", paste(head(pheno$sample_id)), collapse=", ")
  stop("Too few matching samples between expression and sample sheet.")
}
expr_mat <- expr_mat[, common, drop=FALSE]
pheno <- pheno %>% filter(sample_id %in% common) %>% arrange(match(sample_id, common))
stopifnot(all(pheno$sample_id == colnames(expr_mat)))

if(!("probe_id" %in% names(tab_int))) stop("interaction table must include probe_id.")
if(!("t" %in% names(tab_int))) stop("interaction table must include t-stat column 't'.")

probe2entrez <- AnnotationDbi::select(hgu133plus2.db,
                                      keys = unique(tab_int$probe_id),
                                      keytype = "PROBEID",
                                      columns = c("ENTREZID","SYMBOL")) %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::distinct(PROBEID, ENTREZID, .keep_all = TRUE)

tab_int2 <- tab_int %>%
  left_join(probe2entrez, by = c("probe_id"="PROBEID")) %>%
  filter(!is.na(ENTREZID)) %>%
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  ungroup()

ranks <- tab_int2$t
names(ranks) <- tab_int2$ENTREZID
ranks <- sort(ranks, decreasing = TRUE)

hall <- msigdbr(species="Homo sapiens", category="H") %>% select(gs_name, entrez_gene) %>% distinct()
path_h <- split(hall$entrez_gene, hall$gs_name)
fg_h <- fgsea(pathways = path_h, stats = ranks, nperm = 10000, minSize = 10, maxSize = 500) %>% arrange(padj)
readr::write_csv(fg_h, "results/mechanism/GSEA_Hallmark_interaction_t.csv")

reac <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME") %>% select(gs_name, entrez_gene) %>% distinct()
path_r <- split(reac$entrez_gene, reac$gs_name)
fg_r <- fgsea(pathways = path_r, stats = ranks, nperm = 10000, minSize = 10, maxSize = 500) %>% arrange(padj)
readr::write_csv(fg_r, "results/mechanism/GSEA_Reactome_interaction_t.csv")

probe2sym <- AnnotationDbi::select(hgu133plus2.db,
                                  keys = rownames(expr_mat),
                                  keytype="PROBEID",
                                  columns=c("SYMBOL")) %>%
  filter(!is.na(SYMBOL)) %>% distinct(PROBEID, SYMBOL)

expr_sym <- expr_mat[probe2sym$PROBEID, , drop=FALSE]
rownames(expr_sym) <- probe2sym$SYMBOL
expr_sym <- rowsum(expr_sym, group=rownames(expr_sym)) / as.vector(table(rownames(expr_sym)))

pheno$arm <- factor(pheno$arm, levels=c("Control","BLCC"))
pheno$time <- factor(pheno$time, levels=c("week0","week1"))
design <- model.matrix(~ arm * time, pheno)
coef_name <- colnames(design)[grep("armBLCC:timeweek1", colnames(design))]
if(length(coef_name)!=1) stop("Cannot find interaction coefficient armBLCC:timeweek1 in design.")

hall_sym <- msigdbr(species="Homo sapiens", category="H") %>% select(gs_name, gene_symbol) %>% distinct()
gs_h <- split(hall_sym$gene_symbol, hall_sym$gs_name)
gsva_scores <- GSVA::gsva(expr_sym, gs_h, method="gsva", kcdf="Gaussian", mx.diff=TRUE, verbose=FALSE)

corfit <- duplicateCorrelation(gsva_scores, design, block = pheno$patient)
fit <- lmFit(gsva_scores, design, block = pheno$patient, correlation = corfit$consensus.correlation) |> eBayes()
gsva_tab <- topTable(fit, coef = coef_name, number = Inf, sort.by="P") %>% rownames_to_column("pathway")
readr::write_csv(gsva_tab, "results/mechanism/GSVA_Hallmark_DiD_interaction.csv")

prog <- progeny::progeny(expr_sym, scale=TRUE, organism="Human", top=100, perm=1)
prog_mat <- t(as.matrix(prog))
corfit2 <- duplicateCorrelation(prog_mat, design, block = pheno$patient)
fit2 <- lmFit(prog_mat, design, block = pheno$patient, correlation = corfit2$consensus.correlation) |> eBayes()
prog_tab <- topTable(fit2, coef = coef_name, number = Inf, sort.by="P") %>% rownames_to_column("pathway")
readr::write_csv(prog_tab, "results/mechanism/PROGENy_DiD_interaction.csv")

data(dorothea_hs, package="dorothea")
regulon <- dorothea_hs %>% filter(confidence %in% c("A","B","C"))
tf_activity <- viper::viper(expr_sym, regulon, method="scale", minsize=10, verbose=FALSE)

corfit3 <- duplicateCorrelation(tf_activity, design, block = pheno$patient)
fit3 <- lmFit(tf_activity, design, block = pheno$patient, correlation = corfit3$consensus.correlation) |> eBayes()
tf_tab <- topTable(fit3, coef = coef_name, number = Inf, sort.by="P") %>% rownames_to_column("TF")
readr::write_csv(tf_tab, "results/mechanism/DoRothEA_VIPER_TF_DiD_interaction.csv")

message("DONE. CSV saved to results/mechanism/")
