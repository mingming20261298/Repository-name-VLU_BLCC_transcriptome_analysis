suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
})

if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
for(p in c("GEOquery","limma","hgu133plus2.db","AnnotationDbi")){
  if(!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)
}
suppressPackageStartupMessages({
  library(GEOquery)
  library(AnnotationDbi)
  library(hgu133plus2.db)
})

dir.create("results/signature", showWarnings = FALSE, recursive = TRUE)

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
if(is.na(ss_path)) stop("Cannot find sample sheet in data/.")
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

pheno <- pheno %>% mutate(
  sample_id = as.character(sample_id),
  arm = ifelse(grepl("blcc|apligraf", arm, ignore.case=TRUE), "BLCC",
               ifelse(grepl("control|std|standard", arm, ignore.case=TRUE), "Control", as.character(arm))),
  time = ifelse(grepl("0", time), "week0",
                ifelse(grepl("1", time), "week1", as.character(time)))
)

common <- intersect(colnames(expr_mat), pheno$sample_id)
expr_mat <- expr_mat[, common, drop=FALSE]
pheno <- pheno %>% filter(sample_id %in% common) %>% arrange(match(sample_id, common))
stopifnot(all(pheno$sample_id == colnames(expr_mat)))

# ---- signature genes: FDR<0.10 else top 30 ----
if(!("adj.P.Val" %in% names(tab_int))) stop("Need adj.P.Val.")
if(!("probe_id" %in% names(tab_int))) stop("Need probe_id.")
if(!("t" %in% names(tab_int))) stop("Need t.")

sig <- tab_int %>% filter(adj.P.Val < 0.10) %>% arrange(adj.P.Val)
if(nrow(sig) < 10) sig <- tab_int %>% arrange(adj.P.Val) %>% slice_head(n=30)

# probe->symbol and keep one probe per symbol by max |t|
probe2sym <- AnnotationDbi::select(hgu133plus2.db,
                                  keys = unique(sig$probe_id),
                                  keytype="PROBEID",
                                  columns=c("SYMBOL")) %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::distinct(PROBEID, SYMBOL)

sig2 <- sig %>%
  left_join(probe2sym, by=c("probe_id"="PROBEID")) %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  slice_max(order_by=abs(t), n=1, with_ties=FALSE) %>%
  ungroup()

# build SYMBOL expression for GSE84571
probe2sym_all <- AnnotationDbi::select(hgu133plus2.db,
                                      keys = rownames(expr_mat),
                                      keytype="PROBEID",
                                      columns=c("SYMBOL")) %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::distinct(PROBEID, SYMBOL)

expr_sym <- expr_mat[probe2sym_all$PROBEID, , drop=FALSE]
rownames(expr_sym) <- probe2sym_all$SYMBOL
expr_sym <- rowsum(expr_sym, group=rownames(expr_sym)) / as.vector(table(rownames(expr_sym)))

genes <- intersect(sig2$SYMBOL, rownames(expr_sym))
sig2 <- sig2 %>% filter(SYMBOL %in% genes)
if(length(genes) < 5) stop("Too few signature genes after mapping.")

# weights = interaction t
w <- sig2$t; names(w) <- sig2$SYMBOL

# z-score per gene across samples within dataset
Z <- t(scale(t(expr_sym[genes, , drop=FALSE])))
score84571 <- colSums(Z * w[rownames(Z)], na.rm=TRUE) / sum(abs(w), na.rm=TRUE)

df84571 <- tibble(sample_id=names(score84571), score=as.numeric(score84571)) %>%
  left_join(pheno, by="sample_id")

readr::write_csv(sig2, "results/signature/Signature_genes_weights_from_interaction.csv")
readr::write_csv(df84571, "results/signature/Signature_score_GSE84571.csv")

# ---- External: GSE28914 ----
gse_id <- "GSE28914"
gse <- GEOquery::getGEO(gse_id, GSEMatrix=TRUE, getGPL=FALSE)
eset <- gse[[1]]
expr_ext <- exprs(eset)
pd_ext <- pData(eset)

probe2sym_ext <- AnnotationDbi::select(hgu133plus2.db,
                                      keys=rownames(expr_ext),
                                      keytype="PROBEID",
                                      columns=c("SYMBOL")) %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::distinct(PROBEID, SYMBOL)

expr_ext_sym <- expr_ext[probe2sym_ext$PROBEID, , drop=FALSE]
rownames(expr_ext_sym) <- probe2sym_ext$SYMBOL
expr_ext_sym <- rowsum(expr_ext_sym, group=rownames(expr_ext_sym)) / as.vector(table(rownames(expr_ext_sym)))

genes2 <- intersect(genes, rownames(expr_ext_sym))
if(length(genes2) < 5) stop("Too few overlapping genes with GSE28914.")

w2 <- w[genes2]
Zext <- t(scale(t(expr_ext_sym[genes2, , drop=FALSE])))
score28914 <- colSums(Zext * w2[rownames(Zext)], na.rm=TRUE) / sum(abs(w2), na.rm=TRUE)

df28914 <- tibble(sample_id=names(score28914), score=as.numeric(score28914))
pd_ext2 <- as.data.frame(pd_ext) %>%
  tibble::rownames_to_column("sample_id") %>%
  mutate(label = if("title" %in% names(.)) as.character(title) else sample_id)

df28914 <- df28914 %>% left_join(pd_ext2[,c("sample_id","label")], by="sample_id")

readr::write_csv(df28914, "results/signature/Signature_score_GSE28914.csv")

message("DONE. CSV saved to results/signature/")
