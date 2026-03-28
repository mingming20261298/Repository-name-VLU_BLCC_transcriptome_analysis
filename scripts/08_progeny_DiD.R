suppressPackageStartupMessages({
  library(limma)
  library(progeny)
})

dir.create("results_progeny", showWarnings=FALSE, recursive=TRUE)

expr_gene <- readRDS("derived_rma/expr_gene_rma_log2.rds")
ss <- read.csv("derived_rma/sample_sheet_rma.csv", stringsAsFactors=FALSE)

# progeny expects gene symbols
# pathways x samples
path_act <- progeny(as.matrix(expr_gene), organism="Human", top=100, scale=TRUE)

# DiD on pathway activities
arm <- factor(ss$arm, levels=c("Control","BLCC"))
time <- factor(ss$time, levels=c("Week0","Week1"))
patient <- factor(ss$patient)

design <- model.matrix(~ arm * time)
corfit <- duplicateCorrelation(t(path_act), design, block=patient)  # samples x pathways -> transpose
fit <- lmFit(t(path_act), design, block=patient, correlation=corfit$consensus)
fit <- eBayes(fit)

coef_name <- grep("arm.*:time.*", colnames(design), value=TRUE)
stopifnot(length(coef_name)==1)

tab <- topTable(fit, coef=coef_name, number=Inf, sort.by="P")
write.csv(tab, "results_progeny/PROGENy_pathway_DiD.csv", row.names=TRUE)

# also save activity matrix
write.csv(path_act, "results_progeny/PROGENy_pathway_activity.csv")

cat("PROGENy DiD DONE | coef=", coef_name,
    " | FDR<0.10=", sum(tab$adj.P.Val<0.10, na.rm=TRUE),
    " | FDR<0.05=", sum(tab$adj.P.Val<0.05, na.rm=TRUE), "\n", sep="")
