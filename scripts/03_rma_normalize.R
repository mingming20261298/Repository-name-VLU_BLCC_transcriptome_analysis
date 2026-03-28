suppressPackageStartupMessages({
  library(affy)
  library(data.table)
})

setwd("~/JTM_GSE84571")

cel_dir <- "data/cel"
out_dir <- "data"

cel_files <- c(list.files(cel_dir, pattern="\\.CEL$", full.names=TRUE),
               list.files(cel_dir, pattern="\\.cel$", full.names=TRUE))
if (length(cel_files) == 0) stop("No CEL files found in data/cel.")

message("==> CEL files: ", length(cel_files))
raw <- ReadAffy(filenames = cel_files)

message("==> Running RMA with normalize=FALSE (workaround for pthread issue) ...")
eset <- rma(raw, normalize = FALSE)  # 关键：跳过 quantile normalization，避开 preprocessCore pthread

expr <- exprs(eset)

saveRDS(eset, file = file.path(out_dir, "GSE84571_eset_rma_noNorm.rds"))
fwrite(as.data.table(expr, keep.rownames = "probe_id"),
       file = file.path(out_dir, "GSE84571_expr_rma_noNorm.csv"))

message("==> DONE.")
message("Saved: data/GSE84571_eset_rma_noNorm.rds")
message("Saved: data/GSE84571_expr_rma_noNorm.csv")
