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

message("==> Running MAS5 ...")
eset_mas5 <- mas5(raw)  # 不走 preprocessCore 的 RMA 线程实现

# MAS5 输出在非log尺度，通常对其取 log2 用于后续线性模型
expr <- log2(exprs(eset_mas5))

saveRDS(eset_mas5, file = file.path(out_dir, "GSE84571_eset_mas5.rds"))
fwrite(as.data.table(expr, keep.rownames = "probe_id"),
       file = file.path(out_dir, "GSE84571_expr_mas5_log2.csv"))

message("==> DONE.")
message("Saved: data/GSE84571_eset_mas5.rds")
message("Saved: data/GSE84571_expr_mas5_log2.csv")
