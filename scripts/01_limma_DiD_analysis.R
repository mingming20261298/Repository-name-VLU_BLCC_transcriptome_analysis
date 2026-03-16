# ==========================================
# VLU BLCC Transcriptomic Analysis Pipeline
# Dataset: GSE84571
# ==========================================

# Load libraries
library(limma)
library(GEOquery)
library(ggplot2)

# Download dataset
gset <- getGEO("GSE84571", GSEMatrix = TRUE)
expr <- exprs(gset[[1]])

# Metadata
pdata <- pData(gset[[1]])

# Design matrix
group <- factor(pdata$characteristics_ch1)
time <- factor(pdata$characteristics_ch1.1)

design <- model.matrix(~ group*time)

# Fit model
fit <- lmFit(expr, design)
fit <- eBayes(fit)

# Extract interaction results
results <- topTable(fit, coef=4, number=Inf)

# Save results
write.csv(results, "DEG_interaction_results.csv")

# Volcano plot
results$logP <- -log10(results$P.Value)

ggplot(results, aes(logFC, logP)) +
    geom_point(alpha=0.5) +
    theme_minimal() +
    ggtitle("BLCC-specific transcriptional changes")

ggsave("volcano_plot.png")
