library(fgsea)
library(msigdbr)

# Load DEG results
deg <- read.csv("DEG_interaction_results.csv")

# Rank genes
ranked <- deg$logFC
names(ranked) <- deg$Gene

# Run GSEA
pathways <- msigdbr(species="Homo sapiens", category="H")
fgseaRes <- fgsea(pathways, ranked)

write.csv(fgseaRes,"GSEA_results.csv")
