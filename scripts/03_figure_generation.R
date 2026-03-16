library(ggplot2)

deg <- read.csv("DEG_interaction_results.csv")

ggplot(deg, aes(logFC, -log10(P.Value))) +
 geom_point(color="steelblue") +
 theme_classic()

ggsave("Figure3_volcano.png")
