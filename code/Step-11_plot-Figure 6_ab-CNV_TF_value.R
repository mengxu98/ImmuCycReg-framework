rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

eveluateR2 <- read.csv(paste0(pathSave, "TF_regulatory_value.csv"))
eveluateR2 <- eveluateR2[,c("Gene", "Cluster", "L0_Rsquare")]

eveluateCNVs <- read.csv(paste0(pathSave, "CNV_correlation_value.csv"))
eveluateCNVs <- eveluateCNVs[,c("Gene", "Cluster", "CNVs")]

clusterColor <- c("#2874c5", "#008a00", "#c6524a", "#eabf00")

p1 <- bar.plot(eveluateCNVs, barColor = clusterColor) +
  labs(x = "Gene", y = "CNV correlation value") +
  ylim(-0.5, 0.5) +
  coord_flip()

p2 <- bar.plot(eveluateR2, barColor = clusterColor) +
  labs(x = "", y = "TF regulatory value") +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  ylim(0, 1) +
  coord_flip()

p1 + p2 +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(paste0(pathSave, "Figures/Figure 6_ab-CNV_TF_value.pdf"),
       width = 5.5,
       height = 5)
