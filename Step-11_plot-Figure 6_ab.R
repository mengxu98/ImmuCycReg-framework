

rm(list = ls())
library(reshape)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(patchwork)

genes_list <- read.table(paste("data/genes_list", ".txt", sep = ""),
  header = TRUE
)

genes_list <- arrange(genes_list, desc(gene))
eveluate_L0 <- read.csv("results/evaluate_all.csv")
eveluate_L0 <- eveluate_L0[, -4]
eveluate_L0 <- arrange(eveluate_L0, desc(Cluster))
eveluate_L0$Gene <- factor(eveluate_L0$Gene, levels = genes_list$gene)

eveluate_CNV <- read.csv("results/CNV.csv")
eveluate_CNV <- arrange(eveluate_CNV, desc(Cluster))
eveluate_CNV$Gene <- factor(eveluate_CNV$Gene, levels = genes_list$gene)

p1 <- ggplot() +
  geom_bar(
    data = eveluate_L0,
    aes(x = Gene, y = L0_Rsquare, fill = Cluster),
    stat = "identity", position = "dodge", width = 0.9, size = 0.5
  ) +
  theme_bw() +
  scale_fill_manual(values = c("#2874c5", "#008a00", "#c6524a", "#eabf00")) +
  labs(x = "", y = "TF Regulatory Value", fill = "", color = "") +
  theme(
    axis.text.x = element_text(size = rel(1.5), angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) +
  ylim(0, 1) +
  coord_flip()

p2 <- ggplot() +
  geom_bar(
    data = eveluate_CNV,
    aes(x = Gene, y = CNVs, fill = Cluster),
    stat = "identity", position = "dodge", width = 0.9, size = 0.5
  ) +
  theme_bw() +
  scale_fill_manual(values = c("#2874c5", "#008a00", "#c6524a", "#eabf00")) +
  labs(x = "Gene", y = "CNV Correlation Value", fill = "", color = "") +
  theme(
    axis.text.x = element_text(size = rel(1.5), angle = 45, hjust = 1),
    axis.text.y = element_text(size = rel(1.5), hjust = 1),
    legend.position = "none"
  ) +
  ylim(-0.5, 0.5) +
  coord_flip()

p2 + p1 + plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "a")

ggsave(paste0("figure/", "Fig. 6-CNV_TFs.pdf"), width = 5.5, height = 5)
ggsave(paste0("figure/", "Fig. 6-CNV_TFs.png"), width = 5.5, height = 5, dpi = 600)
