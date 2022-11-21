

rm(list = ls())
library(reshape)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(patchwork)

reg_num_all <- read.csv("results/Edge_reg_num.csv")

reg_num_all_dataset <- melt(reg_num_all[, c("Gene", "ImmuCycReg_framework", "Datasets")])
p1 <- ggplot() +
  geom_bar(
    data = reg_num_all_dataset,
    aes(x = Gene, y = value, fill = variable),
    stat = "identity", position = "dodge", width = 0.9, size = 1
  ) +
  theme_bw() +
  scale_fill_manual(values = c("#363636", "#999999")) +
  labs(x = "Gene", y = "", fill = "", color = "") +
  theme(axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1))

reg_num_all_reg <- melt(reg_num_all[, c("Gene", "Negative", "Positive")])
p2 <- ggplot() +
  geom_bar(
    data = reg_num_all_reg,
    aes(x = Gene, y = value, fill = variable),
    stat = "identity", position = "dodge", width = 0.9, size = 1
  ) +
  theme_bw() +
  scale_fill_manual(values = c("#6666ff", "#ff9933")) +
  labs(x = "Gene", y = "", fill = "", color = "") +
  theme(axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1))

p1 / p2 +
  plot_annotation(tag_levels = "a") +
  plot_layout(ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") +
    theme(text = element_text(size = 12))
theme(text = element_text(family = "Times New Roman"))

ggsave(paste0("../manuscript_review/figure/", "Supplementary Figure 4_bc.png"),
  width = 7.5,
  height = 3,
  dpi = 600
)
ggsave(paste0("../manuscript_review/figure/", "Supplementary Figure 4_bc.pdf"),
  width = 7.5,
  height = 3,
  dpi = 600
)
