

rm(list = ls())
library("reshape")
library("ggplot2")
library("ggpubr")
library("ggthemes")
library("patchwork")
library("reshape2")
library("RColorBrewer")
mycol <- c("gray", "white")

results_10nets <- read.csv(paste0("evaluation_gnw_10_", 1, ".csv"))
results_10nets <- results_10nets[, -1]

df_res10 <- melt(results_10nets, id = "Dataset", variable.name = "Method", value.name = "AUROC")
df_res10$Method <- factor(df_res10$Method,
  levels = c("L0Reg.framework", "GENIE3"),
  labels = c("L0Reg framework", "GENIE3")
)

p2 <- ggplot(df_res10, aes(x = Dataset, y = AUROC, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black", width = .6) +
  scale_fill_manual(values = mycol) +
  # geom_errorbar(aes(ymin=AUROC - Sd, ymax=AUROC + Sd), position = position_dodge(.6), width=.2)
  scale_x_discrete(labels = c("Net1", "Net2", "Net3", "Net4", "Net5", "Net6", "Net7", "Net8", "Net9", "Net10")) +
  theme_bw() +
  ylim(0, 1)

results_5nets <- read.csv(paste0("evaluation_gnw_5_", 2, ".csv"))
results_5nets <- results_5nets[, -1]

df_res5 <- melt(results_5nets, id = "Dataset", variable.name = "Method", value.name = "AUROC")
df_res5$Method <- factor(df_res5$Method,
  levels = c("L0Reg.framework", "GENIE3"),
  labels = c("L0Reg framework", "GENIE3")
)

p1 <- ggplot(df_res5, aes(x = Dataset, y = AUROC, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black", width = .6) +
  scale_fill_manual(values = mycol) +
  # geom_errorbar(aes(ymin=AUROC - Sd, ymax=AUROC + Sd), position = position_dodge(.6), width=.2)
  theme_bw() +
  ylim(0, 1)

p1 + p2 + plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "a") +
  # plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") +
    theme(text = element_text(size = 12)) +
    # theme(legend.title = element_text(color="134", size=16, face="bold"))+
    theme(text = element_text(family = "Times New Roman")) # ,face = "bold"
# ggsave("../Results/figure/Supplementary Figure 7.png",width = 8, height = 4, dpi =600)
