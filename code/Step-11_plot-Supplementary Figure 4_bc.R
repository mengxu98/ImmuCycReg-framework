rm(list = ls())
source("functions/Functions.R")

pathSave <- "../../Results/"

netInfor <- read.csv(paste0(pathSave, "Edge_reg_num.csv"))

netEdgeInfor <- melt(netInfor[, c("Gene", "ImmuCycReg_framework", "Datasets")], id = "Gene")
p1 <- bar.plot(netEdgeInfor, barColor = c("white", "gray"))

netRegInfor <- melt(netInfor[, c("Gene", "Negative", "Positive")], id = "Gene")
p2 <- bar.plot(netRegInfor, barColor = c("#6666ff", "#ff9933"))

p3 <- p1 + p2 +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(paste0(pathSave, "Figure/Supplementary Figure 4_bc.pdf"),
       p3,
       width = 7.5,
       height = 3,
       dpi = 600)
