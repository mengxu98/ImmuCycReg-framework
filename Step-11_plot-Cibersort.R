

library(patchwork)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
theme_set(theme_pubclean())
pkgs <- c("matrixStats", "pheatmap", "RColorBrewer", "tidyverse", "cowplot", "ggpubr", "bslib", "ggthemes")
lapply(pkgs, library, character.only = T)
#------------------------------------------------------------------------------#
data_ALL_group <- read.table('data/data_ALL_group.txt',
                             row.names = 1,
                             header = T,
                             check.names = F)
data_ALL_group_t <- as.data.frame(t(data_ALL_group))
CIBERSORT_Results_all <- read.table('results/CIBERSORT_ALL.txt',
                                    sep = '\t',
                                    header = T,
                                    row.names = 1)
CIBERSORT_Results_all_LM14 <- read.table('../Cibersort_LUAD_TPM_nor_LM14/results/CIBERSORT_ALL.txt',
                                         sep = '\t',
                                         header = T,
                                         row.names = 1)

CIBERSORT_Results_all <- cbind.data.frame(CIBERSORT_Results_all, CIBERSORT_Results_all_LM14)

CIBERSORT_Results_all$NK <- CIBERSORT_Results_all$NK.cells.activated+CIBERSORT_Results_all$NK.cells.resting
data_ALL_group_t <- data_ALL_group_t[row.names(CIBERSORT_Results_all),]

CIBERSORT_Results_all$group <- data_ALL_group_t$group
CIBERSORT_Results <- as.data.frame(CIBERSORT_Results_all)

my_comparisons <- list(c("Cluster1", "GTEx"),
                       c("Cluster2", "GTEx"),
                       c("Cluster3", "GTEx"),
                       c("Cluster4", "GTEx"))
#------------------------------------------------------------------------------#
ann_colors_cluster = list(Group = c(Cluster1 = "#2874c5",
                                    Cluster2 = "#008a00",
                                    Cluster3 = "#c6524a",
                                    Cluster4 = "#eabf00",
                                    GTEx = "#696969"))
annotation_cluster = data.frame(Group = CIBERSORT_Results$group)

p_list <- list()
cells <- c("Plasma.cells","Th.cell", "Macrophages.M2","NK") #,"Plasma"

for (i in 1:length(cells)) {
  CIBERSORT_Results_select <- CIBERSORT_Results_all[,c(cells[i],"group")]
  CIBERSORT_barplot_all <- CIBERSORT_Results_select %>%
    as.data.frame() %>%
    pivot_longer(cols = 1:c(ncol(CIBERSORT_Results_select)-1),
                 names_to = "CellType",
                 values_to = "Composition")
  
  p <- ggplot(CIBERSORT_barplot_all,
              aes(x = group,
                  y = Composition)) + 
    guides(fill=guide_legend(title=NULL)) + 
    stat_compare_means(method = "wilcox.test", #wilcox.test, t.test
                       label = "p.signif", # "p.signif"-*, "p.format"-p
                       # label.y = 0.2,
                       comparisons = my_comparisons, 
                       label.y = c(max(CIBERSORT_barplot_all$Composition)-0.18, 
                                   max(CIBERSORT_barplot_all$Composition)-0.12, 
                                   max(CIBERSORT_barplot_all$Composition)-0.06,
                                   max(CIBERSORT_barplot_all$Composition)),
                       # step.increase = 0.01,
                       bracket.size = 0.6,
                       sizen = 4, 
                       color = "#6699cc") +#steelblue
    ylim(0, max(CIBERSORT_barplot_all$Composition) + 0.07) +
    labs(x="",
         y="Percentage")+
    stat_summary(fun.data = "mean_sdl",
                 fun.args = list(mult = 1),
                 geom = "pointrange",
                 color = "black")+
    geom_violin(aes(fill = group),
                trim = FALSE) + 
    geom_boxplot(width = 0.12) +
    scale_fill_manual(values = c("#2874c5",
                                 "#008a00",
                                 "#c6524a",
                                 "#eabf00",
                                 "#696969")) +
    theme(legend.position = "none") +
    scale_color_manual(values = c("#2874c5",
                                  "#008a00",
                                  "#c6524a",
                                  "#eabf00",
                                  "#696969"))+
    facet_wrap(~CellType) +
    theme_gray() +
    theme_bw()+
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 10
      ) )
  p
  p_list[[i]] <- p
  
}

p_list[[1]] + 
  p_list[[2]] +  
  p_list[[3]] + 
  p_list[[4]] + 
  plot_layout(ncol = 4) +
  plot_annotation(tag_levels = 'a') + 
  plot_layout(guides = 'collect')&
  theme(legend.position='bottom')

ggsave(paste0("../manuscript_review/figure/","Fig4.png"), width = 11, height = 4)
ggsave(paste0("../manuscript_review/figure/","Fig4.pdf"), width = 11, height = 4)

