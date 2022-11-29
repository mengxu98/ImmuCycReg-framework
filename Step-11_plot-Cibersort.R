



library(patchwork)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
theme_set(theme_pubclean())
pkgs <- c("matrixStats", "pheatmap", "RColorBrewer", "tidyverse", "cowplot", "ggpubr", "bslib", "ggthemes")
lapply(pkgs, library, character.only = T)

#------------------------------------------------------------------------------#
# <U+6DFB><U+52A0>P<U+503C><U+548C><U+663E><U+8457><U+6027><U+6C34><U+5E73>:
# i)<U+6BD4><U+8F83><U+4E24><U+4E2A><U+6216><U+591A><U+4E2A><U+7EC4><U+7684><U+5747><U+503C>; ii)<U+5E76><U+5C06>p<U+503C><U+548C><U+663E><U+8457><U+6027><U+6C34><U+5E73><U+81EA><U+52A8><U+6DFB><U+52A0><U+5230>ggplot<U+56FE><U+4E2D>
# <U+51FD><U+6570>:
# compare_means(): <U+8BA1><U+7B97><U+5355><U+6B21><U+6216><U+591A><U+6B21><U+5747><U+503C><U+6BD4><U+8F83><U+7684><U+7ED3><U+679C>
# stat_compare_means: <U+53EF><U+5C06>P<U+503C><U+548C><U+663E><U+8457><U+6027><U+6C34><U+5E73><U+81EA><U+52A8><U+6DFB><U+52A0><U+5230>ggplot<U+56FE><U+4E2D>
# T-test	        t.test()	        <U+6BD4><U+8F83><U+4E24><U+7EC4>(<U+53C2><U+6570><U+68C0><U+9A8C>)
# Wilcoxon test	  wilcox.test()	    <U+6BD4><U+8F83><U+4E24><U+7EC4>(<U+975E><U+53C2><U+6570><U+68C0><U+9A8C>)
# ANOVA	          aov() or anova()	<U+6BD4><U+8F83><U+591A><U+7EC4>(<U+53C2><U+6570><U+68C0><U+9A8C>)
# Kruskal-Wallis	kruskal.test()	  <U+6BD4><U+8F83><U+591A><U+7EC4>(<U+975E><U+53C2><U+6570><U+68C0><U+9A8C>)
#------------------------------------------------------------------------------#

data_ALL_group <- read.table("data/data_ALL_group.txt",
  row.names = 1,
  header = T,
  check.names = F
)
data_ALL_group_t <- as.data.frame(t(data_ALL_group))
CIBERSORT_Results_all <- read.table("results/CIBERSORT_ALL.txt",
  sep = "\t",
  header = T,
  row.names = 1
)
CIBERSORT_Results_all_LM14 <- read.table("../Cibersort_LUAD_TPM_nor_LM14/results/CIBERSORT_ALL.txt",
  sep = "\t",
  header = T,
  row.names = 1
)

CIBERSORT_Results_all <- cbind.data.frame(CIBERSORT_Results_all, CIBERSORT_Results_all_LM14)

CIBERSORT_Results_all$NK <- CIBERSORT_Results_all$NK.cells.activated + CIBERSORT_Results_all$NK.cells.resting
data_ALL_group_t <- data_ALL_group_t[row.names(CIBERSORT_Results_all), ] #<U+7528><U+6765><U+6392><U+5217><U+884C><U+540D>,<U+9632><U+6B62><U+51FA><U+9519>

CIBERSORT_Results_all$group <- data_ALL_group_t$group
CIBERSORT_Results <- as.data.frame(CIBERSORT_Results_all)
# write.csv(CIBERSORT_Results,'CIBERSORT_Results.csv')
#<U+5206><U+7EC4><U+6BD4><U+8F83><U+53EA><U+80FD><U+4F7F><U+7528><U+4E00><U+4E2A><U+7EC6><U+80DE><U+4E4B><U+95F4><U+7684><U+6BD4><U+8F83>
my_comparisons <- list(
  c("Cluster1", "GTEx"),
  c("Cluster2", "GTEx"),
  c("Cluster3", "GTEx"),
  c("Cluster4", "GTEx")
) # <U+5206><U+7EC4><U+8BBE><U+5B9A>
#------------------------------------------------------------------------------#

ann_colors_cluster <- list(Group = c(
  Cluster1 = "#33ff00",
  Cluster2 = "#EFC000FF",
  Cluster3 = "#33ccff",
  Cluster4 = "#3333ff",
  GTEx = "#FF7F24"
))
annotation_cluster <- data.frame(Group = CIBERSORT_Results$group)

p_list <- list()
cells <- c("Th.cell", "Macrophages.M2", "NK.cells.activated") # ,"Plasma" "Plasma.cells",

for (i in 1:length(cells)) {
  if (i <= 2) {
    CIBERSORT_Results_select <- CIBERSORT_Results_all[, c(cells[i], "group")]
    CIBERSORT_barplot_all <- CIBERSORT_Results_select %>%
      as.data.frame() %>%
      # rownames_to_column("sample") %>% #bug
      pivot_longer(
        cols = 1:c(ncol(CIBERSORT_Results_select) - 1),
        names_to = "CellType",
        values_to = "Composition"
      )

    p <- ggplot(
      CIBERSORT_barplot_all,
      aes(
        x = group,
        y = Composition
      )
    ) +
      guides(fill = guide_legend(title = NULL)) +
      stat_compare_means(
        method = "wilcox.test", # wilcox.test	#t.test
        label = "p.signif", # "p.signif"<U+7528>*<U+6765><U+8868><U+793A><U+5DEE><U+5F02>,"p.format"<U+7528>p<U+503C>
        # label.y = 0.2,
        comparisons = my_comparisons,
        # label.y = c(max(CIBERSORT_barplot_all$Composition)-0.09,
        #             max(CIBERSORT_barplot_all$Composition)-0.07,
        #             max(CIBERSORT_barplot_all$Composition)-0.05,
        #             max(CIBERSORT_barplot_all$Composition)),
        label.y = c(
          0.7,
          0.65,
          0.6,
          0.55
        ),
        # step.increase = 0.01,
        bracket.size = 0.3,
        sizen = 4,
        color = "#6699cc"
      ) + # steelblue
      # ylim(0, max(CIBERSORT_barplot_all$Composition) + 0.07) +
      ylim(0, 0.75) +
      # stat_compare_means(label.y=max(CIBERSORT_barplot_all$Composition)+0.1)+

      # stat_summary(fun = mean, #fun.y = mean
      #              geom = "point",
      #              shape = 18,
      #              size = 1,
      #              color = "#FC4E07") +
      # <U+5750><U+6807><U+8F74>
      labs(
        x = "",
        y = "Percentage"
      ) +
      # geom_hline(yintercept = 0.1,
      #            lty=4,
      #            col="black",
      #            lwd=0.8) +
      # geom_violin(trim = FALSE) +
      stat_summary(
        fun.data = "mean_sdl",
        fun.args = list(mult = 1),
        geom = "pointrange",
        color = "black"
      ) +
      geom_violin(aes(fill = group),
        # position = position_dodge(0.9),
        trim = FALSE
      ) +
      geom_boxplot(width = 0.12) +
      # geom_boxplot(aes(color = group),
      #              width = 0.15,
      #              position = position_dodge(0.9)) +
      scale_fill_manual(values = c(
        "#2874c5",
        "#008a00",
        "#c6524a",
        "#eabf00",
        "#696969"
      )) + # FF7F24
      theme(legend.position = "none") +
      scale_color_manual(values = c(
        "#2874c5",
        "#008a00",
        "#c6524a",
        "#eabf00",
        "#696969"
      )) +
      facet_wrap(~CellType) +
      theme_gray() +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 10
        )
      ) #+ scale_y_continuous(expand = c(0,0))
    p
    p_list[[i]] <- p
  } else {
    CIBERSORT_Results_select <- CIBERSORT_Results_all[, c(cells[i], "group")]
    CIBERSORT_barplot_all <- CIBERSORT_Results_select %>%
      as.data.frame() %>%
      # rownames_to_column("sample") %>% #bug
      pivot_longer(
        cols = 1:c(ncol(CIBERSORT_Results_select) - 1),
        names_to = "CellType",
        values_to = "Composition"
      )

    p <- ggplot(
      CIBERSORT_barplot_all,
      aes(
        x = group,
        y = Composition
      )
    ) +
      guides(fill = guide_legend(title = NULL)) +
      stat_compare_means(
        method = "wilcox.test", # wilcox.test	#t.test
        label = "p.signif", # "p.signif"<U+7528>*<U+6765><U+8868><U+793A><U+5DEE><U+5F02>,"p.format"<U+7528>p<U+503C>
        # label.y = 0.2,
        comparisons = my_comparisons,
        # label.y = c(max(CIBERSORT_barplot_all$Composition)-0.09,
        #             max(CIBERSORT_barplot_all$Composition)-0.07,
        #             max(CIBERSORT_barplot_all$Composition)-0.05,
        #             max(CIBERSORT_barplot_all$Composition)),
        label.y = c(
          0.23,
          0.21,
          0.19,
          0.17
        ),
        # step.increase = 0.01,
        bracket.size = 0.3,
        sizen = 4,
        color = "#6699cc"
      ) + # steelblue
      # ylim(0, max(CIBERSORT_barplot_all$Composition) + 0.07) +
      ylim(0, 0.25) +
      # stat_compare_means(label.y=max(CIBERSORT_barplot_all$Composition)+0.1)+

      # stat_summary(fun = mean, #fun.y = mean
      #              geom = "point",
      #              shape = 18,
      #              size = 1,
      #              color = "#FC4E07") +
      # <U+5750><U+6807><U+8F74>
      labs(
        x = "",
        y = "Percentage"
      ) +
      # geom_hline(yintercept = 0.1,
      #            lty=4,
      #            col="black",
      #            lwd=0.8) +
      # geom_violin(trim = FALSE) +
      stat_summary(
        fun.data = "mean_sdl",
        fun.args = list(mult = 1),
        geom = "pointrange",
        color = "black"
      ) +
      geom_violin(aes(fill = group),
        # position = position_dodge(0.9),
        trim = FALSE
      ) +
      geom_boxplot(width = 0.12) +
      # geom_boxplot(aes(color = group),
      #              width = 0.15,
      #              position = position_dodge(0.9)) +
      scale_fill_manual(values = c(
        "#2874c5",
        "#008a00",
        "#c6524a",
        "#eabf00",
        "#696969"
      )) + # FF7F24
      theme(legend.position = "none") +
      scale_color_manual(values = c(
        "#2874c5",
        "#008a00",
        "#c6524a",
        "#eabf00",
        "#696969"
      )) +
      facet_wrap(~CellType) +
      theme_gray() +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 10
        )
      ) #+ scale_y_continuous(expand = c(0,0))
    p
    p_list[[i]] <- p
  }
}

p_list[[1]] +
  p_list[[2]] +
  p_list[[3]] +
  # p_list[[4]] +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") +
    theme(text = element_text(size = 17)) +
    # theme(legend.title = element_text(color="134", size=16, face="bold"))+
    theme(text = element_text(family = "Times New Roman"))

ggsave(paste0("../manuscript_review1/figure/", "Fig7.png"), width = 9, height = 4.2)
ggsave(paste0("../manuscript_review1/figure/", "Fig7.pdf"), width = 9, height = 4.2)
