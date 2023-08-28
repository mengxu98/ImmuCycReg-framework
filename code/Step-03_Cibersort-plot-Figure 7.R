source("functions/Functions.R")
source("functions/Cibersort.R")

pathRead <- "../data/"
pathSave <- "../../Results/Cibersort/"

CIBERSORT(paste0(pathRead, "LM22_2.0.txt"),
          "../../Results/Cibersort_data_ALL.txt",
          pathSave = check.dir(pathSave),
          fileName = "CIBERSORT_Results_LM22.txt",
          perm = 1000)
CIBERSORT(paste0(pathRead, "LM14_2.0.txt"),
          "../../Results/Cibersort_data_ALL.txt",
          pathSave = check.dir(pathSave),
          fileName = "CIBERSORT_Results_LM14.txt",
          perm = 1000)

CibersortLM22 <- read.table(paste0(pathSave, "CIBERSORT_Results_LM22.txt"),
                            sep = "\t",
                            header = T,
                            row.names = 1)

CibersortLM14 <- read.table(paste0(pathSave, "CIBERSORT_Results_LM14.txt"),
                            sep = "\t",
                            header = T,
                            row.names = 1)

CibersortResults <- cbind.data.frame(CibersortLM22, CibersortLM14)
CibersortResults$NK <- CibersortResults$NK.cells.activated + CibersortResults$NK.cells.resting

sampleGroup <- read.table("../../Results/Data_ALL_group.txt",
                          row.names = 1,
                          header = T,
                          check.names = F) %>% t() %>% as.data.frame()

sampleGroup <- sampleGroup[rownames(CibersortResults), ]

CibersortResults$group <- sampleGroup$group

comparisons <- list(c("Cluster1", "GTEx"),
                    c("Cluster2", "GTEx"),
                    c("Cluster3", "GTEx"),
                    c("Cluster4", "GTEx"))

clusterColors <- c("#2874c5", "#008a00", "#c6524a", "#eabf00", "#696969")

cells <- c("Th.cell", "Macrophages.M2", "NK.cells.activated")

hightList <- list(c(0.7, 0.65, 0.6, 0.55),
                  c(0.7, 0.65, 0.6, 0.55),
                  c(0.23, 0.21, 0.19, 0.17))
ylims <- c(0.75, 0.75, 0.25)

pList <- list()
for (i in 1:length(cells)) {
  CibersortResultsSelect <- CibersortResults[, c(cells[i], "group")]
  CibersortBarplot <- CibersortResultsSelect %>%
    as.data.frame() %>%
    pivot_longer(cols = 1:c(ncol(CibersortResultsSelect) - 1),
                 names_to = "CellType",
                 values_to = "Composition")
  
  p <- ggplot(CibersortBarplot, aes(x = group, y = Composition)) +
    guides(fill = guide_legend(title = NULL)) +
    stat_compare_means(method = "wilcox.test",
                       label = "p.signif",
                       comparisons = comparisons,
                       label.y = hightList[[i]],
                       bracket.size = 0.3,
                       sizen = 4,
                       color = "#6699cc") +
    ylim(0, ylims[i]) +
    labs(x = "",
         y = "Percentage") +
    stat_summary(fun.data = "mean_sdl",
                 fun.args = list(mult = 1),
                 geom = "pointrange",
                 color = "black") +
    geom_violin(aes(fill = group),
                trim = FALSE) +
    geom_boxplot(width = 0.12) +
    scale_fill_manual(values = clusterColors) +
    theme(legend.position = "none") +
    scale_color_manual(values = clusterColors) +
    facet_wrap(~CellType) +
    theme_gray() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     vjust = 1,
                                     size = 10))
  
  pList[[i]] <- p
}

pList[[1]] +
  pList[[2]] +
  pList[[3]] +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") 

ggsave(paste0(pathSave, "Figures/Figure 7.pdf"),
       width = 9,
       height = 4)
