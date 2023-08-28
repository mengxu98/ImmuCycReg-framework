rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "GTEx-LUAD.Rdata"))

samples20 <- read.csv("../data/Samples_20.csv") %>% .[, 1]
geneList <- read.table("../data/Genes_17.txt", header = T) %>% .[, 1]

sampleLabel <- read.csv("../data/sample_cluster_4.csv",
                        header = FALSE,
                        col.names = c("sample", "cluster"))

TCGA <- tcga_luad[geneList, ]
GTEx <- gtex_luad[geneList, ]

TCGA <- t(scale(t(TCGA)))
TCGA[TCGA >= 2] <- 2
TCGA[TCGA <= -2] <- -2

Heatmap(TCGA[, samples20],
        name = "Z-score",
        show_column_names = FALSE,
        cluster_rows = FALSE,
        col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
        width = unit(10, "cm"))

GTEx <- t(scale(t(GTEx)))
GTEx[GTEx >= 2] <- 2
GTEx[GTEx <= -2] <- -2

expMatList <- list()
for (i in 1:4) {
  expMatList[[i]] <- TCGA[, sampleLabel %>% dplyr::filter(cluster == i) %>% .[, 1]]
}
expMatList[[5]] <- GTEx

htList <- list()
sampleColors <- c("#2874c5", "#008a00", "#c6524a", "#eabf00", "#696969")
for (i in 1:length(expMatList)) {
  if (i < 5) {
    label <- paste0("TCGA-Cluster", i)
  } else {
    label <- "GTEx"
  }
  sampleColor <- sampleColors[i]
  width <- 30 * ncol(expMatList[[i]]) / (ncol(TCGA) + ncol(GTEx))
  annotation <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = sampleColor),
                                                   labels = label,
                                                   labels_gp = gpar(col = "white",
                                                                    fontsize = 10)))
  
  htList[[i]] <- Heatmap(expMatList[[i]],
                         name = "Z-score",
                         show_column_names = FALSE,
                         cluster_rows = FALSE,
                         col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
                         # width = unit(5, "cm"),
                         width = unit(width, "cm"),
                         border = TRUE,
                         top_annotation = annotation)
}

pdf(paste0(pathSave, "SuppFig5.pdf"), width = 14, height = 4.5)
draw(htList[[1]]+htList[[2]]+htList[[3]]+htList[[4]]+htList[[5]])
dev.off()

# Single gene
for (i in 1:length(geneList)) {
  targetGene <- geneList[i]
  data_tcga <- cbind.data.frame(Group = rep("TCGA", length(TCGA[targetGene, ])),
                                Expression = TCGA[targetGene, ])

  data_GTEx <- cbind.data.frame(Group = rep("GTEx", length(GTEx[targetGene, ])),
                                Expression = GTEx[targetGene, ])
  
  data_plot <- rbind.data.frame(data_tcga, data_GTEx)
  box.plot(data_plot, boxColor = c("#e1181f", "#253870"))
  ggsave(paste0(check.dir("../../Results/Figures/boxplot/"), targetGene, ".pdf"),
         width = 2.5,
         height = 4)
}
