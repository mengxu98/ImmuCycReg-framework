rm(list = ls())
source("functions/Functions.R")

# Load --------------------------------------------------------------------
load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "GTEx-LUAD.Rdata"))

samples_20 <- read.csv("../data/Samples_20.csv")
genes_list <- read.table("../data/Genes_17.txt", header = T)
TCGA <- tcga_luad[genes_list$gene, ]
GTEX <- gtex_luad[genes_list$gene, ]

TCGA <- t(scale(t(TCGA)))
TCGA[TCGA >= 2] <- 2
TCGA[TCGA <= -2] <- -2

Heatmap(TCGA,
  name = "expr",
  show_column_names = F,
  col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
  width = unit(10, "cm")
)

GTEX <- t(scale(t(GTEX)))
GTEX[GTEX >= 2] <- 2
GTEX[GTEX <= -2] <- -2
Heatmap(GTEX,
  name = "expr",
  show_column_names = F,
  col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
  width = unit(10, "cm")
)

sample_label <- read.table("../data/sample_cluster_4.csv",
                           header = F,
                           sep = ",",
                           check.names = FALSE
)

cluster1 <- sample_label[which(sample_label$V2 == 1), 1]
cluster1 <- intersect(cluster1, colnames(TCGA))

cluster2 <- sample_label[which(sample_label$V2 == 2), 1]
cluster2 <- intersect(cluster2, colnames(TCGA))

cluster3 <- sample_label[which(sample_label$V2 == 3), 1]
cluster3 <- intersect(cluster3, colnames(TCGA))

cluster4 <- sample_label[which(sample_label$V2 == 4), 1]
cluster4 <- intersect(cluster4, colnames(TCGA))

expmat1 <- TCGA[, cluster1]
expmat2 <- TCGA[, cluster2]
expmat3 <- TCGA[, cluster3]
expmat4 <- TCGA[, cluster4]
expmat5 <- GTEX

ht_list <- Heatmap(expmat1,
                   name = "Z-score",
                   show_column_names = F,
                   # clustering_distance_rows = "spearman",
                   col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
                   width = unit(5, "cm"),
                   border = T,
                   top_annotation = HeatmapAnnotation(foo = anno_block(
                     gp = gpar(fill = "#2874c5"),
                     labels = c("TCGA-Cluster1"),
                     labels_gp = gpar(
                       col = "white",
                       fontsize = 10
                     )
                   ))
) + Heatmap(expmat2,
            name = "Z-score",
            show_column_names = F,
            col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
            width = unit(5, "cm"),
            border = T,
            top_annotation = HeatmapAnnotation(foo = anno_block(
              gp = gpar(fill = "#008a00"),
              labels = c("TCGA-Cluster2"),
              labels_gp = gpar(
                col = "white",
                fontsize = 10
              )
            ))
) + Heatmap(expmat3,
            name = "Z-score",
            show_column_names = F,
            col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
            width = unit(5, "cm"),
            border = T,
            top_annotation = HeatmapAnnotation(foo = anno_block(
              gp = gpar(fill = "#c6524a"),
              labels = c("TCGA-Cluster3"),
              labels_gp = gpar(
                col = "white",
                fontsize = 10
              )
            ))
) + Heatmap(expmat4,
            name = "Z-score",
            show_column_names = F,
            col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
            width = unit(5, "cm"),
            border = T,
            top_annotation = HeatmapAnnotation(foo = anno_block(
              gp = gpar(fill = "#eabf00"),
              labels = c("TCGA-Cluster4"),
              labels_gp = gpar(
                col = "white",
                fontsize = 10
              )
            ))
) + Heatmap(expmat5,
            name = "Z-score",
            show_column_names = F,
            col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
            width = unit(5, "cm"),
            border = T,
            top_annotation = HeatmapAnnotation(foo = anno_block(
              gp = gpar(fill = "#696969"),
              labels = c("GTEx"),
              labels_gp = gpar(
                col = "white",
                fontsize = 10
              )
            ))
)

pdf("../../Results/figures/SuppFig5.pdf", width = 12, height = 4.5)
draw(ht_list, n = 2)
dev.off()

png("../../Results/figures/SuppFig5.png", width = 7300, height = 3000, res = 600)
draw(ht_list, n = 2)
dev.off()

# Single gene
for (i in 1:nrow(genes_list)) {
  target_gene <- genes_list$gene[i]
  data_tcga <- cbind.data.frame(t(TCGA[target_gene, ]), rep("TCGA", ncol(TCGA[target_gene, ])))
  names(data_tcga) <- c("Expression", "Group")
  data_gtex <- cbind.data.frame(t(GTEX[target_gene, ]), rep("GTEx", ncol(GTEX[target_gene, ])))
  names(data_gtex) <- c("Expression", "Group")
  data_plot <- rbind.data.frame(data_tcga, data_gtex)
  data_plot$Sample <- rownames(data_plot)
  ggplot(data = data_plot, aes(x = Group, y = Expression)) +
    geom_boxplot(aes(fill = Group)) +
    ylab(paste0(target_gene, " Expression")) +
    stat_compare_means() +
    theme_bw()
  geom_jitter()
  ggsave(paste0("../../Results/figures/boxplot/", target_gene, ".pdf"), width = 3.5, height = 3.5)
}
