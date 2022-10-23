

rm(list=ls())
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(psych)
library(ggpubr)
library(ggthemes)
library(rtracklayer)
library(svglite)

# Load --------------------------------------------------------------------
if (T) {
  rm(list = ls())
  load("../data/LUAD_raw/luad-rsem-TPM-tcga-t_normlized.Rdata")
  load("../data/LUAD_raw/lung-rsem-TPM-gtex_normlized.Rdata")
}
if (T) {
  rm(list = ls())
  load("../data/LUAD_raw/luad-rsem-TPM-tcga-t_unnormlized.Rdata")
  load("../data/LUAD_raw/lung-rsem-TPM-gtex_unnormlized.Rdata")
}
if (T) {
  rm(list = ls())
  load("../data/LUAD_raw/luad-rsem-fpkm-tcga-t_normlized.Rdata")
  raw_tcga <- tcga_raw
  load("../data/LUAD_raw/lung-rsem-fpkm-gtex_normlized.Rdata")
  raw_gtex <- gtex_raw
}
if (T) {
  rm(list = ls())
  load("../data/LUAD_raw/luad-rsem-fpkm-tcga-t_unnormlized.Rdata")
  raw_tcga <- tcga_raw
  load("../data/LUAD_raw/lung-rsem-fpkm-gtex_unnormlized.Rdata")
  raw_gtex <- gtex_raw
}

if (T) { #Count
  rm(list = ls())
  #load('../data/LUAD_raw/luad-rsem-count-tcga-t-drop-dupli-15.Rdata')
  load('../data/raw_tcga_gtex_mrna.RData')
  raw_tcga <- raw_tcga_mrna; rm(raw_tcga_mrna)
  raw_gtex <- raw_gtex_mrna; rm(raw_gtex_mrna)
}

samples_20 <- read.csv("../data/samples/sample_20.csv")
genes_list <- read.table("../data CNV/all genes.txt",header = T)
TCGA <- raw_tcga[genes_list$gene, ]
GTEX <- raw_gtex[genes_list$gene, ]

TCGA <- t(scale(t(TCGA)))
TCGA[TCGA>=2]=2
TCGA[TCGA<=-2]=-2

Heatmap(TCGA,
        name = "expr",
        show_column_names = F,
        col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
        width = unit(10, "cm")
)

GTEX <- t(scale(t(GTEX)))
GTEX[GTEX>=2]=2
GTEX[GTEX<=-2]=-2
Heatmap(GTEX,
        #matrix(rnorm(nrow(expmat) * ncol(expmat)), ncol = ncol(expmat)),
        name = "expr",
        show_column_names = F,
        col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
        width = unit(10, "cm")
)

if (T) { #492 samples
  sample_label <- read.table('../results/2191/NMF/cluster-nk-ligands-k=8_raw/sample_cluster_4.csv',
                             header = F,
                             sep=',',
                             check.names=FALSE)
  
  cluster1<-sample_label[which(sample_label$V2==1),1]
  cluster1<-intersect(cluster1, colnames(TCGA))
  
  cluster2<-sample_label[which(sample_label$V2==2),1]
  cluster2<-intersect(cluster2, colnames(TCGA))
  
  cluster3<-sample_label[which(sample_label$V2==3),1]
  cluster3<-intersect(cluster3, colnames(TCGA))
  
  cluster4<-sample_label[which(sample_label$V2==4),1]
  cluster4<-intersect(cluster4, colnames(TCGA))
  
  expmat1 <- TCGA[ ,cluster1]
  expmat2 <- TCGA[ ,cluster2]
  expmat3 <- TCGA[ ,cluster3]
  expmat4 <- TCGA[ ,cluster4]
  expmat5 <- GTEX
  
  
  ht_list <- Heatmap(expmat1,
                     #matrix(rnorm(nrow(expmat) * ncol(expmat)), ncol = ncol(expmat)),
                     name = "Z-score", 
                     show_column_names = F,
                     #clustering_distance_rows = "spearman",
                     col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
                     width = unit(5, "cm"),
                     border = T,
                     top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "#2874c5"),
                                                                         labels = c("TCGA-Cluster1"),
                                                                         labels_gp = gpar(col = "white",
                                                                                          fontsize = 10)))
  ) + Heatmap(expmat2,
              #matrix(rnorm(nrow(expmat) * ncol(expmat)), ncol = ncol(expmat)),
              name = "Z-score",
              show_column_names = F,
              col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
              width = unit(5, "cm"),
              border = T,
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "#008a00"),
                                                                  labels = c("TCGA-Cluster2"),
                                                                  labels_gp = gpar(col = "white",
                                                                                   fontsize = 10)))
  ) + Heatmap(expmat3,
              #matrix(rnorm(nrow(expmat) * ncol(expmat)), ncol = ncol(expmat)),
              name = "Z-score",
              show_column_names = F,
              col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
              width = unit(5, "cm"),
              border = T,
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "#c6524a"),
                                                                  labels = c("TCGA-Cluster3"),
                                                                  labels_gp = gpar(col = "white",
                                                                                   fontsize = 10)))
  ) + Heatmap(expmat4,
              #matrix(rnorm(nrow(expmat) * ncol(expmat)), ncol = ncol(expmat)),
              name = "Z-score",
              show_column_names = F,
              col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
              width = unit(5, "cm"),
              border = T,
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "#eabf00"),
                                                                  labels = c("TCGA-Cluster4"),
                                                                  labels_gp = gpar(col = "white",
                                                                                   fontsize = 10)))
  ) + Heatmap(expmat5,
              #matrix(rnorm(nrow(expmat) * ncol(expmat)), ncol = ncol(expmat)),
              name = "Z-score",
              show_column_names = F,
              col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
              width = unit(5, "cm"),
              border = T,
              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "#696969"),
                                                                  labels = c("GTEx"),
                                                                  labels_gp = gpar(col = "white",
                                                                                   fontsize = 10)))
  )
  
  pdf("../manuscript_review/figure/SuppFig5.pdf", width = 12, height = 4.5)
  draw(ht_list,n=2)
  dev.off()
  
  png("../manuscript_review/figure/SuppFig5.png", width = 7300, height = 3000, res = 600)
  draw(ht_list,n=2)
  dev.off()
  
}

TCGA <- raw_tcga[genes_list$x, samples_20$x]; #write.table(TCGA,"TCGA-20samples.txt", sep = "\t",quote = F)
GTEX <- raw_gtex[genes_list$x, ]
TCGA_20 <- t(scale(t(TCGA)))
TCGA_20[TCGA_20>=2]=2
TCGA_20[TCGA_20<=-2]=-2
Heatmap(TCGA_20,
        name = "Z-score",
        show_column_names = T,
        col = circlize::colorRamp2(c(-2, 0, 2), c("#2c58af", "white", "#d43a35")),
        width = unit(10, "cm")
)
# Load --------------------------------------------------------------------
rm(list = ls())
load("../data/LUAD_raw/luad-rsem-TPM-tcga-t_normlized.Rdata")
load("../data/LUAD_raw/lung-rsem-TPM-gtex_normlized.Rdata")

samples_20 <- read.csv("../data/samples/sample_20.csv")
genes_list <- read.table("../data CNV/all genes.txt",header = T)

TCGA <- raw_tcga[genes_list$x, samples_20$x]
GTEX <- raw_gtex[genes_list$x, ]
TCGA <- raw_tcga[genes_list$x, ]
GTEX <- raw_gtex[genes_list$x, ]

for (i in 1:nrow(genes_list)) {
  target_gene <- genes_list$x[i]
  data_tcga <- cbind.data.frame(t(TCGA[target_gene,]),rep("TCGA",ncol(TCGA[target_gene,])))
  names(data_tcga) <- c("Expression", "Group")
  data_gtex <- cbind.data.frame(t(GTEX[target_gene,]),rep("GTEx",ncol(GTEX[target_gene,])))
  names(data_gtex) <- c("Expression", "Group")
  data_plot <- rbind.data.frame(data_tcga, data_gtex)
  data_plot$Sample <- rownames(data_plot)
  ggplot(data=data_plot, aes(x=Group,y=Expression))+
    geom_boxplot(aes(fill=Group))+
    ylab(paste0(target_gene," Expression"))+ 
    stat_compare_means()+
    theme_bw()
    geom_jitter()
  ggsave(paste0("results/boxplot/",target_gene,".pdf"),width = 3.5,height = 3.5)
}

