rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "GTEx-LUAD.Rdata"))

Genes249 <- read.table(paste0(pathRead, "Genes_249.txt"), header = TRUE) %>% .[, 1]
samplesCluster <- read.csv(paste0(pathRead, "sample_cluster_4.csv"), header = F)

foldChange <- 1
padjThreshold <- 0.05

VolcanoPlotList <- list()
for (j in 1:4) {
  message("Cluster ", j, " DESeq2 start......")
  check.dir(paste0(pathSave, "/DESeq2/cluster", j))

  samples <- samplesCluster[which(samplesCluster$V2 == j), 1]
  tcgaCluster <- tcga_luad[Genes249, samples] %>% na.omit()
  gtexCluster <- gtex_luad[Genes249, ] %>% na.omit()
  
  tcga_gtex <- cbind(tcgaCluster, gtexCluster)
  group1 <- rep(c("C1"), each = ncol(tcgaCluster))
  group2 <- rep(c("C2"), each = ncol(gtexCluster))
  group <- c(group1, group2)
  
  groupFile <- data.frame(groupInfo = group)
  rownames(groupFile) <- colnames(tcga_gtex)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(as.matrix(tcga_gtex), 0),
                                colData = groupFile,
                                design = ~groupInfo)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("groupInfo", "C1", "C2"))
  
  diffDenes <- subset(res, padj < padjThreshold & abs(log2FoldChange) > foldChange)
  
  DESeq2Res <- as.data.frame(res[order(res$padj), ])
  DESeq2Res$gene <- rownames(DESeq2Res)
  DESeq2Res <- DESeq2Res[Genes249, c(7, 1:6)]
  
  DESeq2Res$reg <- "not DE"
  DESeq2Res$reg[which(DESeq2Res$log2FoldChange >= 0.5)] <- "min-up-regulated"
  DESeq2Res$reg[which(DESeq2Res$log2FoldChange <= -0.5)] <- "min-down-regulated"
  DESeq2Res$reg[which(DESeq2Res$log2FoldChange >= foldChange & DESeq2Res$padj <= padjThreshold)] <- "up-regulated"
  DESeq2Res$reg[which(DESeq2Res$log2FoldChange <= -foldChange & DESeq2Res$padj <= padjThreshold)] <- "down-regulated"
  write.csv(DESeq2Res,
            file = paste0(pathSave, "DESeq2_results_cluster", j, ".csv"),
            row.names = FALSE)
  
  dataPlot <- as.data.frame(res)
  dataPlot$change <- ifelse(dataPlot$pvalue < padjThreshold & abs(dataPlot$log2FoldChange) >= foldChange,
                            ifelse(dataPlot$log2FoldChange > foldChange,
                                   "Up (log2FC > 1, P-value < 0.05)",
                                   "Down (log2FC < -1, P-value < 0.05)"),
                            "no sig")
  
  dataPlot$gene <- rownames(dataPlot)
  dataPlot$label <- ifelse(dataPlot$pvalue < padjThreshold & abs(dataPlot$log2FoldChange) >= 1,
                           as.character(dataPlot$gene), "")
  
  DESeq2_volcano <- ggplot(dataPlot, aes(x = log2FoldChange, y = -log10(pvalue), colour = change)) +
    geom_point(alpha = 0.5, size = 3) +
    scale_color_manual(values = c("#006699", "#d2dae2", "#990033")) +
    geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8) +
    geom_hline(yintercept = -log10(padjThreshold), lty = 4, col = "black", lwd = 0.8) +
    labs(x = "Log2FoldChange", y = "-Log10(P-value)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          legend.title = element_blank()) +
    geom_text_repel(data = dataPlot, 
                    aes(x = log2FoldChange, y = -log10(pvalue),label = label),
                    size = 3,
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"),
                    segment.color = "black",
                    show.legend = FALSE)
  
  VolcanoPlotList[[j]] <- DESeq2_volcano
  
  message("Cluster ", j, " DESeq2 done......")
}

VolcanoPlotList <- VolcanoPlotList[[1]] +
  VolcanoPlotList[[2]] +
  VolcanoPlotList[[3]] +
  VolcanoPlotList[[4]] +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
VolcanoPlotList

ggsave(paste0(check.dir(paste0(pathSave, "Figures")), "/Figure 3.pdf"),
       VolcanoPlotList,
       width = 8,
       height = 8)
