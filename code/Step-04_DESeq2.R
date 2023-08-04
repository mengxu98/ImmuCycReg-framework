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

for (j in 1:length(table(samplesCluster[, 2]))) {
  message(paste0("Cluster ", j, " DESeq2 start......"))
  check.dir(paste0(pathSave, "/DESeq2/cluster", j))

  samples <- samplesCluster[which(samplesCluster$V2 == j), 1]
  tcgaCluster <- tcga_luad[Genes249, samples] %>% na.omit()
  gtexCluster <- gtex_luad[Genes249, ] %>% na.omit()
  
  tcga_gtex <- cbind(tcgaCluster, gtexCluster)
  group1Nums <- dim(tcgaCluster)[2]
  group2Nums <- dim(gtexCluster)[2]
  group1 <- rep(c("C1"), each = group1Nums)
  group2 <- rep(c("C2"), each = group2Nums)
  group <- c(group1, group2)
  
  groupFile <- data.frame(groupInfo = group)
  row.names(groupFile) <- colnames(tcga_gtex)
  dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(tcga_gtex), 0),
                                colData = groupFile,
                                design = ~groupInfo)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("groupInfo", "C1", "C2"))
  
  diffDenes <- subset(res, padj < padjThreshold & abs(log2FoldChange) > foldChange)
  # write.csv(diffDenes, file = "DESeq2_sig.csv")
  
  DESeq2Res <- as.data.frame(res[order(res$padj), ])
  DESeq2Res$gene <- rownames(DESeq2Res)
  DESeq2Res <- DESeq2Res[Genes249, c(7, 1:6)]
  
  for (i in 1:nrow(DESeq2Res)) {
    if (DESeq2Res[i, "padj"] > padjThreshold) {
      DESeq2Res[i, "reg"] <- "not DE"
    }
    if (DESeq2Res[i, "log2FoldChange"] < 0) {
      DESeq2Res[i, "reg"] <- "not DE"
    }
    if (DESeq2Res[i, "log2FoldChange"] >= 0) {
      DESeq2Res[i, "reg"] <- "not DE"
    }
    if (DESeq2Res[i, "log2FoldChange"] >= 0.5 & DESeq2Res[i, "padj"] <= padjThreshold) {
      DESeq2Res[i, "reg"] <- "min-up-regulated"
    }
    # else{ DESeq2Res[i,'reg'] <- 'not DE'}
    if (DESeq2Res[i, "log2FoldChange"] <= -0.5 & DESeq2Res[i, "padj"] <= padjThreshold) {
      DESeq2Res[i, "reg"] <- "min-down-regulated"
    }
    if (DESeq2Res[i, "log2FoldChange"] >= foldChange & DESeq2Res[i, "padj"] <= padjThreshold) {
      DESeq2Res[i, "reg"] <- "up-regulated"
    }
    # else{ DESeq2Res[i,'reg'] <- 'not DE'}
    if (DESeq2Res[i, "log2FoldChange"] <= -foldChange & DESeq2Res[i, "padj"] <= padjThreshold) {
      DESeq2Res[i, "reg"] <- "down-regulated"
    }
  }
  
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
  # ggsave(paste0(pathSave, "/DESeq2/cluster", j, "/DESeq2_volcano.png"),
  #   DESeq2_volcano,
  #   width = 6,
  #   height = 3.5
  # )
  
  VolcanoPlotList[[j]] <- DESeq2_volcano
  
  new_DESeq2 <- data.frame(
    geneID = dataPlot$gene,
    log2FC = dataPlot$log2FoldChange,
    p_val = dataPlot$pvalue,
    p_val_adj = dataPlot$padj,
    cluster = j
  )
  
  new_DESeq2$label <- ifelse(new_DESeq2$p_val_adj < 0.01, "adjust P-val<0.01", "adjust P-val>=0.01")
  # write.table(new_DESeq2,
  #   paste0(pathSave, "/DESeq2/cluster", j, "/DESeq2Res_sort.csv"),
  #   row.names = FALSE,
  #   sep = ",",
  #   quote = FALSE
  # )
  
  message(paste0("Cluster ", j, " DESeq2 done......"))
}

VolcanoPlotList <- VolcanoPlotList[[1]] +
  VolcanoPlotList[[2]] +
  VolcanoPlotList[[3]] +
  VolcanoPlotList[[4]] +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 13))
VolcanoPlotList

ggsave(paste0(pathSave, "Figure/Figure 3.pdf"), VolcanoPlotList, width = 8, height = 8, dpi = 600)
