
rm(list = ls())
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(openxlsx)

path_read <- "data/"
path_save <- "results/"
path_samples <- paste0(path_save, "/NMF/cluster-rank=6/")

genes_2230 <- read.table(paste0(path_read, "Genes_249.txt"))
samples_cluster <- read.csv(paste0(path_samples, "data/sample_cluster_4.csv"), header = F)

tcga_luad <- readRDS(paste0(path_save, "tcga_luad.rds"))
gtex_luad <- readRDS(paste0(path_save, "gtex_luad.rds"))

DESeq2_volcano_list <- list()
DESeq2_all <- c()
for (j in 1:length(table(samples_cluster$V2))) {
  message(paste0("----- Cluster ", j, " DESeq2 start! -----"))

  if (dir.exists(paste0(path_save, "/DESeq2/cluster", j)) == F) {
    dir.create(file.path(paste0(path_save, "/DESeq2/cluster", j)), recursive = TRUE)
  }

  samples <- samples_cluster[which(samples_cluster$V2 == j), j]
  tcga_luad_cluster <- tcga_luad[genes_2230$V1, samples] %>% na.omit()
  gtex_luad_cluster <- gtex_luad[genes_2230$V1, samples] %>% na.omit()

  tcga_gtex <- cbind(tcga_luad_cluster, gtex_luad_cluster)
  group1_no <- dim(tcga_luad_cluster)[2]
  group2_no <- dim(gtex_luad_cluster)[2]
  group1 <- rep(c("C1"), each = group1_no)
  group2 <- rep(c("C2"), each = group2_no)
  group <- c(group1, group2)

  group_file <- data.frame(group_info = group)
  row.names(group_file) <- colnames(tcga_gtex)
  otu_file <- round(as.matrix(tcga_gtex), 0)
  dds <- DESeqDataSetFromMatrix(countData = otu_file, colData = group_file, design = ~group_info)
  dds <- DESeq(dds)

  suppressMessages(dds)
  res <- results(dds, contrast = c("group_info", "C1", "C2"))
  summary(res)
  foldChange <- 1
  padj_threshold <- 0.05
  diff_gene_deseq2 <- subset(res, padj < padj_threshold & abs(log2FoldChange) > foldChange)
  write.csv(diff_gene_deseq2, file = "DESeq2_sig.csv")
  diffUp <- diff_gene_deseq2[(diff_gene_deseq2$padj < padj_threshold & (diff_gene_deseq2$log2FoldChange > foldChange)), ]
  # write.csv(diffUp,file="DESeq2_up.csv")
  diffDown <- diff_gene_deseq2[(diff_gene_deseq2$padj < padj_threshold & (diff_gene_deseq2$log2FoldChange < (-foldChange))), ]
  # write.csv(diffDown,file="DESeq2_down.csv")
  deseq_res <- as.data.frame(res[order(res$padj), ])
  deseq_res$otu_id <- rownames(deseq_res)

  for (i in 1:nrow(deseq_res)) {
    # print(i) #for debugging
    if (deseq_res[i, "padj"] > padj_threshold) {
      deseq_res[i, "reg"] <- "not DE"
    }
    if (deseq_res[i, "log2FoldChange"] < 0) {
      deseq_res[i, "reg"] <- "not DE"
    }
    if (deseq_res[i, "log2FoldChange"] >= 0) {
      deseq_res[i, "reg"] <- "not DE"
    }
    if (deseq_res[i, "log2FoldChange"] >= 0.5 & deseq_res[i, "padj"] <= padj_threshold) {
      deseq_res[i, "reg"] <- "min-up-regulated"
    }
    # else{ deseq_res[i,'reg'] <- 'not DE'}
    if (deseq_res[i, "log2FoldChange"] <= -0.5 & deseq_res[i, "padj"] <= padj_threshold) {
      deseq_res[i, "reg"] <- "min-down-regulated"
    }
    if (deseq_res[i, "log2FoldChange"] >= foldChange & deseq_res[i, "padj"] <= padj_threshold) {
      deseq_res[i, "reg"] <- "up-regulated"
    }
    # else{ deseq_res[i,'reg'] <- 'not DE'}
    if (deseq_res[i, "log2FoldChange"] <= -foldChange & deseq_res[i, "padj"] <= padj_threshold) {
      deseq_res[i, "reg"] <- "down-regulated"
    }
  }

  DESeq2_res <- deseq_res[c(7, 1:8)][, -8]
  gene_list <- read.table("../../data/all genes.txt",
    header = T,
    row.names = 1
  )

  loc <- match(row.names(gene_list), row.names(DESeq2_res))
  DESeq2_res <- DESeq2_res[loc, ]
  # write.table(DESeq2_res, 'DESeq2_res_sort.txt', row.names = FALSE, sep = '\t', quote = FALSE)
  # write.table(DESeq2_res, 'DESeq2_res_sort.csv', row.names = FALSE, sep = ',', quote = FALSE)

  dataset <- as.data.frame(res)

  dataset$change <- ifelse(dataset$pvalue < padj_threshold & abs(dataset$log2FoldChange) >= foldChange,
    ifelse(dataset$log2FoldChange > foldChange,
      "Up (log2FC > 1, P-value < 0.05)",
      "Down (log2FC < -1, P-value < 0.05)"
    ),
    "no sig"
  )

  dataset$gene <- rownames(dataset)
  dataset$label <- ifelse(dataset$pvalue < padj_threshold & abs(dataset$log2FoldChange) >= 1,
    as.character(dataset$gene), ""
  )

  DESeq2_volcano <- ggplot(
    dataset,
    aes(
      x = log2FoldChange,
      y = -log10(pvalue),
      colour = change
    )
  ) +
    geom_point(
      alpha = 0.5,
      size = 3
    ) +
    scale_color_manual(values = c(
      "#006699",
      "#d2dae2",
      "#990033"
    )) +
    geom_vline(
      xintercept = c(-1, 1),
      lty = 4,
      col = "black",
      lwd = 0.8
    ) +
    geom_hline(
      yintercept = -log10(padj_threshold),
      lty = 4,
      col = "black",
      lwd = 0.8
    ) +
    labs(
      x = "log2FoldChange",
      y = "-log10(P-value)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      legend.title = element_blank()
    ) +
    geom_text_repel(
      data = dataset, aes(
        x = log2FoldChange,
        y = -log10(pvalue),
        label = label
      ),
      size = 3,
      box.padding = unit(
        0.5,
        "lines"
      ),
      point.padding = unit(
        0.8,
        "lines"
      ),
      segment.color = "black",
      show.legend = FALSE
    )

  DESeq2_volcano_list[[j]] <- DESeq2_volcano

  # ggsave('DESeq2_volcano.pdf',
  #        DESeq2_volcano,
  #        width = 12,
  #        height = 8)

  new_DESeq2 <- data.frame(
    geneID = dataset$gene,
    log2FC = dataset$log2FoldChange,
    p_val = dataset$pvalue,
    p_val_adj = dataset$padj,
    cluster = j
  )
  library(tidyr)
  library(ggplot2)
  library(tidyverse)
  library(ggrepel)
  new_DESeq2$label <- ifelse(new_DESeq2$p_val_adj < 0.01, "adjust P-val<0.01", "adjust P-val>=0.01")
  DESeq2_all <- rbind.data.frame(DESeq2_all, new_DESeq2)

  print(paste("------------cluster", j, "DESeq2 done !------------", sep = " "))
}

library(patchwork)
DESeq2_volcano_list[[1]] +
  DESeq2_volcano_list[[2]] +
  DESeq2_volcano_list[[3]] +
  DESeq2_volcano_list[[4]] +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") +
    theme(text = element_text(size = 13)) +
    theme(text = element_text(family = "Times New Roman"))

ggsave(paste0("figure/Figure3.png"), width = 8, height = 8, dpi = 600)
