

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(openxlsx)
workplace <- paste(getwd(),'/results',sep = '')
setwd(workplace)
for (j in 1:4) {
  if (dir.exists(sprintf('workplace/cluster%s',j)) == F) {
    dir.create(file.path(workplace, paste('cluster', j, sep = '')))
  }
  dir <- paste(workplace,'/cluster',j,sep = '')
  setwd(dir)
  load(paste('LUAD_TCGA_GTEX_cluster',j,'.Rdata',sep = ''))
  group_file <- data.frame(group_info=group)
  row.names(group_file) <- colnames(tcga_gtex)
  otu_file <- round(as.matrix(tcga_gtex),0)
  dds <- DESeqDataSetFromMatrix(countData = otu_file, colData = group_file, design = ~group_info)
  dds <- DESeq(dds)
  suppressMessages(dds)
  res <- results(dds, contrast=c('group_info', 'C1', 'C2'))
  summary(res)
  foldChange <- 1
  padj_threshold <- 0.05
  diff_gene_deseq2 <- subset(res, padj < padj_threshold & abs(log2FoldChange) > foldChange)
  write.csv(diff_gene_deseq2,file= "DESeq2_sig.csv")
  diffUp <- diff_gene_deseq2[(diff_gene_deseq2$padj < padj_threshold & (diff_gene_deseq2$log2FoldChange > foldChange)),]
  write.csv(diffUp,file="DESeq2_up.csv")
  diffDown <- diff_gene_deseq2[(diff_gene_deseq2$padj < padj_threshold & (diff_gene_deseq2$log2FoldChange < (-foldChange))),]
  write.csv(diffDown,file="DESeq2_down.csv")
  deseq_res <- as.data.frame(res[order(res$padj), ])
  deseq_res$otu_id <- rownames(deseq_res)
  for (i in 1:nrow(deseq_res)) {
    if (deseq_res[i,'padj'] > padj_threshold ) {
      deseq_res[i,'reg'] <- 'not DE'
    }
    if ( deseq_res[i,'log2FoldChange'] < 0) {
      deseq_res[i,'reg'] <- 'not DE'
    }
    if ( deseq_res[i,'log2FoldChange'] >= 0) {
      deseq_res[i,'reg'] <- 'not DE'
    }
    if (deseq_res[i,'log2FoldChange'] >= 0.5 & deseq_res[i,'padj'] <= padj_threshold) {
      deseq_res[i,'reg'] <- 'min-up-regulated'
    }
    if (deseq_res[i,'log2FoldChange'] <= -0.5 & deseq_res[i,'padj'] <= padj_threshold) {
      deseq_res[i,'reg'] <- 'min-down-regulated'
    }
    if (deseq_res[i,'log2FoldChange'] >= foldChange & deseq_res[i,'padj'] <= padj_threshold) {
      deseq_res[i,'reg'] <- 'up-regulated'
    }
    if (deseq_res[i,'log2FoldChange'] <= -foldChange & deseq_res[i,'padj'] <= padj_threshold) {
      deseq_res[i,'reg'] <- 'down-regulated'
    }
  }
  DESeq2_res <- deseq_res[c(7, 1:8)][,-8]
  gene_list <- read.table('IRgenes.txt',
                          header = T,
                          row.names = 1)
  DESeq2_res <- DESeq2_res[match(row.names(gene_list),row.names(DESeq2_res)),]
  write.table(DESeq2_res, 'DESeq2_res.txt', row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(DESeq2_res, 'DESeq2_res.csv', row.names = FALSE, sep = ',', quote = FALSE)
  dataset <- as.data.frame(res)
  dataset$change = ifelse(dataset$pvalue < padj_threshold & abs(dataset$log2FoldChange) >= foldChange, 
                          ifelse(dataset$log2FoldChange> foldChange ,
                                 'Up(log2FC >= 1, pvalue < 0.05)',
                                 'Down(log2FC <= -1, pvalue < 0.05)'),
                          'no sig')
  dataset$gene <- rownames(dataset)
  dataset$label = ifelse(dataset$pvalue < padj_threshold & abs(dataset$log2FoldChange) >= 1,
                         as.character(dataset$gene),"")
  DESeq2_volcano <- ggplot(
    dataset, 
    aes(x = log2FoldChange, 
        y = -log10(pvalue), 
        colour=change)) +
    geom_point(alpha=0.5, 
               size=3) +
    scale_color_manual(values=c("#006699",
                                "#d2dae2",
                                "#990033")) +
    geom_vline(xintercept=c(-1,1),
               lty=4,
               col="black",
               lwd=0.8) +
    geom_hline(yintercept = -log10(padj_threshold),
               lty=4,
               col="black",
               lwd=0.8) +
    labs(x="log2FoldChange",
         y="-log10(pvalue)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right", 
          legend.title = element_blank()) +
    geom_text_repel(data = dataset, aes(x = log2FoldChange, 
                                          y = -log10(pvalue), 
                                          label = label),
                      size = 3,
                      box.padding = unit(0.5,
                                         "lines"),
                      point.padding = unit(0.8,
                                           "lines"), 
                      segment.color = "black", 
                      show.legend = FALSE)
  ggsave('DESeq2_volcano.pdf',
         DESeq2_volcano,
         width = 12,
         height = 8)
}

