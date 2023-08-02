# 1: select highest score or nearest peak
# 2: make bed file

rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "Peak-LUAD.Rdata"))
load(paste0(pathSave, "ATAC-LUAD.Rdata"))
load(paste0(pathSave, "Geneinfo_df.Rdata"))
row_names <- row.names(genes_adj_peak)
genes_temp <- as.vector(genes_adj_peak[, 1])
genes_list <- strsplit(genes_temp, ",")
sample_name <- intersect(colnames(tcga_luad), colnames(peak_luad))

geneList <- read.table(paste0(pathRead, "Genes_17.txt"), header = TRUE) %>% .[, 1]

# Run ---------------------------------------------------------------------
peakSeqPosAll <- c()
volcanoPlotList <- list()
for (j in 1:length(geneList)) {
  check.dir(paste0(pathSave, "BED"))
  targetGene <- geneList[j]
  peaksAroundGene <- c()
  for (i in 1:nrow(genes_adj_peak)) {
    if (targetGene %in% genes_list[[i]]) {
      peaksAroundGene <- c(peaksAroundGene, row_names[i])
    }
  }
  
  # choice 1: correlation ---------------------------------------------------
  check.dir(paste0(pathSave, "corrgram/", targetGene))
  
  # select peaks according to distance and score correlation (population level)
  Y <- t(tcga_luad[targetGene, sample_name]) %>% as.data.frame()
  X <- t(peak_luad[peaksAroundGene, sample_name]) %>% as.data.frame()
  peak_name <- colnames(X)
  
  for (i in 1:length(peak_name)) {
    peak <- X[, i]
    peak_data <- cbind.data.frame(Y, peak)
    names(peak_data) <- c(targetGene, peak_name[i])
    p <- ggpairs(peak_data,
                 title = targetGene,
                 upper = list(continuous = wrap("cor",
                                                method = "spearman"
                 ))
    )
    print(p)
    ggsave(paste(pathSave, "corrgram/", targetGene, "/", targetGene, "_", peak_name[i], "_corrgram.png", sep = ""),
           width = 4,
           height = 3
    )
  }
  peak_corr_r <- c()
  peak_corr_p <- c()
  for (i in 1:length(peaksAroundGene)) {
    peak_corr <- corr.test(X[, i],
                           Y,
                           method = "spearman",
                           adjust = "none"
    )
    peak_corr_r <- c(peak_corr_r, peak_corr$r)
    peak_corr_p <- c(peak_corr_p, peak_corr$p)
  }
  
  deg.data <- data.frame(
    Symbol = peaksAroundGene,
    corr = peak_corr_r,
    pval = peak_corr_p
  )
  deg.data$logP <- -log10(deg.data$pval)
  deg.data$Group <- "not-significant"
  deg.data$Group[which((deg.data$pval < 0.05) & (deg.data$corr > 0.1))] <- "up-regulated"
  deg.data$Group[which((deg.data$pval < 0.05) & (deg.data$corr < -0.1))] <- "down-regulated"
  
  deg.data$Label <- ""
  deg.data <- deg.data[order(deg.data$pval), ]
  
  up.peak <- head(deg.data$Symbol[which(deg.data$Group == "up-regulated")], 20)
  down.peak <- head(deg.data$Symbol[which(deg.data$Group == "down-regulated")], 20)
  not.sig.peak <- head(deg.data$Symbol[which(deg.data$Group == "not-significant")], 20)
  deg.top.peak <- c(as.character(up.peak), as.character(down.peak), as.character(not.sig.peak))
  deg.data$Label[match(deg.top.peak, deg.data$Symbol)] <- deg.top.peak
  
  write.csv(deg.data,
            file = paste(pathSave, "corrgram/", targetGene, "/", targetGene, "_peak_corr.csv", sep = ""),
            row.names = F
  )
  
  high_cor_peak <- as.character(up.peak[1])
  down_cor_peak <- as.character(down.peak[1])
  
  # choice2: distance----------------------------------------------------------
  # http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/
  
  pos_col_name <- c("chrom", "chromStart", "chromEnd", "strand")
  temp_res <- geneinfo_df[
    which(geneinfo_df[, "gene_name"] == targetGene),
    c("seqnames", "type", "gene_id", "start", "end", "width")
  ]
  temp_res <- temp_res[which(temp_res[, "type"] == "transcript"), ]
  # select the NEAREST long region
  order_index <- order(temp_res[, "start"], decreasing = FALSE)
  mostStart <- temp_res[order_index[1], "start"]
  mostEnd <- temp_res[order_index[1], "end"]
  most_mid <- (mostStart + mostEnd) / 2
  
  peaks_pos_temp <- genes_adj_peak[peaksAroundGene, pos_col_name]
  peaks_pos_temp$mid <- (peaks_pos_temp$chromStart + peaks_pos_temp$chromEnd) / 2
  peaks_pos_temp$mid_distance <- peaks_pos_temp$mid - most_mid
  peaks_pos_temp$start_distance <- peaks_pos_temp$chromEnd - mostStart
  distance_thres <- -8000
  nearest_peaks <- row.names(peaks_pos_temp[which(peaks_pos_temp$start_distance < 0 &
                                                    peaks_pos_temp$start_distance > distance_thres &
                                                    peaks_pos_temp$chromEnd < mostStart), ])
  no_peaks <- length(nearest_peaks)
  nearest_peak <- nearest_peaks[no_peaks] # the most nearest peak
  
  # choice 3:get top peaks in score BUT SCORE IS SAMPLE LEVEL
  top_k <- min(1, length(peaksAroundGene))
  target_sample <- "TCGA-73-A9RS-01"
  
  order_peak_id <- peaksAroundGene[order(peak_luad[peaksAroundGene, target_sample], decreasing = TRUE)]
  highest_score_peak <- order_peak_id[1:top_k]
  
  # get the sequence position of peak
  if (length(nearest_peak) != 0) {
    if (nearest_peak %in% c(high_cor_peak, down_cor_peak)) {
      candidate_peak_id <- nearest_peak
    }
  } else {
    if (highest_score_peak %in% c(high_cor_peak, down_cor_peak)) {
      candidate_peak_id <- highest_score_peak
    } else {
      if (!is.na(high_cor_peak)) {
        candidate_peak_id <- c(high_cor_peak)
      } else {
        candidate_peak_id <- c(down_cor_peak)
      }
    }
  }
  candidate_peak_id <- na.omit(candidate_peak_id)
  
  peak_seq_pos <- genes_adj_peak[candidate_peak_id, pos_col_name]
  peak_seq_pos[, "chromStart"] <- peak_seq_pos[, "chromStart"] - 1
  
  check.dir(paste0(pathSave, "volcano/"))
  
  if (nrow(peak_seq_pos) > 0) {
    volcanoPlot <- ggscatter(deg.data,
                             x = "corr",
                             y = "logP",
                             color = "Group",
                             palette = c("#006699", "#008B00", "#EE3A8C"),
                             label = deg.data$Label,
                             font.label = 8,
                             repel = T,
                             xlab = paste(targetGene, "Correlation"),
                             ylab = "-Log10(P-value)",
                             size = 3) +
      theme_bw() +
      theme(text = element_text(size = 16), legend.position = "none") +
      geom_hline(yintercept = 1.30, linetype = "dashed") +
      geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed")
    volcanoPlot
    ggsave(paste0(pathSave, "volcano/", targetGene, "_volcano.pdf"),
           volcanoPlot,
           width = 6.5,
           height = 5)
    volcanoPlotList[[j]] <- volcanoPlot
    
    write.table(peak_seq_pos[, 1:3],
                file = paste0(pathSave, "BED/", targetGene, "_around_peak.bed"),
                sep = "\t",
                row.names = F,
                col.names = F,
                quote = F)
    
    peak_seq_pos$gene <- targetGene
    peakSeqPosAll <- rbind(peakSeqPosAll, peak_seq_pos)
  }
}

volcanoPlotList[[1]] +
  volcanoPlotList[[2]] +
  volcanoPlotList[[3]] +
  volcanoPlotList[[4]] +
  volcanoPlotList[[5]] +
  volcanoPlotList[[6]] +
  volcanoPlotList[[7]] +
  volcanoPlotList[[8]] +
  plot_layout(ncol = 4)
ggsave(paste0(pathSave, "Figure/Supplementary Figure 1-peaks_1.pdf"), width = 12, height = 6, dpi = 600)

volcanoPlotList[[9]] +
  volcanoPlotList[[10]] +
  volcanoPlotList[[11]] +
  volcanoPlotList[[12]] +
  volcanoPlotList[[13]] +
  volcanoPlotList[[14]] +
  volcanoPlotList[[15]] +
  volcanoPlotList[[16]] +
  volcanoPlotList[[17]] +
  plot_layout(ncol = 4)
ggsave(paste0(pathSave, "Figure/Supplementary Figure 1-peaks_2.pdf"), width = 12, height = 9, dpi = 600)

peakSeqPosAll$peak <- row.names(peakSeqPosAll)
write.table(peakSeqPosAll[, 1:3],
            paste0(pathSave, "all_BED.bed"),
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F
)
write.table(peakSeqPosAll[c(5:6, 1:4)],
            paste0(pathSave, "all_peaks.csv"),
            sep = ",",
            quote = F,
            row.names = F
)
