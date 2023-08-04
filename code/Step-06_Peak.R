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
peakName <- rownames(genes_adj_peak)
peakGeneList <- strsplit(as.vector(genes_adj_peak[, 1]), ",")
sampleATAC <- intersect(colnames(tcga_luad), colnames(peak_luad))

geneList <- read.table(paste0(pathRead, "Genes_17.txt"), header = TRUE) %>% .[, 1]

# Run ---------------------------------------------------------------------
check.dir(paste0(pathSave, "BED"))

peakSeqPosAll <- c()
volcanoPlotList <- list()
for (j in 1:length(geneList)) {
  targetGene <- geneList[j]
  message("Running for ", targetGene, "......")
  peaksAroundGene <- c()
  for (i in 1:length(peakName)) {
    if (targetGene %in% peakGeneList[[i]]) {
      peaksAroundGene <- c(peaksAroundGene, peakName[i])
    }
  }
  
  # Choice 1: correlation ---------------------------------------------------
  check.dir(paste0(pathSave, "corrgram/", targetGene))
  
  # select peaks according to distance and score correlation (population level)
  mrna <- t(tcga_luad[targetGene, sampleATAC]) %>% as.data.frame()
  atac <- t(peak_luad[peaksAroundGene, sampleATAC]) %>% as.data.frame()
  peakNameSelected <- colnames(atac)
  
  for (i in 1:length(peakNameSelected)) {
    peak <- atac[, i]
    peak_data <- cbind.data.frame(mrna, peak)
    names(peak_data) <- c(targetGene, peakNameSelected[i])
    p <- ggpairs(peak_data,
                 title = targetGene,
                 upper = list(continuous = wrap("cor",
                                                method = "spearman"
                 ))
    )
    print(p)
    ggsave(paste(pathSave, "corrgram/", targetGene, "/", targetGene, "_", peakNameSelected[i], "_corrgram.png", sep = ""),
           width = 4,
           height = 3
    )
  }
  peakGeneCorR <- c()
  peakGeneCorP <- c()
  for (i in 1:length(peaksAroundGene)) {
    peakGeneCor <- corr.test(atac[, i],
                             mrna,
                             method = "spearman",
                             adjust = "none")
    peakGeneCorR <- c(peakGeneCorR, peakGeneCor$r)
    peakGeneCorP <- c(peakGeneCorP, peakGeneCor$p)
  }
  
  degData <- data.frame(
    Symbol = peaksAroundGene,
    corr = peakGeneCorR,
    pval = peakGeneCorP)
  
  degData$logP <- -log10(degData$pval)
  degData$Group <- "Not-significant"
  degData$Group[which((degData$pval < 0.05) & (degData$corr > 0.1))] <- "Up-regulated"
  degData$Group[which((degData$pval < 0.05) & (degData$corr < -0.1))] <- "Down-regulated"
  
  degData$Label <- ""
  degData <- degData[order(degData$pval), ]
  
  upPeaks <- head(degData$Symbol[which(degData$Group == "Up-regulated")], 20)
  downPeaks <- head(degData$Symbol[which(degData$Group == "Down-regulated")], 20)
  notSigPeaks <- head(degData$Symbol[which(degData$Group == "Not-significant")], 20)
  labelPeaks <- c(as.character(upPeaks), as.character(downPeaks), as.character(notSigPeaks))
  degData$Label[match(labelPeaks, degData$Symbol)] <- labelPeaks
  
  write.csv(degData,
            file = paste(pathSave, "corrgram/", targetGene, "/", targetGene, "_peakGeneCor.csv", sep = ""),
            row.names = F
  )
  
  highCorPeak <- as.character(upPeaks[1])
  downCorPeak <- as.character(downPeaks[1])
  
  # Choice2: distance ----------------------------------------------------------
  
  posColNames <- c("chrom", "chromStart", "chromEnd", "strand")
  tempRes <- geneinfo_df[which(geneinfo_df[, "gene_name"] == targetGene),
                         c("seqnames", "type", "gene_id", "start", "end", "width")]
  tempRes <- tempRes[which(tempRes[, "type"] == "transcript"), ]
  # Select the NEAREST long region
  orderIndex <- order(tempRes[, "start"], decreasing = FALSE)
  mostStart <- tempRes[orderIndex[1], "start"]
  mostEnd <- tempRes[orderIndex[1], "end"]
  mostMid <- (mostStart + mostEnd) / 2
  
  peaksPosTemp <- genes_adj_peak[peaksAroundGene, posColNames]
  peaksPosTemp$mid <- (peaksPosTemp$chromStart + peaksPosTemp$chromEnd) / 2
  peaksPosTemp$mid_distance <- peaksPosTemp$mid - mostMid
  peaksPosTemp$start_distance <- peaksPosTemp$chromEnd - mostStart
  distanceThres <- -8000
  nearestPeaks <- row.names(peaksPosTemp[which(peaksPosTemp$start_distance < 0 &
                                                 peaksPosTemp$start_distance > distanceThres &
                                                 peaksPosTemp$chromEnd < mostStart), ])
  
  nearestPeak <- nearestPeaks[length(nearestPeaks)] # The most nearest peak
  
  # Choice 3: get top peaks in score but score is sample level
  
  orderPeaks <- peaksAroundGene[order(peak_luad[peaksAroundGene, sampleATAC[1]], decreasing = TRUE)]
  highestScorePeak <- orderPeaks[1 : min(1, length(peaksAroundGene))]
  
  # Get the sequence position of peak
  if (length(nearestPeak) != 0) {
    if (nearestPeak %in% c(highCorPeak, downCorPeak)) {
      candidatePeak <- nearestPeak
    }
  } else {
    if (highestScorePeak %in% c(highCorPeak, downCorPeak)) {
      candidatePeak <- highestScorePeak
    } else {
      if (!is.na(highCorPeak)) {
        candidatePeak <- c(highCorPeak)
      } else {
        candidatePeak <- c(downCorPeak)
      }
    }
  }
  candidatePeak <- na.omit(candidatePeak)
  
  peakSeqPos <- genes_adj_peak[candidatePeak, posColNames]
  peakSeqPos[, "chromStart"] <- peakSeqPos[, "chromStart"] - 1
  
  if (nrow(peakSeqPos) > 0) {
    palette = c(`Down-regulated`="#006699",
                `Not-significant`="gray",
                `Up-regulated`="#EE3A8C")
    
    volcanoPlot <- ggscatter(degData,
                             x = "corr",
                             y = "logP",
                             color = "Group",
                             # palette = palette,
                             label = degData$Label,
                             font.label = 8,
                             repel = TRUE,
                             xlab = paste(targetGene, "correlation"),
                             ylab = "-Log10(P-value)",
                             size = 3) +
      scale_color_manual(values = palette) +
      theme_bw() +
      theme(text = element_text(size = 16), legend.position = "bottom") +
      geom_hline(yintercept = 1.5, linetype = "dashed") +
      geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed")
    volcanoPlot
    
    volcanoPlotList[[j]] <- volcanoPlot
    
    peakSeqPos$gene <- targetGene
    peakSeqPosAll <- rbind(peakSeqPosAll, peakSeqPos)
  }
}

check.dir(paste0(pathSave, "Figure/"))
volcanoPlotList[[1]] +
  volcanoPlotList[[2]] +
  volcanoPlotList[[3]] +
  volcanoPlotList[[4]] +
  volcanoPlotList[[5]] +
  volcanoPlotList[[6]] +
  volcanoPlotList[[7]] +
  volcanoPlotList[[8]] +
  plot_layout(ncol = 4)
ggsave(paste0(pathSave, "Figure/Supplementary Figure 1-peaks_1.pdf"),
       width = 12,
       height = 6,
       dpi = 600)

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
ggsave(paste0(pathSave, "Figure/Supplementary Figure 1-peaks_2.pdf"),
       width = 12,
       height = 9,
       dpi = 600)

peakSeqPosAll$peak <- rownames(peakSeqPosAll)
write.table(peakSeqPosAll[, 1:3],
            paste0(pathSave, "all_BED.bed"),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

write.table(peakSeqPosAll[c(5:6, 1:4)],
            paste0(pathSave, "all_peaks.csv"),
            sep = ",",
            quote = FALSE,
            row.names = FALSE)
