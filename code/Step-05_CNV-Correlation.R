rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "CNV-LUAD.Rdata"))

sampleCluster <- read.csv(paste0(pathRead, "sample_cluster_4.csv"), header = FALSE)
geneList <- read.table(paste0(pathRead, "Genes_17.txt"), header = TRUE) %>% .[, 1]

correctionResultsAll <- c()
for (k in 1:4) {
  sampleInter <- intersect(sampleCluster[which(sampleCluster[, 2] == k), 1],
                           colnames(tcga_cnv))
  cnvData <- tcga_cnv[geneList, sampleInter] %>% t()
  mrnaData <- tcga_luad[geneList, sampleInter] %>% t()
  mrnaData <- log((mrnaData + 1), 2)
  
  corResults <- c()
  corResultsPval <- c()
  for (i in 1:length(geneList)) {
    corLigands <- corr.test(cnvData[, i],
                            mrnaData[, i],
                            method = "spearman",
                            adjust = "fdr")
    corResults <- c(corResults, corLigands$r)
    corResultsPval <- c(corResultsPval, corLigands$p)
  }
  
  interCorMatrix <- c()
  for (i in 1:length(geneList)) {
    interCorResults <- c()
    for (j in 1:length(geneList)) {
      interCorResults <- c(interCorResults, cor(cnvData[, i],
                                                mrnaData[, j],
                                                method = "spearman",
                                                use = "pairwise.complete.obs"))
    }
    interCorMatrix <- rbind(interCorMatrix, interCorResults)
  }
  rownames(interCorMatrix) <- colnames(cnvData)
  colnames(interCorMatrix) <- colnames(cnvData)
  
  correctionResults <- c()
  for (i in 1:length(geneList)) {
    selectGene <- geneList[i]
    correctionValue <- interCorMatrix[selectGene, selectGene]
    message("The correlation coefficient of ", selectGene, ": ", correctionValue)
    correctionResults <- rbind(correctionResults,
                               c(selectGene, correctionValue)) %>% as.data.frame()
  }
  
  correctionResults$Clsuer <- paste0("Cluster", k)
  correctionResultsAll <- rbind.data.frame(correctionResultsAll, correctionResults)
}
names(correctionResultsAll) <- c("Gene", "CNVs", "Cluster")
write.csv(correctionResultsAll,
          paste0(pathSave, "CNV_correlation_value.csv"),
          row.names = FALSE)
