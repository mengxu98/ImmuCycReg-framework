rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/L0Reg-framework/"
penalty <- "L0"

# Load data
load("../../Results/TCGA-LUAD-tpm-tcga-t.Rdata")
tcgaScale <- as.data.frame(scale(tcga_luad))

samplesCluster <- read.csv(paste0(pathRead, "sample_cluster_4.csv"), header = FALSE, row.names = 1)
geneList <- read.table(paste0(pathRead, "Genes_17.txt"), header = TRUE) %>% .[, 1]

TRNClusterFilterAll <- c()
for (j in 1:4) {
  message(paste("Running for cluster", j, "......"))
  check.dir(paste0(pathSave, "cluster", j))
  
  cluster <- rownames(samplesCluster %>% filter(., .[, 1] == j))
  
  selectedTFs <- c()
  TRNClusterFilter <- c()
  for (i in 1:length(geneList)) {
    targetGene <- geneList[i]
    if (file.exists(paste0(pathRead, "TFs/", targetGene, "_TFs_list.txt"))) {
      message(paste0("Running for cluster: ", j, ", ", "NO.", i, " gene: ", targetGene, "......"))
      TFsList <- read.table(paste0(pathRead, "TFs/", targetGene, "_TFs_list.txt"),
                            header = T) %>% .[, 1]
      
      testDataX <- t(tcgaScale[TFsList, cluster])
      testDataY <- t(tcgaScale[targetGene, cluster])
      
      ModelL0 <- L0Learn::L0Learn.fit(testDataX,
                                      testDataY,
                                      penalty = penalty,
                                      maxSuppSize = ncol(testDataX))
      
      ModelL0Infor <- as.data.frame(print(ModelL0)) %>%
        .[order(.$suppSize, decreasing = TRUE), ]
      lambdaL0 <- ModelL0Infor$lambda[1]
      gammaL0 <- ModelL0Infor$gamma[1]
      
      dataYPreL0 <- predict(ModelL0,
                            newx = testDataX,
                            lambda = lambdaL0,
                            gamma = gammaL0) %>% as.vector()
      
      lmFit <- lm(testDataY ~ ., data = data.frame(testDataX))
      lmFitSumm <- as.data.frame(summary(lmFit)$coefficient)
      if (rownames(lmFitSumm)[1] == "(Intercept)") lmFitSumm <- lmFitSumm[-1, ]
      lmFitSumm$Estimate <- coef(ModelL0,
                                 lambda = lambdaL0,
                                 gamma = gammaL0) %>% as.vector() %>% .[-1]
      lmFitSumm$Weight <- abs(lmFitSumm$Estimate)
      lmFitSumm$Weight <- as.numeric(lmFitSumm$Weight / sum(lmFitSumm$Weight))
      lmFitSumm$Regulation <- 1
      lmFitSumm$Regulation[lmFitSumm$Estimate < 0] <- "-1"
      
      lmFitSummFilter <- lmFitSumm[which(lmFitSumm$`Pr(>|t|)` < 0.05), ]
      
      evaluateResultL0 <- evaluate.model(testDataY, dataYPreL0)
      
      if (evaluateResultL0[[4]] < 0.05 & 
          evaluateResultL0[[3]] > 0 &
          nrow(lmFitSummFilter) > 0) {
        singleGeneTRNFilter <- cbind.data.frame("TF" = rownames(lmFitSummFilter),
                                                "Gene" = targetGene,
                                                "Weight" = lmFitSummFilter$Weight,
                                                "Regulation" = lmFitSummFilter$Regulation,
                                                "Cluster" = paste0("Cluster", j))
        
        TRNClusterFilter <- rbind(TRNClusterFilter, singleGeneTRNFilter)
        selectedTFs <- c(selectedTFs, rownames(lmFitSummFilter))
      }
      
    } else {
      print(paste0("No ", targetGene, " TF file......"))
      next
    }
  }
  
  TRNClusterFilterAll <- rbind(TRNClusterFilterAll, TRNClusterFilter)
  
  selectedTFs <- as.data.frame(unique(selectedTFs))
  write.csv(selectedTFs, paste0(pathSave, "cluster", j, "/TFs_list.csv"))
  
  message(paste("Cluster", j, "done......"))
}

write.csv(TRNClusterFilterAll,
          paste0(pathSave, "TRN_cluster_filter_all.csv"))

