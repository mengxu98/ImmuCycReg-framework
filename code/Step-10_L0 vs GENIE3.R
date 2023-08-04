rm(list = ls())

source("functions/Functions.R")
source("functions/Functions-L0REG.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

evaluationGNWsList <- list()
for (n in 1:2) {
  evaluationGNWs <- c()
  if (n == 1) {
    datasetSeclect <- "TRUE"
    datasetNum <- 5
  } else {
    datasetSeclect <- "FALSE"
    datasetNum <- 10
  }
  for (i in 1:datasetNum) {
    if (datasetSeclect) {
      pathway <- paste0(pathRead, "GRN_5_gold_networks/", i)
      dataway <- "/measure/measure.tsv"
      goldway <- "/gold/gold.tsv"
    } else {
      pathway <- paste0(pathRead, "GRN_10_gold_networks/", i)
      dataway <- "/Ecoli-20_nonoise_multifactorial.tsv"
      goldway <- "/Ecoli-20_goldstandard.tsv"
    }
    
    expressionData <- read.table(paste0(pathway, dataway), header = TRUE) %>%
      as.matrix() %>%
      t()
    
    L0REG(t(expressionData),
          L0 = TRUE,
          SVM = TRUE,
          outputFileName = "../output.txt",
          L0PredSampleMin = 20, L0PredSampleMax = 80,
          L0ExpSampleMin = 20, L0ExpSampleMax = 80,
          L0RankThreshold = 5, L0EnsembleSize = 100)
    
    evaluationObject <- prepareEval("../output.txt", paste0(pathway, goldway))
    AUROC_L0 <- calcAUROC(evaluationObject)
    AUPRC_L0 <- calcAUPR(evaluationObject)
    AUROC_L0
    
    weightMat <- GENIE3::GENIE3(expressionData,
                                nCores = 6,
                                verbose = TRUE)
    
    links <- GENIE3::getLinkList(weightMat)
    evaluationObject <- caclEval(links, paste0(pathway, goldway))
    AUROC_GENIE3 <- calcAUROC(evaluationObject)
    AUPRC_GENIE3 <- calcAUPR(evaluationObject)
    AUROC_GENIE3
    
    evaluationGNW <- data.frame(Dataset = paste0("Net-", i),
                                `L0Reg framework` = AUROC_L0,
                                GENIE3 = AUROC_GENIE3)
    evaluationGNWs <- rbind.data.frame(evaluationGNWs, evaluationGNW)
  }
  write.csv(evaluationGNWs, paste0(pathSave, "evaluationGNW_", n, ".csv"))
  evaluationGNWsList[[n]] <- evaluationGNWs
}
