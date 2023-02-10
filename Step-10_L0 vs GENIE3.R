

rm(list = ls())

library("tidyr")
library("tidyverse")
library("L0Learn")
library("glmnet")
library("GENIE3")
source("Functions3.R")
source("Functions.R")

pathRead <- "data/"
pathSave <- "../Results/"

evaluation_gnws_list <- list()
for (n in 1:2) {
  evaluation_gnws <- c()
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

    expressionData <- read.table(paste0(pathway, dataway), header = T) %>%
      as.matrix() %>%
      t()
    # expressionData <- constructExpressionMatrixFromFile(paste0(pathway, dataway))

    # ptm <- proc.time()
    # L0GRN <- L0DWGRN(expressionData,
    #   cores = 6,
    #   penalty = "L0",
    #   crossInteraction = F
    # )
    # running_time <- proc.time() - ptm
    # print(running_time)
    # 
    # evaluationObject <- caclEval(L0GRN, paste0(pathway, goldway))
    # AUROC_L0 <- calcAUROC(evaluationObject)
    # AUPR_L0 <- calcAUPR(evaluationObject)
    # AUROC_L0

    ptm <- proc.time()
    NIMEFI(t(expressionData),
      GENIE = F,
      SVM = T,
      EL = TRUE,
      outputFileName = "../output.txt",
      ELPredSampleMin = 20, ELPredSampleMax = 80,
      ELExpSampleMin = 20, ELExpSampleMax = 80,
      ELRankThreshold = 5, ELEnsembleSize = 100
    )
    running_time_n <- proc.time() - ptm
    print(running_time_n)
    evaluationObject <- prepareEval("../output.txt", paste0(pathway, goldway))
    AUROC_L0_N <- calcAUROC(evaluationObject)
    AUPR_L0_N <- calcAUPR(evaluationObject)
    AUROC_L0_N

    ptm <- proc.time()
    weightMat <- GENIE3(expressionData,
      nCores = 6,
      verbose = TRUE
    )
    running_time_g <- proc.time() - ptm
    print(running_time_g)

    links <- getLinkList(weightMat)
    evaluationObject <- caclEval(links, paste0(pathway, goldway))
    AUROC_GENIE3 <- calcAUROC(evaluationObject)
    AUPR_GENIE3 <- calcAUPR(evaluationObject)
    AUROC_GENIE3

    evaluation_gnw <- data.frame(
      Dataset = paste0("Net-", i),
      `L0Reg framework` = AUROC_L0_N,
      GENIE3 = AUROC_GENIE3
    )
    evaluation_gnws <- rbind.data.frame(evaluation_gnws, evaluation_gnw)
  }
  write.csv(evaluation_gnws, paste0(pathSave, "evaluation_gnw_", n, ".csv"))
  evaluation_gnws_list[[n]] <- evaluation_gnws
}
