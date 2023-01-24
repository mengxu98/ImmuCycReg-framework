

rm(list = ls())

library("tidyr")
library("tidyverse")
library("L0Learn")
library("glmnet")
library("GENIE3")
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
    }else{
      pathway <- paste0(pathRead,"GRN_10_gold_networks/", i)
      dataway <- "/Ecoli-20_nonoise_multifactorial.tsv"
      goldway <- "/Ecoli-20_goldstandard.tsv"
    }
    
    expression_dataset_test <- read.table(paste0(pathway, dataway), header = T) %>% as.matrix() %>% t()
    row.names(expression_dataset_test) <- row.names(expression_dataset_test)
    
    L0GRN <- L0DWGRN(expression_dataset_test)
    evaluationObject <- prepareEval(L0GRN, paste0(pathway, goldway))
    L0_AUROC <- calcAUROC(evaluationObject)
    L0_AUPR <- calcAUPR(evaluationObject)
    L0_AUROC
    
    weightMat <- GENIE3(expression_dataset_test,
                        nCores = 6,
                        verbose = TRUE
    )
    links <- getLinkList(weightMat)
    
    evaluationObject <- prepareEval(links, paste0(pathway, goldway))
    GENIE3_AUROC <- calcAUROC(evaluationObject)
    GENIE3_AUPR <- calcAUPR(evaluationObject)
    GENIE3_AUROC
    
    evaluation_gnw <- data.frame(Dataset = paste0("Net-", i),
                                 `L0Reg framework` = L0_AUROC,
                                 GENIE3 = GENIE3_AUROC)
    evaluation_gnws <- rbind.data.frame(evaluation_gnws, evaluation_gnw)
    
  }
  write.csv(evaluation_gnws, paste0(pathSave, "evaluation_gnw_5_", n,".csv"))
  evaluation_gnws_list[[n]] <- evaluation_gnws
}
