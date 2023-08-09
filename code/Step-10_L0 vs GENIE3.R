rm(list = ls())

source("functions/Functions.R")
source("functions/Functions-L0REG.R")

pathRead <- "../data/"
pathSave <- "../../Results/GRN/"
check.dir(pathSave)

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
      pathWay <- paste0(pathRead, "GRN_5_gold_networks/", i)
      matrixWay <- "/measure/measure.tsv"
      goldWay <- "/gold/gold.tsv"
    } else {
      pathWay <- paste0(pathRead, "GRN_10_gold_networks/", i)
      matrixWay <- "/Ecoli-20_nonoise_multifactorial.tsv"
      goldWay <- "/Ecoli-20_goldstandard.tsv"
    }
    
    matrixData <- read.table(paste0(pathWay, matrixWay),
                             header = TRUE) %>% as.matrix() %>% t()
    
    L0REG(t(matrixData),
          L0 = TRUE,
          L0PredSampleMin = 20, L0PredSampleMax = 80,
          L0ExpSampleMin = 20, L0ExpSampleMax = 80,
          L0RankThreshold = 5, L0EnsembleSize = 100,
          SVM = TRUE,
          SVMPredSampleMin = 20, SVMPredSampleMax = 80,
          SVMExpSampleMin = 20, SVMExpSampleMax = 80,
          SVMRankThreshold = 5, SVMEnsembleSize = 100,
          outputFileName = paste0(pathSave, "GNW", datasetNum, "_L0_Net", i,".txt"))
    
    evaluationObject <- prepareEval(paste0(pathSave, "GNW", datasetNum, "_L0_Net", i,".txt"),
                                    paste0(pathWay, goldWay))
    AUROC_L0 <- calcAUROC(evaluationObject)
    AUPRC_L0 <- calcAUPR(evaluationObject)
    AUROC_L0
    
    weightMat <- GENIE3::GENIE3(matrixData,
                                nCores = 3,
                                verbose = TRUE)
    
    weight <- GENIE3::getLinkList(weightMat)
    write.table(weight, file = paste0(pathSave, "GNW", datasetNum, "_GENIE3_Net", i,".txt"))
    evaluationObject <- prepareEval(paste0(pathSave, "GNW", datasetNum, "_GENIE3_Net", i,".txt"),
                                    paste0(pathWay, goldWay))
    AUROC_GENIE3 <- calcAUROC(evaluationObject)
    AUPRC_GENIE3 <- calcAUPR(evaluationObject)
    AUROC_GENIE3
    
    evaluationGNW <- data.frame(Dataset = paste0("Net", i),
                                L0Reg_framework = AUROC_L0,
                                GENIE3 = AUROC_GENIE3)
    evaluationGNWs <- rbind(evaluationGNWs, evaluationGNW)
  }
  write.csv(evaluationGNWs, paste0(pathRead, "AUROC_Net", datasetNum, ".csv"))
  evaluationGNWsList[[n]] <- evaluationGNWs
}

net5AUROC <- melt(evaluationGNWsList[[1]], id = "Dataset", variable.name = "Method", value.name = "AUROC")
net10AUROC <- melt(evaluationGNWsList[[2]], id = "Dataset", variable.name = "Method", value.name = "AUROC")
allNetData <- rbind(net5AUROC, net10AUROC) %>% .[, -1]

color <- c("white", "gray")
p1 <- bar.plot(net5AUROC, barColor = color, yTitle = "AUROC")
p2 <- bar.plot(net10AUROC, barColor = color)
p3 <- box.plot(allNetData, boxColor = color)

p4 <- p1 + p2 + p3 +
  plot_layout(widths = c(1, 2, 1)) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(paste0(pathSave, "../Figure/Supplementary Figure 7.pdf"),
       p4,
       width = 9,
       height = 4,
       dpi =600)
