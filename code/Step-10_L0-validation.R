rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

load("../../Results/TCGA-LUAD.Rdata")
tcgaScale <- scale(tcga_luad)

samplesCluster <- read.csv(paste0(pathRead, "sample_cluster_4.csv"), header = F) %>% .[, 1]
geneList <- read.table(paste0(pathRead, "Genes_17.txt"), header = TRUE) %>% .[, 1]

evaluateResultAll <- c()
conditions <- c("randomize", "random", "TCGA-LUAD")

for (c in 1:length(conditions)) {
  condition <- conditions[c]
  check.dir(paste0(pathSave, "Vaildation/", condition))
  
  resultTestCluster <- c()
  resultTrainCluster <- c()
  correctionPlotList <- list()
  for (i in 1:length(geneList)) {
    targetGene <- geneList[i]
    if (file.exists(paste0(pathRead, "TFs/", targetGene, "_TFs_list.txt")) == T) {
      TFsList <- read.table(paste0(pathRead, "TFs/", targetGene, "_TFs_list.txt"),
                            header = T) %>% .[, 1]
      
      X <- tcgaScale[TFsList, samplesCluster] %>% t()
      Y <- tcgaScale[targetGene, samplesCluster] %>% t()
      
      # Data -------------------------------------------------------------------
      if (condition == "TCGA-LUAD") {
        message("Using '", condition, "' data......")
      } else if (condition == "randomize") {
        message("Using '", condition, "' data......")
        X <- NMF::randomize(X)
        Y <- NMF::randomize(Y) %>% as.numeric()
      } else if (condition == "random") {
        message("Using '", condition, "' data......")
        Y <- rnorm(nrow(X), 0, 1)
        for (r in 1:ncol(X)) {
          x <- rnorm(nrow(X), 0, 1)
          X[, r] <- x
        }
      }
      
      # Validation -------------------------------------------------------------
      set.seed(2022)
      trainIdx <- sample(nrow(X), 0.7 * nrow(X))
      trainDataX <- X[trainIdx, ]
      trainDataY <- Y[trainIdx]
      testDataX <- X[-trainIdx, ]
      testDataY <- Y[-trainIdx]
      
      L0Model <- L0Learn::L0Learn.fit(trainDataX,
                                      trainDataY,
                                      penalty = "L0",
                                      maxSuppSize = ncol(trainDataX))
      
      # Extract coefficient at middle lambda
      L0ModelInfor <- as.data.frame(print(L0Model))
      L0ModelInfor <- L0ModelInfor[order(L0ModelInfor$suppSize, decreasing = TRUE), ]
      
      trainDataYPre <- predict(L0Model,
                               newx = trainDataX,
                               lambda = L0ModelInfor$lambda[1],
                               gamma = L0ModelInfor$gamma[1]) %>% as.vector()
      
      RMSE_L0Train <- RMSE(trainDataY, trainDataYPre)
      dataFrameTrain <- cbind(trainDataY, trainDataYPre) %>% as.data.frame()
      colnames(dataFrameTrain) <- c("Raw", "Pre")
      
      corDataTrain <- psych::corr.test(dataFrameTrain$Raw,
                                       dataFrameTrain$Pre)
      
      testDataYPre <- predict(L0Model,
                              newx = testDataX,
                              lambda = L0ModelInfor$lambda[1],
                              gamma = L0ModelInfor$gamma[1]) %>% as.vector()
      
      RMSE_L0Test <- RMSE(testDataY, testDataYPre)
      dataFrameTest <- cbind(testDataY, testDataYPre) %>% as.data.frame()
      colnames(dataFrameTest) <- c("Raw", "Pre")
      
      corDataTest <- psych::corr.test(dataFrameTest$Raw,
                                      dataFrameTest$Pre)
      
      correctionPlot <- ggplot(dataFrameTest, aes(x = Raw, y = Pre)) +
        geom_point() +
        theme_bw() +
        stat_cor(data = dataFrameTest) +
        geom_smooth(formula = 'y ~ x',
                    method = "loess", # lm
                    color = "#006699") +
        labs(x = paste0("Expression of ", targetGene), y = "L0Reg framework")
      correctionPlot
      correctionPlotList[[i]] <- correctionPlot
      
      resultTrain <- data.frame(Gene = targetGene,
                                Corr = corDataTrain$r,
                                RMSD = RMSE_L0Train,
                                pval = corDataTrain$p)
      
      resultTrainCluster <- rbind(resultTrainCluster, resultTrain)
      
      resultTest <- data.frame(Gene = targetGene,
                               Corr = corDataTest$r,
                               RMSD = RMSE_L0Test,
                               pval = corDataTest$p)
      
      resultTestCluster <- rbind(resultTestCluster, resultTest)
    } else {
      message("The TFs file of ", targetGene, " not found......")
    }
  }
  
  p <- multiple.plot(correctionPlotList)
  ggsave(paste0(pathSave, "Vaildation/", condition, "/Cor_", condition,".pdf"),
         p,
         width = 11,
         height = 12,
         dpi = 600)
  
  resultTrainClusterFilter <- c()
  for (i in 1:nrow(resultTrainCluster)) {
    if (resultTrainCluster[i, ]$Corr > 0.7 & resultTrainCluster[i, ]$RMSD < 2) {
      resultTrainClusterFilter <- rbind(resultTrainClusterFilter, resultTrainCluster[i, ])
    }
  }
  
  resultTestClusterFilter <- c()
  for (i in 1:nrow(resultTestCluster)) {
    if (resultTestCluster[i, ]$Corr > 0.7 & resultTestCluster[i, ]$RMSD < 2) {
      resultTestClusterFilter <- rbind(resultTestClusterFilter, resultTestCluster[i, ])
    }
  }
  write.csv(resultTestCluster, paste0(pathSave, "Vaildation/", condition, "/", condition, "_Cor-RMSD.csv"))
  write.csv(resultTestClusterFilter, paste0(pathSave, "Vaildation/", condition, "/", condition, "_Cor-RMSD_filter.csv"))
  
  pRMSD <- ggplot() +
    geom_bar(data = resultTestCluster,
             aes(x = Gene, y = RMSD),
             stat = "identity", position = "dodge", width = 0.6) +
    theme_bw() +
    labs(x = "", y = "RMSD", fill = "", color = "") +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 45, hjust = 1))
  pRMSD
  ggsave(paste0(pathSave, "Vaildation/RMSD.pdf"),
         pRMSD,
         width = 5,
         height = 3,
         dpi = 600)
  
  evaluateResult <- data.frame("Condition" = condition,
                               "train_Corr" = mean(resultTrainCluster$Corr),
                               "train_RMSD" = mean(resultTrainCluster$RMSD),
                               "test_Corr" = mean(resultTestCluster$Corr),
                               "test_RMSD" = mean(resultTestCluster$RMSD),
                               "train_0.7_Corr" = mean(resultTrainClusterFilter$Corr),
                               "train_0.7_RMSD" = mean(resultTrainClusterFilter$RMSD),
                               "test_0.7_Corr" = mean(resultTestClusterFilter$Corr),
                               "test_0.7_RMSD" = mean(resultTestClusterFilter$RMSD))
  evaluateResultAll <- rbind(evaluateResultAll, evaluateResult)
}
write.csv(evaluateResultAll, paste0(pathSave, "Vaildation/Cor-RMSD.csv"))
