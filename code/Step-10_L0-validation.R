rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

load("../../Results/TCGA-LUAD.Rdata")
tcgaScale <- scale(tcga_luad)

samplesCluster <- read.csv(paste0(pathRead, "sample_cluster_4.csv"), header = F) %>% .[, 1]
geneList <- read.table(paste0(pathRead, "Genes_17.txt"), header = TRUE) %>% .[, 1]

conditions <- c("Random", "Randomize", "TCGA_LUAD")

evaluateResultAll <- c()
resultTrainClusterList <- list()
resultTestClusterList <- list()
for (c in 1:length(conditions)) {
  condition <- conditions[c]
  check.dir(paste0(pathSave, "Vaildation/", condition))
  
  resultTestCluster <- c()
  resultTrainCluster <- c()
  correctionPlotListTest <- list()
  correctionPlotListTrain <- list()
  for (i in 1:length(geneList)) {
    targetGene <- geneList[i]
    if (file.exists(paste0(pathRead, "TFs/", targetGene, "_TFs_list.txt")) == T) {
      TFsList <- read.table(paste0(pathRead, "TFs/", targetGene, "_TFs_list.txt"),
                            header = T) %>% .[, 1]
      
      X <- t(tcgaScale[TFsList, samplesCluster])
      Y <- t(tcgaScale[targetGene, samplesCluster])
      
      # Condition set
      if (condition == "TCGA_LUAD") {
        message("Using '", condition, "' data......")
      } else if (condition == "Randomize") {
        message("Using '", condition, "' data......")
        X <- NMF::randomize(X)
        Y <- NMF::randomize(Y) %>% as.numeric()
      } else if (condition == "Random") {
        message("Using '", condition, "' data......")
        Y <- rnorm(nrow(X), 0, 1)
        for (r in 1:ncol(X)) {
          x <- rnorm(nrow(X), 0, 1)
          X[, r] <- x
        }
      }
      
      set.seed(2022)
      
      # Split train and test datasets
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
      L0ModelInfor <- as.data.frame(print(L0Model)) %>% .[order(.$suppSize, decreasing = TRUE), ]
      
      trainDataYPre <- predict(L0Model,
                               newx = trainDataX,
                               lambda = L0ModelInfor$lambda[1],
                               gamma = L0ModelInfor$gamma[1]) %>% as.vector()
      
      RMSE_L0Train <- RMSE(trainDataY, trainDataYPre)
      dataFrameTrain <- cbind(trainDataY, trainDataYPre) %>% as.data.frame()
      colnames(dataFrameTrain) <- c("Raw", "Pre")
      
      corDataTrain <- psych::corr.test(dataFrameTrain$Raw,
                                       dataFrameTrain$Pre)
      
      correctionPlotListTrain[[i]] <- scatter.plot(dataFrameTrain,
                                                   title = paste("Gene:", targetGene),
                                                   xTitle = paste0("True expression"),
                                                   yTitle = "Prediction expression")
      
      testDataYPre <- predict(L0Model,
                              newx = testDataX,
                              lambda = L0ModelInfor$lambda[1],
                              gamma = L0ModelInfor$gamma[1]) %>% as.vector()
      
      RMSE_L0Test <- RMSE(testDataY, testDataYPre)
      dataFrameTest <- cbind(testDataY, testDataYPre) %>% as.data.frame()
      colnames(dataFrameTest) <- c("Raw", "Pre")
      
      corDataTest <- psych::corr.test(dataFrameTest$Raw,
                                      dataFrameTest$Pre)
      
      correctionPlotListTest[[i]] <- scatter.plot(dataFrameTest,
                                                  title = paste("Gene:", targetGene),
                                                  xTitle = "True expression",
                                                  yTitle = "Prediction expression")
      
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
  
  p1 <- combine.multiple.plot(correctionPlotListTest)
  ggsave(paste0(pathSave, "Vaildation/", condition, "/Cor_", condition,"_test.pdf"),
         p1,
         width = 11,
         height = 13,
         dpi = 600)
  
  p2 <- combine.multiple.plot(correctionPlotListTrain)
  ggsave(paste0(pathSave, "Vaildation/", condition, "/Cor_", condition,"_train.pdf"),
         p2,
         width = 11,
         height = 13,
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
  
  resultTrainClusterList[[c]] <- resultTrainCluster
  resultTestClusterList[[c]] <- resultTestCluster
}
write.csv(evaluateResultAll, paste0(pathSave, "Vaildation/Cor-RMSD.csv"))

corPlotDataTrain <- data.frame(Gene = resultTrainClusterList[[1]]$Gene,
                               Random = resultTrainClusterList[[1]]$Corr,
                               Randomize = resultTrainClusterList[[2]]$Corr,
                               TCGA_LUAD = resultTrainClusterList[[3]]$Corr) %>% 
  melt(., id = "Gene", variable.name = "Condition", value.name = "Correlation")

corPlotDataTest <- data.frame(Gene = resultTestClusterList[[1]]$Gene,
                              Random = resultTestClusterList[[1]]$Corr,
                              Randomize = resultTestClusterList[[2]]$Corr,
                              TCGA_LUAD = resultTestClusterList[[3]]$Corr) %>% 
  melt(., id = "Gene", variable.name = "Condition", value.name = "Correlation")

RMSDPlotDataTrain <- data.frame(Gene = resultTrainClusterList[[1]]$Gene,
                                Random = resultTrainClusterList[[1]]$RMSD,
                                Randomize = resultTrainClusterList[[2]]$RMSD,
                                TCGA_LUAD = resultTrainClusterList[[3]]$RMSD) %>% 
  melt(., id = "Gene", variable.name = "Condition", value.name = "RMSD")

RMSDPlotDataTest <- data.frame(Gene = resultTestClusterList[[1]]$Gene,
                               Random = resultTestClusterList[[1]]$RMSD,
                               Randomize = resultTestClusterList[[2]]$RMSD,
                               TCGA_LUAD = resultTestClusterList[[3]]$RMSD) %>% 
  melt(., id = "Gene", variable.name = "Condition", value.name = "RMSD")

p3 <- bar.plot(corPlotDataTrain, title = "Test", yTitle = "Correlation")
p4 <- bar.plot(corPlotDataTest, title = "Test", yTitle = "Correlation")
p6 <- bar.plot(RMSDPlotDataTrain, title = "Test", yTitle = "RMSD")
p7 <- bar.plot(RMSDPlotDataTest, title = "Test", yTitle = "RMSD")

p5 <- p3 + p4 +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p5
ggsave(paste0(pathSave, "Vaildation/Contrast_correlation.pdf"),
       p5,
       width = 7,
       height = 3.5,
       dpi = 600)

p8 <- p6 + p7 +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p8
ggsave(paste0(pathSave, "Vaildation/Contrast_RMSD.pdf"),
       p8,
       width = 7,
       height = 3.5,
       dpi = 600)

p9 <- p3 + p4 + p6 + p7 +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p9
ggsave(paste0(pathSave, "Vaildation/Contrast_correlation&RMSD.pdf"),
       p9,
       width = 7,
       height = 6,
       dpi = 600)
