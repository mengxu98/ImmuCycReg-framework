rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

stringData <- read.table(paste0(pathRead, "Databases/STRING_9606.protein.v11.5.txt.gz"))

databases <- c("ASPAR", "ENCODE", "CHEA", "MotifMap", "TRANSFAC_Curated", "TRANSFAC_Predicted")
databasesList <- list()
for (database in databases) {
  databasesList[[database]] <- read.table(paste0(pathRead, "Databases/", database, ".txt.gz"),
                                          header = TRUE,
                                          sep = "\t",
                                          check.names = FALSE)
}

load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "GTEx-LUAD.Rdata"))
load(paste0(pathSave, "CNV-LUAD.Rdata"))
load(paste0(pathSave, "Peak-LUAD.Rdata"))
load(paste0(pathSave, "genes_adj_peak.Rdata"))
load(paste0(pathSave, "geneinfo_df.Rdata"))

candidatePeaks <- read.csv(paste0(pathSave, "all_peaks.csv"))
negativeGenes <- read.table(paste0(pathRead, "Genes_immune_cycle.txt"),
                            header = TRUE,
                            sep = "\t") %>% filter(., Direction == "negative") %>% .[, 1]
targetGenesList <- read.table(paste0(pathRead, "Genes_17.txt"),
                              header = TRUE) %>% .[, 1]

cnvTest <- FALSE

check.dir(paste0(pathSave, "Regulatory_table"))
check.dir(paste0(pathSave, "TFs_list_selected"))
check.dir(paste0(pathSave, "T-test"))

formatResult_all <- c()
tTestResInforAll <- c()
numRegSamplesAll <- c()
edgeNumAll <- c()
geneTFsNodeNumAll <- c()
genesSamplesStatistic <- c()
for (k in seq_along(targetGenesList)) {
  targetGene <- targetGenesList[k]
  message("Running for gene: ", targetGene, " ......")
  candidatePeakID <- candidatePeaks %>% filter(., gene == targetGene) %>% .[, "peak"]
  
  TFsList <- read.table(paste0(pathRead, "TFs/", targetGene, "_TFs_list.txt"),
                        header = T) %>% .[, 1]
  sampleList <- read.table(paste0(pathRead, "peak_luad_sample.csv"),
                            header = T) %>% .[, 1]
  
  gtex_temp <- gtex_luad[TFsList, ]
  tcga_temp <- tcga_luad[TFsList, sampleList]
  
  highExpSamples <- c()
  downExpSamples <- c()
  geneSampleStatistic <- c()
  for (sample_i in seq_along(sampleList)) {
    targetSample <- sampleList[sample_i]
    tRes <- t.test(gtex_luad[targetGene, ], mu = tcga_luad[targetGene, targetSample])
    if (tRes$statistic < 0 & tRes$p.value < 0.001) {
      message("The expression of ", targetGene, " in sample: '", targetSample,
              "' is higher expression than GTEx samples......")
      highExpSamples <- c(highExpSamples, targetSample)
    }
    
    if (tRes$statistic > 0 & tRes$p.value < 0.001) {
      message("The expression of ", targetGene, " in sample: '", targetSample,
              "' is lower expression than GTEx samples......")
      downExpSamples <- c(downExpSamples, targetSample)
    }
    
    geneSampleStatistic <- rbind(geneSampleStatistic,
                                 data.frame(Gene = targetGene,
                                            Sample = targetSample,
                                            `T` = tRes$statistic,
                                            `P-value` = tRes$p.value))
  }
  genesSamplesStatistic <- rbind(genesSamplesStatistic,
                                 geneSampleStatistic)

  if (targetGene %in% negativeGenes) {
    if (length(highExpSamples) <= 0) next
    
    tTestResInfor <- c()
    numRegSamples <- c()
    
    # Positive regulation: 
    #   TF CNV amplified -> Peak score high + TF mRNA high -> target gene high expression
    resultsTFsUp <- list()
    resultsRPsdown <- list()
    
    for (sample_i in seq_along(highExpSamples)) {
      targetSample <- highExpSamples[sample_i]
      
      # Condition 1: PEAK SCORE IN QUANTILE
      upThres <- quantile(peak_luad[, targetSample])[3]
      
      # Condition 2: TF MRNA UP
      upRegTFs <- c()
      downRegTFs <- c()
      tTestResSamples <- c()
      
      peakValue <- peak_luad[candidatePeakID, targetSample]
      message(peakValue, " <------> ", upThres)
      if (peakValue > upThres) {
        for (i in seq_along(TFsList)) {
          tRes_temp <- t.test(gtex_temp[TFsList[i], ],
                              mu = tcga_temp[TFsList[i], targetSample])
          
          if (tRes_temp$statistic < 0 & tRes_temp$p.value < 0.001) {
            upRegTFs <- c(upRegTFs, TFsList[i])
          }
          if (tRes_temp$statistic > 0 & tRes_temp$p.value < 0.001) {
            downRegTFs <- c(downRegTFs, TFsList[i])
          }
          
          if (tRes_temp$p.value < 0.001 & targetGene %in% stringData[, 2]) {
            score <- stringData %>% filter(., .[, 2] == targetGene & .[, 1] == TFsList[i]) %>% .[, 3]
            if (is_empty(score)) score <- 0
          } else {
            score <- 0
          }
          
          tTestResSamples <- rbind(tTestResSamples,
                                   data.frame(TF = rownames(tcga_temp[TFsList[i], ]),
                                              Gene = targetGene,
                                              Sample = targetSample,
                                              `T-statistic` = tRes_temp$statistic,
                                              `P-value` = tRes_temp$p.value,
                                              Score = score))
          
        }
        
      }
      
      if (cnvTest) {
        upRegTFs <- upRegTFs[which(tcga_cnv[upRegTFs, targetSample] != 0)]
        downRegTFs <- downRegTFs[which(tcga_cnv[downRegTFs, targetSample] != 0)]
      }
      
      resultsTFsUp[[targetSample]] <- upRegTFs
      resultsRPsdown[[targetSample]] <- downRegTFs
      
      if (length(tTestResSamples) != 0) tTestResInfor <- rbind(tTestResInfor, tTestResSamples)
      
      num_high <- length(resultsTFsUp[[targetSample]]) + length(resultsRPsdown[[targetSample]])
      numRegSamples <- rbind(numRegSamples,
                             data.frame("Gene" = targetGene,
                                        "Sample" = targetSample,
                                        "Number" = num_high))
    }
    
    resultsRegulatoryFrame <- rbind(formatPositiveResult(resultsTFsUp, targetGene),
                                    formatNegativeResult(resultsRPsdown, targetGene))
  } else {
    
    tTestResInfor <- c()
    numRegSamples <- c()
    
    # Negative regulatory: RP CNV amplified -> peak score high + RP mRNA high -> target gene low expression
    resultsTFdown <- list()
    resultsRPup <- list()
    if (length(downExpSamples) > 0) {
      for (sample_i in seq_along(downExpSamples)) {
        targetSample <- downExpSamples[sample_i]
        upRegTFs <- c()
        downRegTFs <- c()
        tTestResSamples <- c()
        
        # CONDITION 1: PEAK SCORE IN QUANTILE
        upThres <- quantile(peak_luad[, targetSample])[3]
        peakValue <- peak_luad[candidatePeakID, targetSample]
        message(peakValue, " <------> ", upThres)
        # CONDITION 2: RP mRNA high
        if (peakValue < upThres) {
          for (i in seq_along(TFsList)) {
            tRes_temp <- t.test(gtex_temp[TFsList[i], ],
                                mu = tcga_temp[TFsList[i], targetSample])
            
            if (tRes_temp$statistic < 0 & tRes_temp$p.value < 0.001) {
              upRegTFs <- c(upRegTFs, TFsList[i])
            }
            if (tRes_temp$statistic > 0 & tRes_temp$p.value < 0.001) {
              downRegTFs <- c(downRegTFs, TFsList[i])
            }
            
            if (tRes_temp$p.value < 0.001 & targetGene %in% stringData[, 2]) {
              score <- stringData %>% filter(., .[, 2] == targetGene & .[, 1] == TFsList[i]) %>% .[, 3]
              if (is_empty(score)) score <- 0
            } else {
              score <- 0
            }
            
            tTestResSamples <- rbind(tTestResSamples,
                                     data.frame(TF = rownames(tcga_temp[TFsList[i], ]),
                                                Gene = targetGene,
                                                Sample = targetSample,
                                                `T-statistic` = tRes_temp$statistic,
                                                `P-value` = tRes_temp$p.value,
                                                Score = score))
            
          }
        }
        
        if (cnvTest) {
          upRegTFs <- upRegTFs[which(tcga_cnv[upRegTFs, targetSample] != 0)]
          downRegTFs <- downRegTFs[which(tcga_cnv[downRegTFs, targetSample] != 0)]
        }
        
        resultsRPup[[targetSample]] <- upRegTFs
        resultsTFdown[[targetSample]] <- downRegTFs
        
        if (length(tTestResSamples) > 0) tTestResInfor <- rbind(tTestResInfor, tTestResSamples)
        
        num_down <- length(resultsTFdown[[targetSample]]) + length(resultsRPup[[targetSample]])
        numRegSamples <- rbind(numRegSamples,
                               data.frame(Gene = targetGene,
                                          Sample = targetSample,
                                          Number = num_down))
      }
    }
    
    resultsRegulatoryFrame <- rbind(formatPositiveResult(resultsTFdown, targetGene),
                                    formatNegativeResult(resultsRPup, targetGene))
  }
  
  #----------------------------------------------------------------------------#
  
  if (ncol(resultsRegulatoryFrame) == 4) {
    freq_threshold <- 0
    
    resultsRegulatoryFreq <- data.frame(table(resultsRegulatoryFrame$TF))
    resultsRegulatoryFreq$Freq <- resultsRegulatoryFreq$Freq / length(downExpSamples) # length(highExpSamples)
    resultsRegulatoryFreq <- resultsRegulatoryFreq[which(resultsRegulatoryFreq$Freq >= freq_threshold), ]
    #
    resultsRegulatoryTable <- data.frame(TF = resultsRegulatoryFreq$Var1)
    resultsRegulatoryTable <- frame.regulatory.table(resultsRegulatoryTable,
                                                     targetGene,
                                                     databasesList)
    
    write.csv(resultsRegulatoryTable,
              file = paste0(pathSave, "Regulatory_table/", targetGene, "_regulatory_table.csv"))
    
    hitTFs <- resultsRegulatoryTable[which(apply(subset(resultsRegulatoryTable,
                                                        select = names(databasesList)),
                                                 1,
                                                 function (x) {
                                                   sum(as.numeric(x))
                                                 }) != 0),
                                     1] %>% as.vector()
    
    dataset_num <- 0
    for (t in 1:length(hitTFs)) {
      dataset_num <- dataset_num + nrow(resultsRegulatoryFrame[which(resultsRegulatoryFrame$TF == hitTFs[t]), ])
    }
    
    edgeNumAll <- rbind(edgeNumAll,
                        data.frame("Gene" = targetGene,
                                   "Positive" = table(resultsRegulatoryFrame$type)[1],
                                   "Negative" = table(resultsRegulatoryFrame$type)[2],
                                   "ImmuCycReg_framework" = nrow(resultsRegulatoryFrame),
                                   "Database" = dataset_num))
    
    # Update TFs list for L0
    write.table(hitTFs,
                paste0(pathSave, "TFs_list_selected/", targetGene, "_TFs_list.txt"),
                quote = F,
                row.names = F,
                col.names = T,
                sep = "\t")
    
    geneTFsNodeNumAll <- rbind(geneTFsNodeNumAll,
                               data.frame("Gene" = targetGene,
                                          "ImmuCyc_framework" = nrow(resultsRegulatoryFreq),
                                          "Database" = length(hitTFs)))
    
    if (length(tTestResInfor) != 0) {
      tTestResInforAll <- rbind(tTestResInforAll, tTestResInfor)
    }
    numRegSamplesAll <- rbind(numRegSamplesAll, numRegSamples)
    
    formatResult <- c()
    for (i in seq_along(hitTFs)) {
      formatResult <- rbind(formatResult,
                            resultsRegulatoryFrame[which(resultsRegulatoryFrame$TF == hitTFs[i]), ])
    }
    formatResult_all <- rbind(formatResult_all, formatResult)
  }
}

write.csv(genesSamplesStatistic,
          paste0(pathSave, "ImmuCycRegFramework_Gene-Sample.csv"),
          row.names = FALSE)

write.csv(tTestResInforAll,
          paste0(pathSave, "T-test/T-test_results_information.csv"),
          row.names = FALSE)

write.csv(edgeNumAll,
          paste0(pathSave, "Gene_TFs_edge_num.csv"),
          row.names = FALSE)

write.csv(geneTFsNodeNumAll,
          paste0(pathSave, "Gene_TFs_node_num.csv"),
          row.names = FALSE)

