

rm(list = ls())

library("psych")
library("ggpubr")
library("ggthemes")
library("rtracklayer")
library("svglite")
source("Functions.R")

pathRead <- "data/"
pathSave <- "../Results/"

String_names <- read.table(paste0(pathSave, "Datasets/STRING_9606.protein.info.v11.5.txt"))
String_score <- read.table(paste0(pathSave, "Datasets/STRING_9606.protein.links.v11.5.txt.gz"), header = T)

ASPAR_Predicted <- read.table(paste0(pathSave, "Datasets/ASPAR Predicted Transcription Factor Targets.txt.gz"),
                              header = T,
                              check.names = FALSE
)
encode_tf <- read.table(paste0(pathSave, "Datasets/ENCODE Transcription Factor Targets.txt.gz"),
                        header = T,
                        check.names = FALSE
)
CHEA_tf <- read.table(paste0(pathSave, "Datasets/CHEA Transcription Factor Targets.txt.gz"),
                      header = T,
                      check.names = FALSE
)
MotifMap_Predicted <- read.table(paste0(pathSave, "Datasets/MotifMap Predicted Transcription Factor Targets.txt.gz"),
                                 header = T,
                                 sep = "\t",
                                 check.names = FALSE
)
TRANSFAC_Curated <- read.table(paste0(pathSave, "Datasets/TRANSFAC Curated Transcription Factor Targets.txt.gz"),
                               header = T,
                               check.names = FALSE
)
TRANSFAC_Predicted <- read.table(paste0(pathSave, "Datasets/TRANSFAC Predicted Transcription Factor Targets.txt.gz"),
                                 header = T,
                                 check.names = FALSE
)

load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "GTEx-LUAD.Rdata"))
load(paste0(pathSave, "CNV-LUAD.Rdata"))
load(paste0(pathSave, "Peak-LUAD.Rdata"))
load(paste0(pathSave, "genes_adj_peak.Rdata"))
load(paste0(pathSave, "geneinfo_df.Rdata"))

candidate_peaks <- read.csv(paste0(pathRead, "all_peaks.csv"))
negative_genes <- read.table(paste0(pathRead, "Genes_immune_cycle.txt"),header = TRUE, sep = "\t")
negative_genes <- negative_genes$GeneSymbol[negative_genes[, c("Direction")] == "negative"]
target_genes_list <- read.table(paste0(pathRead, "Genes_17.txt"), header = T)
target_genes_list <- target_genes_list$gene

format_result_all <- c()
format_result_all_4 <- c()
JustPositiveGene <- c()
JustNegativeGene <- c()
immune_score_results <- c()
immune_score_results1 <- c()
res_frame_all <- c()
num_TFs_all <- c()
t.test_res_inroemation_all <- c()
num_reg_samples_all <- c()
reg_num_all <- c()
gene_samples <- c()
# Run ---------------------------------------------------------------------

for (k in seq_along(target_genes_list)) {
  target_gene <- target_genes_list[k]
  candidate_peak_id <- candidate_peaks[which(candidate_peaks[, "gene"] == target_gene), "peak"]
  gene_list <- read.table(paste0(pathRead, "TFs/", target_gene, "_TFs_list.txt"),
                          row.names = 1
  )
  sample_list <- read.table(paste0(pathRead, "peak_luad_sample.csv"),
                            row.names = 1,
                            header = T
  )
  row.names(sample_list) <- substr(row.names(sample_list), 1, 15)
  col_matrix <- colnames(tcga_luad) %in% row.names(sample_list)
  allgene_tcga_mrna <- tcga_luad[, which(col_matrix)]
  row_matrix_gtex <- row.names(gtex_luad) %in% row.names(gene_list)
  row_matrix_tcga <- row.names(allgene_tcga_mrna) %in% row.names(gene_list)
  gtex_temp <- gtex_luad[which(row_matrix_gtex), ] # all gtex samples
  tcga_temp <- allgene_tcga_mrna[which(row_matrix_tcga), ] # all luad samples both in ATAC and TCGA
  exsit_gene_list <- row.names(tcga_temp)
  
  if (dir.exists(paste0(pathSave, target_gene)) == F) {
    dir.create(paste0(pathSave, target_gene), recursive = TRUE)
  }
  if (!dir.exists(paste0(pathSave, "TFs_list_selected"))) {
    dir.create(paste0(pathSave, "TFs_list_selected"))
  }
  if (!dir.exists(paste0(pathSave, "T-test"))) {
    dir.create(paste0(pathSave, "T-test"))
  }
  
  high_exp_sample <- c()
  down_exp_sample <- c()
  gene_sample <- c()
  for (sample_i in 1:ncol(allgene_tcga_mrna)) {
    target_sample <- colnames(allgene_tcga_mrna)[sample_i]
    t_res <- t.test(gtex_luad[target_gene, ], mu = allgene_tcga_mrna[target_gene, target_sample])
    if (t_res$statistic < 0 & t_res$p.value < 0.001) {
      print("target gene higher expression than gtex samples !")
      high_exp_sample <- c(high_exp_sample, target_sample)
      gene_sample_con <- data.frame(
        Gene = target_gene,
        Sample = target_sample,
        `T` = t_res$statistic,
        `P-value` = t_res$p.value
      )
    }
    
    if (t_res$statistic > 0 & t_res$p.value < 0.001) {
      print("target gene lower expression than gtex samples !")
      down_exp_sample <- c(down_exp_sample, target_sample)
      gene_sample_con <- data.frame(
        Gene = target_gene,
        Sample = target_sample,
        `T` = t_res$statistic,
        `P-value` = t_res$p.value
      )
    }
    
    gene_sample <- rbind.data.frame(gene_sample, gene_sample_con)
  }
  gene_samples <- rbind.data.frame(gene_samples, gene_sample)
  
  if (target_gene %in% negative_genes) {
    sample_num <- length(high_exp_sample)
    t.test_res_inroemation <- c()
    num_reg_samples <- c()
    
    # Positive regulation: TF CNV amplified -> Peak score high + TF mRNA high -> target gene high expression
    results_TF_up <- list()
    results_RP_down <- list()
    if (length(high_exp_sample) > 0) {
      for (sample_i in seq_along(high_exp_sample)) {
        target_sample <- high_exp_sample[sample_i]
        # CONDITION 1: PEAK SCORE IN QUANTILE
        up_thres <- quantile(peak_luad[, target_sample])[3]
        # CONDITION 2: TF MRNA UP
        up_reg_gene <- c()
        down_reg_gene <- c()
        t.test_res_sample <- c()
        
        print(paste(peak_luad[candidate_peak_id, target_sample], up_thres, sep = "<------>"))
        if (peak_luad[candidate_peak_id, target_sample] > up_thres) {
          for (i in seq_along(exsit_gene_list)) {
            t_res_temp <- t.test(gtex_temp[i, ], mu = tcga_temp[i, target_sample])
            
            if (t_res_temp$statistic < 0 & t_res_temp$p.value < 0.001) {
              up_reg_gene <- c(up_reg_gene, exsit_gene_list[i])
            }
            if (t_res_temp$statistic > 0 & t_res_temp$p.value < 0.001) {
              down_reg_gene <- c(down_reg_gene, exsit_gene_list[i])
            }
            
            if (t_res_temp$p.value < 0.001 & target_gene %in% String_names$V2) {
              target_gene_en <- String_names[which(String_names$V2 == target_gene), 1]
              TF_en <- String_names[which(String_names$V2 == exsit_gene_list[i]), 1]
              
              score1 <- String_score[which(String_score$protein1 == target_gene_en), ]
              score2 <- score1[which(score1$protein2 == TF_en), ]
              
              if (length(t(score2)) > 0) {
                score <- as.numeric(score2$combined_score[1]) / 1000
              } else {
                score <- 0
              }
              t.test_res <- data.frame(target_gene, target_sample, rownames(tcga_temp[i, ]), t_res_temp$statistic, t_res_temp$p.value, score)
              t.test_res_sample <- rbind.data.frame(t.test_res_sample, t.test_res)
            } else {
              score <- 0
              t.test_res <- data.frame(target_gene, target_sample, rownames(tcga_temp[i, ]), t_res_temp$statistic, t_res_temp$p.value, score)
              t.test_res_sample <- rbind.data.frame(t.test_res_sample, t.test_res)
            }
          }
        }
        
        CNV_test <- F
        if (CNV_test) {
          up_reg_gene <- up_reg_gene[which(raw_tcga_cnv[up_reg_gene, target_sample] != 0)] # 现在CNV大于0的转录因子是否合理？选择只要是寻找CNV的即作为突变
          down_reg_gene <- down_reg_gene[which(raw_tcga_cnv[down_reg_gene, target_sample] != 0)]
        }
        TF_up_genes <- up_reg_gene
        results_TF_up[[target_sample]] <- TF_up_genes
        TF_down_genes <- down_reg_gene
        results_RP_down[[target_sample]] <- TF_down_genes
        
        if (length(t.test_res_sample) != 0) {
          t.test_res_inroemation <- rbind.data.frame(t.test_res_inroemation, t.test_res_sample)
        }
        num_high <- length(results_TF_up[[target_sample]]) + length(results_RP_down[[target_sample]])
        num_high_sample <- c(target_gene, target_sample, num_high)
        num_reg_samples <- rbind.data.frame(num_reg_samples, num_high_sample)
      }
      names(num_reg_samples) <- c("Gene", "Sample", "Number")
      names(t.test_res_inroemation) <- c("Gene", "Sample", "TF", "T-statistic", "P-value", "Score")
    }
    
    # write.csv(t.test_res_inroemation, paste0(pathSave, "T-test/", target_gene, "_t.test_res_inroemation.csv"), row.names = F)
    
    res_positive_TF_up <- formatPositiveResult(results_TF_up)
    res_negative_RP_down <- formatNegativeResult(results_RP_down)
    
    results_regulatory_frame <- rbind(res_positive_TF_up, res_negative_RP_down)
  } else {
    sample_num <- length(down_exp_sample)
    
    t.test_res_inroemation <- c()
    num_reg_samples <- c()
    
    # negative regulatory: RP CNV amplified -> peak score high + RP mRNA high -> target gene low expression
    results_TF_down <- list()
    results_RP_up <- list()
    if (length(down_exp_sample) > 0) {
      for (sample_i in seq_along(down_exp_sample)) {
        target_sample <- down_exp_sample[sample_i]
        up_reg_gene <- c()
        down_reg_gene <- c()
        t.test_res_sample <- c()
        # CONDITION 1: PEAK SCORE IN QUANTILE
        up_thres <- quantile(peak_luad[, target_sample])[3]
        print(paste(peak_luad[candidate_peak_id, target_sample], up_thres, sep = "<------>"))
        # CONDITION 2: RP mRNA high
        if (peak_luad[candidate_peak_id, target_sample] < up_thres) {
          for (i in seq_along(exsit_gene_list)) {
            t_res_temp <- t.test(gtex_temp[i, ], mu = tcga_temp[i, target_sample])
            
            if (t_res_temp$statistic < 0 & t_res_temp$p.value < 0.001) {
              up_reg_gene <- c(up_reg_gene, exsit_gene_list[i])
            }
            if (t_res_temp$statistic > 0 & t_res_temp$p.value < 0.001) {
              down_reg_gene <- c(down_reg_gene, exsit_gene_list[i])
            }
            if (t_res_temp$p.value < 0.001 & target_gene %in% String_names$V2) {
              target_gene_en <- String_names[which(String_names$V2 == target_gene), 1]
              TF_en <- String_names[which(String_names$V2 == exsit_gene_list[i]), 1]
              
              score1 <- String_score[which(String_score$protein1 == target_gene_en), ]
              score2 <- score1[which(score1$protein2 == TF_en), ]
              
              if (length(t(score2)) > 0) {
                score <- as.numeric(score2$combined_score[1]) / 1000
              } else {
                score <- 0
              }
              t.test_res <- data.frame(target_gene, target_sample, rownames(tcga_temp[i, ]), t_res_temp$statistic, t_res_temp$p.value, score)
              t.test_res_sample <- rbind.data.frame(t.test_res_sample, t.test_res)
            } else {
              score <- 0
              t.test_res <- data.frame(target_gene, target_sample, rownames(tcga_temp[i, ]), t_res_temp$statistic, t_res_temp$p.value, score)
              t.test_res_sample <- rbind.data.frame(t.test_res_sample, t.test_res)
            }
          }
        }
        
        CNV_test <- F
        if (CNV_test) {
          up_reg_gene <- up_reg_gene[which(raw_tcga_cnv[up_reg_gene, target_sample] != 0)]
          down_reg_gene <- down_reg_gene[which(raw_tcga_cnv[down_reg_gene, target_sample] != 0)]
        }
        
        TF_up_genes <- up_reg_gene
        results_RP_up[[target_sample]] <- TF_up_genes
        TF_down_genes <- down_reg_gene
        results_TF_down[[target_sample]] <- TF_down_genes
        
        if (length(t.test_res_sample) != 0) {
          t.test_res_inroemation <- rbind.data.frame(t.test_res_inroemation, t.test_res_sample)
        }
        num_down <- length(results_TF_down[[target_sample]]) + length(results_RP_up[[target_sample]])
        num_down_sample <- c(target_gene, target_sample, num_down)
        num_reg_samples <- rbind.data.frame(num_reg_samples, num_down_sample)
      }
      names(num_reg_samples) <- c("Gene", "Sample", "Number")
      if (length(t.test_res_inroemation) != 0) {
        names(t.test_res_inroemation) <- c("Gene", "Sample", "TF", "T-statistic", "P-value", "Score")
      }
    }
    
    # write.csv(t.test_res_inroemation, paste0(pathSave, "T-test/", target_gene, "_t.test_res_inroemation.csv"), row.names = F)
    
    res_positive_TF_down <- formatPositiveResult(results_TF_down)
    res_negative_RP_up <- formatNegativeResult(results_RP_up)
    results_regulatory_frame <- rbind(res_positive_TF_down, res_negative_RP_up)
  }
  
  #----------------------------------------------------------------------------#
  
  if (length(results_regulatory_frame) == 4) {
    freq_threshold <- 0
    
    results_regulatory_freq <- data.frame(table(results_regulatory_frame$TF))
    results_regulatory_freq$Freq <- results_regulatory_freq$Freq / sample_num
    results_regulatory_freq <- results_regulatory_freq[which(results_regulatory_freq$Freq >= freq_threshold), ]
    #
    sel_col <- c("ASPAR", "ENCODE", "CHEA", "MotifMap", "TRANSFAC_Curated", "TRANSFAC_Predicted")
    results_regulatory_table <- data.frame(Var1 = results_regulatory_freq$Var1) # 建立联立表
    results_regulatory_table <- FrameRegulatoryTable(results_regulatory_table)
    
    write.table(results_regulatory_table,
                file = paste(pathSave, target_gene, "/", target_gene, "_regulatory_table.csv", sep = ""),
                sep = ",",
                row.names = F,
                col.names = T
    )
    
    hitmaxtirx <- subset(results_regulatory_table, select = sel_col)
    hit_tf <- as.vector(results_regulatory_table[which(apply(hitmaxtirx, 1, sum) != 0), 1])
    dataset_num <- 0
    for (t in 1:length(hit_tf)) {
      dataset_num_tf <- nrow(results_regulatory_frame[which(results_regulatory_frame$TF == hit_tf[t]), ])
      dataset_num <- dataset_num + dataset_num_tf
    }
    
    reg_num <- data.frame(
      "Gene" = target_gene,
      "Positive" = table(results_regulatory_frame$type)[1],
      "Negative" = table(results_regulatory_frame$type)[2],
      "ImmuCycReg_framework" = nrow(results_regulatory_frame),
      "Datasets" = dataset_num
    )
    reg_num_all <- rbind.data.frame(reg_num_all, reg_num)
    
    # Update immune score and L0 input TFs
    write.table(hit_tf,
                paste0(pathSave, target_gene, "_TFs_list.txt"),
                quote = F,
                row.names = F,
                col.names = T,
                sep = "\t"
    )
    
    #
    num_reg_TFs <- length(results_regulatory_freq$Var1)
    num_datasets_TFs <- length(hit_tf)
    num_TFs <- c(target_gene, num_reg_TFs, num_datasets_TFs)
    num_TFs_all <- rbind.data.frame(num_TFs_all, num_TFs)
    # format_result_all
    if (length(t.test_res_inroemation) != 0) {
      t.test_res_inroemation_all <- rbind.data.frame(t.test_res_inroemation_all, t.test_res_inroemation)
    }
    num_reg_samples_all <- rbind.data.frame(num_reg_samples_all, num_reg_samples)
    
    format_result_c <- c()
    for (i in seq_along(hit_tf)) {
      format_result <- results_regulatory_frame[which(results_regulatory_frame$TF == hit_tf[i]), ]
      format_result_c <- rbind.data.frame(format_result_c, format_result)
    }
    format_result_all <- rbind(format_result_all, format_result_c)
  }
}

write.csv(reg_num_all, paste0(pathSave, "Edge_reg_num.csv"), row.names = F)

names(num_TFs_all) <- c("Gene", "RegImmuF", "Dataset")
write.table(num_TFs_all,
            paste0(pathSave, "num_TFs_all.csv"),
            sep = ",",
            row.names = F,
            col.names = T
)

write.table(t.test_res_inroemation_all,
            paste0(pathSave, "t.test_res_inroemation_all.csv"),
            sep = ",",
            row.names = F,
            col.names = T
)
