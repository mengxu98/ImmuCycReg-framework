

rm(list = ls())
library(psych)
library(ggpubr)
library(ggthemes)
library(rtracklayer)
library(svglite)
source("Functions.R")
#------------------------------------------------------------------------------#
String_names <- read.table("../  data/String/9606.protein.info.v11.5.txt")
String_score <- read.table("../  data/String/9606.protein.links.v11.5.txt.gz", header = T)

ASPAR_Predicted <- read.table("../  data/TF-targets/ASPAR Predicted Transcription Factor Targets.txt.gz",
  header = T,
  check.names = FALSE
)
encode_tf <- read.table("../  data/TF-targets/ENCODE Transcription Factor Targets.txt.gz",
  header = T,
  check.names = FALSE
)
CHEA_tf <- read.table("../  data/TF-targets/CHEA Transcription Factor Targets.txt.gz",
  header = T,
  check.names = FALSE
)
MotifMap_Predicted <- read.table("../  data/TF-targets/MotifMap Predicted Transcription Factor Targets.txt.gz",
  header = T,
  sep = "\t",
  check.names = FALSE
)
TRANSFAC_Curated <- read.table("../  data/TF-targets/TRANSFAC Curated Transcription Factor Targets.txt.gz",
  header = T,
  check.names = FALSE
)
TRANSFAC_Predicted <- read.table("../  data/TF-targets/TRANSFAC Predicted Transcription Factor Targets.txt.gz",
  header = T,
  check.names = FALSE
)

load("data/luad-rsem-count-tcga-t.Rdata")
load("data/raw_tcga_gtex_mrna.RData")
load("data/raw_tcga_cnv.Rdata")
load("data/peak_luad.Rdata")
load("data/genes_adj_peak.Rdata")
load("data/allgene_tcga_mrna.Rdata")
load("data/geneinfo_df.Rdata")

candidate_peaks <- read.csv("all_peaks.csv")
negative_gene <- c("CCL22", "CCL20", "CCL28", "EZH2", "NT5E", "LAG3", "RAET1G")

target_gene_list <- read.table("genes.txt", header = T)
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

for (k in seq_along(t(target_gene_list))) {
  target_gene <- t(target_gene_list)[, k]
  candidate_peak_id <- candidate_peaks[which(candidate_peaks[, "gene"] == target_gene), "peak"]
  gene_list <- read.table(paste("data/TFs_list/", target_gene, ".txt", sep = ""),
    row.names = 1
  )
  sample_list <- read.table("peak_luad_sample.csv",
    row.names = 1,
    header = T
  )
  row.names(sample_list) <- substr(row.names(sample_list), 1, 15)
  col_matrix <- colnames(raw_tcga_mrna) %in% row.names(sample_list)
  allgene_tcga_mrna <- raw_tcga_mrna[, which(col_matrix)]
  row_matrix_gtex <- row.names(raw_gtex_mrna) %in% row.names(gene_list)
  row_matrix_tcga <- row.names(allgene_tcga_mrna) %in% row.names(gene_list)
  gtex_temp <- raw_gtex_mrna[which(row_matrix_gtex), ] # all gtex samples
  tcga_temp <- allgene_tcga_mrna[which(row_matrix_tcga), ] # all luad samples both in ATAC and TCGA
  exsit_gene_list <- row.names(tcga_temp)

  if (dir.exists(paste("results/", target_gene, sep = "")) == F) {
    dir.create(paste("results/", target_gene, sep = ""), recursive = TRUE)
  }
  if (!dir.exists("../results ATAC-seq/TFs_list_selected")) {
    dir.create("../results ATAC-seq/TFs_list_selected")
  }
  if (!dir.exists("../results ATAC-seq/T-test")) {
    dir.create("../results ATAC-seq/T-test")
  }

  high_exp_sample <- c()
  down_exp_sample <- c()
  gene_sample <- c()
  for (sample_i in 1:ncol(allgene_tcga_mrna)) {
    target_sample <- colnames(allgene_tcga_mrna)[sample_i]
    t_res <- t.test(raw_gtex_mrna[target_gene, ], mu = allgene_tcga_mrna[target_gene, target_sample])
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

  if (target_gene %in% negative_gene) {
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

    write.csv(t.test_res_inroemation, paste0("T-test/", target_gene, "_t.test_res_inroemation.csv"), row.names = F)

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

    write.csv(t.test_res_inroemation, paste0("T-test/", target_gene, "_t.test_res_inroemation.csv"), row.names = F)

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
      file = paste("results/", target_gene, "/", target_gene, "_regulatory_table.csv", sep = ""),
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
      paste("../L0/data/", target_gene, "_TFs_list.txt", sep = ""),
      quote = F,
      row.names = F,
      col.names = T,
      sep = "\t"
    )
    write.table(hit_tf,
      paste("TFs_list_selected/", target_gene, "_TFs_list.txt", sep = ""),
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

write.csv(reg_num_all, paste0("results/", "Edge_reg_num.csv"), row.names = F)

names(num_TFs_all) <- c("Gene", "RegImmuF", "Dataset")
write.table(num_TFs_all,
  "results/num_TFs_all.csv",
  sep = ",",
  row.names = F,
  col.names = T
)

write.table(t.test_res_inroemation_all,
  "results/t.test_res_inroemation_all.csv",
  sep = ",",
  row.names = F,
  col.names = T
)

# TRN ---------------------------------------------------------------------
format_result_all <- na.omit(format_result_all)
for (i in 1:nrow(format_result_all)) {
  if (format_result_all[i, "target_gene"] == "B2M") {
    format_result_all[i, "immuneCycle"] <- "Recognition of cancer cells by T cells"
  }
  if (format_result_all[i, "target_gene"] == "CXCL17") {
    format_result_all[i, "immuneCycle"] <- "MDSC"
  }
  if (format_result_all[i, "target_gene"] == "CCL20") {
    format_result_all[i, "immuneCycle"] <- "TH cell"
  }
  if (format_result_all[i, "target_gene"] == "CCL22") {
    format_result_all[i, "immuneCycle"] <- "Macrophage.M2"
  }
  if (format_result_all[i, "target_gene"] == "CCL28") {
    format_result_all[i, "immuneCycle"] <- "Treg cell"
  }
  # if (format_result_all[i,'target_gene'] == 'CCR5') {
  #   format_result_all[i,'immuneCycle'] <- 'NKT cell'
  # }
  if (format_result_all[i, "target_gene"] == "CCR6") {
    format_result_all[i, "immuneCycle"] <- "TH cell"
  }
  if (format_result_all[i, "target_gene"] == "CD1C") {
    format_result_all[i, "immuneCycle"] <- "T cell"
  }
  if (format_result_all[i, "target_gene"] == "CD2") {
    format_result_all[i, "immuneCycle"] <- "Priming and activation"
  }
  if (format_result_all[i, "target_gene"] == "CD40") {
    format_result_all[i, "immuneCycle"] <- "Cancer antigen presentation"
  }
  if (format_result_all[i, "target_gene"] == "CXCR5") {
    format_result_all[i, "immuneCycle"] <- "B cell"
  }
  if (format_result_all[i, "target_gene"] == "EZH2") {
    format_result_all[i, "immuneCycle"] <- "Infiltration of immune cells into tumors"
  }
  if (format_result_all[i, "target_gene"] == "HLA-A") {
    format_result_all[i, "immuneCycle"] <- "Cancer antigen presentation"
  }
  if (format_result_all[i, "target_gene"] == "HSP90B1") {
    format_result_all[i, "immuneCycle"] <- "Release of cancer cell antigens"
  }
  if (format_result_all[i, "target_gene"] == "HSPA4") {
    format_result_all[i, "immuneCycle"] <- "Release of cancer cell antigens"
  }
  if (format_result_all[i, "target_gene"] == "KLRK1") {
    format_result_all[i, "immuneCycle"] <- "NK cell"
  }
  if (format_result_all[i, "target_gene"] == "LAG3") {
    format_result_all[i, "immuneCycle"] <- "Killing of cancer cells"
  }
  if (format_result_all[i, "target_gene"] == "NT5E") {
    format_result_all[i, "immuneCycle"] <- "Infiltration of immune cells into tumors"
  }
  if (format_result_all[i, "target_gene"] == "RAET1G") {
    format_result_all[i, "immuneCycle"] <- "Priming and activation"
  }
  if (format_result_all[i, "target_gene"] == "TAP1") {
    format_result_all[i, "immuneCycle"] <- "Recognition of cancer cells by T cells"
  }
}

format_result_all$Score <- ""
format_result_all$Pvalue <- ""
for (t in 1:nrow(format_result_all)) {
  gene_t <- format_result_all$target_gene[t]
  tf_t <- format_result_all$TF[t]
  sample_t <- format_result_all$sample[t]

  score_t1 <- t.test_res_inroemation_all[which(t.test_res_inroemation_all$Gene == gene_t), ]
  score_t2 <- score_t1[which(score_t1$TF == tf_t), ]
  score_t3 <- score_t2[which(score_t2$Sample == sample_t), ]
  if (nrow(score_t3) > 0) {
    score_t <- score_t3$Score
    format_result_all$Score[t] <- score_t
  }
  pvalue <- score_t3$`P-value`
  format_result_all$Pvalue[t] <- pvalue
}

samples_score <- cbind.data.frame(format_result_all$sample, format_result_all$Score)
samples_score <- aggregate(as.numeric(samples_score[, 2]),
  by = list(SampleSelected = samples_score$`format_result_all$Score`), FUN = sum, simplify = TRUE
)

Score_genes <- c()
Score_genes_samples <- c()
for (k in 1:nrow(target_gene_list)) {
  target_gene <- t(target_gene_list)[, k]
  score_gene <- format_result_all[which(format_result_all$target_gene == target_gene), ]
  if (nrow(score_gene) > 0) {
    score_gene$FDR <- p.adjust(score_gene$Pvalue,
      method = "BH"
    )
    score_gene <- score_gene[which(score_gene$FDR < 0.001), ]

    freq_threshold <- 0
    samples_table <- data.frame(table(score_gene$sample))
    samples_table$Freq_com <- samples_table$Freq / nrow(samples_table) / 10

    TFs_table <- data.frame(table(score_gene$TF))
    TFs_table$Freq_com <- TFs_table$Freq / nrow(TFs_table)

    for (i in 1:nrow(score_gene)) {
      TF_sel <- score_gene$TF[i]
      samples_sels <- score_gene[which(score_gene$TF == TF_sel), "sample"]
      for (j in length(samples_sels)) {
        samples_sel <- samples_sels[j]
        sample_sel_score <- samples_table[which(samples_table$Var1 == samples_sel), "Freq_com"]
      }
      TF_sel_score <- TFs_table[which(TFs_table$Var1 == TF_sel), "Freq_com"]
      score_tf_sample <- sample_sel_score * TF_sel_score
      score_gene$Sample_score[i] <- as.numeric(score_gene$Score[i]) + score_tf_sample
    }

    write.csv(score_gene, paste0("results/", target_gene, "_TRN.csv"), row.names = F, quote = F)

    Score_genes_samples <- rbind.data.frame(Score_genes_samples, score_gene)


    Score <- sum(as.numeric(score_gene$Score))
    score_genes <- cbind.data.frame(target_gene, Score)
    Score_genes <- rbind.data.frame(Score_genes, score_genes)
  }
}

write.table(Score_genes,
  "results/Score_genes.csv",
  sep = ",",
  row.names = F,
  col.names = T
)
