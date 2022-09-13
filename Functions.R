

# Check whether the file exists!
check_file.exists <- function(data, format_save = ".rds") {
  if (exists("format_save")) {
    format_save <- format_save
  }
  data_name <- deparse(substitute(data))
  if (file.exists(paste0(path_save, data_name, format_save))) {
    message("----- The file: ", data_name, format_save, " has exist! -----")
  } else {
    if (format_save %in% c(".txt", ".csv")) {
      if (format_save == ".csv") {
        write.csv(data, paste0(path_save, data_name, format_save), quote = FALSE)
      }
      if (format_save == ".txt") {
        write.table(data, paste0(path_save, data_name, format_save), quote = FALSE, sep = "\t")
      }
    } else {
      if (dir.exists(path_save)) {
        saveRDS(data, paste0(path_save, data_name, format_save))
      } else {
        dir.create(path_save)
        saveRDS(data, paste0(path_save, data_name, format_save))
      }
    }
  }
}


# Transformation of COUNT, TPM, FPKM
countToTpm <- function(counts, effLen) {
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen) {
  n <- sum(counts)
  exp(log(counts) + log(1e9) - log(effLen) - log(n))
}

fpkmToTpm <- function(fpkm) {
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

countToEffCounts <- function(counts, len, effLen) {
  counts * (len / effLen)
}



Peak_is_open <- function(candidate_peak_id_input, target_sample_input) {
  up_thres <- quantile(peak_luad[, target_sample_input])[3]
  # print(up_thres)
  # if(TRUE %in% c(peak_luad[candidate_peak_id_input,target_sample_input]>up_thres))
  if (peak_luad[candidate_peak_id_input, target_sample_input] > up_thres) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

formatPositiveResult <- function(results_summary_input) {
  TFs <- c()
  type_inter <- c()
  target_gene_inter <- c()
  sample_level <- c()
  if (length(results_summary_input)) {
    i <- 1
    while (i <= length(results_summary_input)) {
      j <- 1
      while (j <= length(results_summary_input[[i]])) {
        TFs <- c(TFs, results_summary_input[[i]][j])
        type_inter <- c(type_inter, "1") # type 1
        sample_level <- c(sample_level, names(results_summary_input)[i])
        target_gene_inter <- c(target_gene_inter, target_gene)
        j <- j + 1
      }
      i <- i + 1
    }
    format_result <- data.frame(cbind(TFs, type_inter, target_gene_inter, sample_level))
    colnames(format_result) <- c("TF", "type", "target_gene", "sample")
    return(format_result)
  } else {
    print("NULL list")
  }
}

formatNegativeResult <- function(results_summary_input) {
  TFs <- c()
  type_inter <- c()
  target_gene_inter <- c()
  sample_level <- c()
  if (length(results_summary_input)) {
    i <- 1
    while (i <= length(results_summary_input)) {
      j <- 1
      while (j <= length(results_summary_input[[i]])) {
        TFs <- c(TFs, results_summary_input[[i]][j])
        type_inter <- c(type_inter, "2") # type 1
        sample_level <- c(sample_level, names(results_summary_input)[i])
        target_gene_inter <- c(target_gene_inter, target_gene)
        j <- j + 1
      }
      i <- i + 1
    }
    format_result <- data.frame(cbind(TFs, type_inter, target_gene_inter, sample_level))
    colnames(format_result) <- c("TF", "type", "target_gene", "sample")
    return(format_result)
  } else {
    print("NULL list")
  }
}

SampleHit <- function(results_summary_input, hit_tf) {
  hit_count <- c()
  sample_level <- c()
  for (i in 1:length(results_summary_input)) { # each sample
    j <- 1
    count_temp <- 0
    while (j <= length(results_summary_input[[i]])) {
      if (results_summary_input[[i]][j] %in% hit_tf & names(results_summary_input)[i] %in% colnames(raw_tcga_cnv)) {
        if (raw_tcga_cnv[results_summary_input[[i]][j], names(results_summary_input)[i]] > 0) {
          count_temp <- count_temp + 1
        }
      }
      j <- j + 1
    }
    sample_level <- c(sample_level, names(results_summary_input)[i])
    hit_count <- c(hit_count, count_temp)
  }
  format_result <- data.frame(cbind(sample_level, hit_count))
  colnames(format_result) <- c("sample", paste("hitsby", target_gene, sep = ""))
  return(format_result)
}

SampleHitOnlyCNV <- function(high_exp_sample_input, hit_tf) {
  sample_level <- c()
  # tem_ind=0
  for (i in high_exp_sample_input) { # each sample
    if (hit_tf %in% row.names(raw_tcga_cnv) & i %in% colnames(raw_tcga_cnv)) {
      # tem_ind=tem_ind+1
      # print(tem_ind)
      if (raw_tcga_cnv[hit_tf, i] > 0) {
        sample_level <- c(sample_level, i)
      }
    }
  }
  return(sample_level)
}


FrameRegulatoryTable <- function(res_frame_input) {
  ASPAR_hit <- c()
  encode_hit <- c()
  CHEA_hit <- c()
  MotifMap_hit <- c()
  TRANSFAC_Curated_hit <- c()
  TRANSFAC_Predicted_hit <- c()
  for (i in res_frame_input$Var1) {
    ASPAR_res <- ASPAR_Predicted[which(ASPAR_Predicted$source == target_gene & ASPAR_Predicted$target == i), ]
    if (nrow(ASPAR_res)) {
      ASPAR_hit <- c(ASPAR_hit, 1)
    } else {
      ASPAR_hit <- c(ASPAR_hit, 0)
    }
    encode_tf_res <- encode_tf[which(encode_tf$source == target_gene & encode_tf$target == i), ]
    if (nrow(encode_tf_res)) {
      encode_hit <- c(encode_hit, 1)
    } else {
      encode_hit <- c(encode_hit, 0)
    }
    CHEA_tf_res <- CHEA_tf[which(CHEA_tf$source == target_gene & CHEA_tf$target == i), ]
    if (nrow(CHEA_tf_res)) {
      CHEA_hit <- c(CHEA_hit, 1)
    } else {
      CHEA_hit <- c(CHEA_hit, 0)
    }
    MotifMap_Predicted_res <- MotifMap_Predicted[which(MotifMap_Predicted$source == target_gene & MotifMap_Predicted$target == i), ]
    if (nrow(MotifMap_Predicted_res)) {
      MotifMap_hit <- c(MotifMap_hit, 1)
    } else {
      MotifMap_hit <- c(MotifMap_hit, 0)
    }
    TRANSFAC_Curated_res <- TRANSFAC_Curated[which(TRANSFAC_Curated$source == target_gene & TRANSFAC_Curated$target == i), ]
    if (nrow(TRANSFAC_Curated_res)) {
      TRANSFAC_Curated_hit <- c(TRANSFAC_Curated_hit, 1)
    } else {
      TRANSFAC_Curated_hit <- c(TRANSFAC_Curated_hit, 0)
    }
    TRANSFAC_Predicted_res <- TRANSFAC_Predicted[which(TRANSFAC_Predicted$source == target_gene & TRANSFAC_Predicted$target == i), ]
    if (nrow(TRANSFAC_Predicted_res)) {
      TRANSFAC_Predicted_hit <- c(TRANSFAC_Predicted_hit, 1)
    } else {
      TRANSFAC_Predicted_hit <- c(TRANSFAC_Predicted_hit, 0)
    }
  }
  results_regulatory_table$ASPAR <- ASPAR_hit
  results_regulatory_table$ENCODE <- encode_hit
  results_regulatory_table$CHEA <- CHEA_hit
  results_regulatory_table$MotifMap <- MotifMap_hit
  results_regulatory_table$TRANSFAC_Curated <- TRANSFAC_Curated_hit
  results_regulatory_table$TRANSFAC_Predicted <- TRANSFAC_Predicted_hit
  return(results_regulatory_table)
}


FramePositive <- function(res_frame_input) {
  ASPAR_hit <- c()
  encode_hit <- c()
  CHEA_hit <- c()
  MotifMap_hit <- c()
  TRANSFAC_Curated_hit <- c()
  TRANSFAC_Predicted_hit <- c()
  for (i in res_frame_input$Var1) {
    ASPAR_res <- ASPAR_Predicted[which(ASPAR_Predicted$source == target_gene & ASPAR_Predicted$target == i), ]
    if (nrow(ASPAR_res)) {
      ASPAR_hit <- c(ASPAR_hit, 1)
    } else {
      ASPAR_hit <- c(ASPAR_hit, 0)
    }
    encode_tf_res <- encode_tf[which(encode_tf$source == target_gene & encode_tf$target == i), ]
    if (nrow(encode_tf_res)) {
      encode_hit <- c(encode_hit, 1)
    } else {
      encode_hit <- c(encode_hit, 0)
    }
    CHEA_tf_res <- CHEA_tf[which(CHEA_tf$source == target_gene & CHEA_tf$target == i), ]
    if (nrow(CHEA_tf_res)) {
      CHEA_hit <- c(CHEA_hit, 1)
    } else {
      CHEA_hit <- c(CHEA_hit, 0)
    }
    MotifMap_Predicted_res <- MotifMap_Predicted[which(MotifMap_Predicted$source == target_gene & MotifMap_Predicted$target == i), ]
    if (nrow(MotifMap_Predicted_res)) {
      MotifMap_hit <- c(MotifMap_hit, 1)
    } else {
      MotifMap_hit <- c(MotifMap_hit, 0)
    }
    TRANSFAC_Curated_res <- TRANSFAC_Curated[which(TRANSFAC_Curated$source == target_gene & TRANSFAC_Curated$target == i), ]
    if (nrow(TRANSFAC_Curated_res)) {
      TRANSFAC_Curated_hit <- c(TRANSFAC_Curated_hit, 1)
    } else {
      TRANSFAC_Curated_hit <- c(TRANSFAC_Curated_hit, 0)
    }
    TRANSFAC_Predicted_res <- TRANSFAC_Predicted[which(TRANSFAC_Predicted$source == target_gene & TRANSFAC_Predicted$target == i), ]
    if (nrow(TRANSFAC_Predicted_res)) {
      TRANSFAC_Predicted_hit <- c(TRANSFAC_Predicted_hit, 1)
    } else {
      TRANSFAC_Predicted_hit <- c(TRANSFAC_Predicted_hit, 0)
    }
  }
  res_frame_positive$ASPAR <- ASPAR_hit
  res_frame_positive$ENCODE <- encode_hit
  res_frame_positive$CHEA <- CHEA_hit
  res_frame_positive$MotifMap <- MotifMap_hit
  res_frame_positive$TRANSFAC_Curated <- TRANSFAC_Curated_hit
  res_frame_positive$TRANSFAC_Predicted <- TRANSFAC_Predicted_hit
  return(res_frame_positive)
}

FrameNegative <- function(res_frame_input) {
  ASPAR_hit <- c()
  encode_hit <- c()
  CHEA_hit <- c()
  MotifMap_hit <- c()
  TRANSFAC_Curated_hit <- c()
  TRANSFAC_Predicted_hit <- c()
  for (i in res_frame_input$Var1) {
    ASPAR_res <- ASPAR_Predicted[which(ASPAR_Predicted$source == target_gene & ASPAR_Predicted$target == i), ]
    if (nrow(ASPAR_res)) {
      ASPAR_hit <- c(ASPAR_hit, 1)
    } else {
      ASPAR_hit <- c(ASPAR_hit, 0)
    }
    encode_tf_res <- encode_tf[which(encode_tf$source == target_gene & encode_tf$target == i), ]
    if (nrow(encode_tf_res)) {
      encode_hit <- c(encode_hit, 1)
    } else {
      encode_hit <- c(encode_hit, 0)
    }
    CHEA_tf_res <- CHEA_tf[which(CHEA_tf$source == target_gene & CHEA_tf$target == i), ]
    if (nrow(CHEA_tf_res)) {
      CHEA_hit <- c(CHEA_hit, 1)
    } else {
      CHEA_hit <- c(CHEA_hit, 0)
    }
    MotifMap_Predicted_res <- MotifMap_Predicted[which(MotifMap_Predicted$source == target_gene & MotifMap_Predicted$target == i), ]
    if (nrow(MotifMap_Predicted_res)) {
      MotifMap_hit <- c(MotifMap_hit, 1)
    } else {
      MotifMap_hit <- c(MotifMap_hit, 0)
    }
    TRANSFAC_Curated_res <- TRANSFAC_Curated[which(TRANSFAC_Curated$source == target_gene & TRANSFAC_Curated$target == i), ]
    if (nrow(TRANSFAC_Curated_res)) {
      TRANSFAC_Curated_hit <- c(TRANSFAC_Curated_hit, 1)
    } else {
      TRANSFAC_Curated_hit <- c(TRANSFAC_Curated_hit, 0)
    }
    TRANSFAC_Predicted_res <- TRANSFAC_Predicted[which(TRANSFAC_Predicted$source == target_gene & TRANSFAC_Predicted$target == i), ]
    if (nrow(TRANSFAC_Predicted_res)) {
      TRANSFAC_Predicted_hit <- c(TRANSFAC_Predicted_hit, 1)
    } else {
      TRANSFAC_Predicted_hit <- c(TRANSFAC_Predicted_hit, 0)
    }
  }
  res_frame_negative$ASPAR <- ASPAR_hit
  res_frame_negative$ENCODE <- encode_hit
  res_frame_negative$CHEA <- CHEA_hit
  res_frame_negative$MotifMap <- MotifMap_hit
  res_frame_negative$TRANSFAC_Curated <- TRANSFAC_Curated_hit
  res_frame_negative$TRANSFAC_Predicted <- TRANSFAC_Predicted_hit
  return(res_frame_negative)
}

regulation_score_calculation <- function(results_summary_input,
                                         target_gene,
                                         target_sample,
                                         common_samples,
                                         raw_tcga,
                                         immune_cell,
                                         CIBERSORT_Results_tcga) {
  immune_score_all <- c()
  immune_score_res <- 0
  if (length(results_summary_input) > 1) {
    q <- 1
    while (q <= length(results_summary_input)) {
      p <- 1
      while (p <= length(immune_cell)) {
        cor_res <- corr.test(t(raw_tcga[target_gene, common_samples]),
          CIBERSORT_Results_tcga[common_samples, immune_cell[p]],
          method = "spearman",
          adjust = "none"
        )
        # cor_res$r # r值表示在样本中变量间的相关系数，表示相关性的大小
        # cor_res$p # p值是检验值，是检验两变量在样本来自的总体中是否存在和样本一样的相关性
        if (i %in% c("1", "2", "4", "5", "6", "7", "8", "10", "11", "12", "13", "14", "15", "17")) {
          if (abs(cor_res$p) < 0.05) { # abs(cor_res$r) > 0.6 &
            immune_score_res <- immune_score_res + CIBERSORT_Results_tcga[target_sample, p]
          }
        } else {
          if (abs(cor_res$p) < 0.05) {
            immune_score_res <- immune_score_res - CIBERSORT_Results_tcga[target_sample, p]
          }
        }
        p <- p + 1
      }
      ImmuneScore <- length(results_summary_input[[q]]) / nrow(gene_list) / 20 + immune_score_res
      SampleSelected <- names(results_summary_input)[q]
      immune_score_sample <- cbind.data.frame(SampleSelected, ImmuneScore)
      immune_score_all <- rbind.data.frame(immune_score_all, immune_score_sample)
      q <- q + 1
    }
  }
  immune_score_all <- as.data.frame(immune_score_all)
  return(immune_score_all)
}

immune_score_calculation <- function(results_summary_input) {
  immune_score_all <- c()
  if (length(results_summary_input) > 1) {
    q <- 1
    while (q <= length(results_summary_input)) {
      ImmuneScore <- length(results_summary_input[[q]]) / nrow(gene_list) / 10
      SampleSelected <- names(results_summary_input)[q]
      immune_score_sample <- cbind.data.frame(SampleSelected, ImmuneScore)
      immune_score_all <- rbind.data.frame(immune_score_all, immune_score_sample)
      q <- q + 1
    }
  }
  immune_score_all <- as.data.frame(immune_score_all)
  return(immune_score_all)
}
