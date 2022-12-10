

# Packages check, download and library --------------------------------------------------
package.check <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      if (!requireNamespace("dplyr", quietly = TRUE)) {
        install.packages("dplyr")
      }
      library("dplyr")
      if (!requireNamespace("rvest", quietly = TRUE)) {
        install.packages("rvest")
      }
      library("rvest")
      message("[", Sys.time(), "] -----: No package: ", package, " in R environment!")
      CRANpackages <- available.packages() %>%
        as.data.frame() %>%
        select(Package) %>%
        mutate(source = "CRAN")
      url <- "https://www.bioconductor.org/packages/release/bioc/"
      biocPackages <- url %>%
        read_html() %>%
        html_table() %>%
        .[[1]] %>%
        select(Package) %>%
        mutate(source = "BioConductor")
      if (package %in% CRANpackages$Package) {
        message("[", Sys.time(), "] -----: Now install package: ", package, " from CRAN!")
        install.packages(package)
        library(package, character.only = TRUE)
      } else if (package %in% biocPackages$Package) {
        message("[", Sys.time(), "] -----: Now install package: ", package, " from BioConductor!")
        BiocManager::install(package)
        library(package, character.only = TRUE)
      } else { # Bug
        if (!requireNamespace("githubinstall", quietly = TRUE)) {
          install.packages("githubinstall")
        }
        library("githubinstall")
        # githubinstall(package)
        gh_suggest(package)
      }
    } else {
      library(package, character.only = TRUE)
    }
  }
}

# Save R object
# Example: save.file(data, data2, fileName = "test.Rdata")
save.file <- function(..., fileName, pathWay = NULL) {
  if (is.null(pathWay)) {
    pathWay <- ""
  }
  if (as.numeric(...length()) > 1) {
    if (grepl(fileName, pattern = ".Rdata$") | grepl(fileName, pattern = ".rdata$")) {
      save(..., file = fileName)
    } else {
      newFileName <- sub(".Rdata$", ".Rdata", fileName)
      save(..., file = newFileName)
    }
  } else {
    save(..., file = fileName)
  }
}

# Check whether the file exists! ---------------
check.file.exists <- function(data, format_save = ".rds") {
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

# Transformation of COUNT, TPM, FPKM ---------------
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

# To obtain survival data of TCAGA samples ---------------
# It is required to specify the single gene or genes list to obtain survival data
survival.data <- function(cancerType = NULL, immuneGene = NULL) {
  if (is.null(cancerType)) {
    message("----- Pleasure ensure the cancer type! -----")
  } else {
    message(paste0("----- Choose ", cancerType, "and prepare the data! -----"))
    package.check("cgdsr")
    package.check("DT")
    mycgds <- CGDS("http://www.cbioportal.org/")
    message(test(mycgds))
    all <- getCancerStudies(mycgds)
    DT::datatable(all)
    getCaseLists(mycgds, cancerType)[, c(1, 2)]
    getGeneticProfiles(mycgds, cancerType)[, 1]
    getCaseLists(mycgds, cancerType)[, 1]
    getGeneticProfiles(mycgds, cancerType)[, 1]
    mycaselist <- paste0(cancerType, "_tcga_rna_seq_v2_mrna")
    mygeneticprofile <- paste0(cancerType, "_rna_seq_v2_mrna")
    if (is.null(immuneGene)) {
      message("Pleasure select genes!")
    } else {
      # Get expression data
      expr <- getProfileData(mycgds, immuneGene, mygeneticprofile, mycaselist)
      # Get mutation data
      mut_df <- getProfileData(mycgds, caseList = "luad_tcga_sequenced", geneticProfile = "luad_tcga_mutations", genes = immuneGene)
      mut_df <- apply(mut_df, 2, as.factor)
      mut_df[mut_df == "NaN"] <- ""
      mut_df[is.na(mut_df)] <- ""
      mut_df[mut_df != ""] <- "MUT"
      # Get copy number data
      cna <- getProfileData(mycgds,
        caseList = paste0(cancerType, "_sequenced"),
        geneticProfile = paste0(cancerType, "_gistic"),
        genes = immuneGene
      )
    }
    rn <- rownames(cna)
    cna <- apply(cna, 2, function(x) {
      as.character(factor(x,
        levels = c(-2:2),
        labels = c("HOMDEL", "HETLOSS", "DIPLOID", "GAIN", "AMP")
      ))
    })
    cna[is.na(cna)] <- ""
    cna[cna == "DIPLOID"] <- ""
    rownames(cna) <- rn
    myClinicalData <- getClinicalData(mycgds, mycaselist)
    save(expr, myClinicalData, cna, mut_df, file = paste0(path_save, "survival_input.Rdata"))
  }
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
  for (i in 1:length(results_summary_input)) {
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
  for (i in high_exp_sample_input) {
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
