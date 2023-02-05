

# Packages check, download and library --------------------------------------------------
#' Title
#'
#' @param packages
#'
#' @return
#' @export
#'
#' @examples
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

#' save.file Save R object
#'
#' @param ...
#' @param fileName
#' @param pathWay
#'
#' @return
#' @export
#'
#' @examples save.file(data, data2, fileName = "test.Rdata")
save.file <- function(..., fileName, pathWay = NULL) {
  if (is.null(pathWay)) {
    pathWay <- ""
  } else {
    if (!dir.exists(pathWay)) {
      dir.create(pathWay, recursive = TRUE)
    }
  }
  if (as.numeric(...length()) > 1) {
    if (grepl(fileName, pattern = ".Rdata$") | grepl(fileName, pattern = ".rdata$")) {
      save(..., file = paste0(pathWay, fileName))
    } else {
      newFileName <- sub("$", ".Rdata", fileName)
      save(..., file = paste0(pathWay, newFileName))
    }
  } else {
    save(..., file = paste0(pathWay, fileName))
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
survival.data <- function(cancerType = NULL, genes = NULL, pathWay = NULL) {
  if (is.null(pathWay)) {
    pathWay <- ""
  } else {
    if (!dir.exists(pathWay)) {
      dir.create(pathWay, recursive = TRUE)
    }
  }
  if (is.null(cancerType)) {
    message("----- Pleasure ensure the cancer type! -----")
  } else {
    message(paste0("----- Choose ", cancerType, "and prepare the data! -----"))
    package.check("cgdsr")
    package.check("DT")
    mycgds <- CGDS("http://www.cbioportal.org/")
    message(test(mycgds))
    all <- getCancerStudies(mycgds)
    getCaseLists(mycgds, cancerType)[, c(1, 2)]
    getGeneticProfiles(mycgds, cancerType)[, 1]
    getCaseLists(mycgds, cancerType)[, 1]
    mycaselist <- paste0(cancerType, "_rna_seq_v2_mrna")
    mygeneticprofile <- paste0(cancerType, "_rna_seq_v2_mrna")
    if (is.null(genes)) {
      message("----- Pleasure input a single gene or gene list! -----")
    } else {
      # Get expression data
      expr <- getProfileData(mycgds, genes, mygeneticprofile, mycaselist)
      # Get mutation data
      mut_df <- getProfileData(mycgds, caseList = "luad_tcga_sequenced", geneticProfile = "luad_tcga_mutations", genes = genes)
      mut_df <- apply(mut_df, 2, as.factor)
      mut_df[mut_df == "NaN"] <- ""
      mut_df[is.na(mut_df)] <- ""
      mut_df[mut_df != ""] <- "MUT"
      # Get copy number data
      cna <- getProfileData(mycgds,
        caseList = paste0(cancerType, "_sequenced"),
        geneticProfile = paste0(cancerType, "_gistic"),
        genes = genes
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
    save(expr, myClinicalData, cna, mut_df, file = paste0(pathWay, "survival_input.Rdata"))
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

#' LO_fit
#'
#' @param X The rows are samples and the columns are genes of the matrix
#' @param Y
#' @param penalty
#' @param nFolds
#' @param seed
#' @param maxSuppSize
#' @param nGamma
#' @param gammaMin
#' @param gammaMax
#'
#' @return
#' @export
#'
#' @examples
LO_fit <- function(X, Y,
                   penalty = penalty,
                   nFolds = 10,
                   seed = 1,
                   maxSuppSize = maxSuppSize,
                   nGamma = 5,
                   gammaMin = 0.0001,
                   gammaMax = 10) {
  tryCatch(
    {
      fit <- L0Learn::L0Learn.cvfit(X, Y,
        penalty = penalty,
        maxSuppSize = maxSuppSize,
        nFolds = 10,
        seed = 1,
        nGamma = 5,
        gammaMin = 0.0001,
        gammaMax = 10
      )
      fit_inf <- print(fit)
      optimalGammaIndex <- which(unlist(lapply(fit$cvMeans, min)) == min(unlist(lapply(fit$cvMeans, min))))
      gamma <- fit$fit$gamma[optimalGammaIndex]
      lambda_list <- fit_inf[which(fit_inf$gamma == gamma), ]
      if (is.null(maxSuppSize)) {
        lambda <- min(lambda_list$lambda)
      } else {
        if (maxSuppSize %in% lambda_list$maxSuppSize) {
          lambda <- lambda_list$maxSuppSize[which(lambda_list$maxSuppSize == maxSuppSize)]
        } else {
          lambda <- min(lambda_list$lambda)
        }
      }
      temp <- coef(fit,
        lambda = lambda,
        gamma = gamma
      )
    },
    error = function(e) {
      fit <- L0Learn::L0Learn.fit(X, Y,
        penalty = penalty,
        maxSuppSize = maxSuppSize,
        nGamma = 5,
        gammaMin = 0.0001,
        gammaMax = 10
      )
      fit_inf <- print(fit)
      fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
      lambda <- fit_inf$lambda[1]
      gamma <- fit_inf$gamma[1]
      temp <- coef(fit,
        lambda = lambda,
        gamma = gamma
      )
    }
  )
}

LO_fit <- function(X, Y,
                   penalty = penalty,
                   nFolds = 10,
                   seed = 1,
                   maxSuppSize = maxSuppSize,
                   nGamma = 5,
                   gammaMin = 0.0001,
                   gammaMax = 10) {

      fit <- L0Learn::L0Learn.fit(X, Y,
                                  penalty = penalty,
                                  maxSuppSize = maxSuppSize,
                                  nGamma = 5,
                                  gammaMin = 0.0001,
                                  gammaMax = 10
      )
      fit_inf <- print(fit)
      fit_inf <- fit_inf[order(fit_inf$suppSize, decreasing = TRUE), ]
      lambda <- fit_inf$lambda[1]
      gamma <- fit_inf$gamma[1]
      temp <- coef(fit,
                   lambda = lambda,
                   gamma = gamma
     
  )
}

#' Title
#'
#' @param matrix The rows are samples and the columns are genes of the matrix
#' @param penalty
#' @param regulators
#' @param targets
#' @param maxSuppSize
#'
#' @return
#' @export
#'
#' @examples
L0DWGRN <- function(matrix,
                    penalty = NULL,
                    regulators = NULL,
                    targets = NULL,
                    maxSuppSize = NULL,
                    cores = 1) {
  matrix <- as.data.frame(t(matrix))
  weightdf <- c()
  if (is.null(penalty)) {
    penalty <- "L0"
  }
  if (is.null(maxSuppSize)) {
    maxSuppSize <- dim(matrix)[2]
  }
  if (is.null(targets)) {
    targets <- colnames(matrix)
  }
  if (!is.null(regulators)) {
    matrix <- matrix[, regulators]
  } else {
    regulators <- colnames(matrix)
  }
  if (cores == 1) {
    for (i in 1:length(regulators)) {
      X <- as.matrix(matrix[, -which(colnames(matrix) == regulators[i])])
      Y <- matrix[, regulators[i]]
      temp <- LO_fit(X, Y,
        penalty = penalty,
        nFolds = 10,
        seed = 1,
        maxSuppSize = maxSuppSize,
        nGamma = 5,
        gammaMin = 0.0001,
        gammaMax = 10
      )
      temp <- as.vector(temp)
      wghts <- temp[-1]
      wghts <- abs(wghts)
      wghts <- wghts / sum(wghts)
      if (F) {
        wghts <- wghts / max(wghts)
        indices <- sort.list(wghts, decreasing = TRUE)
        zeros <- which(wghts <= 0.8)
        # wghts[1:length(wghts)] <- 1
        wghts[zeros] <- 0
      }
      if (length(wghts) != dim(X)[2]) {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = 0)
      } else {
        weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulators[i], weight = wghts)
      }
      weightdf <- rbind.data.frame(weightdf, weightd)
      if (i == length(regulators)) {
        weightdf <- weightdf[order(weightdf$weight, decreasing = TRUE), ]
      }
    }
  } else {
    cores <- min(parallel::detectCores(logical = F), cores)
    # cl <- parallel::makeCluster(cores)
    # doParallel::registerDoParallel(cl)
    doParallel::registerDoParallel(cores = cores)
    message(paste("\nUsing", foreach::getDoParWorkers(), "cores."))
    "%dopar%" <- foreach::"%dopar%"
    suppressPackageStartupMessages(
      weightdf <- doRNG::"%dorng%"(foreach::foreach(regulator = regulators, .combine = "rbind", .export = "LO_fit"), {
        X <- as.matrix(matrix[, -which(colnames(matrix) == regulator)])
        Y <- matrix[, regulator]
        temp <- LO_fit(X, Y,
          penalty = penalty,
          nFolds = 10,
          seed = 1,
          maxSuppSize = maxSuppSize,
          nGamma = 5,
          gammaMin = 0.0001,
          gammaMax = 10
        )
        temp <- as.vector(temp)
        wghts <- temp[-1]
        wghts <- abs(wghts)
        wghts <- wghts / sum(wghts)
        
        if (F) {
          # indices <- sort.list(wghts, decreasing = TRUE)
          # zeros <- which(wghts <= 0.8)
          # # wghts[1:length(wghts)] <- 1
          # wghts[zeros] <- 0
          
          # Now sort the wghts
          indices <- sort.list(wghts, decreasing = TRUE)
          # Check for zero entries
          zeros <- which(wghts == 0)
          # Now replace by ones that are in the top and are non-zero
          wghts[1:length(wghts)] <- 0
          rankThreshold <- 5
          wghts[indices[1:rankThreshold]] <- 1
          # Set the ones that were zero to zero anyway
          wghts[zeros] <- 0
        }
        
        if (length(wghts) != dim(X)[2]) {
          weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulator, weight = 0)
        } else {
          weightd <- data.frame(regulatoryGene = colnames(X), targetGene = regulator, weight = wghts)
        }
      })
    )
    # parallel::stopCluster(cl)
  }
  weightdf <- weightdf[order(weightdf$weight, decreasing = TRUE), ]
  return(weightdf)
}

# Defines class
setClass("evalClass", representation(
  tpr = "vector",
  fpr = "vector",
  prec = "vector",
  rec = "vector",
  p = "integer",
  n = "integer",
  tpk = "vector",
  fpk = "vector",
  predictors = "vector",
  targets = "vector",
  rank = "vector"
))

caclEval <- function(pred, goldTSV, totalPredictionsAccepted = 100000) {
  library("ROCR")
  library("caTools")
  # pred <- read.table(predictionTSV)
  gold <- convertSortedRankTSVToAdjMatrix(goldTSV)

  # Check arguments
  if ((!isWholeNumber(totalPredictionsAccepted)) || totalPredictionsAccepted < 1) {
    stop("totalPredictionsAccepted should be an integer <0.")
  }

  # Next reduce the predictionTSV to the max allowed predictions if necessary
  if (length(pred[, 1]) > totalPredictionsAccepted) {
    pred <- pred[1:totalPredictionsAccepted, ]
  }

  # Get some metrics from the gold standard
  # The amount of positive links
  p <- sum(gold)
  # The amount of negative links (these could be valid links and are not positive) [self regulating links are not allowed and are consired to be neither positive or negative]
  n <- (dim(gold)[1] * dim(gold)[2]) - p - sum(rownames(gold) %in% colnames(gold))
  # The total amount of valid links (positive +negative)
  t <- p + n

  # Now, remove all the links that are in the prediction file that are not valid predictions (meaning self-regulation or genes not eligble to be predictor or target)
  # 1) Remove all the edges that are in the prediction file but are not recorded in the gold file as an edge or non-edge
  # 2) For those that are in thet network a) Replace value with a one if edge is present in gold or zero otherwise
  # Remove all the non-valid predictors
  pred <- pred[which(as.vector(pred[, 1]) != as.vector(pred[, 2])), , drop = FALSE]

  # Declare some variables for clarity
  firstRow <- as.vector(pred[, 1])
  secondRow <- as.vector(pred[, 2])
  # Will indicate whether this prediction was present in gold standard
  thirdRow <- vector(mode = "integer", length(firstRow))
  # Will indicate the rank of this prediction = (-1 if it is not present and [rank] otherwise]
  rank <- vector(mode = "integer", length(firstRow))
  # Now, loop over all predictions and determine if they are present in gold matri
  correct <- 0
  incorrect <- 0

  for (i in 1:length(firstRow)) {
    if (gold[firstRow[i], secondRow[i]] == 1) {
      correct <- correct + 1
      thirdRow[i] <- 1
      rank[i] <- correct
    } else {
      incorrect <- incorrect - 1
      thirdRow[i] <- 0
      rank[i] <- incorrect
    }
  }
  # Check how many of the gold standard edges we predicted. Any other that still remain are discovered at a uniform rate by definition. (remaining_gold_edges/remaining_edges)
  # Gold links predicted so far

  # If some gold links have not been predicted, calculate the random discovery chance, else set to zero
  if (length(firstRow) < t) {
    odds <- (p - correct) / (t - length(firstRow))
  } else {
    odds <- 0
  }

  # Each guess you have 'odds' chance of getting one right and '1-odds' chance of getting it wrong , now construct a vector till the end
  random_positive <- vector("double", (t - length(firstRow)))
  random_negative <- vector("double", (t - length(firstRow)))
  random_positive[] <- odds
  random_negative[] <- 1 - odds
  # Calculate the amount of true positives and false positives at 'k' guesses depth
  positive <- c(thirdRow, random_positive)
  negative <- c(abs(thirdRow - 1), random_negative)
  tpk <- cumsum(positive)
  fpk <- cumsum(negative)
  # Depth k
  k <- 1:t
  # Calculate true positive rate, false positive rate, precision and recall
  tpr <- tpk / p
  fpr <- fpk / n
  rec <- tpr
  prec <- tpk / k
  predictors <- firstRow
  targets <- secondRow
  # Create and return the object
  return(new("evalClass", tpr = tpr, fpr = fpr, rec = rec, prec = prec, p = as.integer(p), n = as.integer(n), tpk = tpk, fpk = fpk, predictors = predictors, targets = targets, rank = rank))
}

convertSortedRankTSVToAdjMatrix <- function(inputFilename = NULL, input = NULL, outputFilename = NULL) {
  # Read the table from file
  if (!is.null(inputFilename)) {
    tbl <- read.table(inputFilename)
  } else {
    tbl <- input
  }
  # Sort by second column
  tbl <- tbl[sort.list(tbl[, 2]), ]
  # Get the number of unique elements in the predictors
  pred <- unique(tbl[, 1])
  # Get the number of unique elements in the targets
  targ <- unique(tbl[, 2])
  # Pre allocate return matrix
  m <- matrix(0.0, length(pred), length(targ))
  # Set the col/row-names
  colnames(m) <- targ
  rownames(m) <- pred
  # Get the duplicates
  dups <- duplicated(tbl[, 2])
  # Get the starIndices of another column
  startIndices <- which(FALSE == dups)
  for (i in 1:(length(startIndices) - 1)) {
    targetToAdd <- tbl[startIndices[i], 2]
    predToAdd <- tbl[startIndices[i]:(startIndices[i + 1] - 1), 1]
    valuesToAdd <- tbl[startIndices[i]:(startIndices[i + 1] - 1), 3]
    colIndex <- which(colnames(m) %in% targetToAdd)
    tmp <- c()
    for (i in predToAdd) {
      tmp <- c(tmp, which(rownames(m) == i))
    }
    rowIndexes <- tmp
    m[rowIndexes, colIndex] <- valuesToAdd
  }
  targetToAdd <- tbl[startIndices[length(startIndices)], 2]
  predToAdd <- tbl[startIndices[length(startIndices)]:length(tbl[, 1]), 1]
  valuesToAdd <- tbl[startIndices[length(startIndices)]:length(tbl[, 1]), 3]
  tmp <- c()
  for (i in predToAdd) {
    tmp <- c(tmp, which(rownames(m) == i))
  }
  rowIndexes <- tmp
  colIndex <- which(colnames(m) %in% targetToAdd)
  m[rowIndexes, colIndex] <- valuesToAdd

  if (!is.null(outputFilename)) {
    write.table(m, outputFilename)
  } else {
    return(m)
  }
}

areAllNullOrNonNull <- function(...) {
  e <- ((length((which(c(...) == FALSE)))))
  return(e != 1)
}

# 	isWholeNumber - check for integer number
isWholeNumber <- function(x, tol = .Machine$double.eps^0.5) {
  return(is.numeric(x) && abs(x - round(x)) < tol)
}

calcAUROC <- function(prediction) {
  return(trapz(prediction@fpr, prediction@tpr))
}

calcAUPR <- function(prediction) {
  return(trapz(prediction@rec, prediction@prec) / (1 - 1 / prediction@p))
}
