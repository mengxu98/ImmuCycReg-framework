

rm(list = ls())
library(tidyr)
library(tidyverse)
library(L0Learn)
library(glmnet)
library(GENIE3)

TRN_L0_all_list <- list()
links_list <- list()
for (d in c(1, 3, 4)) {
  dataset <- read.table(
    paste0(
      "DREAM5_network_inference_challenge/Network",
      d,
      "/input data/net",
      d,
      "_expression_data.tsv"
    ),
    header = T
  )
  if (T) {
    genie3_data <- as.matrix(dataset) %>% t()
    row.names(genie3_data) <- row.names(genie3_data)
    pre_time <- proc.time()
    weightMat <- GENIE3(genie3_data,
      nCores = 4,
      verbose = TRUE
    )
    time_GRN <- proc.time() - pre_time
    links <- getLinkList(weightMat)
    write.table(links,
      paste0(
        "results/GENIE3/DREAM5_NetworkInference_myteam_Network",
        d,
        ".txt"
      ),
      col.names = F,
      row.names = F,
      sep = "\t",
      quote = F
    )
    links_list[[d]] <- links
  }

  if (T) {
    maxSNVSize <- ncol(dataset)
    TRN_L0_all <- c()
    evaluate <- c()
    for (i in 1:ncol(dataset)) {
      target_gene <- colnames(dataset)[i]

      message(paste0("Computing gene ", i, "/", maxSNVSize, ": ", target_gene))

      Y <- dataset[, i]
      X <- dataset[, -i] %>% as.matrix()
      row.names(X) <- row.names(X)

      nfolds <- round(sqrt(dim(X)[1]))

      # Cross validate
      cross <- glmnet::cv.glmnet(X, Y, nfolds = nfolds, family = "gaussian", standardize = FALSE)
      # Retrain optimal
      bestModel <- glmnet::glmnet(X, Y, family = "gaussian", lambda = cross$lambda.min, standardize = FALSE)
      wghts <- abs(as.vector(bestModel$beta))

      # L0
      cvfit_L0 <- L0Learn.cvfit(X, Y, nFolds = 10)
      fit_L0 <- L0Learn.fit(X, Y,
        penalty = "L0",
        maxSuppSize = maxSNVSize
      )
      # Extract coefficient at middle lambda
      fit_L0_information <- as.data.frame(print(fit_L0))
      fit_L0_information <- fit_L0_information[order(fit_L0_information$suppSize, decreasing = TRUE), ]
      lambda_L0 <- fit_L0_information$lambda[1]
      gamma_L0 <- fit_L0_information$gamma[1]

      y_cat <- predict(fit_L0,
        newx = X,
        lambda = lambda_L0,
        gamma = gamma_L0
      )

      y_hat <- as.vector(y_cat)

      temp <- coef(fit_L0,
        lambda = lambda_L0,
        gamma = gamma_L0
      )

      temp <- as.vector(temp)
      temp <- temp[-1]
      temp <- which(temp != 0)
      temp <- colnames(X)[temp]
      temp <- na.omit(temp)

      if (length(temp) == 1) {
        X_Y <- cbind(X[, temp], Y)
        colnames(X_Y)[1] <- temp
        X_Y_frame <- as.data.frame(X_Y)
      } else {
        X_Y <- cbind(X[, temp], Y)
        X_Y_frame <- as.data.frame(X_Y)
      }

      lmfit <- lm(Y ~ ., data = X_Y_frame)
      fit_temp <- summary(lmfit)
      res_data <- as.matrix(fit_temp$coefficients)
      res_data_f <- as.data.frame(res_data)
      res_data_f <- res_data_f[which(res_data_f$`Pr(>|t|)` <= 0.05), ]

      if (nrow(res_data_f) > 0) {
        if (rownames(res_data_f)[1] == "(Intercept)") {
          res_data_f <- res_data_f[-1, ]
        }
        res_data_f$strength <- abs(res_data_f$Estimate)
        for (k in 1:nrow(res_data_f)) {
          if (res_data_f$Estimate[k] < 0) {
            res_data_f$reg[k] <- "2"
          } else {
            res_data_f$reg[k] <- "1"
          }
        }
      }

      TRN_L0 <- cbind(
        "TF" = row.names(res_data_f),
        "Gene" = target_gene,
        "Strength" = res_data_f$strength
      ) %>% as.data.frame()

      if (i > 2) {
        if (ncol(TRN_L0_all) == ncol(TRN_L0)) {
          TRN_L0_all <- rbind.data.frame(TRN_L0_all, TRN_L0)
        }
      } else {
        TRN_L0_all <- rbind.data.frame(TRN_L0_all, TRN_L0)
      }
    }
  }
  TRN_L0_all$Strength <- as.numeric(TRN_L0_all$Strength)
  TRN_L0_all <- TRN_L0_all[order(-TRN_L0_all$Strength), ]
  TRN_L0_all_list[[d]] <- TRN_L0_all
  write.table(TRN_L0_all,
    paste0(
      "results/L0/DREAM5_NetworkInference_myteam_Network",
      d,
      ".txt"
    ),
    col.names = F,
    row.names = F,
    sep = "\t",
    quote = F
  )
}
