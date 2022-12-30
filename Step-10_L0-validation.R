

rm(list = ls())
library(broom)
library(olsrr)
library(car)
library(HH)
library(Metrics)
library(plotmo)
library(L0Learn)
library(caret)
library(glmnet)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(ggthemes)
library(psych)
library(ggpubr)

load("data/luad-rsem-count-tcga-t.Rdata")
load("data/luad-rsem-fpkm-tcga-t_normlized.Rdata")
raw_tcga <- scale(raw_tcga)

sample_label <- read.table("../results/2191/NMF/cluster-nk-ligands-k=8/sample_cluster_4.csv",
  header = FALSE,
  sep = ",",
  check.names = FALSE
)

if (T) {
  genes_list <- read.table(paste("data/genes_list", ".txt", sep = ""),
    header = TRUE,
    row.names = 1
  )
} else {
  genes_list <- read.table(paste("data/genes_list_cluster", j, ".txt", sep = ""),
    header = TRUE,
    row.names = 1
  )
}

maxSNVSize <- 100
correction_plots_list_all <- list()
evaluating_indicator_all <- c()
conditions <- c("randomize", "random", "TCGA-LUAD")

for (c in 1:length(conditions)) {
  evaluate_cluster <- c()
  deg.data_cluster <- c()
  deg.data_cluster_train <- c()
  condition <- conditions[c]
  if (dir.exists(paste0("results/validation_", condition)) == F) {
    dir.create(paste0("results/validation_", condition))
  }
  correction_plots_list <- list()
  for (i in 1:nrow(genes_list)) {
    target_gene <- row.names(genes_list)[i]
    if (file.exists(paste0("data/", target_gene, "_TFs_list.txt")) == T) {
      TFs_list <- read.table(paste0("data/", target_gene, "_TFs_list.txt"),
        header = T,
        row.names = 1
      )

      X <- raw_tcga[row.names(TFs_list), sample_label$V1]
      X <- t(X) %>% as.matrix()
      row.names(X) <- row.names(X)

      Y <- raw_tcga[target_gene, sample_label$V1]
      Y <- t(Y) %>% as.numeric()

      # Data -------------------------------------------------------------------
      if (condition == "TCGA-LUAD") {
        message("Use ", condition, " data")
      } else if (condition == "randomize") {
        message("Use ", condition, " data")
        X <- NMF::randomize(X)
        Y <- as.data.frame(Y)
        Y <- NMF::randomize(Y) %>% as.numeric()
      } else if (condition == "random") {
        message("Use ", condition, " data")
        samples_num <- nrow(X)
        for (i in seq_len(ncol(X))) {
          x <- rnorm(samples_num, 0, 1)
          X[, i] <- x
        }
        Y <- rnorm(samples_num, 0, 1)
      }

      # Validation -------------------------------------------------------------
      set.seed(2022)
      train_idx <- sample(nrow(X), 0.7 * nrow(X))
      train_data_X <- X[train_idx, ]
      train_data_Y <- Y[train_idx]
      test_data_X <- X[-train_idx, ]
      test_data_Y <- Y[-train_idx]

      cvfit_L0_validation <- L0Learn::L0Learn.cvfit(train_data_X, train_data_Y, nFolds = 10)
      plot(cvfit_L0_validation)
      fit_L0_validation <- L0Learn::L0Learn.fit(train_data_X, train_data_Y,
        penalty = "L0",
        maxSuppSize = maxSNVSize
      )

      # Extract coefficient at middle lambda
      fit_L0_information <- as.data.frame(print(fit_L0_validation))
      fit_L0_information <- fit_L0_information[order(fit_L0_information$suppSize, decreasing = TRUE), ]
      lambda_L0 <- fit_L0_information$lambda[1]
      gamma_L0 <- fit_L0_information$gamma[1]

      train_data_y_cat <- predict(fit_L0_validation,
        newx = train_data_X,
        lambda = lambda_L0,
        gamma = gamma_L0
      )
      train_data_y_hat <- as.vector(train_data_y_cat)
      L0_RMSE_train <- RMSE(train_data_Y, train_data_y_hat)
      re_status_train <- cbind(train_data_Y, train_data_y_hat) %>% as.data.frame()
      cor(re_status_train)
      head(re_status_train)
      colnames(re_status_train) <- c("raw", "pre")

      peak_corr <- psych::corr.test(re_status_train$raw,
        re_status_train$pre,
        # method = "spearman",
        adjust = "none"
      )

      deg.data <- data.frame(
        Gene = target_gene,
        corr = peak_corr$r,
        RMSD = L0_RMSE_train,
        pval = peak_corr$p
      )

      deg.data_cluster_train <- rbind.data.frame(deg.data_cluster_train, deg.data)

      test_data_y_cat <- predict(fit_L0_validation,
        newx = test_data_X,
        lambda = lambda_L0,
        gamma = gamma_L0
      )
      test_data_y_hat <- as.vector(test_data_y_cat)
      ###
      L0_RMSE <- RMSE(test_data_Y, test_data_y_hat)
      re_status <- cbind(test_data_Y, test_data_y_hat) %>% as.data.frame()
      head(re_status)
      colnames(re_status) <- c("raw", "pre")

      peak_corr <- psych::corr.test(re_status$raw,
        re_status$pre,
        # method = "spearman",
        adjust = "none"
      )

      deg.data <- data.frame(
        Gene = target_gene,
        corr = peak_corr$r,
        RMSD = L0_RMSE,
        pval = peak_corr$p
      )

      deg.data_cluster <- rbind.data.frame(deg.data_cluster, deg.data)

      cor(re_status)
      correction_plot <- ggplot(data = re_status, mapping = aes(x = raw, y = pre)) +
        geom_point() +
        geom_smooth(method = "lm", se = F, color = "red") +
        theme_bw() +
        stat_cor(data = re_status) + # , method = "spearman"
        geom_smooth(method = "lm", se = F, color = "#006699", size = 1) +
        labs(x = paste0("Expression of ", target_gene), y = "L0Reg framework")
      correction_plot
      correction_plots_list[[i]] <- correction_plot
      ggsave(paste0("results/validation_", condition, "/", target_gene, "_corrgram.png"),
        width = 3,
        height = 3
      )
    } else {
      print(paste("---------- No", target_gene, "TFs file ! ----------", sep = " "))
    }
  }

  ###
  deg.data_cluster_filter_train <- c()
  for (i in 1:nrow(deg.data_cluster_train)) {
    if (deg.data_cluster_train[i, ]$corr > 0.7 & deg.data_cluster_train[i, ]$RMSD < 2) {
      deg.data_cluster_filter_train <- rbind.data.frame(deg.data_cluster_filter_train, deg.data_cluster_train[i, ])
    }
  }

  deg.data_cluster_filter <- c()
  for (i in 1:nrow(deg.data_cluster)) {
    if (deg.data_cluster[i, ]$corr > 0.7 & deg.data_cluster[i, ]$RMSD < 2) {
      deg.data_cluster_filter <- rbind.data.frame(deg.data_cluster_filter, deg.data_cluster[i, ])
    }
  }

  write.csv(deg.data_cluster, paste0("results/", condition, "_Cor-RMSD.csv"))

  p1 <- ggplot() +
    geom_bar(
      data = deg.data_cluster,
      aes(x = Gene, y = RMSD),
      stat = "identity", position = "dodge", width = 0.9, size = 0.5
    ) +
    theme_bw() +
    labs(x = "", y = "RMSD", fill = "", color = "") +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 45, hjust = 1))

  evaluating_indicator <- data.frame(
    "Condition" = condition,
    "train_Corr" = mean(deg.data_cluster_train$corr),
    "train_RMSD" = mean(deg.data_cluster_train$RMSD),
    "test_Corr" = mean(deg.data_cluster$corr),
    "test_RMSD" = mean(deg.data_cluster$RMSD),
    "train_0.7_Corr" = mean(deg.data_cluster_filter_train$corr),
    "train_0.7_RMSD" = mean(deg.data_cluster_filter_train$RMSD),
    "test_0.7_Corr" = mean(deg.data_cluster_filter$corr),
    "test_0.7_RMSD" = mean(deg.data_cluster_filter$RMSD)
  )
  evaluating_indicator_all <- rbind.data.frame(evaluating_indicator_all, evaluating_indicator)
  correction_plots_list_all[[c]] <- correction_plots_list
}
