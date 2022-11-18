

rm(list = ls())
library("igraph")
library("RColorBrewer")
library("Metrics")
library("L0Learn")
library("glmnet")
library("tidyr")
library("tidyverse")
library("ggthemes")
library("psych")
library("ggpubr")
library("GENIE3")

load("data/luad-rsem-fpkm-tcga-t_normlized.Rdata")
raw_tcga <- scale(raw_tcga)

sample_label <- read.table("NMF/cluster-rank=8/sample_cluster_4.csv",
  header = FALSE,
  sep = ",",
  check.names = FALSE
)

deg.data_clusters <- c()
evaluate_all <- c()
res_TRN_all <- c()
Contrast_all <- c()

for (j in 1:4) {
  survival_tfs <- c()
  res_TRN_cluster <- c()
  deg.data_cluster <- c()
  Contrast_cluster <- c()
  if (dir.exists(paste0("results/cluster", j)) == F) {
    dir.create(file.path(paste0("results/cluster", j)))
  }

  message(paste0("------------cluster", j, "start------------"))

  genes_list <- read.table(paste("data/genes_list", ".txt", sep = ""),
    header = TRUE,
    row.names = 1
  )

  cluster <- sample_label[which(sample_label$V2 == j), 1]
  Y_raw <- raw_tcga[row.names(genes_list), cluster]
  maxSNVSize <- 100
  evaluate_cluster <- c()
  for (i in 1:nrow(genes_list)) {
    target_gene <- row.names(genes_list)[i]
    if (file.exists(paste("data/", target_gene, "_TFs_list.txt", sep = "")) == T) {
      message(paste0("---------- the cluster", j, ",", "number", i, "gene:", target_gene, "computation is start ! ----------"))

      TFs_list <- read.table(paste("data/", target_gene, "_TFs_list.txt", sep = ""),
        header = T,
        row.names = 1
      )

      X <- t(raw_tcga[row.names(TFs_list), cluster]) %>% as.matrix()
      row.names(X) <- row.names(X)
      Y <- t(raw_tcga[target_gene, cluster]) %>% as.numeric()

      # Validation --------------------------------------------------------------
      set.seed(2022)
      train_idx <- sample(nrow(X), 0.7 * nrow(X))
      train_data_X <- X[train_idx, ]
      train_data_Y <- Y[train_idx]
      test_data_X <- X[-train_idx, ]
      test_data_Y <- Y[-train_idx]

      cvfit_L0_validation <- L0Learn.cvfit(train_data_X, train_data_Y, nFolds = 10)
      fit_L0_validation <- L0Learn.fit(train_data_X, train_data_Y,
        penalty = "L0",
        maxSuppSize = maxSNVSize
      )
      # Extract coefficient at middle lambda
      fit_L0_information <- as.data.frame(print(fit_L0_validation))
      fit_L0_information <- fit_L0_information[order(fit_L0_information$suppSize, decreasing = TRUE), ]
      lambda_L0 <- fit_L0_information$lambda[1]
      gamma_L0 <- fit_L0_information$gamma[1]

      test_data_y_cat <- predict(fit_L0_validation,
        newx = test_data_X,
        lambda = lambda_L0,
        gamma = gamma_L0
      )
      test_data_y_hat <- as.vector(test_data_y_cat)

      re_status <- cbind(test_data_Y, test_data_y_hat)
      head(re_status)
      re_status <- as.data.frame(re_status)
      colnames(re_status) <- c("raw", "pre")

      peak_corr <- corr.test(re_status$raw,
        re_status$pre,
        # method = "spearman",
        adjust = "none"
      )

      deg.data <- data.frame(
        Cluster = paste("Cluster", j),
        Gene = target_gene,
        corr = peak_corr$r,
        pval = peak_corr$p
      )
      deg.data_cluster <- rbind.data.frame(deg.data_cluster, deg.data)

      cor(re_status)
      ggplot(data = re_status, mapping = aes(x = raw, y = pre)) +
        geom_point() +
        geom_smooth(method = "lm", se = F, color = "red") +
        theme_bw() +
        stat_cor(data = re_status) + # , method = "spearman"
        geom_smooth(method = "lm", se = F, color = "#006699", size = 1) +
        labs(
          x = "True Expression", y = "L0Reg framework",
          title = paste0("Gene: ", target_gene)
        )

      ggsave(paste("results/cluster", j, "/", target_gene, "_corrgram.png", sep = ""),
        width = 3,
        height = 3
      )

      if (peak_corr$p < 1 && peak_corr$r > 0) {
        #L0
        fit_L0 <- L0Learn.fit(X, Y,
          penalty = "L0",
          maxSuppSize = maxSNVSize
        )
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

        L0_Rsquare <- 1 - mse(Y, y_hat) / var(Y)
        evaluate_L0 <- data.frame(
          Cluster = paste0("Cluster", j),
          Gene = target_gene,
          L0_Rsquare = 1 - mse(Y, y_hat) / var(Y),
          L0_RMSE = RMSE(Y, y_hat)
        )

        evaluate <- evaluate_L0
        evaluate_cluster <- rbind(evaluate_cluster, evaluate)

        res_data_f <- as.data.frame(res_data)
        res_data_f <- res_data_f[which(res_data_f$`Pr(>|t|)` <= 0.05), ]

        if (nrow(res_data_f) > 0) {
          if (rownames(res_data_f)[1] == "(Intercept)") {
            res_data_f <- res_data_f[-1, ]
          }

          for (k in 1:nrow(res_data_f)) {
            res_data_f$strength[k] <- abs(res_data_f$Estimate[k])
            if (res_data_f$Estimate[k] < 0) {
              res_data_f$reg[k] <- "2"
            } else {
              res_data_f$reg[k] <- "1"
            }
          }

          write.csv(res_data_f, sprintf(paste0("results/cluster", j, "/%s_TFs_Selected_L0.csv"), target_gene))

          res_TRN <- cbind(
            "TF" = row.names(res_data_f),
            "Reg" = res_data_f$reg,
            "Gene" = target_gene,
            "Cluster" = paste("Cluster", j, sep = ""),
            "Strength" = res_data_f$strength
          )
          res_TRN_cluster <- rbind(res_TRN_cluster, res_TRN)
          survival_tfs_gene <- rownames(res_data_f) %>% as.data.frame()
          survival_tfs <- rbind.data.frame(survival_tfs, survival_tfs_gene, target_gene)

          X <- raw_tcga[row.names(TFs_list), cluster]
          X <- rbind(X, Gene = Y_raw[target_gene, ]) %>% as.matrix()
          row.names(X) <- row.names(X)

          L0_matrix <- 0
          for (l in 1:nrow(res_data_f)) {
            gene <- rownames(res_data_f)[l]

            L0_gene <- X[gene, ] * res_data_f$Estimate[l]
            L0_matrix <- L0_matrix + L0_gene
          }

          # GENIE3
          if (T) {
            Gene <- target_gene
            pre_time <- proc.time()
            weightMat <- GENIE3(X, verbose = TRUE)
            time_GRN <- proc.time() - pre_time

            if (exists("weightMat")) {
              links <- getLinkList(weightMat)
              rm(weightMat)
              links <- links[which(links$targetGene == "Gene"), ]
              links$targetGene <- Gene

              GENIE3_matrix <- 0
              for (g in 1:nrow(links)) {
                gene <- links$regulatoryGene[g]
                if (gene != target_gene) {
                  GENIE3_gene <- X[gene, ] * links$weight[g]
                  GENIE3_matrix <- GENIE3_matrix + GENIE3_gene
                }
              }
            } else {
              GENIE3_matrix <- 0
            }
          }

          # TRN_Plot
          # Data_L0
          links_L0 <- res_TRN %>% as.data.frame()
          links_L0$weight <- as.numeric(links_L0$Strength)
          links_L0 <- links_L0[, c("TF", "Gene", "weight")]
          nodes_L0 <- data.frame(
            name = c(target_gene, links_L0$TF),
            carac = c(rep("Gene", 1), rep("TF", nrow(links_L0)))
          )

          network_L0 <- graph_from_data_frame(d = links_L0, vertices = nodes_L0, directed = F)
          coul <- brewer.pal(3, "Set1")
          # Create a vector of color
          color_L0 <- coul[as.numeric(as.factor(V(network_L0)$carac))]
          # Data_GENIE3
          links <- as.data.frame(links)
          links$regulatoryGene <- as.character(links$regulatoryGene)
          nodes <- data.frame(
            name = c(target_gene, links$regulatoryGene),
            carac = c(rep("Gene", 1), rep("TF", nrow(links)))
          )

          network <- graph_from_data_frame(d = links, vertices = nodes, directed = F)
          # Make a palette of 3 colors
          coul <- brewer.pal(3, "Set1")
          # Create a vector of color
          my_color <- coul[as.numeric(as.factor(V(network)$carac))]
          my_color

          png(paste0("results/cluster", j, "/", target_gene, "_GENIE3-L0.png"),
            width = 6000, height = 3000, res = 600
          )

          par(mfrow = c(1, 2))
          plot(network,
            vertex.color = my_color,
            vertex.label.color = c("black"),
            vertex.label.font = c(1),
            vertex.size = c(35),
            edge.width = E(network)$weight * 20
          )
          plot(network_L0,
            vertex.color = color_L0,
            vertex.label.color = c("black"),
            vertex.label.font = c(1),
            vertex.size = c(35),
            edge.width = E(network_L0)$weight * 20
          )
          legend("bottomleft",
            legend = levels(as.factor(V(network_L0)$carac)),
            col = coul, bty = "n", pch = 20, pt.cex = 3,
            cex = 1, text.col = coul, horiz = FALSE,
            inset = c(0, 0)
          )
          dev.off()

          # L0 VS GENIE3
          Contrast <- data.frame(
            Cluster = paste0("Cluster", j),
            Gene = target_gene,
            L0 = cor(X["Gene", ], L0_matrix),
            GENIE3 = cor(X["Gene", ], GENIE3_matrix)
          )
          Contrast_cluster <- rbind.data.frame(Contrast_cluster, Contrast)
        }
      }
      print(paste("---------- The cluster", j, ",", "number", i, "gene:", target_gene, "computation is completed ! ----------", sep = " "))
    } else {
      print(paste("---------- No", target_gene, "TFs file ! ----------", sep = " "))
    }
  }

  Contrast_all <- rbind.data.frame(Contrast_all, Contrast_cluster)

  deg.data_clusters <- rbind.data.frame(deg.data_clusters, deg.data_cluster)

  write.csv(evaluate_cluster, sprintf("results/cluster%s_evaluate.csv", j))

  evaluate_all <- rbind(evaluate_all, evaluate_cluster)
  res_TRN_cluster <- as.data.frame(res_TRN_cluster)
  write.csv(res_TRN_cluster, sprintf("results/res_TRN_cluster%s.csv", j))
  res_TRN_all <- rbind(res_TRN_all, res_TRN_cluster)
  survival_tfs <- survival_tfs[!duplicated(survival_tfs), ] %>% as.data.frame()
  write.csv(survival_tfs, paste0("../SurvivalAnalysis/genes_list_cluster", j, ".csv"))
  print(paste("---------- cluster", j, "done ! ----------", sep = " "))
}

# res_TRN_all <- res_TRN_all[-which(res_TRN_all$TFs=='(Intercept)'),]
write.csv(res_TRN_all, "results/res_TRN_all.csv", row.names = T)
write.csv(evaluate_all, "results/evaluate_all.csv", row.names = F)
write.csv(Contrast_all, "results/Contrast_all.csv", row.names = T)

Contrast_all_box_data <- Contrast_all[, -c(1, 2)]
names(Contrast_all_box_data) <- c("L0Reg framework", "GENIE3")
mean(Contrast_all_box_data$`L0Reg framework`)
mean(Contrast_all_box_data$GENIE3)

theme_set(theme_pubclean())
pkgs <- c("matrixStats", "pheatmap", "RColorBrewer", "tidyverse", "cowplot", "ggpubr", "bslib", "ggthemes")
lapply(pkgs, library, character.only = T)
Contrast_all_box <- Contrast_all_box_data %>%
  as.data.frame() %>%
  # rownames_to_column("sample") %>% #bug
  pivot_longer(
    cols = 1:2,
    names_to = "Method",
    values_to = "Correction"
  )
my_comparisons <- list(c("L0Reg framework", "GENIE3"))

ggplot(data = Contrast_all_box, aes(x = Method, y = Correction)) +
  geom_boxplot(aes(fill = Method)) +
  # ylab(paste0(target_gene," Expression"))+
  stat_compare_means() +
  theme_bw() +
  # geom_jitter(color="gray")+
  geom_jitter()

