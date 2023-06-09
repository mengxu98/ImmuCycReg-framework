rm(list = ls())
source("functions/Functions.R")

# Load data
load("../../Results/TCGA-LUAD.Rdata")
tcga_luad <- scale(tcga_luad)

sampleLabel <- read.table("../data/sample_cluster_4.csv",
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
  if (dir.exists(paste0("../../Results/L0Reg-framework/cluster", j)) == FALSE) {
    dir.create(file.path(paste0("../../Results/L0Reg-framework/cluster", j)))
  }

  message(paste("Running for cluster", j, "......"))

  genes_list <- read.table(paste("../data/Genes_17", ".txt", sep = ""),
    header = TRUE,
    row.names = 1
  )

  cluster <- sampleLabel[which(sampleLabel$V2 == j), 1]
  Y_raw <- tcga_luad[row.names(genes_list), cluster]
  maxSNVSize <- 100
  evaluate_cluster <- c()
  for (i in 1:nrow(genes_list)) {
    targetGene <- row.names(genes_list)[i]
    if (file.exists(paste0("../data/TFs/", targetGene, "_TFs_list.txt")) == T) {
      message(paste0("Running for cluster: ", j, ", ", "NO.", i, " gene: ", targetGene, "......"))

      TFs_list <- read.table(paste0("../data/TFs/", targetGene, "_TFs_list.txt"),
        header = T,
        row.names = 1
      )

      X <- t(tcga_luad[row.names(TFs_list), cluster]) %>% as.matrix()
      row.names(X) <- row.names(X)
      Y <- t(tcga_luad[targetGene, cluster]) %>% as.numeric()

      # Validation -------------------------------------------------------------
      set.seed(2022)

      # Split train dataset and test dataset
      train_idx <- sample(nrow(X), 0.7 * nrow(X))
      train_data_X <- X[train_idx, ]
      train_data_Y <- Y[train_idx]
      test_data_X <- X[-train_idx, ]
      test_data_Y <- Y[-train_idx]

      fit_L0_validation <- L0Learn::L0Learn.fit(
        train_data_X,
        train_data_Y,
        penalty = "L0",
        maxSuppSize = maxSNVSize
      )

      # Extract coefficient at middle lambda
      fitInformation <- print(fit_L0_validation) %>% as.data.frame()
      fitInformation <- fitInformation[order(fitInformation$suppSize, decreasing = TRUE), ]
      lambda_L0 <- fitInformation$lambda[1]
      gamma_L0 <- fitInformation$gamma[1]

      test_data_y <- predict(fit_L0_validation,
        newx = test_data_X,
        lambda = lambda_L0,
        gamma = gamma_L0
      ) %>% as.vector()

      re_status <- cbind(test_data_Y, test_data_y) %>% as.data.frame()

      colnames(re_status) <- c("raw", "pre")

      peak_corr <- psych::corr.test(re_status$raw,
        re_status$pre,
        # method = "spearman",
        adjust = "none"
      )

      deg.data <- data.frame(
        Cluster = paste("Cluster", j),
        Gene = targetGene,
        corr = peak_corr$r,
        pval = peak_corr$p
      )
      deg.data_cluster <- rbind.data.frame(deg.data_cluster, deg.data)

      ggplot(data = re_status, mapping = aes(x = raw, y = pre)) +
        geom_point() +
        geom_smooth(method = "lm", se = F, color = "red") +
        theme_bw() +
        stat_cor(data = re_status) + # , method = "spearman"
        geom_smooth(method = "lm", se = F, color = "#006699", size = 1) +
        labs(
          x = "True Expression", y = "L0Reg framework",
          title = paste0("Gene: ", targetGene)
        )

      ggsave(paste("../../Results/L0Reg-framework/cluster", j, "/", targetGene, "_corrgram.png", sep = ""),
        width = 3,
        height = 3
      )

      if (peak_corr$p < 0.05 && peak_corr$r > 0) {
        fit_L0 <- L0Learn::L0Learn.fit(X, Y,
          penalty = "L0",
          maxSuppSize = maxSNVSize
        )
        fitInformation <- as.data.frame(print(fit_L0))
        fitInformation <- fitInformation[order(fitInformation$suppSize, decreasing = TRUE), ]
        lambda_L0 <- fitInformation$lambda[1]
        gamma_L0 <- fitInformation$gamma[1]

        y_hat <- predict(fit_L0,
          newx = X,
          lambda = lambda_L0,
          gamma = gamma_L0
        ) %>% as.vector()
        temp <- coef(fit_L0,
          lambda = lambda_L0,
          gamma = gamma_L0
        ) %>% as.vector()
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
          Gene = targetGene,
          L0_Rsquare = 1 - mse(Y, y_hat) / var(Y),
          L0_RMSE = Metrics::rmse(Y, y_hat)
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

          write.csv(res_data_f, sprintf(paste0("../../Results/L0Reg-framework/cluster", j, "/%s_TFs_Selected_L0.csv"), targetGene))

          res_TRN <- cbind(
            "TF" = row.names(res_data_f),
            "Reg" = res_data_f$reg,
            "Gene" = targetGene,
            "Cluster" = paste("Cluster", j, sep = ""),
            "Strength" = res_data_f$strength
          )
          res_TRN_cluster <- rbind(res_TRN_cluster, res_TRN)
          survival_tfs_gene <- rownames(res_data_f) %>% as.data.frame()
          survival_tfs <- rbind.data.frame(survival_tfs, survival_tfs_gene, targetGene)

          X <- tcga_luad[row.names(TFs_list), cluster]
          X <- rbind(X, Gene = Y_raw[targetGene, ]) %>% as.matrix()
          row.names(X) <- row.names(X)

          L0_matrix <- 0
          for (l in 1:nrow(res_data_f)) {
            gene <- rownames(res_data_f)[l]

            L0_gene <- X[gene, ] * res_data_f$Estimate[l]
            L0_matrix <- L0_matrix + L0_gene
          }

          # GENIE3
          if (T) {
            Gene <- targetGene
            pre_time <- proc.time()
            weightMat <- GENIE3::GENIE3(X, verbose = TRUE)
            time_GRN <- proc.time() - pre_time

            if (exists("weightMat")) {
              links <- GENIE3::getLinkList(weightMat)
              rm(weightMat)
              links <- links[which(links$targetGene == "Gene"), ]
              links$targetGene <- Gene

              GENIE3_matrix <- 0
              for (g in 1:nrow(links)) {
                gene <- links$regulatoryGene[g]
                if (gene != targetGene) {
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
            name = c(targetGene, links_L0$TF),
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
            name = c(targetGene, links$regulatoryGene),
            carac = c(rep("Gene", 1), rep("TF", nrow(links)))
          )

          network <- graph_from_data_frame(d = links, vertices = nodes, directed = F)
          # Make a palette of 3 colors
          coul <- brewer.pal(3, "Set1")
          # Create a vector of color
          my_color <- coul[as.numeric(as.factor(V(network)$carac))]
          my_color

          png(paste0("../../Results/L0Reg-framework/cluster", j, "/", targetGene, "_GENIE3-L0.png"),
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
            Gene = targetGene,
            L0 = cor(X["Gene", ], L0_matrix),
            GENIE3 = cor(X["Gene", ], GENIE3_matrix)
          )
          Contrast_cluster <- rbind.data.frame(Contrast_cluster, Contrast)
        }
      }
    } else {
      print(paste0("No ", targetGene, " TF file......"))
      break
    }
  }

  Contrast_all <- rbind.data.frame(Contrast_all, Contrast_cluster)

  deg.data_clusters <- rbind.data.frame(deg.data_clusters, deg.data_cluster)

  write.csv(evaluate_cluster, sprintf("../../Results/L0Reg-framework/cluster%s_evaluate.csv", j))

  evaluate_all <- rbind(evaluate_all, evaluate_cluster)
  res_TRN_cluster <- as.data.frame(res_TRN_cluster)
  write.csv(res_TRN_cluster, sprintf("../../Results/L0Reg-framework/res_TRN_cluster%s.csv", j))
  res_TRN_all <- rbind(res_TRN_all, res_TRN_cluster)
  survival_tfs <- survival_tfs[!duplicated(survival_tfs), ] %>% as.data.frame()
  write.csv(survival_tfs, paste0("../../Results/L0Reg-framework/genes_list_cluster", j, ".csv"))
  
  message(paste0("Cluster", j, "done......"))
}

write.csv(res_TRN_all, "../../Results/L0Reg-framework/res_TRN_all.csv", row.names = T)
write.csv(evaluate_all, "../../Results/L0Reg-framework/evaluate_all.csv", row.names = F)
write.csv(Contrast_all, "../../Results/L0Reg-framework/Contrast_all.csv", row.names = T)

# Plot
Contrast_all_box <- Contrast_all[, -c(1, 2)]
names(Contrast_all_box) <- c("L0Reg framework", "GENIE3")

Contrast_all_box <- Contrast_all_box %>%
  as.data.frame() %>%
  pivot_longer(
    cols = 1:2,
    names_to = "Method",
    values_to = "Correction"
  )

ggplot(data = Contrast_all_box, aes(x = Method, y = Correction)) +
  geom_boxplot(aes(fill = Method)) +
  stat_compare_means() +
  theme_bw() +
  scale_fill_manual(values = c("white", "gray")) +
  geom_jitter()

ggsave(paste("../../Results/L0Reg-framework/", "boxplot.png", sep = ""),
  width = 4,
  height = 3,
  dpi = 600
)
