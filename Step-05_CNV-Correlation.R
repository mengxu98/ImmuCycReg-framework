

rm(list = ls())
library(psych)
library(ggcorrplot)

load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "CNV-LUAD.Rdata"))

# Split data in the target genes and samples with mRNA data
sample_label_sur <- read.table("../results/2191/NMF/cluster-nk-ligands-k=8/sample_cluster_4.csv",
  header = FALSE,
  sep = ",",
  check.names = FALSE
)

correction_results_all <- c()
for (k in 1:4) {
  cluster <- sample_label_sur[which(sample_label_sur$V2 == k), 1]
  if (TRUE) {
    gene_list <- read.table(paste("../data CNV/all genes", ".txt", sep = ""),
      header = TRUE,
      row.names = 1
    )
  } else {
    gene_list <- read.table(sprintf("../L0/data/genes_list_cluster%s.txt", k),
      header = TRUE,
      row.names = 1
    )
  }

  booble <- colnames(raw_tcga_cnv) %in% cluster
  CNV_cluster <- raw_tcga_cnv[row.names(gene_list), booble]
  mRNA_cluster <- raw_tcga[row.names(gene_list), colnames(CNV_cluster)]

  cnv_temp <- t(CNV_cluster)
  mrna_temp <- t(mRNA_cluster)
  mrna_temp <- log((mrna_temp + 1), 2)
  num_ligands <- ncol(cnv_temp)

  cor_results <- c()
  cor_results_p <- c()
  for (i in 1:num_ligands) {
    corr_ligands <- corr.test(
      cnv_temp[, i],
      mrna_temp[, i],
      method = "spearman",
      adjust = "fdr"
    )
    cor_results <- c(cor_results, corr_ligands$r)
    cor_results_p <- c(cor_results_p, corr_ligands$p)
  }
  strong_cor_self <- order(cor_results, decreasing = TRUE)
  colnames(cnv_temp)[strong_cor_self[1:3]]
  colnames(cnv_temp)[strong_cor_self[length(strong_cor_self)]]
  negative_value <- cor_results[strong_cor_self[length(strong_cor_self)]]

  intercor_matrix <- c()
  for (i in 1:num_ligands) {
    intercor_results <- c()
    for (j in 1:num_ligands) {
      corr_ligands <- cor(
        cnv_temp[, i],
        mrna_temp[, j],
        method = "spearman",
        use = "pairwise.complete.obs"
      )
      intercor_results <- c(intercor_results, corr_ligands)
    }
    intercor_matrix <- rbind(intercor_matrix, intercor_results)
  }
  rownames(intercor_matrix) <- colnames(cnv_temp)
  colnames(intercor_matrix) <- colnames(cnv_temp)

  correction_results <- c()
  for (i in 1:nrow(mRNA_cluster)) {
    selectgene <- row.names(mRNA_cluster[i, ])
    intercor_matrix[selectgene, selectgene]
    if (all(intercor_matrix[selectgene, selectgene] >= 0.4)) {
      print(selectgene)
    } else {
      print("The correlation coefficient of selected genes was less than 0.4")
    }
    correction_results1 <- cbind(selectgene, intercor_matrix[selectgene, selectgene])
    correction_results <- rbind(correction_results, correction_results1)
  }
  write.table(correction_results,
    paste0(pathSave, "correction_results_cluster", k, ".txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  correction_results <- as.data.frame(correction_results)
  correction_results$Clsuer <- paste0("Cluster", k)
  correction_results_all <- rbind.data.frame(correction_results_all, correction_results)
}
names(correction_results_all) <- c("Gene", "CNVs", "Cluster")
write.csv(correction_results_all, paste0(pathSave, "CNV values.csv"), row.names = FALSE)
