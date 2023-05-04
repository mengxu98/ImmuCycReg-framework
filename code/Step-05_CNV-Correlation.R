rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "CNV-LUAD.Rdata"))

# Split data in the target genes and samples with mRNA data
samples_cluster <- read.csv(paste0(pathRead, "sample_cluster_4.csv"), header = F)
correction_results_all <- c()
for (k in 1:4) {
  cluster <- samples_cluster[which(samples_cluster$V2 == k), 1]
  gene_list <- read.table(paste(paste0(pathRead, "Genes_17.txt"), sep = ""),
    header = TRUE
  )
  samples_inter <- intersect(cluster, colnames(tcga_cnv))
  CNV_cluster <- tcga_cnv[gene_list$gene, samples_inter]
  mRNA_cluster <- tcga_luad[gene_list$gene, samples_inter]

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
    # if (all(intercor_matrix[selectgene, selectgene] >= 0.4)) {
    #   print(selectgene)
    # } else {
    #   print("The correlation coefficient of selected genes was less than 0.4")
    # }
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
