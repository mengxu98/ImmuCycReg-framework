

rm(list = ls())
source("Functions.R")
library("tidyr")
pathRead <- "data/"
pathSave <- "../Results/"

# TCGA --------------------------------------------------------------------
tcga_raw <- read.table("data/luad-rsem-fpkm-tcga-t_normlized.txt.gz",
                       header = T,
                       sep = "\t",
                       row.names = 1,
                       check.names = FALSE
)
tcga_raw <- tcga_raw[, -1]
names(tcga_raw) <- substr(colnames(tcga_raw), 1, 15)
raw_tcga <- apply(tcga_raw, 2, fpkmToTpm) %>% as.data.frame()
write.table(raw_tcga, "../Results/TCGA-LUAD-TPM_normlized.txt", sep = "\t", quote = F)
save(raw_tcga, file = "../Results/TCGA-LUAD-TPM_normlized.Rdata")

# GTEx --------------------------------------------------------------------
gtex_raw <- read.table("data/lung-rsem-fpkm-gtex_normlized.txt.gz",
                       header = T,
                       sep = "\t",
                       row.names = 1,
                       check.names = FALSE
)
gtex_raw <- gtex_raw[, -1]
names(gtex_raw) <- substr(colnames(gtex_raw), 1, 15)
raw_gtex <- apply(gtex_raw, 2, fpkmToTpm)
raw_gtex[1:3, ]
colSums(raw_gtex)
raw_gtex <- as.data.frame(raw_gtex)
write.table(raw_gtex, "../Results/GTEx-LUAD-TPM_normlized.txt", sep = "\t", quote = F)
save(raw_gtex, file = "../Results/GTEx-LUAD-TPM_normlized.Rdata")
# ------------------------------------------------------------------------------#
load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "GTEx-LUAD.Rdata"))
LM22_2.0_gene <- read.table("data/LM22_2.0.txt",
                            header = T,
                            row.names = 1,
                            sep = "\t"
)

raw_tcga <- tcga_luad[row.names(LM22_2.0_gene), ]
raw_tcga <- na.omit(raw_tcga)
raw_gtex <- gtex_luad[row.names(LM22_2.0_gene), ]
raw_gtex <- na.omit(raw_gtex)
sample_label <- read.table("data/sample_cluster_4.csv",
                           header = F,
                           sep = ",",
                           check.names = FALSE
)
cluster1 <- sample_label[which(sample_label$V2 == 1), 1]
cluster2 <- sample_label[which(sample_label$V2 == 2), 1]
cluster3 <- sample_label[which(sample_label$V2 == 3), 1]
cluster4 <- sample_label[which(sample_label$V2 == 4), 1]

tcga_cluster1 <- raw_tcga[, cluster1]
tcga_cluster2 <- raw_tcga[, cluster2]
tcga_cluster3 <- raw_tcga[, cluster3]
tcga_cluster4 <- raw_tcga[, cluster4]
Cibersort_data_tcga_gtex <- cbind.data.frame(
  tcga_cluster1,
  tcga_cluster2,
  tcga_cluster3,
  tcga_cluster4,
  raw_gtex
)
write.table(Cibersort_data_tcga_gtex, "../Results/Cibersort_data_ALL.txt",
            sep = "\t",
            quote = F,
            row.names = T
)
# -----------------------------------------------------------------------------
tcga_cluster1_t <- as.data.frame(t(tcga_cluster1))
tcga_cluster2_t <- as.data.frame(t(tcga_cluster2))
tcga_cluster3_t <- as.data.frame(t(tcga_cluster3))
tcga_cluster4_t <- as.data.frame(t(tcga_cluster4))
raw_gtex_t <- as.data.frame(t(raw_gtex))

tcga_cluster1_t$group <- "Cluster1"
tcga_cluster2_t$group <- "Cluster2"
tcga_cluster3_t$group <- "Cluster3"
tcga_cluster4_t$group <- "Cluster4"
raw_gtex_t$group <- "GTEx"

tcga_cluster1_g <- as.data.frame(t(tcga_cluster1_t))
tcga_cluster2_g <- as.data.frame(t(tcga_cluster2_t))
tcga_cluster3_g <- as.data.frame(t(tcga_cluster3_t))
tcga_cluster4_g <- as.data.frame(t(tcga_cluster4_t))
raw_gtex_g <- as.data.frame(t(raw_gtex_t))
# -----------------------------------------------------------------------------
data_ALL_group <- cbind.data.frame(tcga_cluster1_g,
                                   tcga_cluster2_g,
                                   tcga_cluster3_g,
                                   tcga_cluster4_g,
                                   raw_gtex_g
)
write.table(data_ALL_group, 
            "../Results/data_ALL_group.txt",
            sep = "\t",
            quote = F,
            row.names = T
)
