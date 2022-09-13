

rm(list = ls())

library(NMF)
library(ConsensusClusterPlus)
library(limma)
library(doMPI)
library(pheatmap)
library(tidyr)
library(tidyverse)
library(survival)
library(survminer)

source("Functions.R")

path_read <- "data/"
path_save <- "results/"
seed <- "20210101"

tcga_luad <- readRDS(paste0(path_save, "tcga_luad.rds"))
genes_2230 <- read.csv(paste0(path_read, "Genes_2230.csv"))
data_nmf <- tcga_luad[genes_2230$Gene, ] %>% as.matrix()

data_nmf <- log(data_nmf + 1, 10)
gene_no <- 30
mads <- apply(data_nmf, 1, mad)
data_nmf <- data_nmf[rev(order(mads)), ]
dataset <- data_nmf[1:gene_no, ]


max_cluster_num <- 6
title <- paste0(path_save, "/NMF/cluster-nk-ligands-k=6")
results <- ConsensusClusterPlus(dataset,
  maxK = max_cluster_num,
  reps = 50,
  pItem = 0.8,
  pFeature = 0.8,
  clusterAlg = "pam",
  distance = "euclidean",
  seed = seed,
  title = title,
  corUse = "complete.obs",
  plot = "pdf"
)
icl <- calcICL(results, title = title, plot = "pdf")

icl[["itemConsensus"]][4:6, ]

for (i in 2:max_cluster_num) {
  write.table(results[[i]][["consensusClass"]],
    file = paste(title, paste(paste("sample_cluster_", as.character(i), sep = ""), ".csv", sep = ""), sep = "//"),
    na = "",
    col.names = FALSE,
    sep = ","
  )
}

Kvec <- 2:6
x1 <- 0.1
x2 <- 0.9 # threshold defining the intermediate sub-interval
PAC <- rep(NA, length(Kvec))
names(PAC) <- paste("K=", Kvec, sep = "") # from 2 to maxK

for (i in Kvec) {
  M <- results[[i]]$consensusMatrix
  Fn <- ecdf(M[lower.tri(M)])
  PAC[i - 1] <- Fn(x2) - Fn(x1)
} # end for i

# The optimal K
optK <- Kvec[which.min(PAC)]
optK

table(results[[4]]$consensusClass)

Cluster <- predict(res_2) %>% as.data.frame()
Cluster <- predict(res_3) %>% as.data.frame()
Cluster <- predict(res_4) %>% as.data.frame()

Cluster <- results[[4]]$consensusClass %>% as.data.frame()

Cluster$sample <- rownames(Cluster)
load("../SurvivalAnalysis/survival_LUAD.Rdata")
rownames(Cluster) <- gsub("-", ".", rownames(Cluster))
samples <- intersect(rownames(Cluster), rownames(myclinicaldata))
Cluster <- Cluster[samples, ] %>% as.data.frame()
myclinicaldata <- myclinicaldata[samples, ]
identical(rownames(Cluster), rownames(myclinicaldata))
choose_columns <- c(
  "AJCC_METASTASIS_PATHOLOGIC_PM",
  "AJCC_NODES_PATHOLOGIC_PN",
  "AJCC_PATHOLOGIC_TUMOR_STAGE",
  "AJCC_TUMOR_PATHOLOGIC_PT",
  "AGE",
  "SEX",
  "OS_STATUS",
  "OS_MONTHS",
  "DFS_MONTHS",
  "DFS_STATUS"
)
myclinicaldata <- myclinicaldata[, choose_columns]
meta <- myclinicaldata
meta$OS_STATUS <- gsub("1:DECEASED", "1", meta$OS_STATUS)
meta$OS_STATUS <- gsub("0:LIVING", "0", meta$OS_STATUS)

meta$Cluster <- as.character(Cluster$.)

sfit <- survfit(Surv(OS_MONTHS, OS_STATUS == "1") ~ Cluster,
  data = meta
)
ggsurvplot(sfit,
  pval = T,
  palette = c("#0050ef", "#008a00", "#fa6800", "#ffff88")
)

###

exp <- dataset %>% as.data.frame()
colnames(exp) <- gsub("-", ".", colnames(exp))
exp <- exp[, samples]
#
# PCA
library(tinyarray)
draw_pca(exp, Cluster$., addEllipses = F)
# Tsne
library(Rtsne)
tsne_out <- Rtsne(t(exp), perplexity = 30)
pdat <- data.frame(tsne_out$Y, factor(Cluster$.))
colnames(pdat) <- c("Y1", "Y2", "group")
head(pdat)

library(ggplot2)
library(paletteer)
ggplot(pdat, aes(Y1, Y2)) +
  geom_point(aes(Y1, Y2, fill = group), shape = 21, color = "black") +
  stat_ellipse(aes(color = group, fill = group),
    geom = "polygon",
    alpha = 0.3,
    linetype = 2
  ) +
  scale_color_paletteer_d("RColorBrewer::Set3") +
  scale_fill_paletteer_d("RColorBrewer::Set3") +
  theme_classic() +
  theme(legend.position = "top")

##########
sample_label <- read.table(paste(title, "/sample_cluster_4.csv", sep = ""),
  header = F,
  sep = ",",
  check.names = FALSE
)
cluster1 <- sample_label[which(sample_label$V2 == 1), 1]
cluster2 <- sample_label[which(sample_label$V2 == 2), 1]
cluster3 <- sample_label[which(sample_label$V2 == 3), 1]
cluster4 <- sample_label[which(sample_label$V2 == 4), 1]

tcga_cluster1 <- data_nmf[, cluster1]
tcga_cluster2 <- data_nmf[, cluster2]
tcga_cluster3 <- data_nmf[, cluster3]
tcga_cluster4 <- data_nmf[, cluster4]

group1_no <- dim(tcga_cluster1)[2]
group2_no <- dim(tcga_cluster2)[2]
group3_no <- dim(tcga_cluster3)[2]
group4_no <- dim(tcga_cluster4)[2]

group1 <- rep(c("C1"), each = group1_no)
group2 <- rep(c("C2"), each = group2_no)
group3 <- rep(c("C3"), each = group3_no)
group4 <- rep(c("C4"), each = group4_no)

group_label <- c(group1, group2, group3, group4)
annotation_c <- data.frame(group_label)

sorted_samples <- cbind(tcga_cluster1, tcga_cluster2)
sorted_samples <- cbind(sorted_samples, tcga_cluster3)
sorted_samples <- cbind(sorted_samples, tcga_cluster4)
rownames(annotation_c) <- colnames(sorted_samples)

pheatmap(sorted_samples,
  # scale='row',
  cluster_row = TRUE, cluster_col = FALSE,
  show_colnames = F,
  annotation_col = annotation_c
)

pheatmap(sorted_samples,
  cluster_row = TRUE,
  cluster_col = FALSE,
  show_colnames = F,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)
# validate data whether identical
smallset <- sorted_samples
bigset <- data_nmf
identical_flag <- TRUE
for (i in 1:nrow(smallset)) {
  for (j in 1:ncol(smallset)) {
    row_name <- rownames(smallset)[i]
    col_name <- colnames(smallset)[j]
    if (smallset[i, j] != bigset[row_name, col_name]) {
      print("not equal!! wrong!")
      identical_flag <- FALSE
    }
    if (smallset[i, j] == bigset[row_name, col_name]) {
      print(i)
    }
  }
}
if (identical_flag) {
  print("data is OK")
}

#------------------------------------------------------------------------------#
#### NMF####
dataset <- data_nmf
dataset <- as.matrix(data_derived)
dataset <- log(dataset + 1, 2)
gene_no <- dim(dataset)[1] * 0.1
gene_no <- 30
gene_no <- 100
mads <- apply(dataset, 1, mad)
dataset <- dataset[rev(order(mads))[1:gene_no], ]
dataset <- na.omit(dataset)

nrun <- 10
seed <- 135792468 # 2022

res <- nmf(dataset, 2:6, method = "KL", nrun = nrun, seed = seed)
res_ls <- nmf(dataset, 2:6, method = "ls-nmf", nrun = nrun, seed = seed)
plot(res)

summary(res)

seed <- 20220621
seed <- 2021
nruns <- c(10, 20, 30, 40, 50)
cophenetic_table_all <- c()

res_list <- list()

cophenetic_table <- c()
for (i in 1:length(nruns)) {
  nrun <- nruns[i]
  res <- nmf(dataset, 2:6, nrun = nrun, seed = seed)
  cophenetic <- cbind.data.frame(
    seed, nrun,
    res$measures$cophenetic[1],
    res$measures$cophenetic[2],
    res$measures$cophenetic[3]
  )
  cophenetic_table <- rbind.data.frame(cophenetic_table, cophenetic)
  res_list[[i]] <- res
}
cophenetic_table_all <- rbind.data.frame(cophenetic_table_all, cophenetic_table)

# generate a synthetic dataset with known classes: 50 features, 18 samples (5+5+8)
n <- 50
counts <- c(5, 5, 8)
V <- syntheticNMF(n, counts)
cl <- unlist(mapply(rep, 1:3, counts))

purity(res_4, predict(res_2))
entropy(res_4, predict(res_2))

plot(res_list[[1]])
plot(res_list[[2]])
plot(res_list[[3]])
plot(res_list[[4]])
plot(res_list[[5]])

plot(res_list[[1]], res_list[[3]])

res10 <- nmf(dataset, 2:6, nrun = 10, seed = seed)
plot(res, res10)

res$measures$cophenetic
res10$measures$cophenetic

V.random <- randomize(dataset)
estim.r.random <- nmf(V.random, 2:6, nrun = 10, seed = seed)
# plot measures on same graph
plot(res, estim.r.random)
plot(2:6, res$measures$cophenetic, type = "b", col = "purple")

aheatmap(dataset)
consensusmap(res)

res_2 <- nmf(dataset, 2, nrun = nrun, seed = seed)
res_3 <- nmf(dataset, 3, nrun = nrun, seed = seed)
res_4 <- nmf(dataset, 4, nrun = nrun, seed = seed)

jco <- c("#2874C5", "#EABF00", "#C6524A", "#868686")
index <- extractFeatures(res_4, "max")
sig.order <- unlist(index)
NMF.Exp.rank4 <- dataset[sig.order, ]
NMF.Exp.rank4 <- na.omit(NMF.Exp.rank4)
group <- predict(res_4)
consensusmap(res_4,
  labRow = NA,
  labCol = NA,
  annCol = data.frame("cluster" = group[colnames(NMF.Exp.rank4)]),
  annColors = list(cluster = c("1" = jco[1], "2" = jco[2], "3" = jco[3], "4" = jco[4]))
)
###

cluster2 <- predict(res_2) %>% as.data.frame()
table(cluster2$.)
cluster3 <- predict(res_3) %>% as.data.frame()
table(cluster3$.)
cluster4 <- predict(res_4) %>% as.data.frame()
table(cluster4$.)

cluster42 <- predict(res_42) %>% as.data.frame()
table(cluster42$.)

consensusmap(res_2)
consensusmap(res_3)
consensusmap(res_4)


arg <- "samples" # “columns”, “rows”, “samples”, “features”, “consensus”

si2 <- silhouette(res_2, what = arg)
plot(si2)

si3 <- silhouette(res_3, what = arg)
plot(si3)

si4 <- silhouette(res_4, what = arg)
plot(si4)

###############
w <- basis(res_4)
dim(w)
h <- coef(res_4)
dim(h)
opar <- par(mfrow = c(1, 2))
# basis components
basismap(res_4, subsetRow = TRUE) # subsetRow=TRUE, select features
basismap(res_4)
# mixture coefficients
coefmap(res_4)
par(opar)
consensusmap(res_4)
consensusmap(res_4, annCol = dataset, tracks = NA)

res_cluter4 <- read.csv("../results/2191/NMF/cluster-nk-ligands-k=8/sample_cluster_4.csv", header = F)
table(res_cluter4$V2)
boxplot(table(res_cluter4$V2))
hist(res_cluter4$V2)

res_cluter3 <- read.csv("../results/2191/NMF/cluster-nk-ligands-k=8/sample_cluster_3.csv", header = F)
table(res_cluter3$V2)
boxplot(table(res_cluter4$V2))
hist(res_cluter3$V2)
hist(res_cluter3$V2,
  breaks = c(0, 1, 2, 3),
  col = rainbow(4),
  border = NA
)

res_cluter2 <- read.csv("../results/2191/NMF/cluster-nk-ligands-k=8/sample_cluster_2.csv", header = F)
table(res_cluter2$V2)

