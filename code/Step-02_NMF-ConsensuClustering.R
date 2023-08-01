rm(list = ls())

source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

load(paste0(pathSave, "TCGA-LUAD.Rdata"))
genes2230 <- read.csv(paste0(pathRead, "genes2230.csv"))
dataNMF <- tcga_luad[genes2230$Gene, ] %>% as.matrix()

# --------------------------------------------------
seed <- "20210101"
dataNMF <- log(dataNMF + 1, 10)
genesNums <- 30
mads <- apply(dataNMF, 1, mad)
dataNMF <- dataNMF[rev(order(mads)), ]
dataset <- dataNMF[1:genesNums, ]

maxClusterNums <- 6
title <- paste0(pathSave, "NMF/")
results <- ConsensusClusterPlus(
  dataset,
  maxK = maxClusterNums,
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
icl[["itemConsensus"]][1:6, ]

for (i in 2:maxClusterNums) {
  write.table(results[[i]][["consensusClass"]],
              file = paste0(title, "sample_cluster_", i, ".csv"),
              na = "",
              col.names = FALSE,
              sep = ","
  )
}

Kvec <- 2:6
x1 <- 0.1
x2 <- 0.9 # threshold defining the intermediate sub-interval
PAC <- rep(NA, length(Kvec))
names(PAC) <- paste0("K=", Kvec) # from 2 to maxK

for (i in Kvec) {
  M <- results[[i]]$consensusMatrix
  Fn <- ecdf(M[lower.tri(M)])
  PAC[i - 1] <- Fn(x2) - Fn(x1)
}

# The optimal K
optK <- Kvec[which.min(PAC)]
optK

table(results[[4]]$consensusClass)

Cluster <- results[[4]]$consensusClass %>% as.data.frame()
Cluster$sample <- rownames(Cluster)

survival.data(cancerType = "luad_tcga", genes = "IL2", pathWay = pathSave)
load(paste0(pathSave, "survival_input.Rdata"))
rownames(Cluster) <- gsub("-", ".", rownames(Cluster))
samples <- intersect(rownames(Cluster), rownames(ClinicalData))
Cluster <- Cluster[samples, ] %>% as.data.frame()
ClinicalData <- ClinicalData[samples, ]
identical(rownames(Cluster), rownames(ClinicalData))
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
ClinicalData <- ClinicalData[, choose_columns]
meta <- ClinicalData
meta$OS_STATUS <- gsub("1:DECEASED", "1", meta$OS_STATUS)
meta$OS_STATUS <- gsub("0:LIVING", "0", meta$OS_STATUS)

meta$Cluster <- as.character(Cluster$.)

sfit <- survfit(Surv(OS_MONTHS, OS_STATUS == "1") ~ Cluster,
                data = meta
)
ggsurvplot(sfit,
           pval = TRUE,
           palette = c("#0050ef", "#008a00", "#fa6800", "#ffff88")
)

###
exp <- dataset %>% as.data.frame()
colnames(exp) <- gsub("-", ".", colnames(exp))
exp <- exp[, samples]

# PCA
draw_pca(exp, Cluster$., addEllipses = F)
# Tsne
tsne_out <- Rtsne(t(exp), perplexity = 30)
pdat <- data.frame(tsne_out$Y, factor(Cluster$.))
colnames(pdat) <- c("Y1", "Y2", "group")

ggplot(pdat, aes(Y1, Y2)) +
  geom_point(aes(Y1, Y2, fill = group), shape = 21, color = "black") +
  stat_ellipse(aes(color = group, fill = group),
               geom = "polygon",
               alpha = 0.3,
               linetype = 2
  ) +
  scale_color_paletteer_d("RColorBrewer::Set3") +
  scale_fill_paletteer_d("RColorBrewer::Set3") +
  theme_bw() +
  theme(legend.position = "top")
