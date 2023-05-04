rm(list = ls())

pathRead <- "../data/"
pathSave <- "../../results/"

load(paste0(pathSave, "TCGA-LUAD.Rdata"))
feature_genes <- read.csv(paste0(pathRead, "Genes_2230.csv"))
dataset <- tcga_luad[feature_genes$Gene, ]

# --------------------------------------------------
# Determine the rank value
if (F) {
  seeds <- c(20220101, 2022)
  percent_genes <- c(0.05, 0.1, 0.15, 0.2, 30, 50, 100)
  nruns <- c(10, 20, 30, 40, 50)
  cophenetic_table_s <- c()
  res_list_s <- list()
  for (s in 1:length(seeds)) {
    seed <- seeds[s]
    res_list_p <- list()
    cophenetic_table_p <- c()
    for (p in 1:length(percent_genes)) {
      percent_gene <- percent_genes[p]
      dataset <- as.matrix(data_derived)
      dataset <- log(dataset + 1, 2)
      if (percent_gene < 1) {
        gene_no <- dim(dataset)[1] * percent_gene
      } else {
        gene_no <- percent_gene
      }
      mads <- apply(dataset, 1, mad)
      dataset <- dataset[rev(order(mads))[1:gene_no], ]
      dataset <- na.omit(dataset)

      res_list_n <- list()
      cophenetic_table_n <- c()
      for (n in 1:length(nruns)) {
        nrun <- nruns[n]
        process.check <- tryCatch(
          {
            res <- nmf(dataset, 2:6, nrun = nrun, seed = seed)
          },
          error = function(e) {
            next
            print(1)
          }
        )
        cophenetic <- cbind.data.frame(
          percent_gene, seed, nrun,
          res$measures$cophenetic[1],
          res$measures$cophenetic[2],
          res$measures$cophenetic[3],
          res$measures$cophenetic[4],
          res$measures$cophenetic[5]
        )
        cophenetic_table_n <- rbind.data.frame(cophenetic_table_n, cophenetic)
        res_list_n[[n]] <- res
      }
      cophenetic_table_p <- rbind.data.frame(cophenetic_table_p, cophenetic_table_n)
      res_list_p[[p]] <- res_list_n
    }
    cophenetic_table_s <- rbind.data.frame(cophenetic_table_s, cophenetic_table_p)
    res_list_p[[s]] <- res_list_p
  }
}

plot(res_list[[1]])

nrun <- 50
seed <- "20210101"
data_nmf <- log(data_nmf + 1, 10)
gene_no <- 30
mads <- apply(data_nmf, 1, mad)
data_nmf <- data_nmf[rev(order(mads)), ]
dataset <- data_nmf[1:gene_no, ]

rand <- nmfEstimateRank(dataset, 2:6, nrun = nrun, stop = FALSE, seed = seed, verbose = T)
res <- nmf(dataset, 2:6, nrun = nrun, seed = seed)
plot(res)
plot(res, c("cophenetic"))
plot(res, rand)
plot(res, rand, c("cophenetic", "silhouette"))

res$measures$cophenetic
V.random <- randomize(dataset)
estim.r.random <- nmf(V.random, 2:6, nrun = nrun, seed = seed)
plot(res, estim.r.random)

# --------------------------------------------------
Cluster <- predict(res_4) %>% as.data.frame()
table(Cluster)
write.table(Cluster,
  file = paste0(title, "sample_cluster_4.csv"),
  na = "",
  col.names = FALSE,
  sep = ","
)

Cluster$sample <- rownames(Cluster)
survival.data(cancerType = "luad_tcga", genes = "IL2", pathWay = pathSave)
load(paste0(pathSave, "survival_input.Rdata"))
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
  # palette = "jco",
  palette = c("#2874C5", "#008a00", "#C6524A", "#EABF00"), #
  xlab = "Survival time (month)"
)

dataset <- as.matrix(data_derived)
dataset <- log(dataset + 1, 2)
gene_no <- dim(dataset)[1] * 0.1
gene_no <- 30
mads <- apply(dataset, 1, mad) # mad(x)
dataset <- dataset[rev(order(mads))[1:gene_no], ]
dataset <- na.omit(dataset)
res_4 <- nmf(dataset, 4, nrun = 20, seed = 20220101)
group <- predict(res_4)
table(group)
jco <- c("#2874C5", "#EABF00", "#C6524A", "#868686")
index <- extractFeatures(res_4, "max")
sig.order <- unlist(index)
NMF.Exp.rank4 <- dataset[sig.order, ]
NMF.Exp.rank4 <- na.omit(NMF.Exp.rank4)

consensusmap(res_4,
  labRow = NA,
  labCol = NA,
  annCol = data.frame("cluster" = group[colnames(NMF.Exp.rank4)]),
  annColors = list(cluster = c("1" = jco[1], "2" = jco[2], "3" = jco[3], "4" = jco[4]))
)

dataset <- as.matrix(data_derived)
dataset <- log(dataset + 1, 2)
gene_no <- dim(dataset)[1] * 0.15
gene_no <- 30
mads <- apply(dataset, 1, mad) # mad(x)
dataset <- dataset[rev(order(mads))[1:gene_no], ]
dataset <- na.omit(dataset)
res_4 <- nmf(dataset, 4, nrun = 60, seed = 2022)
group <- predict(res_4)
###
arg <- "samples" # “columns”, “rows”, “samples”, “features”, “consensus”
si2 <- silhouette(res_2, what = arg)
plot(si2)
si3 <- silhouette(res_3, what = arg)
plot(si3)
si4 <- silhouette(res_4, what = arg)
plot(si4)

w <- basis(res_4)
dim(w)
h <- coef(res_4)
dim(h)
opar <- par(mfrow = c(1, 2))
# basis components
basismap(res_4, subsetRow = TRUE)
basismap(res_4)
# mixture coefficients
coefmap(res_4)
par(opar)
consensusmap(res_4)
consensusmap(res_4, annCol = dataset, tracks = NA)
