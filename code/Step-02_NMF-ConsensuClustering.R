#!/usr/bin/env Rscript

if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
  stop("Please install ConsensusClusterPlus before running this script.")
}

path_save <- "results"
expr_file <- file.path(path_save, "TCGA-LUAD.Rdata")
gene_file <- file.path("data", "Genes_2230.csv")
out_dir <- file.path(path_save, "NMF")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(expr_file)
if (!exists("tcga_luad")) {
  stop("The input file must contain an object named 'tcga_luad'.")
}

genes2230 <- read.csv(gene_file, stringsAsFactors = FALSE)
if (!"Gene" %in% colnames(genes2230)) {
  stop("Genes_2230.csv must contain a 'Gene' column.")
}

tcga_luad <- tcga_luad[, !duplicated(colnames(tcga_luad))]
genes <- genes2230$Gene
missing_genes <- setdiff(genes, rownames(tcga_luad))
if (length(missing_genes) > 0) {
  stop(
    "The following genes are missing from tcga_luad: ",
    paste(missing_genes, collapse = ", ")
  )
}

data_nmf <- as.matrix(tcga_luad[genes, ])
data_nmf <- log(data_nmf + 1, 10)

genes_nums <- 30
mads <- apply(data_nmf, 1, mad)
dataset <- data_nmf[rev(order(mads))[seq_len(genes_nums)], ]

max_cluster_nums <- 6
seed <- 135792468

results <- ConsensusClusterPlus::ConsensusClusterPlus(
  dataset,
  maxK = max_cluster_nums,
  reps = 50,
  pItem = 0.8,
  pFeature = 0.8,
  clusterAlg = "pam",
  distance = "euclidean",
  seed = seed,
  title = out_dir,
  corUse = "complete.obs",
  plot = "pdf"
)

icl <- ConsensusClusterPlus::calcICL(results, title = out_dir, plot = "pdf")

for (k in 2:max_cluster_nums) {
  cluster <- data.frame(
    sample = names(results[[k]][["consensusClass"]]),
    cluster = as.integer(results[[k]][["consensusClass"]]),
    stringsAsFactors = FALSE
  )
  write.table(
    cluster,
    file = file.path(out_dir, paste0("sample_cluster_", k, ".csv")),
    row.names = FALSE,
    col.names = FALSE,
    sep = ","
  )
}

Kvec <- 2:max_cluster_nums
x1 <- 0.1
x2 <- 0.9
PAC <- rep(NA_real_, length(Kvec))
names(PAC) <- paste0("K=", Kvec)

for (k in Kvec) {
  M <- results[[k]]$consensusMatrix
  Fn <- ecdf(M[lower.tri(M)])
  PAC[k - 1] <- Fn(x2) - Fn(x1)
}

write.csv(
  data.frame(K = Kvec, PAC = PAC),
  file.path(out_dir, "pac_values.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "ConsensusClusterPlus/PAM clustering",
    paste("Input expression file:", expr_file),
    paste("Input gene file:", gene_file),
    paste("Samples after duplicate removal:", ncol(tcga_luad)),
    paste("Feature genes:", length(genes)),
    paste("MAD-selected genes:", paste(rownames(dataset), collapse = ", ")),
    paste("Parameters: maxK=6, reps=50, pItem=0.8, pFeature=0.8, clusterAlg=pam, distance=euclidean, seed=", seed, sep = ""),
    "K=4 cluster sizes:",
    paste(capture.output(print(table(results[[4]][["consensusClass"]]))), collapse = "\n"),
    paste("Optimal K by PAC:", Kvec[which.min(PAC)])
  ),
  file.path(out_dir, "run_summary.txt")
)

print(icl[["itemConsensus"]][1:6, ])
print(table(results[[4]][["consensusClass"]]))
