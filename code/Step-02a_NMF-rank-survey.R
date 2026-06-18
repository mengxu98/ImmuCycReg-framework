#!/usr/bin/env Rscript

if (!requireNamespace("NMF", quietly = TRUE)) {
  stop("Please install the NMF package before running this script.")
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
data_nmf <- log(data_nmf + 1, 2)

genes_nums <- 60
mads <- apply(data_nmf, 1, mad)
dataset <- data_nmf[rev(order(mads))[seq_len(genes_nums)], ]
dataset <- na.omit(dataset)

nrun <- 50
seed <- 20220101
ranks <- 2:6

rank_survey <- NMF::nmfEstimateRank(
  dataset,
  ranks,
  nrun = nrun,
  stop = FALSE,
  seed = seed,
  verbose = TRUE
)

rank_fit <- NMF::nmf(
  dataset,
  ranks,
  nrun = nrun,
  seed = seed
)

pdf(file.path(out_dir, "NMF.pdf"), width = 8, height = 6)
plot(rank_fit, rank_survey)
dev.off()

pdf(file.path(out_dir, "NMF_cophenetic.pdf"), width = 6, height = 5)
plot(rank_fit, c("cophenetic"))
dev.off()

rank4_fit <- NMF::nmf(
  dataset,
  4,
  nrun = nrun,
  seed = seed
)

pdf(file.path(out_dir, "NMF_consensus_rank4.pdf"), width = 9, height = 8)
NMF::consensusmap(
  rank4_fit,
  labRow = NA,
  labCol = NA
)
dev.off()

rank4_cluster <- predict(rank4_fit)
cluster4 <- data.frame(
  sample = names(rank4_cluster),
  cluster = as.integer(rank4_cluster),
  stringsAsFactors = FALSE
)

write.table(
  cluster4,
  file = file.path(out_dir, "nmf_rank4_predicted_clusters.csv"),
  row.names = FALSE,
  col.names = FALSE,
  sep = ","
)

writeLines(
  c(
    "NMF rank survey and rank-4 consensus map",
    paste("Input expression file:", expr_file),
    paste("Input gene file:", gene_file),
    paste("Samples after duplicate removal:", ncol(tcga_luad)),
    paste("Feature genes:", length(genes)),
    paste("MAD-selected genes:", paste(rownames(dataset), collapse = ", ")),
    paste("Parameters: ranks=2:6, nrun=", nrun, ", seed=", seed, sep = ""),
    "Rank-4 NMF cluster sizes:",
    paste(capture.output(print(table(cluster4$cluster))), collapse = "\n")
  ),
  file.path(out_dir, "NMF_run_summary.txt")
)

print(rank_fit$measures)
print(table(cluster4$cluster))
