

rm(list = ls())

source("Functions.R")

path_read <- "data/"
path_save <- "results/"

tcga_cnv <- read.table(paste0(path_read, "all_thresholded.by_genes_whitelisted.tsv.gz"),
     header = TRUE,
     row.names = 1,
     sep = "\t",
     check.names = FALSE
)
colnames(tcga_cnv) <- substr(colnames(tcga_cnv), 1, 15)

check_file.exists(tcga_cnv)

# tcga_cnv <- readRDS(paste0(path_save, "tcga_cnv.rds"))

tcga_luad <- read.table(paste0(path_read, "luad-rsem-count-tcga-t.txt.gz"),
     header = TRUE,
     row.names = 1,
     sep = "\t",
     check.names = FALSE
)
colnames(tcga_luad) <- substr(colnames(tcga_luad), 1, 15)

check_file.exists(tcga_luad)


gtex_luad <- read.table(paste0(path_read, "lung-rsem-count-gtex.txt.gz"),
     header = TRUE,
     row.names = 1,
     sep = "\t",
     check.names = FALSE
)
colnames(gtex_luad) <- substr(colnames(gtex_luad), 1, 15)

check_file.exists(gtex_luad)

peak_all <- read.table(paste0(path_read, "TCGA_ATAC_peak_Log2Counts_dedup_sample.gz"),
     header = TRUE,
     row.names = 1,
     check.names = FALSE
)
colnames(peak_all) <- substr(colnames(peak_all), 1, 15)

peak_luad <- peak_all[, which(colnames(peak_all) %in% colnames(tcga_luad))]
peak_luad <- peak_luad[grep(pattern = "LUAD", row.names(peak_luad)), ]

check_file.exists(peak_luad)

genes_adj_peak <- read.table(paste0(path_read, "TCGA_ATAC_peak.all.probeMap"),
     header = TRUE,
     row.names = 1,
     check.names = FALSE,
     sep = "\t"
)
genes_adj_peak <- genes_adj_peak[grep(pattern = "LUAD", row.names(genes_adj_peak)), ]
check_file.exists(genes_adj_peak)

tcga_luad_peak_samples <- tcga_luad[, which(colnames(tcga_luad) %in% colnames(peak_luad))]
check_file.exists(tcga_luad_peak_samples)

gtf <- rtracklayer::import(paste0(path_read, "hg38.ncbiRefSeq.gtf"))
gtf_df <- as.data.frame(gtf)
geneinfo_df <- dplyr::select(gtf_df, c(
     "seqnames", "start", "end", "width", "strand",
     "type", "gene_id", "transcript_id", "gene_name"
))
check_file.exists(geneinfo_df)
