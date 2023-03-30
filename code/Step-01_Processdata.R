

rm(list = ls())
source("Functions.R")
pathRead <- "data/"
pathSave <- "../Results/"
# --------------------------------------------------
tcga_luad <- read.table(paste0(pathRead, "luad-rsem-count-tcga-t.txt.gz"),
     header = TRUE,
     row.names = 1,
     sep = "\t",
     check.names = FALSE
)
colnames(tcga_luad) <- substr(colnames(tcga_luad), 1, 15)
tcga_luad <- tcga_luad[, -1]
save.file(tcga_luad, fileName = "TCGA-LUAD.Rdata", pathWay = pathSave)
# tcga_luad <- readRDS(paste0(pathSave, "tcga_luad.rds"))

# --------------------------------------------------
gtex_luad <- read.table(paste0(pathRead, "lung-rsem-count-gtex.txt.gz"),
     header = TRUE,
     row.names = 1,
     sep = "\t",
     check.names = FALSE
)
colnames(gtex_luad) <- substr(colnames(gtex_luad), 1, 15)
gtex_luad <- gtex_luad[, -1]
save.file(gtex_luad, fileName = "GTEx-LUAD.Rdata", pathWay = pathSave)

# --------------------------------------------------
tcga_cnv <- read.table(paste0(pathRead, "all_thresholded.by_genes_whitelisted.tsv.gz"),
     header = TRUE,
     row.names = 1,
     sep = "\t",
     check.names = FALSE
)
colnames(tcga_cnv) <- substr(colnames(tcga_cnv), 1, 15)

tcga_cnv <- tcga_cnv[
     intersect(rownames(tcga_luad), rownames(tcga_cnv)),
     intersect(colnames(tcga_luad), colnames(tcga_cnv))
]
save.file(tcga_cnv, fileName = "CNV-LUAD.Rdata", pathWay = pathSave)

# --------------------------------------------------
peak_all <- read.table(paste0(pathRead, "TCGA_ATAC_peak_Log2Counts_dedup_sample.gz"),
     header = TRUE,
     row.names = 1,
     check.names = FALSE
)
colnames(peak_all) <- substr(colnames(peak_all), 1, 15)

peak_luad <- peak_all[
     grep(pattern = "LUAD", rownames(peak_all)),
     intersect(colnames(tcga_luad), colnames(peak_all))
]
tcga_luad_peak_samples <- tcga_luad[, which(colnames(tcga_luad) %in% colnames(peak_luad))]
save.file(tcga_luad_peak_samples, fileName = "Peak-LUAD.Rdata", pathWay = pathSave)

# --------------------------------------------------
genes_adj_peak <- read.table(paste0(pathRead, "TCGA_ATAC_peak.all.probeMap.gz"),
     header = TRUE,
     row.names = 1,
     check.names = FALSE,
     sep = "\t"
)
genes_adj_peak <- genes_adj_peak[grep(pattern = "LUAD", rownames(genes_adj_peak)), ]
save.file(genes_adj_peak, fileName = "ATAC-LUAD.Rdata", pathWay = pathSave)

# --------------------------------------------------
gtf <- rtracklayer::import(paste0(pathRead, "hg38.ncbiRefSeq.gtf.gz"))
gtf_df <- as.data.frame(gtf)
geneinfo_df <- dplyr::select(gtf_df, c(
     "seqnames",
     "start",
     "end",
     "width",
     "strand",
     "type",
     "gene_id",
     "transcript_id",
     "gene_name"
))
check.file.exists()
save.file(geneinfo_df, fileName = "Geneinfo_df.Rdata", pathWay = pathSave)

# survival.data(cancerType = "luad_tcga", immuneGene = "IL2") # Immune genes (accept genes list) need to be selected according to analysis!