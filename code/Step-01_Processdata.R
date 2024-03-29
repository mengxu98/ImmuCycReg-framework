rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

# --------------------------------------------------
dataClass <- c("count-tcga-t", "fpkm-tcga-t", "fpkm-tcga-t_normlized", "fpkm-tcga-t_unnormlized")
for (i in 1:length(dataClass)) {
  tcga_luad <- read.table(paste0(pathRead, paste0("luad-rsem-", dataClass[i],".txt.gz")),
                          header = TRUE,
                          row.names = 1,
                          sep = "\t",
                          check.names = FALSE) %>% .[, -1]
  colnames(tcga_luad) <- substr(colnames(tcga_luad), 1, 15)
  
  if (grepl("fpkm", dataClass[i])) {
    save.file(tcga_luad,
              fileName = paste0("TCGA-LUAD-", dataClass[i], ".Rdata"),
              pathWay = pathSave)
    
    tcga_luad <- fpkm.to.tpm(tcga_luad)
    save.file(tcga_luad,
              fileName = paste0("TCGA-LUAD-", gsub("fpkm", "tpm", dataClass[i]), ".Rdata"),
              pathWay = pathSave)
  } else {
    save.file(tcga_luad,
              fileName = "TCGA-LUAD.Rdata",
              pathWay = pathSave)
  }
}

dataClass <- c("count-gtex", "fpkm-gtex_normlized", "fpkm-gtex_unnormlized")
for (i in 1:length(dataClass)) {
  gtex_luad <- read.table(paste0(pathRead, paste0("lung-rsem-", dataClass[i],".txt.gz")),
                          header = TRUE,
                          row.names = 1,
                          sep = "\t",
                          check.names = FALSE) %>% .[, -1]
  colnames(gtex_luad) <- substr(colnames(gtex_luad), 1, 15)
  
  if (grepl("fpkm", dataClass[i])) {
    save.file(gtex_luad,
              fileName = paste0("GTEx-LUAD-", dataClass[i], ".Rdata"),
              pathWay = pathSave)
    
    gtex_luad <- fpkm.to.tpm(gtex_luad)
    save.file(gtex_luad,
              fileName = paste0("GTEx-LUAD-", gsub("fpkm", "tpm", dataClass[i]), ".Rdata"),
              pathWay = pathSave)
  } else {
    save.file(gtex_luad,
              fileName = "GTEx-LUAD.Rdata",
              pathWay = pathSave)
  }
}

# --------------------------------------------------
if (!file.exists(paste0(pathSave, "all_thresholded.by_genes_whitelisted.tsv"))) {
  download.file(
    "https://api.gdc.cancer.gov/data/7d64377f-2cea-4ee3-917f-8fcfbcd999e7",
    paste0(pathSave, "all_thresholded.by_genes_whitelisted.tsv"))
}
tcga_cnv <- read.table(paste0(pathSave, "all_thresholded.by_genes_whitelisted.tsv"),
                       header = TRUE,
                       row.names = 1,
                       sep = "\t",
                       check.names = FALSE)
colnames(tcga_cnv) <- substr(colnames(tcga_cnv), 1, 15)

tcga_cnv <- tcga_cnv[intersect(rownames(tcga_luad), rownames(tcga_cnv)),
                     intersect(colnames(tcga_luad), colnames(tcga_cnv))]
save.file(tcga_cnv,
          fileName = "CNV-LUAD.Rdata",
          pathWay = pathSave)

# --------------------------------------------------
if (!file.exists(paste0(pathSave, "TCGA_ATAC_peak_Log2Counts_dedup_sample.gz"))) {
  download.file(
    "https://tcgaatacseq.s3.us-east-1.amazonaws.com/download/TCGA_ATAC_peak_Log2Counts_dedup_sample.gz",
    paste0(pathSave, "TCGA_ATAC_peak_Log2Counts_dedup_sample.gz"))
}
peak_all <- read.table(paste0(pathSave, "TCGA_ATAC_peak_Log2Counts_dedup_sample.gz"),
                       header = TRUE,
                       row.names = 1,
                       check.names = FALSE)
colnames(peak_all) <- substr(colnames(peak_all), 1, 15)

peak_luad <- peak_all[grep(pattern = "LUAD", rownames(peak_all)),
                      intersect(colnames(tcga_luad), colnames(peak_all))]
tcga_luad_peak_samples <- tcga_luad[, which(colnames(tcga_luad) %in% colnames(peak_luad))]
save.file(tcga_luad_peak_samples,
          fileName = "Peak-LUAD.Rdata",
          pathWay = pathSave)

# --------------------------------------------------
if (!file.exists(paste0(pathSave, "TCGA_ATAC_peak.all.probeMap"))) {
  download.file(
    "https://tcgaatacseq.s3.us-east-1.amazonaws.com/download/TCGA_ATAC_peak.all.probeMap",
    paste0(pathSave, "TCGA_ATAC_peak.all.probeMap"))
}
genes_adj_peak <- read.table(paste0(pathSave, "TCGA_ATAC_peak.all.probeMap"),
                             header = TRUE,
                             row.names = 1,
                             check.names = FALSE,
                             sep = "\t")
genes_adj_peak <- genes_adj_peak[grep(pattern = "LUAD", rownames(genes_adj_peak)), ]
save.file(genes_adj_peak,
          fileName = "ATAC-LUAD.Rdata",
          pathWay = pathSave)

# --------------------------------------------------
if (!file.exists(paste0(pathSave, "hg38.ncbiRefSeq.gtf.gz"))) {
  download.file(
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz",
    paste0(pathSave, "hg38.ncbiRefSeq.gtf.gz"))
}
gtf_df <- rtracklayer::import(paste0(pathSave, "hg38.ncbiRefSeq.gtf.gz")) %>% as.data.frame()
geneinfo_df <- dplyr::select(gtf_df,
                             c("seqnames",
                               "start",
                               "end",
                               "width",
                               "strand",
                               "type",
                               "gene_id",
                               "transcript_id",
                               "gene_name"))
save.file(geneinfo_df,
          fileName = "Geneinfo_df.Rdata",
          pathWay = pathSave)
