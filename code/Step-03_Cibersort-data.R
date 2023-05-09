rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

# TCGA --------------------------------------------------------------------
tcga_raw <- read.table(paste0(pathRead, "luad-rsem-fpkm-tcga-t_normlized.txt.gz"),
                       header = T,
                       sep = "\t",
                       row.names = 1,
                       check.names = FALSE
) %>% .[, -1]
names(tcga_raw) <- substr(colnames(tcga_raw), 1, 15)
raw_tcga <- apply(tcga_raw, 2, fpkmToTpm) %>% as.data.frame()
write.table(raw_tcga, paste0(pathSave, "TCGA-LUAD-TPM_normlized.txt"), sep = "\t", quote = F)

# GTEx --------------------------------------------------------------------
gtex_raw <- read.table(paste0(pathRead, "lung-rsem-fpkm-gtex_normlized.txt.gz"),
                       header = T,
                       sep = "\t",
                       row.names = 1,
                       check.names = FALSE
) %>% .[, -1]
names(gtex_raw) <- substr(colnames(gtex_raw), 1, 15)
raw_gtex <- apply(gtex_raw, 2, fpkmToTpm)
raw_gtex <- as.data.frame(raw_gtex)
write.table(raw_gtex, paste0(pathSave, "GTEx-LUAD-TPM_normlized.txt"), sep = "\t", quote = F)

# ------------------------------------------------------------------------------#
load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "GTEx-LUAD.Rdata"))
LM22_2.0_genes <- read.table("../data/LM22_2.0.txt",
                            header = T,
                            # row.names = 1,
                            sep = "\t"
) %>% .[, 1]

sample_label <- read.table("../data/sample_cluster_4.csv",
                           header = F,
                           sep = ",",
                           check.names = FALSE,
                           col.names = c("sample", "cluster")
)

raw_tcga <- tcga_luad[LM22_2.0_genes, ] %>% na.omit()
raw_gtex <- gtex_luad[LM22_2.0_genes, ] %>% na.omit()

dataList <- list()
dataList[[1]] <- raw_tcga[, sample_label %>% dplyr::filter(cluster == 1) %>% .[, 1]]
dataList[[2]] <- raw_tcga[, sample_label %>% dplyr::filter(cluster == 2) %>% .[, 1]]
dataList[[3]] <- raw_tcga[, sample_label %>% dplyr::filter(cluster == 3) %>% .[, 1]]
dataList[[4]] <- raw_tcga[, sample_label %>% dplyr::filter(cluster == 4) %>% .[, 1]]
dataList[[5]] <- raw_gtex

Cibersort_data_tcga_gtex <- map_dfc(1:5, function(x){
  dataList[[x]]
})

write.table(Cibersort_data_tcga_gtex, "Cibersort_data_ALL.txt",
            sep = "\t",
            quote = F,
            row.names = T
)

# -----------------------------------------------------------------------------

data_ALL_group1 <- map_dfc(1:5, function(x){
  tcga_cluster1_t <- as.data.frame(t(dataList[[x]]))
  tcga_cluster1_t$group <- "Cluster1"
  as.data.frame(t(tcga_cluster1_t))
})

write.table(data_ALL_group, 
            "data_ALL_group.txt",
            sep = "\t",
            quote = F,
            row.names = T
)
