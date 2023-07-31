rm(list = ls())
source("functions/Functions.R")

pathRead <- "../data/"
pathSave <- "../../Results/"

# ------------------------------------------------------------------------------#
load(paste0(pathSave, "TCGA-LUAD.Rdata"))
load(paste0(pathSave, "GTEx-LUAD.Rdata"))
LM22Genes <- read.table("../data/LM22_2.0.txt",
                        header = T,
                        # row.names = 1,
                        sep = "\t"
) %>% .[, 1]

sampleLabel <- read.table("../data/sample_cluster_4.csv",
                          header = F,
                          sep = ",",
                          check.names = FALSE,
                          col.names = c("sample", "cluster")
)

rawTCGA <- tcga_luad[LM22Genes, ] %>% na.omit()
rawGTEx <- gtex_luad[LM22Genes, ] %>% na.omit()

dataList <- list()
dataList[[1]] <- rawTCGA[, sampleLabel %>% dplyr::filter(cluster == 1) %>% .[, 1]]
dataList[[2]] <- rawTCGA[, sampleLabel %>% dplyr::filter(cluster == 2) %>% .[, 1]]
dataList[[3]] <- rawTCGA[, sampleLabel %>% dplyr::filter(cluster == 3) %>% .[, 1]]
dataList[[4]] <- rawTCGA[, sampleLabel %>% dplyr::filter(cluster == 4) %>% .[, 1]]
dataList[[5]] <- rawGTEx

CibersortData <- map_dfc(1:5, function(x){
  dataList[[x]]
})

write.table(CibersortData, "Cibersort_data_ALL.txt",
            sep = "\t",
            quote = F,
            row.names = T
)

# -----------------------------------------------------------------------------

data_ALL_group <- map_dfc(1:5, function(x){
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
