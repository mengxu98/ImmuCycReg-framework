rm(list = ls())

pathRead <- "../data/"
pathSave <- "../../Results/"

load(paste0(pathSave, "TCGA-LUAD-tpm-tcga-t_unnormlized.Rdata"))
load(paste0(pathSave, "GTEx-LUAD-tpm-gtex_unnormlized.Rdata"))

LM22Genes <- read.table("../data/LM22_2.0.txt",
                        header = TRUE,
                        sep = "\t") %>% .[, 1]

sampleLabel <- read.table("../data/sample_cluster_4.csv",
                          header = FALSE,
                          sep = ",",
                          check.names = FALSE,
                          col.names = c("sample", "cluster"))

rawTCGA <- tcga_luad[LM22Genes, ] %>% na.omit()
rawGTEx <- gtex_luad[LM22Genes, ] %>% na.omit()

dataList <- list()
for (i in 1:4) {
  dataList[[i]] <- rawTCGA[, sampleLabel %>% dplyr::filter(cluster == i) %>% .[, 1]]
}
dataList[[5]] <- rawGTEx

CibersortData <- purrr::map_dfc(1:5, function(x){
  dataList[[x]]
})

write.table(CibersortData, "Cibersort_data_ALL.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE
)

CibersortDataLable <- purrr::map_dfc(1:5, function(x){
  cluster <- as.data.frame(t(dataList[[x]]))
  if (x < 5) {
    clusterLable <- paste0("Cluster", x)
  } else {
    clusterLable <- "GTEx"
  }
  cluster$group <- clusterLable
  as.data.frame(t(cluster))
})

write.table(CibersortDataLable, 
            "data_ALL_group.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE
)
