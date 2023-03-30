

source("Cibersort.R")
cibersort_all <- CIBERSORT("data/LM22_2.0.txt",
                           "../Results/Cibersort_data_ALL.txt",
                           pathSave = "../Results/Cibersort/",
                           perm = 1000,
                           QN = F)
