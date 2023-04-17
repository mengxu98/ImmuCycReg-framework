

Sys.setenv(LANG = "en_US.UTF-8")
source("functions/Functions.R")
packages <- c(
  "backports", "BiocManager", "usethis", "devtools", "dplyr",
  "DT", "feather", "ggthemes", "gplots", "grImport2",
  "gWidgetsRGtk2", "HH", "igraph", "markdown", "Metrics",
  "networkD3", "plotly", "plotmo", "processx", "RColorBrewer",
  "rmarkdown", "tidymodels", "tidyverse", "bgsmtr", "BiocVersion",
  "cgdsr", "clusterProfiler", "ComplexHeatmap", "DESeq2", "pheatmap",
  "pathview", "rtracklayer", "org.Hs.eg.db", "L0Learn"
)
package.check(packages)
