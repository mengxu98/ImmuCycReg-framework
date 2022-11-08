

Sys.setenv(LANG="en_US.UTF-8")
# Packages check, download and library --------------------------------------------------
package.check <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      if (!requireNamespace("dplyr", quietly = TRUE)) {
        install.packages("dplyr")
      }
      library(dplyr)
      if (!requireNamespace("rvest", quietly = TRUE)) {
        install.packages("rvest")
      }
      library(rvest)
      message("[", Sys.time(), "] -----: No package: ", package, " in R environment!")
      CRANpackages <- available.packages() %>%
        as.data.frame() %>%
        select(Package) %>%
        mutate(source = "CRAN")
      url <- "https://www.bioconductor.org/packages/release/bioc/"
      biocPackages <- url %>%
        read_html() %>%
        html_table() %>%
        .[[1]] %>%
        select(Package) %>%
        mutate(source = "BioConductor")
      if (package %in% CRANpackages$Package) {
        message("[", Sys.time(), "] -----: Now install package: ", package, " from CRAN!")
        install.packages(package)
        library(package, character.only = T)
      } else if (package %in% biocPackages$Package) {
        message("[", Sys.time(), "] -----: Now install package: ", package, " from BioConductor!")
        BiocManager::install(package)
        library(package, character.only = T)
      } else { # Bug
        if (!requireNamespace("githubinstall", quietly = TRUE)) {
          install.packages("githubinstall")
        }
        library(githubinstall)
        # githubinstall(package)
        gh_suggest(package)
      }
    } else {
      library(package, character.only = T)
    }
  }
}

packages <- c("backports", "BiocManager", "usethis", "devtools", "dplyr",
             "DT", "feather", "ggthemes", "gplots", "grImport2",
             "gWidgetsRGtk2", "HH", "igraph", "markdown", "Metrics",
             "networkD3", "plotly", "plotmo", "processx", "RColorBrewer",
             "rmarkdown", "tidymodels", "tidyverse", "bgsmtr", "BiocVersion",
             "cgdsr", "clusterProfiler", "ComplexHeatmap", "DESeq", "pheatmap",
             "pathview", "rtracklayer", "org.Hs.eg.db", "L0Learn")
package.check(packages)
