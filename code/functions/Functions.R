if (!requireNamespace("magrittr", quietly = TRUE)) install.packages("magrittr")
library(magrittr)

#' @title Check, install and library packages
#'
#' @param packages Packages list
#'
package.check <- function(packages) {
  for (i in 1:length(packages)) {
    package <- packages[i]
    message("Checking: '", package, "' in R environment......")
    if (grepl("/", package)) {
      if (!requireNamespace(strsplit(package, "/")[[1]][2], quietly = TRUE)) {
        if (!requireNamespace("devtools")) install.packages("devtools")
        message("Now install package: '", package, "' from Github......")
        tryCatch({
          devtools::install_github(package)
        }, error = function(e) {
          message("Unsuccessful installation package: '", package, "'......")
        })
        
        insRes <- try(library(strsplit(package, "/")[[1]][2], character.only = TRUE))
        
        if (class(insRes) == "try-error") {
          message("Unsuccessful install package: '", package, "'......")
          next
        }
      }
    } else {
      if (!requireNamespace(package)) {
        if (!requireNamespace("dplyr")) install.packages("dplyr")
        library(dplyr)
        if (!requireNamespace("rvest")) install.packages("rvest")
        library(rvest)
        CRANpackages <- available.packages() %>%
          as.data.frame() %>%
          # mutate(source = "CRAN") %>%
          # select(Package) %>%
          .[, 1]
        url <- "https://www.bioconductor.org/packages/release/bioc/"
        biocPackages <- url %>%
          read_html() %>%
          html_table() %>%
          .[[1]]  %>%
          as.data.frame() %>%
          # mutate(source = "BioConductor")  %>%
          # select(Package) %>% 
          .[, 1]
        
        tryCatch({
          if (package %in% CRANpackages) {
            message("Now install package: '", package, "' from CRAN......")
            install.packages(package)
          } else if (package %in% biocPackages) {
            message("Now install package: '", package, "' from BioConductor......")
            BiocManager::install(package)
          }
        }, error = function(e) {
          message("Unsuccessful installation package: '", package, "'......")
        })

        insRes <- try(library(package, character.only = TRUE))
        
        if (class(insRes) == "try-error") {
          message("Unsuccessful install package: '", package, "'......")
          next
        }

      }  else {
        library(package, character.only = TRUE)
      }
    }
  }
}

packages <- read.table("requirements.txt")
package.check(packages[, 1])

# Transformation of COUNT, TPM, FPKM
Rcpp::sourceCpp("functions/DataFormatConversion.cpp")
fpkm.to.tpm <- function(fpkmData) {
  fpkmMatrix <- as.matrix(fpkmData)
  tpmMatrix <- fpkmMatrixToTpmMatrixCpp(fpkmMatrix)
  rownames(tpmMatrix) <- rownames(fpkmData)
  colnames(tpmMatrix) <- colnames(fpkmData)
  return(as.data.frame(tpmMatrix))
}

#' Check dir exist
#'
#' @param dirPath 
#'
#' @return
#' @export
#'
check.dir <- function(dirPath) {
  if (dir.exists(dirPath)) {
    message("'", dirPath, "'", " existed......")
  } else {
    message("'", dirPath, "' not exist, creat it......")
    dir.create(dirPath, recursive = TRUE)
  }
  return(dirPath)
}

#' Save R object
#'
#' @param ...
#' @param fileName
#' @param pathWay
#'
#' @return
#' @export
#'
save.file <- function(..., fileName, pathWay = NULL) {
  if (is.null(pathWay)) {
    pathWay <- ""
  } else {
    check.dir(pathWay)
  }
  if (as.numeric(...length()) > 1) {
    if (grepl(fileName, pattern = ".Rdata$") | grepl(fileName, pattern = ".rdata$")) {
      save(..., file = paste0(pathWay, fileName))
    } else {
      newFileName <- sub("$", ".Rdata", fileName)
      save(..., file = paste0(pathWay, newFileName))
    }
  } else {
    save(..., file = paste0(pathWay, fileName))
  }
}

#' To obtain survival data of TCGA samples
#'
#' @param cancerType
#' @param genes It is required specify gene or genes list to obtain survival data
#' @param pathWay
#'
#' @return
#' @export
#'
survival.data <- function(cancerType = NULL,
                          genes = NULL,
                          pathWay = NULL) {
  if (is.null(cancerType)) stop("Pleasure ensure the cancer type......")
  if (is.null(genes)) stop("Pleasure input a single gene or gene list......")
  message("Choose '", cancerType, "' and preparing data......")
  package.check("cgdsr")
  package.check("DT")
  cgdsLoc <- cgdsr::CGDS("http://www.cbioportal.org/")
  
  # Get expression data
  expr <- getProfileData(cgdsLoc,
                         genes = genes,
                         caseList = paste0(cancerType, "_rna_seq_v2_mrna"),
                         geneticProfile = paste0(cancerType, "_rna_seq_v2_mrna"))
  
  # Get mutation data
  mut <- getProfileData(cgdsLoc,
                        genes = genes,
                        caseList = paste0(cancerType, "_sequenced"),
                        geneticProfile = paste0(cancerType, "_mutations"))
  mut <- apply(mut, 2, as.factor)
  mut[mut == "NaN"] <- ""
  mut[is.na(mut)] <- ""
  mut[mut != ""] <- "MUT"
  
  # Get copy number data
  cna <- getProfileData(cgdsLoc,
                        genes = genes,
                        caseList = paste0(cancerType, "_sequenced"),
                        geneticProfile = paste0(cancerType, "_gistic"))
  cnaNames <- rownames(cna)
  cna <- apply(cna, 2, function(x) {
    as.character(factor(x,
                        levels = c(-2:2),
                        labels = c("HOMDEL", "HETLOSS", "DIPLOID", "GAIN", "AMP")
    ))
  })
  cna[is.na(cna)] <- ""
  cna[cna == "DIPLOID"] <- ""
  rownames(cna) <- cnaNames
  
  # Get clinical data
  clinicalData <- getClinicalData(cgdsLoc,
                                  caseList = paste0(cancerType, "_rna_seq_v2_mrna"))
  
  if (is.null(pathWay)) {
    pathWay <- ""
  } else {
    check.dir(pathWay)
  }
  save(expr, clinicalData, cna, mut, file = paste0(pathWay, "survival_input.Rdata"))
  message("Successed make survuval data......")
}

#' formatPositiveResult
#'
#' @param inputData
#'
#' @return
#' @export
#'
format.regulation <- function(inputData,
                              targetGene,
                              regulation) {
  if (!(regulation %in% c(1, "positive", "-1", "negative"))) {
    stop("Please set 'regulation' as 1 / 'positive', or '-1' / 'negative'......")
  }
  TFs <- c()
  typeInter <- c()
  targetGeneInter <- c()
  sampleLevel <- c()
  if (length(inputData)) {
    i <- 1
    while (i <= length(inputData)) {
      j <- 1
      while (j <= length(inputData[[i]])) {
        TFs <- c(TFs, inputData[[i]][j])
        typeInter <- c(typeInter, regulation)
        sampleLevel <- c(sampleLevel, names(inputData)[i])
        targetGeneInter <- c(targetGeneInter, targetGene)
        j <- j + 1
      }
      i <- i + 1
    }
    formatResult <- cbind.data.frame(TFs, typeInter, targetGeneInter, sampleLevel)
    colnames(formatResult) <- c("TF", "Gene", "Sample", "Type")
    return(formatResult)
  } else {
    print("NULL list")
  }
}

#' FrameRegulatoryTable
#'
#' @param regulatoryTable
#'
#' @return
#' @export
#'
frame.regulatory.table <- function(regulatoryTable,
                                   gene,
                                   databasesLsit) {

  for (database in names(databasesLsit)) {
    weights <- c()
    for (TF in regulatoryTable[, 1]) {
      data <- databasesLsit[[database]]
      result <- data %>% filter(., source == targetGene & target == TF)
      if (nrow(result) > 0) {
        weights <- c(weights, result$weight)
      } else {
        weights <- c(weights, 0)
      }
    }
    weights <- as.data.frame(weights)
    names(weights) <- database
    regulatoryTable <- cbind(regulatoryTable, weights)
  }

  return(regulatoryTable)
}


#' combine.multiple.plot
#'
#' @param ggplotObj 
#' @param ncol 
#' @param legend 
#'
#' @return
#' @export
#'
combine.multiple.plot <- function(ggplotObject,
                                  ncol = 4,
                                  legend = "none") {
  if (!is.list(ggplotObject)) stop("Please provide a list of ggplot objects......")
  g <- 0
  for (i in 1:length(ggplotObject)) {
    x <- ggplotObject[[i]]
    if (is.ggplot(x)) {
      g = g + 1
      if (g == 1) {
        p <- x
      } else {
        p <- p + x
      }
    }
  }
  
  p <- p +
    plot_layout(ncol = ncol) +
    plot_layout(guides = "collect") &
    theme(legend.position = legend)
  
  return(p)
}

#' scatter.plot
#'
#' @param data 
#' @param method 
#' @param lineColor 
#' @param titleColor 
#' @param title 
#' @param xTitle 
#' @param yTitle 
#' @param legendTitle 
#' @param legend 
#'
#' @return
#' @export
#'
scatter.plot <- function(data,
                         method = NULL,
                         lineColor = NULL,
                         titleColor = NULL,
                         title = NULL,
                         xTitle = NULL,
                         yTitle = NULL,
                         legendTitle = NULL,
                         legend = "bottom") {
  if (is.null(method)) {
    method = "loess"
  } else if (!(method %in% c("lm", "loess"))) {
    stop("Select method in 'lm' or 'loess'")
  }
  if (ncol(data) == 3) {
    colnames(data) <- c("Raw", "Pre", "Group")
    if (is.null(lineColor)) lineColor = c("#e1181f", "#253870", "#09474f", "gray", "white", "black")
    p <- ggplot(data = data, mapping = aes(x = Raw, y = Pre, group = Group, col = Group)) +
      geom_smooth(formula = 'y ~ x',
                  method = method)
  } else if (ncol(data) == 2) {
    colnames(data) <- c("Raw", "Pre")
    if (is.null(lineColor)) lineColor = "#09474f"
    p <- ggplot(data = data, mapping = aes(x = Raw, y = Pre)) + 
      geom_smooth(formula = 'y ~ x',
                  color = lineColor,
                  method = method)
  }
  if (is.null(titleColor)) titleColor = "#006699"
  
  p <- p +
    geom_point() +
    theme_bw() +
    stat_cor(data = data, method = "pearson") +
    labs(title = title, x = xTitle, y = yTitle, fill = legendTitle) + 
    theme(plot.title = element_text(color = titleColor),
          legend.position = legend) +
  scale_color_manual(values = lineColor)
  
  return(p)
}

#' bar.plot
#'
#' @param data 
#' @param barColor 
#' @param titleColor 
#' @param title 
#' @param xTitle 
#' @param yTitle 
#' @param legendTitle 
#' @param legend 
#'
#' @return
#' @export
#'
bar.plot <- function(data,
                     barColor = NULL,
                     titleColor = NULL,
                     title = NULL,
                     xTitle = NULL,
                     yTitle = NULL,
                     legendTitle = NULL,
                     legend = "bottom") {
  if (ncol(data) == 3) colnames(data) <- c("id", "Variable", "Value")
  if (is.null(barColor)) barColor <- c("#09474f", "#99bac7", "gray", "white", "black")
  if (is.null(titleColor)) titleColor = "#006699"
  
  data$Variable <- factor(data$Variable,
                          levels = gtools::mixedsort(unique(data$Variable)),
                          labels = gtools::mixedsort(unique(data$Variable)))
  data$id <- factor(data$id,
                    levels = gtools::mixedsort(unique(data$id)),
                    labels = gtools::mixedsort(unique(data$id)))
  
  p <- ggplot() +
    geom_bar(data = data,
             mapping = aes(x = id, y = Value, group = Variable, fill = Variable),
             stat = "identity",
             position = "dodge", # "jitter"
             color = "black",
             width = 0.8) +
    theme_bw() +
    scale_fill_manual(values = barColor) +
    labs(title = title, x = xTitle, y = yTitle, fill = legendTitle) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.8),
          plot.title = element_text(color = titleColor),
          legend.position = legend)
  
  return(p)
}

#' box.plot
#'
#' @param data 
#' @param boxColor 
#' @param titleColor 
#' @param title 
#' @param xTitle 
#' @param yTitle 
#' @param legendTitle 
#' @param legend 
#'
#' @return
#' @export
#'
box.plot <- function(data,
                     boxColor = NULL,
                     titleColor = NULL,
                     title = NULL,
                     xTitle = NULL,
                     yTitle = NULL,
                     legendTitle = NULL,
                     legend = "bottom") {
  if (ncol(data) == 2) colnames(data) <- c("Variable", "Value")
  if (is.null(boxColor)) boxColor <- c("#09474f", "#99bac7", "gray", "white", "black")
  p <- ggplot(data = data, aes(x = Variable, y = Value)) +
    geom_boxplot(aes(fill = Variable)) +
    geom_jitter() +
    stat_compare_means() +
    theme_bw() +
    scale_fill_manual(values = boxColor) +
    labs(title = title, x = xTitle, y = yTitle, fill = legendTitle)  +
    theme(plot.title = element_text(color = titleColor),
          legend.position = legend)
  
  return(p)
}

#' network.plot
#'
#' @param network 
#' @param title 
#' @param legend 
#'
#' @return
#' @export
#'
network.plot <- function(network,
                         title = NULL,
                         legend = TRUE,
                         legendPosition = "bottomright") {
  # Make a palette of 3 colors
  # color <- brewer.pal(5, "Set1")
  color <- c("#e5244f", "#99bac7")
  pointSize = 20
  
  # Prepare network data
  colnames(network) <- c("TF", "Gene", "weight")
  network$TF <- as.vector(network$TF)
  network$Gene <- as.vector(network$Gene)
  network$weight <- as.numeric(network$weight / sum(network$weight))
  nodes <- data.frame(name = c(unique(network$Gene), network$TF),
                        carac = c(rep("Gene", 1), rep("TF", nrow(network))))
  
  networkData <- graph_from_data_frame(d = network, vertices = nodes)
  # Create a vector of color
  plot(networkData,
       vertex.color = color[as.numeric(as.factor(V(networkData)$carac))],
       vertex.label.color = "black",
       vertex.label.font = c(1),
       vertex.size = c(pointSize),
       edge.width = E(networkData)$weight * 20)
  title(main = title)
  if (legend) {
    legend(legendPosition,
           legend = levels(as.factor(V(networkData)$carac)),
           col = color, text.col = color, 
           bty = "n", pch = 20, pt.cex = 2,
           cex = 1, horiz = TRUE)
  }
}

Rcpp::sourceCpp("functions/compute_expression_vector.cpp")

#' evaluate.model
#'
#' @param rawData 
#' @param preData 
#'
#' @return
#' @export
#'
evaluate.model <- function(rawData,
                           preData) {
  evaluateResult <- list()
  
  R2 <- caret::postResample(rawData, preData)[2]
  evaluateResult[[1]] <- as.numeric(R2)
  
  # RMSE <- Metrics::rmse(rawData, preData)
  RMSE <- caret::postResample(rawData, preData)[1]
  evaluateResult[[2]] <- as.numeric(RMSE)
  
  corResult <- psych::corr.test(rawData, preData)
  evaluateResult[[3]] <- as.numeric(corResult$r)
  evaluateResult[[4]] <- as.numeric(corResult$p)
  
  names(evaluateResult) <- c("R^2", "RMSE", "Cor_R", "Cor_P")
  return(evaluateResult)
}
