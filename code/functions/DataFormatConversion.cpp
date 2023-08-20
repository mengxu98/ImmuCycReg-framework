#include <Rcpp.h>
using namespace Rcpp;

// Function to convert FPKM matrix to TPM matrix
// [[Rcpp::export]]
NumericMatrix fpkmMatrixToTpmMatrixCpp(NumericMatrix fpkm) {
  int n_genes = fpkm.nrow();
  int n_samples = fpkm.ncol();
  NumericMatrix tpm(n_genes, n_samples);

  for (int j = 0; j < n_samples; ++j) {
    double sum_fpkm = sum(fpkm(_, j)); // Calculate the sum of FPKM for each sample
    for (int i = 0; i < n_genes; ++i) {
      tpm(i, j) = exp(log(fpkm(i, j)) - log(sum_fpkm) + log(1e6)); // Calculate TPM
    }
  }

  return tpm;
}

// Function to convert counts vector to TPM vector
// [[Rcpp::export]]
NumericVector countToTpmCpp(NumericVector counts, NumericVector effLen) {
  int n = counts.size();
  NumericVector rate = log(counts) - log(effLen);
  double denom = log(sum(exp(rate)));
  
  NumericVector tpmVector(n);
  for (int i = 0; i < n; ++i) {
    tpmVector[i] = exp(rate[i] - denom + log(1e6));
  }
  
  return tpmVector;
}

// Function to convert counts matrix to TPM matrix
// [[Rcpp::export]]
NumericMatrix countMatrixToTpmMatrixCpp(NumericMatrix counts, NumericVector effLen) {
  int nrow = counts.nrow();
  int ncol = counts.ncol();
  NumericMatrix tpmMatrix(nrow, ncol);

  for (int j = 0; j < ncol; ++j) {
    // NumericVector rate = log(counts(_, j)) - log(effLen(j));
    NumericVector rate = log(counts(_, j)) - log(effLen);
    double denom = log(sum(exp(rate)));
    for (int i = 0; i < nrow; ++i) {
      tpmMatrix(i, j) = exp(rate[i] - denom + log(1e6));
    }
  }
  
  return tpmMatrix;
}


// Function to convert counts to FPKM
// [[Rcpp::export]]
NumericVector countToFpkmCpp(NumericVector counts, NumericVector effLen) {
  int n = counts.size();
  double N = sum(counts);
  
  NumericVector fpkm(n);
  for (int i = 0; i < n; ++i) {
    fpkm[i] = exp(log(counts[i]) + log(1e9) - log(effLen[i]) - log(N));
  }
  
  return fpkm;
}

// Function to convert counts matrix to FPKM matrix
// [[Rcpp::export]]
NumericMatrix countMatrixToFpkmMatrixCpp(NumericMatrix counts, NumericVector effLen) {
  int nrow = counts.nrow();
  int ncol = counts.ncol();
  NumericMatrix fpkmMatrix(nrow, ncol);
  
  for (int j = 0; j < ncol; ++j) {
    double N = sum(counts(_, j));
    NumericVector count = counts(_, j);
    
    for (int i = 0; i < nrow; ++i) {
      fpkmMatrix(i, j) = exp(log(count[i]) + log(1e9) - log(effLen[i]) - log(N));
    }
    
  }
  
  return fpkmMatrix;
}


// Function to convert counts matrix to effective counts matrix
// [[Rcpp::export]]
NumericMatrix countMatrixToEffCountsMatrixCpp(NumericMatrix counts, NumericVector len, NumericVector effLen) {
  int nrow = counts.nrow();
  int ncol = counts.ncol();
  NumericMatrix effCountsMatrix(nrow, ncol);

  for (int j = 0; j < ncol; ++j) {
    for (int i = 0; i < nrow; ++i) {
      effCountsMatrix(i, j) = counts(i, j) * (len(i) / effLen(j));
    }
  }

  return effCountsMatrix;
}

// Function to convert FPKM to TPM
// [[Rcpp::export]]
NumericVector fpkmToTpmCpp(NumericVector fpkm) {
  int n = fpkm.size();
  double total = 0.0;

  for (int i = 0; i < n; ++i) {
    total += fpkm[i];
  }

  NumericVector tpm(n);
  for (int i = 0; i < n; ++i) {
    tpm[i] = exp(log(fpkm[i]) - log(total) + log(1e6));
  }

  return tpm;
}


// Function to convert counts to effective counts
// [[Rcpp::export]]
NumericVector countToEffCountsCpp(NumericVector counts, NumericVector len, NumericVector effLen) {
  int n = counts.size();

  NumericVector effCounts(n);
  for (int i = 0; i < n; ++i) {
    effCounts[i] = counts[i] * (len[i] / effLen[i]);
  }

  return effCounts;
}


// Reference: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
// R code and examples
// countToTpm <- function(counts,
//                        effLen) {
//   rate <- log(counts) - log(effLen)
//   denom <- log(sum(exp(rate)))
//   exp(rate - denom + log(1e6))
// }
// 
// countToFpkm <- function(counts, effLen) {
//   N <- sum(counts)
//   exp( log(counts) + log(1e9) - log(effLen) - log(N) )
// }
// 
// fpkmToTpm <- function(fpkm) {
//   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
// }
// 
// countToEffCounts <- function(counts, len, effLen) {
//   counts * (len / effLen)
// }
// 
//
// # Examples
// 
// Rcpp::sourceCpp("DataFormatConversion.cpp")
// 
// # Testing for vector
// 
// cnts <- c(4250, 3300, 200, 1750, 50, 0)
// lens <- c(900, 1020, 2000, 770, 3000, 1777)
// 
// # Assume a mean(FLD) = 203.7
// meanFLD = 203.7
// countDf <- data.frame(count = cnts, length = lens)
// 
// countDf$effLength <- countDf$length - meanFLD + 1
// countDf$tpm <- with(countDf, countToTpm(count, effLength))
// countDf$fpkm <- with(countDf, countToFpkm(count, effLength))
// with(countDf, all.equal(tpm, fpkmToTpm(fpkm)))
// countDf$effCounts <- with(countDf, countToEffCounts(count, length, effLength))
// 
// countDfCpp <- data.frame(count = cnts, length = lens)
// countDfCpp$effLength <- countDfCpp$length - meanFLD + 1
// countDfCpp$tpm <- with(countDfCpp, countToTpmCpp(count, effLength))
// countDfCpp$fpkm <- with(countDfCpp, countToFpkmCpp(count, effLength))
// with(countDfCpp, all.equal(tpm, fpkmToTpmCpp(fpkm)))
// countDfCpp$effCounts <- with(countDfCpp, countToEffCountsCpp(count, length, effLength))
// 
// 
// 
// # Testing for matrix
// 
// lens <- c(900, 1020, 2000)
// countMatrix <- matrix(c(4250, 3300, 200, 1750, 20, 175, 510, 850, 10), nrow = 3)
// 
// fpkmMatrix <- apply(countMatrix, 2, countToFpkm, lens)
// tpmMatrix <- apply(fpkmMatrix, 2, fpkmToTpm)
// tpmMatrix2 <- apply(countMatrix, 2, countToTpm, lens)
// 
// fpkmMatrixCpp <- apply(countMatrix, 2, countToFpkmCpp, lens)
// tpmMatrixCpp <- apply(fpkmMatrix, 2, fpkmToTpmCpp)
// tpmMatrixCpp2 <- apply(countMatrix, 2, countToTpmCpp, lens)
// 
// fpkmMatrixCppM <- countMatrixToFpkmMatrixCpp(countMatrix, lens) # W
// tpmMatrixCppM <- fpkmMatrixToTpmMatrixCpp(fpkmMatrix)
// tpmMatrixCppM2 <- countMatrixToTpmMatrixCpp(countMatrix, lens)
// 
// for (i in 1:ncol(fpkmMatrix)) {
//   print(all.equal(fpkmMatrix[, i], fpkmMatrixCpp[, i]))
// }
// for (i in 1:ncol(fpkmMatrix)) {
//   print(all.equal(fpkmMatrix[, i], fpkmMatrixCppM[, i]))
// }
// 
// for (i in 1:ncol(fpkmMatrix)) {
//   print(all.equal(tpmMatrix[, i], tpmMatrixCpp[, i]))
// }
// for (i in 1:ncol(fpkmMatrix)) {
//   print(all.equal(tpmMatrix[, i], tpmMatrixCpp2[, i]))
// }
// for (i in 1:ncol(fpkmMatrix)) {
//   print(all.equal(tpmMatrix[, i], tpmMatrixCppM[, i]))
// }
// for (i in 1:ncol(fpkmMatrix)) {
//   print(all.equal(tpmMatrix[, i], tpmMatrixCppM2[, i]))
// }
// 
// 
// # Testing for tcga matrix
// pathRead <- "../data/"
// pathSave <- "../../Results/"
// 
// tcga_luad <- read.table(paste0(pathRead, paste0("luad-rsem-fpkm-tcga-t.txt.gz")),
//                         header = TRUE,
//                         row.names = 1,
//                         sep = "\t",
//                         check.names = FALSE) %>% .[, -1]
// colnames(tcga_luad) <- substr(colnames(tcga_luad), 1, 15)
// 
// tpmMatrixTcga <- apply(tcga_luad, 2, fpkmToTpm) %>% as.data.frame()
// tpmMatrixTcgaCpp <- apply(tcga_luad, 2, fpkmToTpmCpp) %>% as.data.frame()
// tcga_luad_m <- as.matrix(tcga_luad)
// tpmMatrixTcgaCppM <- fpkmMatrixToTpmMatrixCpp(tcga_luad_m) %>% as.data.frame()
// 
// for (i in 1:ncol(tpmMatrixTcga)) {
//   print(all.equal(tpmMatrixTcga[, i], tpmMatrixTcgaCpp[, i]))
// }
// for (i in 1:ncol(tpmMatrixTcga)) {
//   print(all.equal(tpmMatrixTcga[, i], tpmMatrixTcgaCppM[, i]))
// }
// 
// colSums(tpmMatrixTcga)
// colSums(tpmMatrixTcgaCpp)
// colSums(tpmMatrixTcgaCppM)
