#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector compute_expression_vector(DataFrame weightDT, NumericMatrix rawMatrix) {
  // Get column names from rawMatrix
  CharacterVector rawMatrixColNames = colnames(rawMatrix);
  
  CharacterVector colnames = weightDT.names();
  colnames[0] = "regulatoryGene";
  colnames[1] = "targetGene";
  colnames[2] = "Weight";
  weightDT.attr("names") = colnames;
  
  CharacterVector regulatoryGene = weightDT["regulatoryGene"];
  NumericVector weightVector = weightDT["Weight"];
  NumericVector expressionVector(rawMatrix.nrow(), 0.0);
  
  IntegerVector sortingIndex = match(regulatoryGene, rawMatrixColNames) - 1; // Adjusting for 0-based indexing
  
  for (int i = 0; i < weightDT.nrows(); ++i) {
    NumericVector column = rawMatrix(_, sortingIndex[i]);
    double weight = weightVector[i];
    expressionVector = expressionVector + column * weight;
  }
  
  return expressionVector;
}

// R code
// #' compute.expression.vector
// #'
// #' @param weightDT 
// #' @param rawMatrix 
// #'
// #' @return
// #' @export
// #'
// compute.expression.vector <- function(weightDT,
//                                       rawMatrix) {
//   colnames(weightDT) <- c("regulatoryGene", "targetGene", "Weight")
//   weightDT$regulatoryGene <- as.vector(weightDT$regulatoryGene)
//   expressionVector <- 0
//   for (i in 1:nrow(weightDT)) {
//     gene <- weightDT$regulatoryGene[i]
//     expressionVector <- expressionVector + rawMatrix[, gene] * weightDT$Weight[i]
//   }
//   return(expressionVector)
// }