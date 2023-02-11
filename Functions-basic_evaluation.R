

library("caTools")
library("ROCR")

# Defines class
#' evalClass
#'
#' @slot tpr vector.
#' @slot fpr vector.
#' @slot prec vector.
#' @slot rec vector.
#' @slot p integer.
#' @slot n integer.
#' @slot tpk vector.
#' @slot fpk vector.
#' @slot predictors vector.
#' @slot targets vector.
#' @slot rank vector.
#'
#' @return
#' @export
#'
#' @examples
setClass(
    "evalClass",
    representation(
        tpr = "vector",
        fpr = "vector",
        prec = "vector",
        rec = "vector",
        p = "integer",
        n = "integer",
        tpk = "vector",
        fpk = "vector",
        predictors = "vector",
        targets = "vector",
        rank = "vector"
    )
)

#' prepareEval
#'  Creates an evalClass object from a prediction file in TSV format and a given gold standard in tsv format
#' @param predictionTSV prediction file in tsv format (format: TF[TAB]TARGET_GENE[TAB]RANK -  with most confident gene at top)
#' @param goldTSV gold standard in tsv format ((format: TF[TAB]TARGET_GENE[TAB]1/0] - 1 for true link, 0 for other links
#' @param totalPredictionsAccepted [DEFAULT=100000] maximum amount of predictions allowed
#'
#' @return an EvalClass object, which can be queried for AUROC and AUPR scores
#' @export
#'
#' @examples
prepareEval <- function(predictionTSV,
                        goldTSV,
                        totalPredictionsAccepted = 100000) {
    pred <- read.table(predictionTSV)
    gold <- convertSortedRankTSVToAdjMatrix(goldTSV)

    # Check arguments
    if ((!isWholeNumber(totalPredictionsAccepted)) ||
        totalPredictionsAccepted < 1) {
        stop("totalPredictionsAccepted should be an integer <0.")
    }

    # Next reduce the predictionTSV to the max allowed predictions if necessary
    if (length(pred[, 1]) > totalPredictionsAccepted) {
        pred <- pred[1:totalPredictionsAccepted, ]
    }

    # Get some metrics from the gold standard
    # The amount of positive links
    p <- sum(gold)
    # The amount of negative links (these could be valid links and are not positive) [self regulating links are not allowed and are consired to be neither positive or negative]
    n <- (dim(gold)[1] * dim(gold)[2]) - p - sum(rownames(gold) %in% colnames(gold))
    # The total amount of valid links (positive +negative)
    t <- p + n

    # Now, remove all the links that are in the prediction file that are not valid predictions (meaning self-regulation or genes not eligble to be predictor or target)
    # 1) Remove all the edges that are in the prediction file but are not recorded in the gold file as an edge or non-edge
    # 2) For those that are in thet network a) Replace value with a one if edge is present in gold or zero otherwise
    # Remove all the non-valid predictors
    pred <- pred[which(as.vector(pred[, 1]) != as.vector(pred[, 2])), , drop = FALSE]

    # Declare some variables for clarity
    firstRow <- as.vector(pred[, 1])
    secondRow <- as.vector(pred[, 2])

    # Will indicate whether this prediction was present in gold standard
    thirdRow <- vector(mode = "integer", length(firstRow))

    # Will indicate the rank of this prediction = (-1 if it is not present and [rank] otherwise]
    rank <- vector(mode = "integer", length(firstRow))

    # Now, loop over all predictions and determine if they are present in gold matri
    correct <- 0
    incorrect <- 0

    for (i in 1:length(firstRow)) {
        if (gold[firstRow[i], secondRow[i]] == 1) {
            correct <- correct + 1
            thirdRow[i] <- 1
            rank[i] <- correct
        } else {
            incorrect <- incorrect - 1
            thirdRow[i] <- 0
            rank[i] <- incorrect
        }
    }

    # Check how many of the gold standard edges we predicted. Any other that still remain are discovered at a uniform rate by definition. (remaining_gold_edges/remaining_edges)
    # Gold links predicted so far
    # If some gold links have not been predicted, calculate the random discovery chance, else set to zero
    if (length(firstRow) < t) {
        odds <- (p - correct) / (t - length(firstRow))
    } else {
        odds <- 0
    }

    # Each guess you have 'odds' chance of getting one right and '1-odds' chance of getting it wrong , now construct a vector till the end
    random_positive <- vector("double", (t - length(firstRow)))
    random_negative <- vector("double", (t - length(firstRow)))
    random_positive[] <- odds
    random_negative[] <- 1 - odds

    # Calculate the amount of true positives and false positives at 'k' guesses depth
    positive <- c(thirdRow, random_positive)
    negative <- c(abs(thirdRow - 1), random_negative)
    tpk <- cumsum(positive)
    fpk <- cumsum(negative)

    # Depth k
    k <- 1:t

    # Calculate true positive rate, false positive rate, precision and recall
    tpr <- tpk / p
    fpr <- fpk / n
    rec <- tpr
    prec <- tpk / k

    predictors <- firstRow
    targets <- secondRow

    # Create and return the object
    return(new("evalClass",
        tpr = tpr,
        fpr = fpr,
        rec = rec,
        prec = prec,
        p = as.integer(p),
        n = as.integer(n),
        tpk = tpk,
        fpk = fpk,
        predictors = predictors,
        targets = targets,
        rank = rank
    ))
}

#' calcAUROC
#'  Returns the Area under the ROC curve given an evalClass object
#' @param prediction An evalClass object
#'
#' @return AUROC score
#' @export
#'
#' @examples
calcAUROC <- function(prediction) {
    return(trapz(prediction@fpr, prediction@tpr))
}

#' calcAUPR
#'  Returns the Area under the Precision Recall curve given an evalClass object
#' @param prediction An evalClass object
#'
#' @return AUPR score
#' @export
#'
#' @examples
calcAUPR <- function(prediction) {
    return(trapz(prediction@rec, prediction@prec) / (1 - 1 / prediction@p))
}
