

# Author: Joeri Ruyssinck (joeri.ruyssinck@intec.ugent.be)

###############################################################################
# 	genieSolve
#
#	
# 	Solves the feature selection problem using GENIE3
#
#	Parameters	(required):
#		-- expressionMatrix: expression matrix as returned by constructExpressionMatrixFromFile.
#							 Format: samples x genes ; colnames are unique gene labels ; rownames are ignored
#		-- predictorIndices: indicates the rows of the matrix that should be considered as predictors for a regression problem.
#		-- targetIndices   : indicates the colums of the matrix that should be considered as targets for a regression problem.
#									
#	
#	Parameters (optional):
# 
#		-- candidateSplitCount[DEFAULT="SQRT"]: number of variables randomly selected as candidates at each split
#									format: - "SQRT": square root of total number of input genes
#											- "ALL": all genes
#											- an integer [1,all]
#		-- treesInEnsembleCount[DEFAULT=1000]: amount of trees in ensemble for each target gene 
#
#		-- importanceMeasure[DEFAULT="IncNodePurity"]: 
#									format: - "IncNodePurity": for importance measure based on decrease of residual sum of squares
#											- "IncMSE": for importance measure obtained by permutation of OOB data
#		-- randomSeed[DEFAULT=NULL]: seed off the random generator, default is not reset
#
#		-- trace[DEFAULT=TRUE]: index of target currently computed is reported to stdout
#		
#		-- ... : the other parameters are passed to the randomForest function.
#
#
#	Returns:								
#		
#		A matrix where each cell(row i, column j) specifies the importance of variable i to target j.
#									
#

genieSolve <- function(expressionMatrix,
                       predictorIndices,
                       targetIndices,
                       candidateSplitCount="SQRT",
                       treesInEnsembleCount=1000,
                       importanceMeasure="IncNodePurity",
                       randomSeed=NULL,
                       trace=TRUE,
                       ...){
	
	# Load libary
	library("randomForest")
	
	# Check the  the optional arguments, the required parameters have been checked in previous methods and should be correct
	if(!(candidateSplitCount=="SQRT" || candidateSplitCount=="ALL" || (isWholeNumber(candidateSplitCount) && candidateSplitCount > 0 && candidateSplitCount < length(predictorIndices)))){
		stop("Stopping the computation. \"candidateSplitCount\" should be either \"SQRT\", \"ALL\" or a positive non zero integer smaller or equal than the amount of predictor variables.")
	}
	if(!(importanceMeasure=="IncNodePurity" || importanceMeasure=="incMSE")){
		stop("Stopping the computation. \"importanceMeasure\" should be either \"IncNodePurity\" or \"incMSE\" ")
	}
	if(!(is.null(randomSeed) || isWholeNumber(randomSeed))){
		stop("Stopping the computation. \"randomSeed\" should be an integer")
	}
	if(!(isWholeNumber(treesInEnsembleCount) && treesInEnsembleCount >0)){
		stop("Stopping the computation. \"treesInEnsembleCount\" should be a non zero positive integer")
	}
	# End check
	
	# Start the algorithm
	# Set the random seed if supplied
	if(!is.null(randomSeed)){ 
		set.seed(randomSeed) 
		print("Setting the random seed")
	}
	# Parse the candidateSplitCount argument if it is specified as a string
	if(candidateSplitCount=="SQRT"){
		candidateSplitCount <- round(sqrt(length(predictorIndices)-1))
	}
	else if(candidateSplitCount=="ALL"){
		candidateSplitCount <- length(predictorIndices)-1
	}
	
	# Extract the predictors
	predictors <- expressionMatrix[,predictorIndices]
	
	# Preallocate the return matrix
	resultMatrix <- matrix(0.0,length(colnames(predictors)),length(targetIndices))
	rownames(resultMatrix) <- colnames(expressionMatrix)[predictorIndices]
	colnames(resultMatrix) <- colnames(expressionMatrix)[targetIndices]
	
	# Main loop over all the targets
	cat(paste("Starting main loop with ",treesInEnsembleCount, "trees/target \nand",candidateSplitCount,"candidate input genes/tree node \n\n"))
	flush.console()
	i <- 0
	for ( targetIndex in targetIndices){
		# Report on progress if requested
		if(trace){ 
			i <- i + 1
			cat(paste("Computing target ",targetIndex," iteration ", i , "\\", length(targetIndices),"\n"))
			flush.console()
		}
		# Check if the target has already been removed from the predictors, if not, remove it
		if(targetIndex %in% predictorIndices){
			nameCol <- colnames(expressionMatrix)[targetIndex]
			indexInPredictorsToRemove <- which(nameCol == colnames(predictors))
			predictorsWithoutTarget <- predictors[,-indexInPredictorsToRemove]
		}
		else{
			predictorsWithoutTarget <- predictors
		}
		# The target considered in this iteration
		target 	 <- expressionMatrix[,targetIndex]
		# Call random forest, pass on any other parameters if needed
		randomFResult <- randomForest(predictorsWithoutTarget,target,mtry=candidateSplitCount,ntree=treesInEnsembleCount,importance=TRUE,...)
		# Extract the feature importance results
		im <- importance(randomFResult)[,importanceMeasure]
		# Paste into the return matrix
		resultMatrix[names(im),colnames(expressionMatrix)[targetIndex]] <- im
	}
	# Divide by number of samples
	return (resultMatrix/ dim(expressionMatrix)[1])
	
}

