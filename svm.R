
# Author: Joeri Ruyssinck (joeri.ruyssinck@intec.ugent.be)

###############################################################################

# SVMVariableEnsembleSolve
#
# Solves the feature selection problem using an ensemble of SVM
#
#
# 
# #	Parameters	(required):
#		-- expressionMatrix: expression matrix as returned by constructExpressionMatrixFromFile.
#							 Format: samples x genes ; colnames are unique gene labels ; rownames are ignored
#		-- predictorIndices: indicates the rows of the matrix that should be considered as predictors for a regression problem.
#		-- targetIndices   : indicates the colums of the matrix that should which genes should be considered as targets for a regression problem.
#									
#	
#	
#	
#	Parameters (optional):
# 		
#		-- ensembleSize[DEFAULT=2000] : the amount of models in the ensemble
#		-- minSampleSize[DEFAULT=round((dim(expressionMatrix)[1])/5)] : mininimum sample size of the experiments
#		-- maxSampleSize[DEFAULT=4*round((dim(expressionMatrix)[1])/5)] : maximum sample size of the experiments
#		-- predictorSampleSizeMin[DEFAULT=round(length(predictorIndices)/5)] : mininimum sample size of predictors
#		-- predictorSampleSizeMax[DEFAULT=4*round(length(predictorIndices)/5)] : maximum sample size of the predictors
#		-- rankThreshold[DEFAULT=round((length(predictorIndices)/20))] : the amount of top ranked features that should be awarded 1 instead of 0 during a feature ranking in an iteration step
#		-- trace[DEFAULT=TRUE]: index of target currently computed is reported to stdout
#       -- traceDetail[DEFAULT=FALSE] : index of subproblen is reported to stdout


#		-- ... : the other parameters are passed to the SVM method (svm)
#
#
#	Returns:								
#		
#		A matrix where each cell(row i, column j) specifies the importance of variable i to target j.
#
#
#
#
#
#

SVMVariableEnsembleSolve <- function(expressionMatrix,predictorIndices=NULL,targetIndices=NULL,ensembleSize=2000,
		minSampleSize=round((dim(expressionMatrix)[1])/5),maxSampleSize=4*round((dim(expressionMatrix)[1])/5),
		predictorSampleSizeMin=round(length(predictorIndices)/5),predictorSampleSizeMax=4*round(length(predictorIndices)/5),
		rankThreshold=round((length(predictorIndices)/20)),trace=TRUE,traceDetail=FALSE,...){
	
	
	
	# Check input
	if( (minSampleSize<=0 || maxSampleSize > dim(expressionMatrix)[1] || maxSampleSize < minSampleSize  )){
		stop("Please specify a valid sample sample size minimum and maximum.")
	}
	
	# Preallocate the return matrix
	resultMatrix <- matrix(0.0,length(predictorIndices),length(targetIndices))
	rownames(resultMatrix) <- colnames(expressionMatrix)[predictorIndices]
	colnames(resultMatrix) <- colnames(expressionMatrix)[targetIndices]
	
	
	# Take a sample this many times
	for (i in 1:ensembleSize){
		if(minSampleSize==maxSampleSize){
			sampleSize <- minSampleSize
		}else{
			sampleSize <- sample(minSampleSize:maxSampleSize,1)
		}
		
		# report on progress
		if(trace){
				print(paste("Sample: ",i, "/",ensembleSize ," of size: ",sampleSize ," with a rank threshold of ", rankThreshold))
		}
		
		# Take sample of specified size 
		ind <- sample(dim(expressionMatrix)[1],sampleSize)
		sampleMatrix <- expressionMatrix[ind,,drop=FALSE]
		
		# Call svm
		resultMatrix <-	resultMatrix + SVMRankedSolve(sampleMatrix,predictorIndices=predictorIndices,targetIndices=targetIndices,traceDetail=traceDetail,rankThreshold=rankThreshold,predictorSampleSizeMin = predictorSampleSizeMin, predictorSampleSizeMax = predictorSampleSizeMax,...)

	}
	
	
	    # return
		return(resultMatrix)

	
}




# 	SVMRankedSolve
#
#	
# 	Solves the feature selection problem using ranked svm.
#
#
#	Parameters	(required):
#		-- expressionMatrix: expression matrix as returned by constructExpressionMatrixFromFile.
#							 Format: samples x genes ; colnames are unique gene labels ; rownames are ignored
#		-- predictorIndices: indicates the rows of the matrix that should be considered as predictors for a regression problem.
#										 It is highly recommened to use the output of getPredictorIndices
#										 If null, getIndicesOfGenesInMatrixs will be called on the matrix and all rows will be considered predictors if needed for the algorithm.
#		-- targetIndices [DEFAULT=NULL]: indicates the colums of the matrix that should which genes should be considered as targets for a regression problem.
#										 It is highly recommened to use the output of getIndicesOfGenesInMatrixs or groupRegressionProblemIndices.
#										 If null, getIndidicesForRegressionProblems will be called on the matrix and all columns will be considered targets if needed for the algorithm.
#		-- rankThreshold : the amount of top ranked features that should be awarded 1 instead of 0 during a feature ranking in an iteration step
#		-- predictorSampleSizeMin : mininimum sample size of predictors
#		-- predictorSampleSizeMax : maximum sample size of the predictors
#		-- traceDetail: index of target currently computed is reported to stdout			
#	
#	Parameters (optional):
# 		
#		-- ... : the other parameters are passed to the SVM method (svm)
#
#
#	Returns:								
#		
#		A matrix where each cell(row i, column j) specifies the importance of variable i to target j.
#									
SVMRankedSolve <- function(expressionMatrix,predictorIndices,targetIndices, rankThreshold,traceDetail,predictorSampleSizeMin,predictorSampleSizeMax, ...){
	
	# Check input
	if( (predictorSampleSizeMin<=0 || predictorSampleSizeMax >= length(predictorIndices) || predictorSampleSizeMin > predictorSampleSizeMax  )){
		stop("Please specify a valid predictor sample size minimum and maximum (between 1 inclusive and the length of the predictorIndices exclusive).")
	}
	# lib should be set to a correct valulie, also load the library
	library("e1071")
	# End parameter check and loading libraries
	
	# Execute algorithm
	# Extract the predictors
	predictors <- expressionMatrix[,predictorIndices]
	# Preallocate the return matrix
	resultMatrix <- matrix(0.0,length(colnames(predictors)),length(targetIndices))
	rownames(resultMatrix) <- colnames(expressionMatrix)[predictorIndices]
	colnames(resultMatrix) <- colnames(expressionMatrix)[targetIndices]
	
	
	i <- 0
	for ( targetIndex in targetIndices){
		# Report on progress if requested
		if(traceDetail){ 
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
		
		# Sample on possible predictors
		if(predictorSampleSizeMin==predictorSampleSizeMax){
			predictorSampleSize <- predictorSampleSizeMin
			pIndices <- sample(1:(dim(predictorsWithoutTarget)[2]),predictorSampleSize)	
			predictorsWithoutTarget <- predictorsWithoutTarget[,pIndices,drop=FALSE]
		}
		else{
			predictorSampleSize <- sample(predictorSampleSizeMin:predictorSampleSizeMax,1)
			pIndices <- sample(1:(dim(predictorsWithoutTarget)[2]),predictorSampleSize)	
			predictorsWithoutTarget <- predictorsWithoutTarget[,pIndices,drop=FALSE]
		}	
		
		
		# The target considered in this iteration
		target 	 <- expressionMatrix[,targetIndex]
		
		# Call
		sv <- svm(predictorsWithoutTarget,target,kernel="linear")

		# Calc weights
		wghts <- t(sv$coefs) %*% sv$SV

		
		# Now sort the wghts
		wghts <- abs(wghts)
		indices <- sort.list(wghts,decreasing=TRUE)
		# Check for zero entries
		zeros <- which(wghts==0)
		
		# Now replace by ones that are in the top and are non-zero
		wghts[1:length(wghts)] <- 0 
		wghts[indices[1:rankThreshold]] <- 1
		
		
		# Set the ones that were zero to zero anyway
		wghts[zeros] <- 0
		
		# write result
		resultMatrix[colnames(predictorsWithoutTarget),colnames(expressionMatrix)[targetIndex]] <- wghts
		
	
	}

	return (resultMatrix)
}
