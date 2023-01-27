
# Author: Joeri Ruyssinck (joeri.ruyssinck@intec.ugent.be)

###############################################################################

# 	constructExpressionMatrixFromFile
#
#	Reads and constructs a normalized expression matrix (samples are rows x genes are cols) from different file formats.
#
#	Parameters	(required):
#		-- filename: name or full path to the expression matrix file

#	Parameters (optional):
#		-- rowFormat: format of the file. Supported:
#													- A row is a series of measurements of one gene    (TRUE) 
#													- A column is a series of measurements of one gene (FALSE)
#		-- sep [DEFAULT=\t]: seperator of the values, default tab
#		-- isLabeled [DEFAULT=TRUE]: first row or column of file are gene labels.
#		-- labels [DEFAULT=NULL] : unique labels for each gene
#									- labels!=NULL : use labels from input argument	
#									- labels==NULL :
#													- isLabeled==FALSE: use default generated labels
#													- isLabeled==TRUE : use labels extracted from the data in the file
#		-- sampleNames [DEFAULT=FALSE]: indicates whether the first row/col is the name of the sample. In this case, the row/col will be dropped
#		-- normalize [DEFAULT=TRUE] : indicates whether the expression matrix should be normalized
#	Returns:								
#		(normalized) expression matrix in format (samples x genes), i.e. every col is a gene
#									
#
constructExpressionMatrixFromFile <- function(filename, rowFormat=FALSE, seperator="\t", isLabeled=TRUE, labels=NULL,sampleNames=FALSE,normalize=TRUE) {
	# Default gene prefix
	DEFAULT.GENE.LABEL = "GENE_"
	# Read data
	m <- read.table(filename,sep=seperator,header=FALSE,as.is=TRUE,colClasses="character")
	# Set the labels of the genes if they are not provided
	if(is.null(labels)) {
		if(isLabeled){ # Labels should be extracted from the file
			if (rowFormat) {labels <- m[,1]}
			else		   {labels <- m[1,]}
		}
		else{ # Labels should be generated 
			labels <- paste(DEFAULT.GENE.LABEL,seq(from=1,to=dim(m)[1]),sep="")
		}
	}
	# Check the labels, quit on not unique
	if(any(duplicated(colnames(m)))){stop("Labels for genes are not unique") }
	# Remove gene and/or sample labels
	if (rowFormat) {
		if(sampleNames) { m <- m[-1,] }
		if(isLabeled) 	{ m <- m[,-1] }
	}
	else{	  
		if(sampleNames) { m <- m[,-1] }
		if(isLabeled) 	{ m <- m[-1,] }
	}
	# Coerce to numeric
	m <- as.matrix(m)
	m <- apply(m,2, as.numeric)
	# Transpose to samples x genes format if needed
	if (rowFormat){ m <- t(m)}
	# Add gene labels
	colnames(m) <- labels
	# normalize if requested
	if (normalize){ m <- scale(m) }
	# return expression matrix
	return(m)
}





# 	getIndicesOfGenesInMatrix
#
#	Returns and checks the column indices to be used as targets or predictors.
#	If not all genes specified are present in the matrix. This method stops executing.
#
#	Parameters	(required):
#		-- expressionMatrix: expression matrix as returned by constructExpressionMatrixFromFile.
#							 Format: samples x genes ; colnames are unique gene labels ; rownames are ignored
#
#	Parameters (optional):
#		-- genes[DEFAULT=NULL]: vector indicating the genes whose column indices need to be returned. 
#								 Parameters genes and inputFile can not be both non-null.
#										- format: 1- A vector of indices
#												  2- A vector of gene names (strings)
#		--inputFile[DEFAULT=NULL]: path to a file containing a tab or new line seperated file with the genes specified as in inputfileFormat.
# 									parameters genes and inputFile cannot be both non-null.
#		--inputFileFormat[DEFAULT=character()]: indicates whether the genes are specified as strings or indices in the file.
#	Returns:								
#		sorted indices list of columns
#									
#
getIndicesOfGenesInMatrix <- function(expressionMatrix,genes=NULL,inputFile=NULL,seperator="\t",inputFileFormat=character()){	
	
	# Check for conflicting arguments
	if ((!is.null(genes)) && (!is.null(inputFile))){ stop("Aborting. Parameters \"genes\" and \"inputFile\" cannot be both non-null")}
	# Read the genes from file if requested
	if (!is.null(inputFile)) { genes <- (scan(inputFile,character(),sep=seperator,what=inputFileFormat))}
	# Check if genes are provided, if not all genes are considered
	if (is.null(genes)) {
		genes= 1:length(colnames(expressionMatrix))
	}
	else{ # genes are provided
		# Check if indices or gene names are given
		if(!isPositiveNonZeroIntegerVector(genes)){ # Do conversions if genes are specified by gene name
			# Check first if all geneTargets are in the expressionMatrix
			tmp <- setdiff(genes,colnames(expressionMatrix)) 
			if (length(tmp) > 0){
				for (i in tmp){
					warning(paste("Gene ", i, " not in expression matrix\n",sep=""))
				}
				warning("Certain genes not present in expression matrix") 	
			}
			# If correct, convert into indices array
			tmp <- mat.or.vec(length(genes),1)
			for (i in 1:length(genes)){
				tmp[i] = match(genes[i],colnames(expressionMatrix))
			}
			
			
			tmp <- tmp[!is.na(tmp)]
			tmp <- tmp[sort.list(as.numeric(tmp))]
		}
		else{ # already an indices array, just sort
			genes <- sort(geneTargets)
			# Check for invalid indexes
			if(genes[length(genes)] > length(colnames(expressionMatrix))){
				warning("Target index exceeds amount of genes/columns.")
			}
		}
		# Remove doubles from the list
		genes <- unique(genes)
	}
	return (genes) 
}
