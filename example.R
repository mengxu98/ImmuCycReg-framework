
# Author: Joeri Ruyssinck (joeri.ruyssinck@intec.ugent.be)

# Following packages should be installed:

# glmnet, randomForest, e1071, ROCR


###############################################################################



# Load all the script files
source("nimefi_main.R")


expression_dataset <- expression_dataset[1:100,1:100]
expression_dataset <- abs(expression_dataset[1:100,1:100])
# Read the expression dataset
expression_dataset <- constructExpressionMatrixFromFile("exampleData/expression_dataset.txt")
min(expression_dataset)
# Optional: Read a list of TF
# tf <- getIndicesOfGenesInMatrix(expression_dataset,inputFile="")

# Run algorithm using default settings, writes result to output.txt
NIMEFI(expression_dataset,GENIE=F,SVM=F,EL=TRUE,outputFileName = "output/output.txt")
evaluationObject <- prepareEval("output/output.txt","exampleData/gold.tsv")
L0_AUROC <- calcAUROC(evaluationObject)
L0_AUPR <- calcAUPR(evaluationObject)

NIMEFI(expression_dataset,GENIE=T,SVM=F,EL=F,outputFileName = "output/output.txt")
# Calculate AUROC/AUPR
evaluationObject <- prepareEval("output/output.txt","exampleData/gold.tsv")
# Print AUROC/AUPR
GENIE3_AUROC <- calcAUROC(evaluationObject)
GENIE3_AUPR <- calcAUPR(evaluationObject)



expression_dataset <- abs(expression_dataset)
# Optional: Read a list of TF
# tf <- getIndicesOfGenesInMatrix(expression_dataset,inputFile="")

# Run algorithm using default settings, writes result to output.txt
NIMEFI(expression_dataset,GENIE=F,SVM=F,EL=TRUE,outputFileName = "output/output.txt")
evaluationObject <- prepareEval("output/output.txt","exampleData/gold.tsv")
L0_AUROC_abs <- calcAUROC(evaluationObject)
L0_AUPR_abs <- calcAUPR(evaluationObject)

NIMEFI(expression_dataset,GENIE=T,SVM=F,EL=F,outputFileName = "output/output.txt")
# Calculate AUROC/AUPR
evaluationObject <- prepareEval("output/output.txt","exampleData/gold.tsv")
# Print AUROC/AUPR
GENIE3_AUROC_abs <- calcAUROC(evaluationObject)
GENIE3_AUPR_abs <- calcAUPR(evaluationObject)

