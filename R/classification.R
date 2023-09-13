# Ensure needed methods are loaded first
# Create Classifier -------------------------------------------------------
#' Create a Classifier using one of the Rank Entropy Algorithms
#'
#' @param expression Gene expression data, can be either a SummarizedExperiment
#'      object, or an array type object which implements [], and also
#'      rownames/colnames if phenotypes or geneNetworks use character vector
#'      indices
#' @param phenotype1,phenotype2 Index vectors for the locations within
#'      expression of the samples within each phenotype, can be a numeric,
#'      character, or logical vector
#' @param geneNetwork Index vector for the gene network of interest,
#'      can be a logical, numeric or character vector
#' @param method Method to use to calculate the rank entropy, can be DIRAC,
#'      RACE, CRANE, or INFER
#' @param phenotype1Name,phenotype2Name Names of the phenotypes, used for
#'      creating default name for new colData column in the SummarizedExperiment
#'      input (if a SummarizedExperiment is used as input)
#' @param geneNetworkName Name of the gene network, used for
#'      creating default name for new colData column in the SummarizedExperiment
#'      input (if a SummarizedExperiment is used as input)
#' @param assayName Name of the assay containing gene expression data in the
#'      SummarizedExperiment object (if a SummarizedExperiment is used as input)
#'
#' @return A classifier function
#'      If the input is not a SummarizedExperiment object:
#'          The classification function will take a gene expression array,
#'          with rows representing different genes, and columns representing
#'          different samples. The function will return a logical vector
#'          the same length as the number of columns in the input array,
#'          with TRUE meaning classified as phenotype1, and FALSE meaning
#'          classified as phenotype2
#'      If the input is a SummarizedExperiment object:
#'          The function will take a SummarizedExperiment as input,
#'          along with optional parameters for assayName and ColDataName.
#'          The assayName argument is used to determine which assay to use
#'          as input for the classifier (defaults to the same as assayName
#'          argument to this function). The colDataName is used to name the
#'          output column in the returned SummarizedExperiment object.
#'          The return will be a SummarizedExperiment object with all the same
#'          data as the input, except with an added column in the colData
#'          DataFrame, named according to colDataName (or
#'          geneNetworkName_phenotype1Name if no argument is provided) which is
#'          logical, with TRUE representing that sample being classified as
#'          phenotype1, and FALSE representing that sample being classified as
#'          phenotype2
#' @export
#'
#' @examples
#' @include utils.R crane.R race.R dirac.R infer.R
createClassifier <- function(expression, phenotype1, phenotype2, geneNetwork,
                             method=c("DIRAC","CRANE","RACE"),
                             phenotypeNames=c("p1","p2"),
                             geneNetworkName="gn",
                             assayName="counts"){
    method <- match.args(method)
    classifierFactoryList <- list(DIRAC=.diracEntropyClassifier,
                                  CRANE=.craneEntropyClassifier,
                                  RACE=.raceEntropyClassifier)
    classifierFactory <- classifierFactoryList[[method]]
    if(is(expression, "SummarizedExperiment")){
        p1Index <- .checkPhenotype(seObject=expression, phenotype=phenotype1)
        p2Index <- .checkPhenotype(seObject=expression, phenotype=phenotype2)
        numRows <- SummarizedExperiment::dims(expression)[1]
        rowNames <- SummarizedExperiment::rownames(expression)
        geneNetwork <- .checkNetwork(network=geneNetwork,
                                   numRows=numRows,
                                   rowNames=rowNames)
        classifier <- classifierFactory(
            expression=SummarizedExperiment::assays(expression)[[assayName]],
            phenotype1=p1Index,
            phenotype2=p2Index,
            geneIndex=geneNetwork)
        return(.wrapSeClassifier(classifier = classifier,
                                 colDataNameDefault =
                                     paste(geneNetworkName,
                                           phenotypeNames[[1]],
                                           sep="_"),
                                 assayDefault = assayName))

    }
    p1Index <- .convertPhenotype(phenotype = phenotype1,
                                 expressionMatrix = expression)
    p2Index <- .convertPhenotype(phenotype = phenotype2,
                                 expressionMatrix = expression)
    geneNetwork <- .convertNetwork(geneNetwork=geneNetwork,
                                   expressionMatrix=expression)
    classifier <- classifierFactory(expression=expression,
                                    phenotype1=p1Index,
                                    phenotype2=p2Index,
                                    geneIndex=geneNetwork)
    classifier
}

# Phenotype Classifier ---------------------------------------------------------

#' Creates a phenotype classifier based on DIRAC entropy
#'
#' @param expression Numeric matrix of gene expression values, rows are genes
#'    columns are samples
#' @param phenotype1,phenotype2 Numeric vectors with indices for each of the
#'    phenotypes
#' @param geneIndex Numeric vector of indices for the genes within the gene
#'    network being used for classification
#'
#' @returns A function with expression argument, which will return a logical
#'      vector, with TRUE for the samples predicted to be phenotype1 and FALSE
#'      for the samples predicted to be phenotype2
#' @examples
#' # example code
#' @export
#' @include dirac.R
.diracEntropyClassifier <- function(expression,
                                   phenotype1,
                                   phenotype2,
                                   geneIndex){
    # get the rank matrices for the expression matrix
    p1RankMatrix <- diracRankFunction(expression[geneIndex, phenotype1])
    p2RankMatrix <- diracRankFunction(expression[geneIndex, phenotype2])
    # Get the rank templates for each of the phenotypes
    p1RankTemplate <- diracRankTemplate(p1RankMatrix)
    p2RankTemplate <- diracRankTemplate(p2RankMatrix)
    function(expression) {
        # Get the rank matrix
        rankMatrix <- diracRankFunction(expression[geneIndex, ])
        # Get the p1 matching score
        p1MatchingScore <- diracMatrixMatchingScore(
            rankMatrix,
            p1RankTemplate
        )
        p2MatchingScore <- diracMatrixMatchingScore(
            rankMatrix,
            p2RankTemplate
        )
        # Return the difference between the two
        rmd <- p1MatchingScore - p2MatchingScore
        (rmd > 0)
    }
}

#' Creates a phenotype classifier based on RACE entropy
#'
#' @param expression Numeric matrix of gene expression values, rows are genes
#'    columns are samples
#' @param phenotype1,phenotype2 Numeric vectors with indices for each of the
#'    phenotypes
#' @param geneIndex Numeric vector of indices for the genes within the gene
#'    network being used for classification
#'
#' @returns A function with expression argument, which will return a logical
#'      vector, with TRUE for the samples predicted to be phenotype1 and FALSE
#'      for the samples predicted to be phenotype2
#' @examples
#' # example code
#' @export
#' @include race.R
.raceEntropyClassifier <- function(expression,
                                  phenotype1,
                                  phenotype2,
                                  geneIndex){
    p1RankMatrix <- expression[geneIndex, phenotype1]
    p2RankMatrix <- expression[geneIndex, phenotype2]
    function(expression){
        expressionRankMatrix <- expression[geneIndex,]
        # The cor of two matrices returns a matrix, with the columns of
        # the first matrix becoming the rows of the return, and the columns
        # of the second matrix becoming the columns of the return.
        p1Cor <- stats::cor(p1RankMatrix, expressionRankMatrix)
        p2Cor <- stats::cor(p2RankMatrix, expressionRankMatrix)
        # Get the mean correlation between each column of the input
        # expression matrix, and the columns of the phenotype rank matrices
        p1MeanCor <- base::colMeans(p1Cor)
        p2MeanCor <- base::colMeans(p2Cor)
        # Find which phenotype each has a higher correlation with
        corDiff <- p1MeanCor-p2MeanCor
        corDiff > 0
    }

}

#' Create a phenotype classifier based on CRANE entropy
#'
#' When the distance to the two rank centroids is equal, this will return
#' TRUE, so the ordering of phenotypes will affect classification
#'
#' @param expression Numeric matrix of gene expression values, rows are genes
#'    columns are samples
#' @param phenotype1,phenotype2 Numeric vectors with indices for each of the
#'    phenotypes
#' @param geneIndex Numeric vector of indices for the genes within the gene
#'    network being used for classification
#'
#' @returns A function with expression argument, which will return a logical
#'      vector, with TRUE for the samples predicted to be phenotype1 and FALSE
#'      for the samples predicted to be phenotype2
#' @examples
#' # example code
#' @export
#' @include crane.R
.craneEntropyClassifier <- function(expression,
                                   phenotype1,
                                   phenotype2,
                                   geneIndex){
    # get the rank matrices for the expression matrix
    p1RankMatrix <- craneRankFunction(expression[geneIndex, phenotype1])
    p2RankMatrix <- craneRankFunction(expression[geneIndex, phenotype2])
    # Get the rank centroids for the two phenotypes
    p1RankCentroid <- rowMeans(p1RankMatrix)
    p2RankCentroid <- rowMeans(p2RankMatrix)
    function(expression){
        expressionRankMatrix <- craneRankFunction(expression[geneIndex,])
        p1Distance <- sqrt(colSums((expressionRankMatrix-p1RankCentroid)^2))
        p2Distance <- sqrt(colSums((expressionRankMatrix-p1RankCentroid)^2))
        distDiff <- p2Distance - p1Distance
        distDiff>=0
    }
}


# Helper Functions --------------------------------------------------------

.wrapSeClassifier <- function(classifier,
                              colDataNameDefault,
                              assayDefault){
    function(seObject, assayName=assayDefault, colDataName=colDataNameDefault){
        res <- classifier(SummarizedExperiment::assays(seObject)[[assayName]])
        SummarizedExperiment::colData(seObject)[colDataName] <- res
        seObject
    }
}
