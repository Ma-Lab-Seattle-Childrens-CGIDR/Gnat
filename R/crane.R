# Bootstrap Score ---------------------------------------------------------
#' Find a p-value for the difference in entropy between two phenotypes according
#' to the CRANE method using a bootstrapped distribution
#'
#' @param expression Numeric matrix of gene expression, with rows representing
#'    genes, and columns representing samples.
#' @param phenotype1,phenotype2 Index vectors, representing the indices of the
#'    phenotypes in the expression matrix
#' @param geneNetwork Index Vector, representing the indices of the genes in the
#'    gene network
#' @param bootstrapIterations Integer, number of iterations to perform when
#'    creating the null distribution
#' @param replace Boolean, whether sampling the phenotypes should be done with
#'    replacement
#'
#' @return Named list of results, specifically p1Score, p2Score,
#'    absoluteDifference, and pValue
#' @export
#'
#' @examples
craneBootstrapScore <- function(expression, geneNetwork, phenotype1,
                                phenotype2, bootstrapIterations=1000,
                                replace=TRUE, BPPARAM=bpparam()){
    bootstrapScore(geneNetwork = geneNetwork, expression = expression,
                   rankFun=craneRankFunction, scoreFun=craneScoreFunction,
                   phenotype1=phenotype1, phenotype2=phenotype2,
                   bootstrapIterations = bootstrapIterations,
                   replace=replace, BPPARAM=BPPARAM)
}



# Compare Phenotypes ------------------------------------------------------

#' Compare the entropy difference between phenotypes for a list of gene
#' networks using the CRANE method
#'
#' @param expression Numeric matrix of gene expression values, rows
#'    representing genes, and columns representing samples
#' @param geneNetworkList List of index vectors, each vector representing the
#'    indices of a different gene networks
#' @param phenotype1,phenotype2 Index vectors, representing the indices of the
#'    phenotypes in the expression matrix
#' @param bootstrapIterations Integer, number of iterations to perform when
#'    creating the null distribution
#' @param replace Boolean, whether sampling the phenotypes should be done with
#'    replacement
#' @param asFrame Boolean, whether the return should be a data.frame (TRUE),
#'    or a named list (FALSE)
#'
#' @return Either a data.frame, or a named list depending on the asFrame
#'    parameter. In the data.frame form, with columns for geneNetwork, p1Score,
#'    p2Score, absoluteDifference, and pValue. The list return is named
#'    according to the names in the geneNetworkList parameter, with each value
#'    being a named list with p1Score, p2Score, absoluteDifference, and pValue.
#' @export
#'
#' @examples
craneComparePhenotypes <- function(expression, geneNetworkList,
                                   phenotype1, phenotype2,
                                   bootstrapIterations=1000,
                                   replace=TRUE, asFrame=TRUE,
                                   BPPARAM=bpparam()){
    comparePhenotypes(expression=expression, geneNetworkList = geneNetworkList,
                      phenotype1=phenotype1, phenotype2=phenotype2,
                      rankFun=craneRankFunction, scoreFun=craneScoreFunction,
                      bootstrapIterations=bootstrapIterations,replace=replace,
                      asFrame=asFrame, BPPARAM=BPPARAM)
}


# Rank Function -----------------------------------------------------------

#' Rank function for the CRANE method, rank each value in a column
#'
#' @param filteredExpression Filtered expression, where each row is a gene
#'      within the gene network, and each column is a sample
#'
#' @return Numeric matrix, with each value representing the rank of the
#'      expression within a sample compared to the other genes in the network
#' @export
#'
#' @examples
#' testMat <- matrix(runif(20), ncol=4, nrow=5)
#' print(craneRankFunction(testMat))
#' print(craneRankFunction(testMat, margin=1))
craneRankFunction <- function(filteredExpression){
    simpleRank(filteredExpression, margin=2)
}


# Score Function ----------------------------------------------------------

#' Function to score a rank matrix using the CRANE method, the average distance
#' between samples and the rank centroid
#'
#' @param rankMatrix Numeric matrix, a matrix where each value represents
#'      the rank of the gene's expression within the gene network
#'
#' @return Numeric score, value representing the mean distance to the rank
#'      centroid
#' @export
#'
#' @examples
craneScoreFunction <- function(rankMatrix){
    centroid <- apply(rankMatrix, MARGIN=1, mean, na.rm=TRUE)
    mean(sqrt(apply((rankMatrix-centroid)^2, MARGIN=2, sum)))
}



# Sample Score ------------------------------------------------------------
#' Calculate the by sample entropy scores using the CRANE method
#'
#' @param filteredExpression Numeric Matrix, rows representing genes in
#'      the gene network, columns representing samples in the phenotype
#'
#' @return Numeric vector, with a score for the entropy of each sample
#' @export
#'
#' @examples
craneSampleScore <- function(filteredExpression){
    rankMatrix <- craneRankFunction(filteredExpression = filteredExpression)
    centroid <- apply(rankMatrix, MARGIN=1, mean, na.rm=TRUE)
    sqrt(apply((rankMatrix-centroid)^2, MARGIN=2, sum))
}
