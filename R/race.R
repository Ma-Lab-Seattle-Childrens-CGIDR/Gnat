# Bootstrap Score ---------------------------------------------------------


#' Find a p-value for the difference in entropy between two phenotypes acording
#' to the RACE method using a bootstrapped distribution
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
raceBootstrapScore <- function(expression, geneNetwork, phenotype1,
                                phenotype2, bootstrapIterations=1000,
                                replace=TRUE, BPPARAM=bpparam()){
    bootstrapScore(geneNetwork = geneNetwork, expression = expression,
                   rankFun=raceRankFunction, scoreFun=raceScoreFunction,
                   phenotype1=phenotype1, phenotype2=phenotype2,
                   bootstrapIterations = bootstrapIterations,
                   replace=replace, BPPARAM=BPPARAM)
}

# Compare Phenotypes ------------------------------------------------------

#' Compare the entropy difference between phenotypes for a list of gene
#' networks using the RACE method
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
raceComparePhenotypes <- function(expression, geneNetworkList, phenotype1,
                                  phenotype2, bootstrapIterations=1000,
                                  replace=TRUE, asFrame=TRUE,
                                  BPPARAM=bpparam()){
    comparePhenotypes(expression=expression, geneNetworkList = geneNetworkList,
                      phenotype1=phenotype1, phenotype2=phenotype2,
                      rankFun=raceRankFunction, scoreFun=raceScoreFunction,
                      bootstrapIterations=bootstrapIterations,replace=replace,
                      asFrame=asFrame, BPPARAM=BPPARAM)
}


# Rank Function -----------------------------------------------------------
#' Rank function for the RACE method, simply the identity function
#'
#' @param filteredExpression Filtered expression, where each row is a gene
#'      within the gene network, and each column is a sample
#'
#' @return Numeric matrix, identical to input
#'
#' @details
#' This function just returns the same matrix that it is provided, it is used
#' for providing compatability with the other methods
#'
#'
#' @export
#'
#' @examples
#' testMat <- matrix(runif(20), ncol=4, nrow=5)
#' print(raceRankFunction(testMat))
#' print(raceRankFunction(testMat, margin=1))
raceRankFunction <- function(filteredExpression){
    filteredExpression
}


# Score Function ----------------------------------------------------------

#' Function to score a rank matrix using RACE method, the average of the
#' Kendall Tau rank correlation coefficient between every pair of samples
#' in a phenotype
#'
#' @param rankMatrix Numeric matrix, for this method this is the same as the
#'      expression matrix
#'
#' @return Numeric score, value representing the mean Kendall Tau between all
#'      pairs of samples within a phenotype
#'
#'
#' @export
#'
#' @examples
raceScoreFunction <- function(rankMatrix){
    correlation <- stats::cor(rankMatrix, method="kendall")
    mean(correlation[lower.tri(correlation)])
}



# Sample Score ------------------------------------------------------------

#' Calculate the by sample entropy scores using the RACE method
#'
#' @param filteredExpression Numeric Matrix, rows representing genes in
#'      the gene network, columns representing samples in the phenotype
#'
#' @return Numeric vector, with a score for the entropy of each sample
#' @export
#'
#' @examples
raceSampleScore <- function(filteredExpression){
    rankMatrix <- raceRankFunction(filteredExpression=filteredExpression)
    correlation <- stats::cor(rankMatrix, method = "kendall")
    # Remove diagonal values from the means
    diag(correlation) <- NA
    rowMeans(correlation, na.rm=TRUE)
}




