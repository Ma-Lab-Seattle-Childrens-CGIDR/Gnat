# Bootstrap Score ---------------------------------------------------------
#' Find a p-value for the difference in entropy between two phenotypes according
#' to the INFER method using a bootstrapped distribution
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
inferBootstrapScore <- function(expression, geneNetwork, phenotype1,
                                phenotype2, bootstrapIterations=1000,
                                replace=TRUE, BPPARAM=bpparam()){
    bootstrapScore(geneNetwork = geneNetwork, expression = expression,
                   rankFun=inferRankFunction, scoreFun=inferScoreFunction,
                   phenotype1=phenotype1, phenotype2=phenotype2,
                   bootstrapIterations = bootstrapIterations,
                   replace=replace, BPPARAM=BPPARAM)
}


# Compare Phenotypes ------------------------------------------------------

#' Compare the entropy difference between phenotypes for a list of gene
#' networks using the INFER method
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
#' @include bootstrap_score.R
inferComparePhenotypes <- function(expression, geneNetworkList,
                                   phenotype1, phenotype2,
                                   bootstrapIterations=1000,
                                   replace=TRUE, asFrame=TRUE,
                                   BPPARAM=bpparam()){
    comparePhenotypes(expression=expression, geneNetworkList = geneNetworkList,
                      phenotype1=phenotype1, phenotype2=phenotype2,
                      rankFun=inferRankFunction, scoreFun=inferScoreFunction,
                      bootstrapIterations=bootstrapIterations,replace=replace,
                      asFrame=asFrame, BPPARAM=BPPARAM)
}

# Rank Function -----------------------------------------------------------
#' Rank function for the INFER method, rank each value in a column
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
inferRankFunction <- function(filteredExpression){
    simpleRank(filteredExpression, margin=2)
}

# Score Function ----------------------------------------------------------

#' Determine Entropy of a Vector
#'
#' `vectorEntropy` finds the information entropy of a vector
#'
#' @param vec Numeric vector (or other class which implements table and
#'    length methods), for which the entropy is calculated
#'
#' @return Float, representing the entropy of the vector
#' @export
#'
#' @examples
#' # Create an example vector
#' example_vector <- c(1, 2, 3, 2, 1, 2, 3, 4, 5, 6, 6, 7, 1)
#' # Print the entropy of the vector
#' print(entropy.vector(example_vector))
vectorEntropy <- function(vec){
    t <- table(vec)/length(vec)
    -sum(t*log2(t))
}

#' Compute entropy for a matrix
#'
#' @param mat Numeric matrix to calculate the entropy of
#' @param margin Integer, determines which axis the calculation is performed
#'    along, 1 for rows, 2 for columns
#'
#' @return Numeric Vector, each element representing the entropy of
#'    corresponding row or column
#' @export
#'
#' @examples
#' # Create an example matrix
#' example_matrix <- matrix(rnorm(4 * 4), ncol = 4)
#' # Print result vector
#' print(entropy.matrix(example_matrix, 2))
matrixEntropy <- function(mat, margin = 1) {
    apply(mat, MARGIN = margin, vectorEntropy)
}

#' Function to score a rank matrix using the INFER method, the average
#' information entropy of genes within the network
#'
#' @param rankMatrix Numeric matrix, a matrix where each value represents
#'      the rank of the gene's expression within the gene network
#'
#' @return Numeric score, value representing the mean information entropy
#'      for the genes within the network in the given phenotype
#' @export
#'
#' @examples
inferScoreFunction <- function(rankMatrix){
    ## Defining these internally to this function so that they will
    ## work with BiocParallel
    vectorEntropy <- function(vec){
        t <- table(vec)/length(vec)
        -sum(t*log2(t))
    }
    matrixEntropy <- function(mat, margin = 1) {
        apply(mat, MARGIN = margin, vectorEntropy)
    }
    mean(matrixEntropy(rankMatrix))
}



# Gene Entropy ------------------------------------------------------------

#' Find the rank entropy for each gene in a network
#'
#' `inferGeneEntropy` Determines the rank entropy for each gene in a network
#'
#' @param filteredExpression Numeric matrix representing gene expression, with
#'    rows representing genes in the network and columns representing samples
#'    in the phenotype.
#'
#' @return Numeric vector representing the entropy for each gene in the network
#'    in the given phenotype.
#' @export
#'
#' @examples
inferGeneEntropy <- function(filteredExpression) {
    matrixEntropy(inferRankFunction(filteredExpression = filteredExpression))
}
