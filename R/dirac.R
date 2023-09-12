# Citation ----------------------------------------------------------------

# Based on: Eddy, James A., et al. “Identifying Tightly Regulated and Variably
# Expressed Networks by Differential Rank Conservation (DIRAC).” PLoS
# Computational Biology, vol. 6, no. 5, May 2010.
# go-gale-com.offcampus.lib.washington.edu,
# https://doi.org/10.1371/journal.pcbi.1000792.

# Bootstrap Score ---------------------------------------------------------
#' Find a p-value for the difference in entropy between two phenotypes according
#' to the DIRAC method using a bootstrapped distribution
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
diracBootstrapScore <- function(expression,
                                geneNetwork,
                                phenotype1,
                                phenotype2,
                                bootstrapIterations=1000,
                                replace=TRUE,
                                asFrame=TRUE,
                                BPPARAM=bpparam()){
    bootstrapScore(geneNetwork = geneNetwork, expression=expression,
                   rankFun=diracRankFunction, scoreFun = diracScoreFunction,
                   phenotype1 = phenotype1, phenotype2 = phenotype2,
                   bootstrapIterations = bootstrapIterations,
                   replace = replace, BPPARAM=BPPARAM)
}


# Compare Phenotypes ------------------------------------------------------

#' Compare the entropy difference between phenotypes for a list of gene
#' networks using the DIRAC method
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
diracComparePhenotypes <- function(expression,
                                   geneNetworkList,
                                   phenotype1,
                                   phenotype2,
                                   bootstrapIterations=1000,
                                   replace=TRUE,
                                   asFrame=TRUE,
                                   BPPARAM=bpparam()){
    comparePhenotypes(expression=expression, geneNetworkList=geneNetworkList,
                      phenotype1=phenotype1, phenotype2=phenotype2,
                      rankFun=diracRankFunction, scoreFun=diracScoreFunction,
                      bootstrapIterations = bootstrapIterations,
                      replace=replace, asFrame=asFrame, BPPARAM=BPPARAM)
}


# Rank Function -----------------------------------------------------------
#' Determine DIRAC ranking vector
#'
#' `dirac.rank_vector` creates a DIRAC ranking vector from expression data.
#'
#' @param expression A numeric vector of gene expression data
#' @param expression.length Optional, length of expression matrix, useful for
#'    repeated calls to this function with same expression length
#' @param matrix.index Optional, a Boolean vector or matrix with TRUE for
#'    entries of gene-by-gene square matrix to be included in the rank_vector.
#'    Mainly used for repeated calls to this function so it doesn't have to be
#'    recalculated each time.
#'
#' @return Numeric rank vector
#'
#' @examples
#' dirac.rank_vector(c(4, 2, 1, 3))
#' dirac.rank_vector(c(4, 2, 1, 3), 4)
#' @export
diracRankVector <- function(expressionVector, expressionLength=NULL,
                            matrixIndex=NULL){
    ## Calculate expression Length if needed
    if(is.null(expressionLength)){
        expressionLength <- length(expressionVector)
    }
    ## Calculate matrixIndex if needed
    if(is.null(matrixIndex)){
        matrixIndex <- lower.tri(matrix(1, nrow = expressionLength,
                                        ncol = expressionLength,))
    }
    if(expressionLength==0){
        stop("Expression vector must not be empty")
    }
    expressionMatrix <- matrix(expressionVector, nrow = expressionLength,
                               ncol = expressionLength, byrow=TRUE)
    ((expressionMatrix-t(expressionMatrix)) < 0)[matrixIndex]
}

#' Determine DIRAC ranking matrix
#'
#' Creates a DIRAC ranking matrix from expression data.
#'
#' @param expression A gene expression matrix, rows are genes, columns are
#'    samples
#'
#' @return matrix whose columns represent the rank_vector for each sample
#'
#' @examples
#' dirac.rank_matrix(matrix(1:16, ncol = 4, nrow = 4))
#'
#' @export
diracRankFunction <- function(filteredExpression){
    # Calculate expression length
    expressionLength <- nrow(filteredExpression)
    matrixIndex <- as.vector(
        lower.tri(matrix(1, nrow = expressionLength, ncol = expressionLength))
    )
    ## Define function here, so that the rank function will work with
    ## BiocParallel using snow cluster
    diracRankVector <- function(expressionVector, expressionLength,
                                matrixIndex){
        expressionMatrix <- matrix(expressionVector, nrow = expressionLength,
                                   ncol = expressionLength, byrow=TRUE)
        ((expressionMatrix-t(expressionMatrix)) < 0)[matrixIndex]
    }
    apply(filteredExpression, 2, diracRankVector,
          expressionLength = expressionLength, matrixIndex=matrixIndex)
}
# Score Function ----------------------------------------------------------
#' Generate rank template
#'
#' Creates a DIRAC ranking template from expression data.
#'
#' Strict inequality is used for finding the rank template.
#'
#' @param rank_matrix Rank matrix, from dirac.rank_matrix
#'
#' @return Numeric vector, rank template
#'
#' @examples
#' dirac.rank_template(matrix(1:16, ncol = 4, nrow = 4))
#'
#' @export
diracRankTemplate <- function(rankMatrix){
    rowMeans(rankMatrix) > 0.5
}

#' Find rank matching score
#'
#' Finds the rank matching score between a rank vector, and a rank template.
#'
#' @param rank_vector Boolean vector representing the rank vector, must be
#'    the same length as rank_template.
#' @param rank_template Boolean vector representing the rank template.
#'
#' @returns Float, representing the rank matching score
#'
#' @examples
#' dirac.rank_matching_score(c(1, 0, 0, 1, 0), c(0, 1, 0, 1, 0))
#' @export
diracRankMatchingScore <- function(rankVector, rankTemplate){
    if(length(rankVector)!=length(rankTemplate)){
        stop("rankVector and rankTemplate must be same length")
    }
    if(length(rankVector)==0){
        stop("Neither rankVector, nor rankTemplate can be empty")
    }
    mean(rankVector==rankTemplate)
}

#' Find rank matching score across samples
#'
#' Finds the rank matching score between each
#'    sample in a rank_matrix, and a rank_template
#'
#' @param rank_matrix Boolean matrix, columns represent each samples rank
#'      vector
#' @param rank_template Boolean vector, represents template to match against
#'
#' @returns Numeric vector of rank matching scores
#'
#' @examples
#' expression <- matrix(1:16, ncol = 4, nrow = 4)
#' rank_matrix <- dirac.rank_matrix(expression)
#' rank_template <- dirac.rank_template(rank_matrix)
#' dirac.rank_matching_score.vector(rank_matrix, rank_template)
#' @export
diracMatrixMatchingScore <- function(rankMatrix, rankTemplate){
    apply(rankMatrix, 2, diracRankMatchingScore, rankTemplate=rankTemplate)
}

#' Score using DIRAC
#'
#' Finds the rank rank conservation index in gene expression data
#'
#' @param rank_matrix Boolean matrix, columns represent each samples rank
#'      vector
#' @param rank_template Boolean vector, represents template to match against
#'
#' @returns Float representing the rank conservation index
#'
#' @examples
#' expression <- matrix(1:16, ncol = 4, nrow = 4)
#' rankMatrix <- diracRankMatrix(expression)
#' rankTemplate <- diracRankTemplate(rankMatrix)
#' diracScoreFunction(rankMatrix, rankTemplate)
#' @export
diracScoreFunction <- function(rankMatrix){
    rankTemplate <- (rowMeans(rankMatrix) > 0.5)
    ## Implemented this way so that it works with BiocParallel on snow
    ## clusters
    mean(
        apply(rankMatrix, 2,
              function(rankVector, rankTemplate) mean(rankVector==rankTemplate),
              rankTemplate=rankTemplate)
    )
}


# Sample Score ------------------------------------------------------------

#' Compute Sample-wise entropy using DIRAC
#'
#' @param filteredExpression Numeric matrix, representing gene expression,
#'      with rows representing genes in the network of interest, and columns
#'      representing samples within a phenotype
#'
#' @return Numeric vector, DIRAC Rank Matching Score for each sample
#' @export
#'
#' @examples
diracSampleScore <- function(filteredExpression){
    rankMatrix <- diracRankFunction(filteredExpression = filteredExpression)
    rankTemplate <- diracRankTemplate(rankMatrix)
    diracMatrixMatchingScore(rankMatrix=rankMatrix, rankTemplate=rankTemplate)
}


# Gene Network Classification Comparison ----------------------------------
#' Find the Rank Difference Score
#'
#' Finds the rank difference score for a network
#' between 2 phenotypes. The score against the rank template for the second
#' phenotype is subtracted from the score against the rank template for the
#' first phenotype.
#'
#' @param expression Numeric matrix with gene expression data, rows are
#'    genes, columns are samples.
#' @param phenotype1,phenotype2 Numeric vectors containing the indexes for the
#'    phenotypes.
#' @param geneInex Numeric vector of indices for the genes in the network
#'
#' @returns A list of two vectors, each containing the rank difference scores
#'    for the samples within the phenotypes.
#'
#' @examples
#' expression <- matrix(1:16, ncol = 4, nrow = 4)
#' dirac.rank_difference_score(expression,
#'     phenotype1 = c(1, 4), phenotype2 =
#'         c(2, 3), gene_index = c(1, 2, 4)
#' )
#' @export
diracRankDifferenceScore <- function(expression, phenotype1, phenotype2,
                                     geneIndex) {
    # Find the rank matrix for each phenotype
    p1RankMatrix <- diracRankFunction(expression[geneIndex, phenotype1])
    p2RankMatrix <- diracRankFunction(expression[geneIndex, phenotype2])
    # Find the rank templates for the gene network in each of the phenotypes
    p1RankTemplate <- diracRankTemplate(p1RankMatrix)
    p2RankTemplate <- diracRankTemplate(p2RankMatrix)
    # Find the rank difference score for each phenotype
    p1RankDifference <- diracMatrixMatchingScore(
        p1RankMatrix,
        p1RankTemplate
    ) -
        diracMatrixMatchingScore(p1RankMatrix, p2RankTemplate)
    p2RankDifference <- diracMatrixMatchingScore(
        p2RankMatrix,
        p1RankTemplate
    ) -
        diracMatrixMatchingScore(p2RankMatrix, p2RankTemplate)
    names(p1RankDifference) <- phenotype1
    names(p2RankDifference) <- phenotype2
    list(p1RankDifference, p2RankDifference)
}


#' Find the classification rate for single network
#'
#' Finds the classification rate between two
#'    phenotypes based on a gene network.
#'
#' Classification rate is defined as the average between sensitivity and
#'   specificity. classification_rate =
#'   Pr(rank_diff>0|phenotype=A)*Pr(phenotype=A) +
#'   Pr(rank_diff<=0|phenotype=B)*Pr(phenotype=B)
#'
#' Due to the way that the classification rate splits 0's, this value depends
#' on the ordering of the gene_index vector. The order of this vector should
#' therefor be kept consistent.
#'
#' @param expression Numeric matrix with gene expression data, rows are genes
#'     columns are samples.
#' @param phenotype1,phenotype2 Numeric vectors containing the indexes for the
#'    two phenotypes to compare.
#' @param gene_index Numeric vector of index for genes in network
#'
#' @returns Float representing the classification rate based on DIRAC
#' @examples
#' # example code
#' @export
diracClassificationRate <- function(geneIndex, expression, phenotype1,
                                    phenotype2) {
    # Get the rank differences
    rankDiff <- diracRankDifferenceScore(expression,
                                         phenotype1 = phenotype1,
                                         phenotype2 = phenotype2,
                                         geneIndex = geneIndex
    )
    # Unpack the list
    p1RankDiff <- rankDiff[[1]]
    p2RankDiff <- rankDiff[[2]]
    # calculate Pr(diff_score>O|phenotype=A), number true positives
    tp <- mean(p1RankDiff > 0)
    # calculate Pr(diff_score<=O|phenotype=B), number true negative
    tn <- mean(p2RankDiff <= 0)
    # Calculate Pr(phenotype=A) and Pr(phenotype=B)
    p1Num <- length(phenotype1)
    p2Num <- length(phenotype2)
    total <- p1Num + p2Num
    p1Prob <- p1Num / total
    p2Prob <- p2Num / total
    # Calculate the classification rate
    (tp * p1Prob + tn * p2Prob)
}

#' Compare classification rate for different gene networks
#'
#' `dirac.classification_rate.compare` compares the classification rate between
#'    multiple gene networks for two phenotypes
#'
#' @param expression Numeric matrix, represents gene expression, rows are genes
#'    columns are samples
#' @param phenotype1,phenotype2 Numeric vectors, contain the indices for the
#'    samples corresponding to each phenotype.
#' @param gene_index List of numeric vectors, each containing the indices for
#'    one of the gene networks.
#' @param parallel Boolean indicating Whether the computation should be run in
#'   parallel. On unix makes use of the R core library parallel, on windows
#'   requires foreach, and doParallel.
#' @param cores Integer, indicates number of cores to use for parallel
#'    calculations. If greater than parallel::detectCores() return value,
#'    will instead be set to said return value.
#' @param as.frame Boolean, whether return value should be data.frame or named
#'    list, see return for more information.
#'
#' @returns
#'    if as.frame is TRUE:
#'      data.frame with gene_network names as one column, and classification
#'      rate as another.
#'    if as.frame is FALSE:
#'      Numeric vector containing the classification rate for each of the
#'      provided gene networks. Named based on names of gene_index argument.
#' @examples
#' # example code
#' @export
diracClassificationRateCompare <- function(expression, phenotype1,
                                           phenotype2, geneIndexList,
                                           asFrame = FALSE) {

    geneIndexList <- ensureNamed(geneIndexList)
    res <- unlist(
        lapply(geneIndexList, diracClassificationRate,
                 expression = expression, phenotype1 = phenotype1,
                 phenotype2 = phenotype2
        )
    )
    if (!asFrame) {
        # Rename the results based on the names of the gene_index list arg
        names(res) <- names(geneIndexList)
        # Return the renamed results
        return(res)
    }
    data.frame(geneNetwork = names(geneIndex), classificationRate = res)
}

#' Compare classification rate for different gene networks
#'
#' `dirac.classification_rate.best` finds the network with the best
#'    classification rate between multiple gene networks for two phenotypes.
#'
#' @param expression Numeric matrix, represents gene expression, rows are genes
#'    columns are samples
#' @param phenotype1,phenotype2 Numeric vectors, contain the indices for the
#'    samples corresponding to each phenotype.
#' @param gene_index List of numeric vectors, each containing the indices for
#'    one of the gene networks.
#' @param parallel Boolean indicating Whether the computation should be run in
#'   parallel. On unix makes use of the R core library parallel, on windows
#'   requires foreach, and doParallel.
#' @param cores Integer, indicates number of cores to use for parallel
#'    calculations. If greater than parallel::detectCores() return value,
#'    will instead be set to said return value.
#'
#' @returns String, name of best gene network
#' @examples
#' # example code
#' @export
diracClassificationRateBest <-
    function(expression,
             phenotype1,
             phenotype2,
             gene_index) {
        names(which.max(
            diracClassificationRateCompare(
                expression,
                phenotype1,
                phenotype2,
                gene_index,
                asFrame = FALSE
            )
        ))
    }
