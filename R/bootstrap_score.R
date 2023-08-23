# Main Function -----------------------------------------------------------
#' Compute the pValue for the difference in entropy between two phenotypes for
#' list of gene networks
#'
#' @param expression Numeric matrix of gene expression values, rows
#'    representing genes, and columns representing samples
#' @param phenotype1,phenotype2 Index vectors, representing the indices of the
#'    phenotypes in the expression matrix
#' @param geneNetworkList List of index vectors, each vector representing the
#'    indices of a different gene networks
#' @param rankFun Function for creating a rank matrix, should take a matrix of
#'    dim (gene, sample) as input, and return a rank matrix, with samples as the
#'    columns.
#' @param scoreFun Function for scoring entropy in a phenotype, should take a
#'    rank matrix as input, and return a score (the rank matrix will only
#'    include samples within 1 phenotype)
#' @param bootstrapIterations Integer, number of iterations to perform when
#'    creating the null distribution
#' @param asFrame Boolean, whether the return should be a data.frame (TRUE),
#'    or a named list (FALSE)
#' @param replace Boolean, whether sampling the phenotypes should be done with
#'    replacement
#'
#' @return Either a data.frame, or a named list depending on the asFrame
#'    parameter. In the data.frame form, with columns for geneNetwork, p1Score,
#'    p2Score, absoluteDifference, and pValue. The list return is named
#'    according to the names in the geneNetworkList parameter, with each value
#'    being a named list with p1Score, p2Score, absoluteDifference, and pValue.
#' @export
#'
#' @examples
comparePhenotypes <-
  function(expression,
           geneNetworkList,
           phenotype1,
           phenotype2,
           rankFun,
           scoreFun,
           bootstrapIterations = 1000,
           replace = TRUE,
           asFrame = TRUE,
           BPPARAM=bpparam()) {
    # Ensure list is named
    geneNetworkList <- ensureNamed(geneNetworkList, prefix = "gn_")
    # Run the bootstrapping
    resList <- lapply(
      geneNetworkList,
      bootstrapScore,
      expression = expression,
      phenotype1 = phenotype1,
      phenotype2 = phenotype2,
      rankFun=rankFun,
      scoreFun=scoreFun,
      bootstrapIterations = bootstrapIterations,
      replace = replace,
      BPPARAM=BPPARAM
    )
    if (!asFrame) {
      names(resList) <- names(geneNetworkList)
      return(resList)
    }
    p1ScoreList <- sapply(resList, function(x)
      x$p1Score)
    p2ScoreList <- sapply(resList, function(x)
      x$p2Score)
    absoluteDifferenceList <- sapply(resList, function(x) x$absoluteDifference)
    pValuesList <- sapply(resList, function(x)
      x$pValue)
    geneNetworks <- names(geneNetworkList)
    data.frame(
      geneNetwork = geneNetworks,
      p1Score = p1ScoreList,
      p2Score = p2ScoreList,
      absoluteDifference = absoluteDifferenceList,
      pValue = pValuesList
    )
  }


# Helper Functions --------------------------------------------------------
#' Compare two phenotypes for a single gene network
#'
#' @param rankMatrix Numeric matrix With samples as the columns. See details
#'    for further information.
#' @param phenotype1,phenotype2 Index vector representing the indices of each
#'    phenotype within the `rankMatrix`.
#' @param scoreFun Function which calculates a score form the
#'    rank matrix.
#'
#' @return `numeric` Absolute difference of the scores for the two phenotypes
#' @export
#'
#' @details
#' Rank matrix is normally generated from a passed rank_fun function, so the
#' columns should represent samples, but what the rows represent depends
#' on the specific rank function. For RACE, CRANE, and INFER the rows simply
#' represent genes, while for DIRAC they represent specific pairwise
#' comparison.
#'
#' The phenotype vectors can contain any values that can be used with the
#' `[` matrix function, and the `c` function.
#'
#' @examples
comparePhenotypesSingle <-
  function(rankMatrix,
           phenotype1,
           phenotype2,
           scoreFun) {
    p1RankMatrix <- rankMatrix[, phenotype1]
    p2RankMatrix <- rankMatrix[, phenotype2]
    p1Score <- scoreFun(p1RankMatrix)
    p2Score <- scoreFun(p2RankMatrix)
    abs(p1Score - p2Score)
  }

#' Compare two shuffled phenotypes
#'
#' @param iteration Value representing the current bootstrap iteration,
#'    used to make signature match the required signature for `bplapply`.
#' @param rankMatrix Numeric matrix with samples as the columns. See details
#'    for further information.
#' @param combined Vector representing the concatenated indices of
#'    phenotype 1 and 2.
#' @param p1Size,p2Size Integer representing the size of the two phenotypes
#' @param scoreFun Function which calculates a score form the
#'    rank matrix.
#' @param replace Boolean representing Whether sampling should be done with
#'    replacement
#'
#' @return
#' @export
#'
#' @details
#' Rank matrix is normally generated from a passed rank_fun function, so the
#' columns should represent samples, but what the rows represent depends
#' on the specific rank function. For RACE, CRANE, and INFER the rows simply
#' represent genes, while for DIRAC they represent specific pairwise
#' comparison.
#'
#' The combined vector can contain any values that can be used with the
#' `[` matrix function, as well as the vector `c` and `sample` functions.
#'
#' The replace parameter dictates whether the sampling to create shuffled
#' phenotype index vectors should be done with replacement. If TRUE, sampling
#' with replacement if performed to create the two shuffled index vectors
#' of the appropriate sizes. If FALSE, random sampling without replacement
#' occurs once to select phenotype1 indices, and the remaining indices are used
#' for phenotype2.
#'
#' @examples
scoreShuffled <-
  function(iteration,
           rankMatrix,
           combined,
           p1Size,
           p2Size,
           scoreFun,
           replace = TRUE) {
      if (replace) {
        p1Idx <- sample(combined, p1Size, replace = replace)
        p2Idx <- sample(combined, p2Size, replace = replace)
      } else {
        p1Idx <- sample(combined, p1Size, replace = replace)
        p2Idx <- setdiff(combined, p1Idx)
      }
      ## Can't call comparePhenotypeSingle due to snow cluster,
      ## instead, just perform the actions here
      p1RankMatrix <- rankMatrix[, p1Idx]
      p2RankMatrix <- rankMatrix[, p2Idx]
      p1Score <- scoreFun(p1RankMatrix)
      p2Score <- scoreFun(p2RankMatrix)
      abs(p1Score-p2Score)
  }



#' Find a p-value for the difference in entropy between two phenotypes using a
#' bootstrapped distribution
#'
#' @param expression Numeric matrix of gene expression, with rows representing
#'    genes, and columns representing samples.
#' @param phenotype1,phenotype2 Index vectors, representing the indices of the
#'    phenotypes in the expression matrix
#' @param geneNetwork Index Vector, representing the indices of the genes in the
#'    gene network
#' @param rankFun Function for creating a rank matrix, should take a matrix of
#'    dim (gene, sample) as input, and return a rank matrix, with samples as the
#'    columns.
#' @param scoreFun Function for scoring entropy in a phenotype, should take a
#'    rank matrix as input, and return a score (the rank matrix will only
#'    include samples within 1 phenotype)
#' @param bootstrapIterations Integer, number of iterations to perform when
#'    creating the null distribution
#' @param replace Boolean, whether sampling the phenotypes should be done with
#'    replacement
#'
#' @return Named list of results, specifically p1Score, p2Score,
#'    absoluteDifference, and pValue
#' @importFrom BiocParallel bplapply SnowParam SerialParam MulticoreParam
#' @export
#'
#' @examples
bootstrapScore <-
  function(geneNetwork,
           expression,
           rankFun,
           scoreFun,
           phenotype1,
           phenotype2,
           bootstrapIterations = 1000,
           replace = TRUE,
           BPPARAM=bpparam()) {
    # Get the lengths of the phenotypes
    p1Size <- length(phenotype1)
    p2Size <- length(phenotype2)
    # Get the combined phenotype index vector
    combinedPhenotypes <- c(phenotype1, phenotype2)
    # Compute the rank matrix using the provided method
    rankMatrix <-
      rankFun(expression[geneNetwork,])
    ## Having to rank the whole dataset to avoid issues with changes in index,
    ## should find a way to only rank the needed phenotypes
    ## Compute the unshuffled score for the two phenotypes
    p1Score <- scoreFun(rankMatrix[, phenotype1])
    p2Score <- scoreFun(rankMatrix[, phenotype2])
    absDiff <- abs(p1Score - p2Score)
    res <-
      unlist(
        bplapply(
          seq_len(bootstrapIterations),
          scoreShuffled,
          rankMatrix = rankMatrix,
          combined = combinedPhenotypes,
          p1Size = p1Size,
          p2Size = p2Size,
          scoreFun=scoreFun,
          replace = replace,
          BPPARAM = BPPARAM
        )
      )
    # Create the empirical CDF
    bootCdf <- stats::ecdf(res)
    pValue <- 1 - bootCdf(absDiff)
    # Return results list
    list(
      p1Score = p1Score,
      p2Score = p2Score,
      absoluteDifference = absDiff,
      pValue = pValue
    )
}
