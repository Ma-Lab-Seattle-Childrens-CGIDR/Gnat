# Citation ----------------------------------------------------------------

# Based on: Eddy, James A., et al. “Identifying Tightly Regulated and Variably
# Expressed Networks by Differential Rank Conservation (DIRAC).” PLoS
# Computational Biology, vol. 6, no. 5, May 2010.
# go-gale-com.offcampus.lib.washington.edu,
# https://doi.org/10.1371/journal.pcbi.1000792.


# Main Functions ----------------------------------------------------------

#' Creates a phenotype classifier based on DIRAC
#'
#' `DIRAC.classifier` Creates a function which takes in a matrix of gene
#' expression values and returns a classification of phenotype1, or phenotype2
#' for each of the samples.
#'
#' @param expression Numeric matrix of gene expression values, rows are genes
#'    columns are samples
#' @param phenotype1,phenotype2 Numeric vectors with indices for each of the
#'    phenotypes
#' @param gene_index Numeric vector of indices for the genes within the gene
#'    network being used for classification
#' @param return String describing desired return value from created function:
#'    "name":  returns a character vector of the names of the predicted
#'      phenotypes
#'    "is phenotype1": returns a boolean vector indicating if the
#'      predicted phenotype is phenotype1
#'    "is phenotype2": returns a boolean vector indicating if the predicted
#'      phenotype is phenotype2
#'    "phenotype number": returns numeric vector, with 1 representing the
#'      predicted phenotype is phenotype 1, and 2 representing that the
#'      predicted phenotype is phenotype2
#'    Any other value will lead the function to just return the rank matching
#'      difference
#' @param names optional Character vector with the desired names for the
#'      phenotypes
#'
#' @returns A function which takes as an argument a gene expression numeric
#'    matrix with genes as the rows and samples as the columns, and returns
#'    a vector predicting which samples are from which phenotype (see return
#'    argument).
#' @examples
#' # example code
#' @export
diracClassifier <- function(expression, phenotype1, phenotype2, geneIndex,
                             return = "is phenotype1",
                             names = c("phenotype1", "phenotype2")) {
    # get the rank matrices for the expression matrix
    p1RankMatrix <- diracRankFunction(expression[geneIndex, phenotype1])
    p2RankMatrix <- diracRankFunction(expression[geneIndex, phenotype2])
    # Get the rank templates for each of the phenotypes
    p1RankTemplate <- diracRankTemplate(p1RankMatrix)
    p2RankTemplate <- diracRankTemplate(p2RankMatrix)
    # Create an intermediate function which returns a vector of the rank
    #   matching score differences from an expression matrix
    rankMatchingDiff <- function(expression) {
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
        p1MatchingScore - p2MatchingScore
    }
    # Based on desired return type, define a function
    if (return == "name") {
        if (is.null(names)) {
            names <- c("phenotype1", "phenotype2")
        }
        return(
            function(expression) {
                # Get the rank matching score difference
                rmd <- rankMatchingDiff(expression)
                ifelse(rmd > 0, names[1], names[2])
            }
        )
    } else if (return == "is phenotype1") {
        return(
            function(expression) {
                # Get the rank matching score difference
                rmd <- rankMatchingDiff(expression)
                (rmd > 0)
            }
        )
    } else if (return == "is phenotype2") {
        return(
            function(expression) {
                # Get the rank matching score difference
                rmd <- rankMatchingDiff(expression)
                (rmd <= 0)
            }
        )
    } else if (return == "phenotype number") {
        return(
            function(expression) {
                # Get the rank matching score difference
                rmd <- rankMatchingDiff(expression)
                ifelse(rmd > 0, 1, 2)
            }
        )
    } else {
        return(
            function(expression) {
                # Get the rank matching score difference
                rmd <- rankMatchingDiff(expression)
                rmd
            }
        )
    }
}





# Bootstrap Score ---------------------------------------------------------

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
  mean(rankVector==rankTemplate)
}

#' Find rank matching score across samples
#'
#' `dirac.rank_matching_score.vector` finds the rank matching score between each
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
#' rank_matrix <- dirac.rank_matrix(expression)
#' rank_template <- dirac.rank_template(rank_matrix)
#' dirac.rank_conservation_index(rank_matrix, rank_template)
#' @export
diracScoreFunction <- function(rankMatrix, rankTemplate){
  mean(diracMatrixMatchingScore(rankMatrix=rankMatrix,
                                rankTemplate=rankTemplate))
}


# Sample Score ------------------------------------------------------------

diracSampleScore <- function(filteredExpression){
  rankMatrix <- diracRankFunction(filteredExpression=filteredExpression)
  rankTemplate <- diracRankTemplate(rankMatrix=rankMatrix)
  diracMatrixMatchingScore(rankMatrix=rankMatrix, rankTemplate=rankTemplate)
}



# Helper Functions --------------------------------------------------------


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
                                              phenotype2, geneIndex,
                                              parallel = TRUE, cores = 4,
                                              asFrame = FALSE) {
    if (parallel) {
        cores <- if (parallel::detectCores() > cores) cores else parallel::detectCores()
    }
    os_type <- .Platform$OS.type
    if (parallel && os_type == "unix") {
        res <- unlist(
            parallel::mclapply(geneIndex, diracClassificationRate,
                expression = expression, phenotype1 = phenotype1,
                phenotype2 = phenotype2,
                mc.cores = cores
            )
        )
    } else if (parallel && os_type == "windows") {
        # make cluster
        cl <- parallel::makeCluster(cores)
        # Export needed functions to cluster
        parallel::clusterExport(cl, list(
            "dirac.rank_matching_score",
            "dirac.rank_vector",
            "dirac.rank_template",
            "dirac.rank_matrix",
            "dirac.rank_matching_score.vector",
            "dirac.rank_difference_score",
            "dirac.classification_rate"
        ))
        # Perform calculation
        res <- tryCatch(
            expr = {
                unlist(parallel::parLapply(cl, geneIndex, diracClassificationRate,
                    expression = expression,
                    phenotype1 = phenotype1,
                    phenotype2 = phenotype2
                ))
            },
            # Whether there is an error or not, stop the cluster
            finally = {
                parallel::stopCluster(cl)
            }
        )
    } else if (parallel) {
        stop("Unsupported OS for parallel operation")
    } else {
        res <- unlist(
            lapply(geneIndex, diracClassificationRate,
                expression = expression,
                phenotype1 = phenotype1, phenotype2 = phenotype2
            )
        )
    }
    if (!asFrame) {
        # Rename the results based on the names of the gene_index list arg
        names(res) <- names(geneIndex)
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
             gene_index,
             parallel = TRUE,
             cores = 4) {
        names(which.max(
            diracClassificationRateCompare(
                expression,
                phenotype1,
                phenotype2,
                gene_index,
                parallel,
                cores,
                as.frame = FALSE
            )
        ))
    }
