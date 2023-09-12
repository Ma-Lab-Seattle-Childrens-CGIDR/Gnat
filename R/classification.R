# Ensure needed methods are loaded first
#' @include crane.R race.R dirac.R infer.R




# Create Classifier -------------------------------------------------------
createClassifier <- function(expression, phenotype1, phenotype2, geneIndex,
                             method=c("DIRAC","CRANE","RACE"),
                             assay="counts"){
    method <- match.args(method)

}






# Classifier Methods ------------------------------------------------------

# Phenotype Classifier ----------------------------------------------------------

#' Creates a phenotype classifier based on DIRAC entropy
#'
#' @param expression Numeric matrix of gene expression values, rows are genes
#'    columns are samples
#' @param phenotype1,phenotype2 Numeric vectors with indices for each of the
#'    phenotypes
#' @param gene_index Numeric vector of indices for the genes within the gene
#'    network being used for classification
#' @param names optional Character vector with the desired names for the
#'      phenotypes
#'
#' @returns A function with expression argument, which will return a logical
#'      vector, with TRUE for the samples predicted to be phenotype1 and FALSE
#'      for the samples predicted to be phenotype2
#' @examples
#' # example code
#' @export
#' @include dirac.R
diracEntropyClassifier <- function(expression,
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
#' @param gene_index Numeric vector of indices for the genes within the gene
#'    network being used for classification
#' @param names optional Character vector with the desired names for the
#'      phenotypes
#'
#' @returns A function with expression argument, which will return a logical
#'      vector, with TRUE for the samples predicted to be phenotype1 and FALSE
#'      for the samples predicted to be phenotype2
#' @examples
#' # example code
#' @export
#' @include race.R
raceEntropyClassifier <- function(expression,
                                  phenotype1,
                                  phenotype2,
                                  geneIndex){

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
#' @param gene_index Numeric vector of indices for the genes within the gene
#'    network being used for classification
#' @param names optional Character vector with the desired names for the
#'      phenotypes
#'
#' @returns A function with expression argument, which will return a logical
#'      vector, with TRUE for the samples predicted to be phenotype1 and FALSE
#'      for the samples predicted to be phenotype2
#' @examples
#' # example code
#' @export
#' @include crane.R
craneEntropyClassifier <- function(expression,
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
        p1Distance <- sqrt(colSums((p1RankMatrix-p1RankCentroid)^2))
        p2Distance <- sqrt(colSums((p2RankMatrix-p1RankCentroid)^2))
        distDiff <- p2Distance - p1Distance
        distDiff>=0
    }
}
