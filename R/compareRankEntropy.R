# File includes main method for comparing two phenotypes using the rank
# entropy methods

# Main Function -----------------------------------------------------------
#' Compare Phenotypes using Rank Entropy Methods
#'
#' Compare the rank entropy between two phenotypes for a given gene network,
#' using either DIRAC, CRANE, RACE, or INFER. The expression data can be an
#' array type object which implements [], if a logical vector input is
#' desired it will also need to implement length and dim, and if a character
#' vector input is desired it will also have to implement rownames, and colnames
#' functions. Index vectors (for phenotypes and gene networks) can be
#' numeric (representing the numeric index of the gene/sample), character
#' (representing the row or column name), or logical (TRUE for samples/genes
#' in the phenotype/network, and FALSE otherwise).
#'
#' @param expression Gene expression data, can be either a SummarizedExperiment
#'      object, or an array type object, see details for more information
#' @param method Method to use to calculate the rank entropy, can be DIRAC,
#'      RACE, CRANE, or INFER
#' @param phenotype1,phenotype2 Index vectors for the locations within
#'      expression of the samples within each phenotype, can be a numeric,
#'      character, or logical vector
#' @param geneNetworkList (Named) list of index vectors, each entry must be
#'      a character, numeric, or logical vector describing the locations
#'      within expression of the genes within the network
#' @param assayName Name of assay within SummarizedExperiment object holding
#'      gene expression data
#' @param bootstrapIterations Number of iterations to perform for boostrapping
#'      the null distribution for computing the p-value
#' @param replace Whether sampling during bootstrapping should be done with
#'      replacement
#' @param asFrame Whether the return value should be a data.frame, otherwise
#'      a named list is returned
#' @param BPPARAM Biocparallel parameter for choosing a parallel backend
#'      to use for bootstrapping computations
#'
#' @return Either a data.frame, or a named list depending on the asFrame
#'    parameter. In the data.frame form, with columns for geneNetwork, p1Score,
#'    p2Score, absoluteDifference, and pValue. The list return is named
#'    according to the names in the geneNetworkList parameter, with each value
#'    being a named list with p1Score, p2Score, absoluteDifference, and pValue.
#' @export
#'
#' @examples
#' @include seIntegration.R
compareRankEntropy <- function(expression,
                              method =
                                  c("DIRAC", "RACE", "CRANE","INFER"),
                              phenotype1,
                              phenotype2,
                              geneNetworkList,
                              assayName = "counts",
                              bootstrapIterations = 1000,
                              replace = TRUE,
                              asFrame = TRUE,
                              BPPARAM = bpparam()) {
    methodFunctionList <- list(
        DIRAC=diracComparePhenotypes,
        CRANE=craneComparePhenotypes,
        INFER=inferComparePhenotypes,
        RACE=raceComparePhenotypes
    )
    method <- match.arg(method)
    methodFunction <- methodFunctionList[[method]]
    if(methods::is(expression, "SummarizedExperiment")){
        return(
            .wrapMethodComparePhenotypes(method=methodFunction,
                                         seObject = expression,
                                         phenotype1 = phenotype1,
                                         phenotype2 = phenotype2,
                                         geneNetworkList = geneNetworkList,
                                         assayName = assayName,
                                         bootstrapIterations =
                                         bootstrapIterations,
                                         replace = replace,
                                         asFrame = asFrame,
                                         BPPARAM=BPPARAM))
    }
    p1Index <- .convertPhenotype(phenotype = phenotype1,
                                 expressionMatrix = expression)
    p2Index <- .convertPhenotype(phenotype = phenotype2,
                                 expressionMatrix = expression)
    geneNetworkList <- .convertNetworkList(geneNetworkList=geneNetworkList,
                                           expressionMatrix=expression)
    methodFunction(expression=expression,
                   geneNetworkList=geneNetworkList,
                   phenotype1=p1Index,
                   phenotype2=p2Index,
                   bootstrapIterations=bootstrapIterations,
                   replace=replace,
                   asFrame=asFrame,
                   BPPARAM=BPPARAM)
}

