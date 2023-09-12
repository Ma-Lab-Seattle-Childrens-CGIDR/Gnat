# File includes main method for comparing two phenotypes using the rank
# entropy methods


# Main Function -----------------------------------------------------------
#' Compare Phenotypes using Rank Entropy Methods
#'
#' Compare the rank entropy between two phenotypes for a given gene network,
#' using either DIRAC, CRANE, RACE, or INFER
#'
#' @param expression Gene expression data, can be either a SummarizedExperiment
#'      object, or an array type object which implements [], and also
#'      rownames/colnames if phenotypes or geneNetworks use character vector
#'      indices
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
                                  c("CRANE", "DIRAC", "INFER", "RACE"),
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
        RACE=raceComparePhenotypes,
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
                                            asFrame = asFrame))
    }
    p1Index <- .convertPhenotype(phenotype = phenotype1,
                                 expressionMatrix = expression)
    p2Index <- .convertPhenotype(phenotype = phenotype2,
                                 expressionMatrix = expression)
    geneNetworkList <- .convertNetworkList(geneNetworkList=geneNetworkList,
                                           expressionMatrix=expressionMatrix)
    methodFunction(expression=expression,
                   geneNetworkList=geneNetworkList,
                   phenotype1=p1Index,
                   phenotype2=p2Index,
                   bootstrapIterations=bootstrapIterations,
                   replace=replace,
                   asFrame=asFrame,
                   BPPARAM=BPPARAM)
}



# Helper Functions --------------------------------------------------------
.convertPhenotype <- function(phenotype, expressionMatrix){
    numCols <- dim(expressionMatrix)[2]
    if(is.logical(phenotype)){
        if(length(phenotype)==numCols){
            pIndex <- seq(numCols)[phenotype]
            return(pIndex)
        } else {
            stop("If phenotype is logical, must be same length as",
                 " columns of expression matrix")
        }
    }
    if(is.numeric(phenotype)){
        if(!(max(phenotype)<=numCols)){
            stop("Phenotype index out of range")
        }
        return(phenotype)
    }
    columnNames <- colnames(expressionMatrix)
    if(is.character(phenotype)){
        if(all(phenotype %in% columnNames)){
            pIndex <- match(phenotype, colNames)
            return(pIndex)
        } else {
            stop(paste("Couldn't find all phenotype entries in expression ",
                       "matrix. Missing: ",
                       setdiff(phenotype, columnNames)))
        }
    }
    stop(paste("Couldn't Parse Phenotype, should be logical, character, ",
               "or numeric vector"))
}

.convertNetwork <- function(geneNetwork, expressionMatrix){
    numRows <- dim(expressionMatrix)[1]
    if(is.logical(geneNetwork)){
        if(length(geneNetwork)==numRows){
            gIndex <- seq(numRows)[geneNetwork]
            return(gIndex)
        } else {
            stop(paste("If Gene Network is logical, must be same length as",
                       " rows of expression matrix"))
        }
    }
    if(is.numeric(geneNetwork)){
        if(!(max(geneNetwork)<=numRows)){
            stop("Gene Network index out of range")
        }
        return(geneNetwork)
    }
    rowNames <- colnames(expressionMatrix)
    if(is.character(geneNetwork)){
        if(all(geneNetwork %in% rowNames)){
            gIndex <- match(geneNetwork, colNames)
            return(gIndex)
        } else {
            stop(paste("Couldn't find all gene network entries in expression ",
                       "matrix. Missing: ",
                       setdiff(geneNetwork, rowNames)))
        }
    }
    stop(paste("Couldn't parse gene network, should be logical, character, ",
               "or numeric vector"))
}

.convertNetworkList <- function(geneNetworkList, expressionMatrix){
    geneNetworkList <- ensureNamed(geneNetworkList)
    lapply(geneNetworkList, .convertNetwork, expressionMatrix=expressionMatrix)
}

