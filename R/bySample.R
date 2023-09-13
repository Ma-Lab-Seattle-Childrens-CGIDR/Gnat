## Contains wrapper function for the bySampleEntropy Calculations


#' Compute the Sample wise entropy using DIRAC, RACE, or CRANE
#'
#' The return value depends on the type of the provided expression data. If
#' the data is a SummarizedExperiment object, that object will be returned
#' with two additional colData columns, one for each phenotype. These columns
#' will include the entropy of the gene network for that sample within its
#' phenotype (samples not in that phenotype will have an NA). If the input
#' is a matrix, the retun will be a list named according to the phenotypeNames
#' argument, with the elements being numeric vectors representing the
#' sample wise entropy of the samples within the two provided phenotypes.
#'
#' @param expression Gene expression data, can be either a SummarizedExperiment
#'      object, or an array type object which implements [], and also
#'      rownames/colnames if phenotypes or geneNetworks use character vector
#'      indices
#' @param method Method to use to calculate the rank entropy, can be DIRAC,
#'      RACE, CRANE, or INFER
#' @param phenotype1,phenotype2 Index vectors for the locations within
#'      expression of the samples in each phenotype, can be a numeric,
#'      character, or logical vector
#' @param geneNetwork Index vector for the locations within expression of
#'      the genes in the network, can be a numeric, character, or a logical
#'      vector
#' @param phenotypeNames Character vector of names of the phenotypes
#' @param geneNetworkName Name of the gene network
#' @param assayName Name of the assay in SummarizedExperiment object containing
#'      the gene expression data, (if SummarizedExperiment object is used as
#'      input)
#'
#' @return Either the sample wise entropy values as a named list of vectors,
#'      or a SummarizedExperiment object with additional colData for the
#'      sample wise entropy depending on the type of the input expression
#' @export
#'
#' @examples
#' @include seIntegration.R
bySampleEntropy <- function(expression,
                            method=c("DIRAC", "CRANE", "RACE"),
                            phenotype1,
                            phenotype2,
                            geneNetwork,
                            phenotypeNames=c("p1", "p2"),
                            geneNetworkName="gn",
                            assayName="counts"){
    # These functions take in filtered expression
    methodFunctionList <- list(
        DIRAC = diracSampleScore,
        CRANE = craneSampleScore,
        RACE = raceSampleScore
    )
    method <- match.arg(method)
    methodFunction <- methodFunctionList[[method]]
    if(methods::is(expression, "SummarizedExperiment")){
        return(
            .wrapMethodBySample(seObject=expression,
                                method=methodFunction,
                                phenotype1=phenotype1,
                                phenotype2=phenotype2,
                                geneNetwork=geneNetwork,
                                phenotypeNames=phenotypeNames,
                                geneNetworkName=geneNetworkName,
                                assayName=assayName)
        )
    }
    p1Index <- .convertPhenotype(phenotype=phenotype1,
                                 expressionMatrix = expression)
    p2Index <- .convertPhenotype(phenotype=phenotype2,
                                 expressionMatrix = expression)
    geneIndex <- .convertNetwork(geneNetwork = geneNetwork,
                                 expressionMatrix = expression)
    p1SampleEntropy <- methodFunction(expression[geneIndex, p1Index])
    p2SampleEntropy <- methodFunction(expression[geneIndex, p2Index])
    retList <- list(p1SampleEntropy, p2SampleEntropy)
    names(retlist) <- phenotypeNames
    retList
}
