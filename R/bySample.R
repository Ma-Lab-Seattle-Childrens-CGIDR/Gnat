## Contains wrapper function for the bySampleEntropy Calculations


#' Compute the Sample wise entropy using DIRAC, RACE, or CRANE
#'
#' The return value depends on the type of the provided expression data. If
#' the data is a SummarizedExperiment object, that object will be returned
#' with an additional colData column. The column will include the entropy of
#' the gene network for that sample within its phenotype (samples not in that
#' phenotype will have an NA). If the input is an matrix, the return will be a
#' numeric vector representing the sample wise entropy of the samples within
#' the provided phenotype.
#'
#' @param expression Gene expression data, can be either a SummarizedExperiment
#'      object, or an array type object which implements [], and also
#'      rownames/colnames if phenotypes or geneNetworks use character vector
#'      indices, as well as length and dim functions if a logical vector
#'      input is desired
#' @param method Method to use to calculate the rank entropy, can be DIRAC,
#'      RACE, CRANE, or INFER
#' @param phenotype Index vector for samples within the phenotype, can be a
#'      numeric, character, or logical vector
#' @param geneNetwork Index vector for the locations within expression of
#'      the genes in the network, can be a numeric, character, or a logical
#'      vector
#' @param colDataName Name to give new column in the colData of the
#'      SummarizedExperiment object which will hold the sample entropy data
#' @param assayName Name of the assay in SummarizedExperiment object containing
#'      the gene expression data, (if SummarizedExperiment object is used as
#'      input)
#'
#' @return Either the sample wise entropy values as a vector,
#'      or a SummarizedExperiment object with additional colData for the
#'      sample wise entropy depending on the type of the input expression
#' @export
#'
#' @examples
#' @include seIntegration.R
bySampleEntropy <- function(expression,
                            method=c("DIRAC", "CRANE", "RACE"),
                            phenotype,
                            geneNetwork,
                            colDataName=paste(method,"sampleEntropy", sep="_"),
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
                                methodName=method,
                                phenotype=phenotype,
                                geneNetwork=geneNetwork,
                                colName=colDataName,
                                assayName=assayName)
        )
    }
    pIndex <- .convertPhenotype(phenotype=phenotype,
                                 expressionMatrix = expression)
    geneIndex <- .convertNetwork(geneNetwork = geneNetwork,
                                 expressionMatrix = expression)
    sampleEntropy <- methodFunction(expression[geneIndex, pIndex])
    sampleEntropy
}
