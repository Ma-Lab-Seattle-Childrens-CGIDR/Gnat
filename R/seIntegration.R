## Include the required packages before loading this one
#' @include race.R crane.R dirac.R infer.R
## Import needed Summarized Experiment Functions
# Phenotype Comparison Function -------------------------------------------


#' Perform Phenotype Comparison using the CRANE method on a Summarized
#' Experiment Object
#'
#' @param seObject SummarizedExperiment object used as input, should contain
#'      an assay with gene expression data, and genes as rownames
#' @param phenotype1,phenotype2 Index vector describing the location of the
#'      two phenotypes in the seObject, can be a numeric, character, or logical
#'      vector
#' @param geneNetworkList (Named) list of gene network index vectors, each
#'      describing the indices of a gene network, these can be integer,
#'      character, or logical vectors
#' @param assayName Name of the assay containing the gene expression data
#' @param bootstrapIterations Number of iterations to use for bootstrapping
#'      the null distribution
#' @param replace Logical determining if the permutation sampling done as
#'      part of the bootstrapping should be done with or without replacement
#'      (TRUE for replacement, FALSE for no replacement)
#' @param asFrame Logical determining if the return value should be in the form
#'      of a data.frame
#'
#' @return Either a data.frame, or a named list depending on the asFrame
#'    parameter. The data.frame form will have columns for geneNetwork, p1Score,
#'    p2Score, absoluteDifference, and pValue. The list return is named
#'    according to the names in the geneNetworkList parameter, with each value
#'    being a named list with p1Score, p2Score, absoluteDifference, and pValue.
#' @export
#'
#' @examples
seCrane <- function(seObject, phenotype1, phenotype2, geneNetworkList,
                    assayName = "count", bootstrapIterations=1000, replace=TRUE,
                    asFrame=TRUE){
    .wrapMethod(method=craneComparePhenotypes,
                seObject = seObject,
                phenotype1 = phenotype1,
                phenotype2 = phenotype2,
                geneNetworkList = geneNetworkList,
                assayName = assayName,
                bootstrapIterations = bootstrapIterations,
                replace = replace,
                asFrame = asFrame)
}

#' Perform Phenotype Comparison using the DIRAC method on a Summarized
#' Experiment Object
#'
#' @param seObject SummarizedExperiment object used as input, should contain
#'      an assay with gene expression data, and genes as rownames
#' @param phenotype1,phenotype2 Index vector describing the location of the
#'      two phenotypes in the seObject, can be a numeric, character, or logical
#'      vector
#' @param geneNetworkList (Named) list of gene network index vectors, each
#'      describing the indices of a gene network, these can be integer,
#'      character, or logical vectors
#' @param assayName Name of the assay containing the gene expression data
#' @param bootstrapIterations Number of iterations to use for bootstrapping
#'      the null distribution
#' @param replace Logical determining if the permutation sampling done as
#'      part of the bootstrapping should be done with or without replacement
#'      (TRUE for replacement, FALSE for no replacement)
#' @param asFrame Logical determining if the return value should be in the form
#'      of a data.frame
#'
#' @return Either a data.frame, or a named list depending on the asFrame
#'    parameter. The data.frame form will have columns for geneNetwork, p1Score,
#'    p2Score, absoluteDifference, and pValue. The list return is named
#'    according to the names in the geneNetworkList parameter, with each value
#'    being a named list with p1Score, p2Score, absoluteDifference, and pValue.
#' @export
#'
#' @examples
seDirac <- function(seObject, phenotype1, phenotype2, geneNetworkList,
                    assayName = "count", bootstrapIterations=1000, replace=TRUE,
                    asFrame=TRUE){
    .wrapMethod(method=diracComparePhenotypes,
                seObject = seObject,
                phenotype1 = phenotype1,
                phenotype2 = phenotype2,
                geneNetworkList = geneNetworkList,
                assayName = assayName,
                bootstrapIterations = bootstrapIterations,
                replace = replace,
                asFrame = asFrame)
}

#' Perform Phenotype Comparison using the RACE method on a Summarized
#' Experiment Object
#'
#' @param seObject SummarizedExperiment object used as input, should contain
#'      an assay with gene expression data, and genes as rownames
#' @param phenotype1,phenotype2 Index vector describing the location of the
#'      two phenotypes in the seObject, can be a numeric, character, or logical
#'      vector
#' @param geneNetworkList (Named) list of gene network index vectors, each
#'      describing the indices of a gene network, these can be integer,
#'      character, or logical vectors
#' @param assayName Name of the assay containing the gene expression data
#' @param bootstrapIterations Number of iterations to use for bootstrapping
#'      the null distribution
#' @param replace Logical determining if the permutation sampling done as
#'      part of the bootstrapping should be done with or without replacement
#'      (TRUE for replacement, FALSE for no replacement)
#' @param asFrame Logical determining if the return value should be in the form
#'      of a data.frame
#'
#' @return Either a data.frame, or a named list depending on the asFrame
#'    parameter. The data.frame form will have columns for geneNetwork, p1Score,
#'    p2Score, absoluteDifference, and pValue. The list return is named
#'    according to the names in the geneNetworkList parameter, with each value
#'    being a named list with p1Score, p2Score, absoluteDifference, and pValue.
#' @export
#'
#' @examples
seRace <- function(seObject, phenotype1, phenotype2, geneNetworkList,
                    assayName = "count", bootstrapIterations=1000, replace=TRUE,
                    asFrame=TRUE){
    .wrapMethod(method=raceComparePhenotypes,
                seObject = seObject,
                phenotype1 = phenotype1,
                phenotype2 = phenotype2,
                geneNetworkList = geneNetworkList,
                assayName = assayName,
                bootstrapIterations = bootstrapIterations,
                replace = replace,
                asFrame = asFrame)
}

#' Perform Phenotype Comparison using the INFER method on a Summarized
#' Experiment Object
#'
#' @param seObject SummarizedExperiment object used as input, should contain
#'      an assay with gene expression data, and genes as rownames
#' @param phenotype1,phenotype2 Index vector describing the location of the
#'      two phenotypes in the seObject, can be a numeric, character, or logical
#'      vector
#' @param geneNetworkList (Named) list of gene network index vectors, each
#'      describing the indices of a gene network, these can be integer,
#'      character, or logical vectors
#' @param assayName Name of the assay containing the gene expression data
#' @param bootstrapIterations Number of iterations to use for bootstrapping
#'      the null distribution
#' @param replace Logical determining if the permutation sampling done as
#'      part of the bootstrapping should be done with or without replacement
#'      (TRUE for replacement, FALSE for no replacement)
#' @param asFrame Logical determining if the return value should be in the form
#'      of a data.frame
#'
#' @return Either a data.frame, or a named list depending on the asFrame
#'    parameter. The data.frame form will have columns for geneNetwork, p1Score,
#'    p2Score, absoluteDifference, and pValue. The list return is named
#'    according to the names in the geneNetworkList parameter, with each value
#'    being a named list with p1Score, p2Score, absoluteDifference, and pValue.
#' @export
#'
#' @examples
seInfer <- function(seObject, phenotype1, phenotype2, geneNetworkList,
                    assayName = "count", bootstrapIterations=1000, replace=TRUE,
                    asFrame=TRUE){
    .wrapMethod(method=craneComparePhenotypes,
                seObject = seObject,
                phenotype1 = phenotype1,
                phenotype2 = phenotype2,
                geneNetworkList = geneNetworkList,
                assayName = assayName,
                bootstrapIterations = bootstrapIterations,
                replace = replace,
                asFrame = asFrame)
}



# By Sample Functions -----------------------------------------------------
seDiracBySample <- function(){

}

seRaceBySample <- function(){

}

seCraneBySample <- function(){

}




# Helper Functions --------------------------------------------------------

.wrapMethod <- function(
        method,
        seObject,
        phenotype1,
        phenotype2,
        geneNetworkList,
        assayName = "count",
        bootstrapIterations = 1000,
        replace = TRUE,
        asFrame = TRUE
){
    args <- .checkArgs(seObject, phenotype1, phenotype2, geneNetworkList)
    # Unpack args
    p1Index <- args$p1Index
    p2Index <- args$p2Index
    geneNetworks <- args$geneNetworks
    # Ensure that assay argument is an assay in the seObject
    if(!(assayName %in% names(assays(seObject)))){
        stop(paste("Assay ", assayName, " not found in Summarized Experiment"))
    }
    # Use method to calculate result
    res <- method(
        SummarizedExperiment::assays(seObject)[[assayName]],
        geneNetworkList = geneNetworks,
        phenotype1 = p1Index,
        phenotype2 = p2Index,
        bootstrapIterations = bootstrapIterations,
        replace = replace,
        asFrame = asFrame
    )
}


.checkArgs <- function(seObject,phenotype1, phenotype2, geneNetworks){
    if(!is(seObject, "SummarizedExperiment")){
        # This will not be called with subclasses (i.e. RangedExperiment)
        stop("First argument must be a SummarizedExperiment object")
    }
    p1Index <- .checkPhenotype(seObject, phenotype1)
    p2Index <- .checkPhenotype(seObject, phenotype2)
    seDims <- SummarizedExperiment::dims(seObject)
    geneNetworks <- ensureNamed(geneNetworks)
    geneNetworks <- lapply(geneNetworks, .checkNetwork, numRows = seDims[1],
                           rownames(seObject))
    list(p1Index=p1Index, p2Index=p2Index, geneNetworks=geneNetworks)
}


.checkPhenotype <- function(seObject, phenotype){
    dimSeObject <- SummarizedExperiment::dim(seObject)
    if(is.logical(phenotype)){
        if(length(phenotype)!=dimSeObject[2]){
            stop(paste("Logical Phenotype Vectors should be ",
                       "the same length as the number of samples, but ",
                       "phenotype has length: ", length(phenotype),
                       " not ", dimSeObject[2]))
        }
        pIndex <- SummarizedExperiment::colnames(seObject)[phenotype]
    } else {
        pIndex <- phenotype
    }
    pIndex
}

.checkNetwork <- function(network, numRows, rowNames){
    if(is.logical(network)){
        if(length(phenotype)!=numRows){
            stop(paste("Logical Gene Network Vectors should be ",
                 "the same length as the number of genes, but ",
                 "network has has length: ", length(network),
                 " not ", numRows))
        }
        gIndex <- rowNames[network]
    } else {
        gIndex <- network
    }
    gIndex
}

