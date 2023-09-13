## File includes the methods for performing the rank entropy phenotype
## comparison, and the sample wise entropy calculation with
## SummarizedExperiment objects


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
#' @include crane.R
seCrane <- function(seObject, phenotype1, phenotype2, geneNetworkList,
                    assayName = "counts", bootstrapIterations=1000,
                    replace=TRUE, asFrame=TRUE, BPPARAM=bpparam()){
    .wrapMethodComparePhenotypes(method=craneComparePhenotypes,
                seObject = seObject,
                phenotype1 = phenotype1,
                phenotype2 = phenotype2,
                geneNetworkList = geneNetworkList,
                assayName = assayName,
                bootstrapIterations = bootstrapIterations,
                replace = replace,
                asFrame = asFrame,
                BPPARAM=BPPARAM)
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
#' @include dirac.R
seDirac <- function(seObject, phenotype1, phenotype2, geneNetworkList,
                    assayName = "counts", bootstrapIterations=1000,
                    replace=TRUE, asFrame=TRUE, BPPARAM=bpparam()){
    .wrapMethodComparePhenotypes(method=diracComparePhenotypes,
                seObject = seObject,
                phenotype1 = phenotype1,
                phenotype2 = phenotype2,
                geneNetworkList = geneNetworkList,
                assayName = assayName,
                bootstrapIterations = bootstrapIterations,
                replace = replace,
                asFrame = asFrame,
                BPPARAM=BPPARAM)
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
#' @include race.R
seRace <- function(seObject, phenotype1, phenotype2, geneNetworkList,
                   assayName = "counts", bootstrapIterations=1000,
                   replace=TRUE, asFrame=TRUE, BPPARAM=bpparam()){
    .wrapMethodComparePhenotypes(method=raceComparePhenotypes,
                seObject = seObject,
                phenotype1 = phenotype1,
                phenotype2 = phenotype2,
                geneNetworkList = geneNetworkList,
                assayName = assayName,
                bootstrapIterations = bootstrapIterations,
                replace = replace,
                asFrame = asFrame,
                BPPARAM=BPPARAM)
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
#' @include infer.R
seInfer <- function(seObject, phenotype1, phenotype2, geneNetworkList,
                    assayName = "counts", bootstrapIterations=1000,
                    replace=TRUE, asFrame=TRUE, BPPARAM=bpparam()){
    .wrapMethodComparePhenotypes(method=craneComparePhenotypes,
                seObject = seObject,
                phenotype1 = phenotype1,
                phenotype2 = phenotype2,
                geneNetworkList = geneNetworkList,
                assayName = assayName,
                bootstrapIterations = bootstrapIterations,
                replace = replace,
                asFrame = asFrame,
                BPPARAM=BPPARAM)
}



# By Sample Functions -----------------------------------------------------

#' Compute the Sample Wise entropy using DIRAC for a SummarizedExperiment
#'
#' @param seObject SummarizedExperiment object used as input, should contain
#'      an assay with gene expression data, and genes as rownames
#' @param phenotype1,phenotype2 Index vector describing the location of the
#'      two phenotypes in the seObject, can be a numeric, character, or logical
#'      vector
#' @param phenotypeNames Character vector with names for the phenotypes,
#'      phenotype1 should be the first entry and phenotype2 should be the
#'      second
#' @param geneNetwork Index vector describing the location of the
#'      gene network in the seObject, can be a numeric, character, or logical
#'      vector
#' @param geneNetworkName Name for the gene network
#' @param assayName Name of the assay containing the gene expression data
#'
#' @return SummarizedExperiment object with two columns added to colData,
#'      one for the entropy of each of the phenotypes
#' @export
#'
#' @examples
seDiracBySample <- function(seObject, phenotype1, phenotype2,
                            phenotypeNames=c("phenotype1", "phenotype2"),
                            geneNetwork, geneNetworkName="gn",
                            assayName="counts"){
    .wrapMethodBySample(seObject = seObject, method=diracSampleScore,
                        methodName = "DIRAC", phenotype1=phenotype1,
                        phenotype2=phenotype2, phenotypeNames = phenotypeNames,
                        geneNetwork=geneNetwork,
                        geneNetworkName=geneNetworkName, assayName=assayName)
}

#' Compute the Sample Wise entropy using RACE for a SummarizedExperiment
#'
#' @param seObject SummarizedExperiment object used as input, should contain
#'      an assay with gene expression data, and genes as rownames
#' @param phenotype1,phenotype2 Index vector describing the location of the
#'      two phenotypes in the seObject, can be a numeric, character, or logical
#'      vector
#' @param phenotypeNames Character vector with names for the phenotypes,
#'      phenotype1 should be the first entry and phenotype2 should be the
#'      second
#' @param geneNetwork Index vector describing the location of the
#'      gene network in the seObject, can be a numeric, character, or logical
#'      vector
#' @param geneNetworkName Name for the gene network
#' @param assayName Name of the assay containing the gene expression data
#'
#' @return SummarizedExperiment object with two columns added to colData,
#'      one for the entropy of each of the phenotypes
#' @export
#'
#' @examples
seRaceBySample <- function(seObject, phenotype1, phenotype2,
                            phenotypeNames=c("phenotype1", "phenotype2"),
                            geneNetwork, geneNetworkName="gn",
                            assayName="counts"){
    .wrapMethodBySample(seObject = seObject, method=raceSampleScore,
                        methodName = "RACE", phenotype1=phenotype1,
                        phenotype2=phenotype2, phenotypeNames = phenotypeNames,
                        geneNetwork=geneNetwork,
                        geneNetworkName=geneNetworkName, assayName=assayName)
}

#' Compute the Sample Wise entropy using CRANE for a SummarizedExperiment
#'
#' @param seObject SummarizedExperiment object used as input, should contain
#'      an assay with gene expression data, and genes as rownames
#' @param phenotype1,phenotype2 Index vector describing the location of the
#'      two phenotypes in the seObject, can be a numeric, character, or logical
#'      vector
#' @param phenotypeNames Character vector with names for the phenotypes,
#'      phenotype1 should be the first entry and phenotype2 should be the
#'      second
#' @param geneNetwork Index vector describing the location of the
#'      gene network in the seObject, can be a numeric, character, or logical
#'      vector
#' @param geneNetworkName Name for the gene network
#' @param assayName Name of the assay containing the gene expression data
#'
#' @return SummarizedExperiment object with two columns added to colData,
#'      one for the entropy of each of the phenotypes
#' @export
#'
#' @examples
seCraneBySample <- function(seObject, phenotype1, phenotype2,
                           phenotypeNames=c("phenotype1", "phenotype2"),
                           geneNetwork, geneNetworkName="gn",
                           assayName="counts"){
    .wrapMethodBySample(seObject = seObject, method=craneSampleScore,
                        methodName = "CRANE", phenotype1=phenotype1,
                        phenotype2=phenotype2, phenotypeNames = phenotypeNames,
                        geneNetwork=geneNetwork,
                        geneNetworkName=geneNetworkName, assayName=assayName)
}




# Helper Functions --------------------------------------------------------

.wrapMethodBySample <- function(seObject, method, methodName, phenotype1,
                                phenotype2, phenotypeNames, geneNetwork,
                                geneNetworkName, assayName="counts"){
    # Ensure that seObject is a SummarizedExperiment object
    if(!methods::is(seObject, "SummarizedExperiment")){
        stop("First argument must be a SummarizedExperiment object or subclass")
    }
    # Copy colData, it will be modified, and then the colData will be set
    # to the modified value
    columnData <- SummarizedExperiment::colData(seObject)
    # Get dims of the seObject
    dims <- SummarizedExperiment::dim(seObject)
    p1Index <- .checkPhenotype(seObject=seObject, phenotype=phenotype1)
    p2Index <- .checkPhenotype(seObject=seObject, phenotype=phenotype2)
    geneIndex <- .checkNetwork(network = geneNetwork,
                               numRows = dims[[1]],
                               rowNames =
                                   SummarizedExperiment::rownames(seObject))
    # Create names for the additional columns being added to the colData
    p1ColName <- paste(methodName, geneNetworkName, phenotypeNames[1],
                       "entropy", sep="_")
    p2ColName <- paste(methodName, geneNetworkName, phenotypeNames[2],
                       "entropy", sep="_")
    # Get expression values for the two phenotypes
    p1FilteredExpression <- assays(seObject)[[assayName]][geneIndex, p1Index]
    p2FilteredExpression <- assays(seObject)[[assayName]][geneIndex, p2Index]
    # Get entropy values
    p1Entropy <- method(p1FilteredExpression)
    p2Entropy <- method(p2FilteredExpression)
    # Update the coldata to include the columns for entropy
    columnData[p1ColName] <- NA
    columnData[p1Index, p1ColName] <- p1Entropy
    columnData[p2ColName] <- NA
    columnData[p2Index, p2ColName] <- p2Entropy
    # Update the SummarizedExperiment object with
    colData(seObject) <- columnData
    seObject
}

.wrapMethodComparePhenotypes <- function(
        seObject,
        method,
        phenotype1,
        phenotype2,
        geneNetworkList,
        assayName = "counts",
        bootstrapIterations = 1000,
        replace = TRUE,
        asFrame = TRUE,
        BPPARAM=bpparam()
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
        asFrame = asFrame,
        BPPARAM = bpparam()
    )
    res
}



