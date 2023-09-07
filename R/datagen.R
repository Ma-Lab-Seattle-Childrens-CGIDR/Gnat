# Main Generating Functions -----------------------------------------------

#' Generates gene expression data with ordered sub network
#'
#' @param ngenesOrdered Integer, number of genes in the ordered sub network.
#' @param ngenesTotal Integer, total number of genes desired.
#' @param nsamplesOrdered Integer, number of samples in the ordered phenotype.
#' @param nsamplesTotal Integer, total number of samples
#' @param dist Function, distribution to use for generating data. It is passed
#'    ngenes (for number of values to generate), and any other keyword
#'    arguments. Must return a numeric vector.
#' @param reorderGenes Boolean, whether the order of genes within the ordered
#'   sub network should be shuffled. If TRUE, the order of the genes in the
#'   network will be randomly reshuffled. If FALSE, the genes within the network
#'   will always have their expression values increase with increasing row
#'   index.
#' @param ... Keyword arguments, passed to the dist function
#'
#' @return Named list with:
#'    expression: Numeric matrix of gene expression values with an ordered
#'      sub network. There will be ngenes.total rows representing genes, and
#'      nsamples.total columns representing samples.
#'    orderedGenes: Integer vector with the indices of the genes within
#'      the ordered sub network. These will be in ascending order, so the gene
#'      with the lowest expression values will be in the row represented by the
#'      first index, and the gene with the highest expression will be
#'      represented by the last index.
#'    orderedSamples: Integer vector containing the indices of the samples
#'      within the ordered phenotype.
#' @export
#'
#' @examples
datagenGenerate <- function(ngenesOrdered, ngenesTotal, nsamplesOrdered,
                             nsamplesTotal, dist, reorderGenes = TRUE, ...) {
    expressionMatrix <- datagenUnorderedMatrix(
        ngenesTotal,
        nsamplesTotal,
        dist, ...
    )
    orderedExpression <- datagenOrderedMatrix(
        ngenesOrdered, nsamplesOrdered,
        dist, ...
    )
    orderedGenes <- sample(ngenesTotal, ngenesOrdered, replace = FALSE)
    ordered.samples <- sample(nsamplesTotal, nsamplesOrdered, replace = FALSE)
    if (!reorderGenes) {
        orderedGenes <- sort(orderedGenes)
    }
    for (i in 1:ngenesOrdered) {
        expressionMatrix[orderedGenes[i], ordered.samples] <-
            orderedExpression[i, ]
    }
    list(
        expression = expressionMatrix,
        orderedGenes = orderedGenes,
        orderedSamples = sort(ordered.samples)
    )
}



# Summarized Experiment Generator -----------------------------------------
#' Generates a summarized experiment object with gene expression data, which
#' has an ordered subnetwork
#'
#' @param ngenesOrdered Integer, number of genes in the ordered sub network.
#' @param ngenesTotal Integer, total number of genes desired.
#' @param nsamplesOrdered Integer, number of samples in the ordered phenotype.
#' @param nsamplesTotal Integer, total number of samples
#' @param dist Function, distribution to use for generating data. It is passed
#'    ngenes (for number of values to generate), and any other keyword
#'    arguments. Must return a numeric vector.
#' @param reorderGenes Boolean, whether the order of genes within the ordered
#'   sub network should be shuffled. If TRUE, the order of the genes in the
#'   network will be randomly reshuffled. If FALSE, the genes within the network
#'   will always have their expression values increase with increasing row
#'   index.
#' @param ... Keyword arguments, passed to the dist function
#'
#' @return Named list with:
#'      seObject: The generated summarized experiment, with
#'          colData:
#'              sampleNames, phenotypeNum, and phenotypeStr
#'          rowData:
#'              geneNames, networkNum, netowrkStr
#'      orderedGenes: A numeric vector with the indices of the ordered genes
#'      orderedSamples: A numeric vector with the indices of the ordered samples
#'      orderedGeneNames: A character vector with the ordered gene names
#'      orderedSampleNames: A character vector with the ordered sample names
#'      unorderedGenes: A numeric vector with the indices of the unordered genes
#'      unorderedSamples: A numeric vector with the indices of the unorrdered
#'          samples
#'      unorderedGeneNames: A character vector with the unordered gene names
#'      unorderedSampleNames: A character vector with the unordered sample
#'          names
#' @export
#'
#' @examples
generateSummarizedExperiment <- function(ngenesOrdered, ngenesTotal,
                                         nsamplesOrdered, nsamplesTotal,
                                         dist, reorderGenes=FALSE,
                                         genePrefix="g_",
                                         samplePrefix="s_", ...){
    expData <- datagenGenerate(ngenesOrdered = ngenesOrdered,
                               ngenesTotal = ngenesTotal,
                               nsamplesOrdered = nsamplesOrdered,
                               nsamplesTotal = nsamplesTotal,
                               dist=dist, reorderGenes = reorderGenes, ...)
    expression <- expData$expression
    orderedGenes <- expData$orderedGenes
    unorderedGenes <- setdiff(1:ngenesTotal, orderedGenes)
    orderedSamples <- expData$orderedSamples
    unorderedSamples <- setdiff(1:nsamplesTotal, orderedSamples)
    geneNames <- vapply(1:ngenesTotal,
                         function(x) paste(genePrefix, x, sep=""),
                         character(1))
    sampleNames <- vapply(1:nsamplesTotal,
                           function(x) paste(samplePrefix, x, sep=""),
                           character(1))
    phenotypeNum <- numeric(nsamplesTotal)
    phenotypeNum[orderedSamples] <-  1
    phenotypeNum[unorderedSamples] <- 2
    phenotypeStr <- character(nsamplesTotal)
    phenotypeStr[orderedSamples] <- "one"
    phenotypeStr[unorderedSamples] <- "two"
    columnDataFrame <- S4Vectors::DataFrame(sampleNames, phenotypeNum,
                                            phenotypeStr,
                                            row.names = sampleNames)
    networkNum <- numeric(ngenesTotal)
    networkNum[orderedGenes] <- 1
    networkNum[unorderedGenes] <- 2
    networkStr <- character(ngenesTotal)
    networkStr[orderedGenes] <- "one"
    networkStr[unorderedGenes] <- "two"
    rowDataFrame <- S4Vectors::DataFrame(geneNames, networkNum, networkStr,
                                         row.names=geneNames)
    summarizedExperimentResult <-
        SummarizedExperiment::SummarizedExperiment(
            assays=list(geneExpression = expression),
            rowData = rowDataFrame, colData = columnDataFrame
        )
    orderedGeneNames <- geneNames[orderedGenes]
    orderedSampleNames <- sampleNames[orderedSamples]
    unorderedGeneNames <- geneNames[unorderedGenes]
    unorderedSampleNames <- sampleNames[unorderedSamples]
    list(seObject=summarizedExperimentResult, orderedGenes=orderedGenes,
         orderedSamples=orderedSamples, orderedGeneNames=orderedGeneNames,
         orderedSampleNames=orderedSampleNames, unorderedGenes=unorderedGenes,
         unorderedSamples=unorderedSamples,
         unorderedGeneNames=unorderedGeneNames,
         unorderedSampleNames=unorderedSampleNames)
}

# Main Noise Functions ----------------------------------------------------

#' Perform nswaps on the expression data, and return the result
#'
#' For each swap, a random sample, and two random genes are chosen. The values
#' for gene expression for those two genes, in that sample are then swapped.
#'
#' @param expression Numeric matrix of gene expression data, rows represent
#'    genes, and columns represent samples
#' @param nswaps Integer, number of swaps to perform
#'
#' @return Numeric matrix, expression data but with nswaps performed
#' @export
#'
#' @examples
datagenSwap <- function(expression, nswaps) {
    nsamples <- ncol(expression)
    ngenes <- nrow(expression)
    sample_list <- sample(nsamples, size = nswaps, replace = TRUE)
    for (i in 1:nswaps) {
        genes <- sample(ngenes, size = 2, replace = FALSE)
        temp <- expression[genes[1], sample_list[i]]
        expression[genes[1], sample_list[i]] <- expression[genes[2], sample_list[i]]
        expression[genes[2], sample_list[i]] <- temp
    }
    expression
}

#' Adds noise from a distribution to the expression data
#'
#' @param expression Numeric matrix, gene expression data. Rows are genes,
#'    columns are samples.
#' @param dist Function, distribution function to generate random numbers to be
#'    added as noise to expression. Must take number of random numbers to
#'    generate as first argument, and is also passed other keyword arguments.
#' @param ... Keyword arguments, passed to dist function.
#'
#' @return Numeric matrix, expression data with noise added.
#' @export
#'
#' @examples
datagenNoise <- function(expression, dist, ...) {
    ngenes <- nrow(expression)
    nsamples <- ncol(expression)
    noiseMatrix <- matrix(dist(ngenes * nsamples, ...), nrow = ngenes,
                          ncol = nsamples)
    expression + noiseMatrix
}

#' Dropout a proportion of expression data
#'
#' Sets dropout.rate proportion of data to 0.
#'
#' @param expression Numeric matrix, represents gene expression data with genes
#'    as rows, and samples as columns.
#' @param dropoutRate Double, representing the proportion of values to set to
#'    0, 0<=dropout.rate<=1.
#'
#' @return Numeric matrix, expression data with dropout.rate proportion of
#'    values set to 0.
#' @export
#'
#' @examples
datagenDropout <- function(expression, droupoutRate) {
    ngenes <- nrow(expression)
    nsamples <- ncol(expression)
    nvaluesDrop <- floor(droupoutRate * (ngenes * nsamples))
    dropIndex <- sample(ngenes * nsamples, size = nvaluesDrop, replace = FALSE)
    dropVector <- logical(ngenes * nsamples)
    dropVector[dropIndex] <- TRUE
    expression[dropVector] <- 0
    expression
}


# Helper Functions --------------------------------------------------------

#' Generate ordered expression vector from a provided distribution function
#'
#' @param ngenes Integer, number of genes to generate expression values for.
#' @param dist Function, distribution function (rnorm, rnbinom, etc.), which
#'    will be passed the ngenes, and other keyword arguments, and should
#'    return a numeric vector of length ngenes.
#' @param ... Keyword arguments passed to distribution function
#'
#' @return Numeric vector of ordered gene expression values
#' @export
#'
#' @examples
datagenOrderedVector <- function(ngenes, dist, ...) {
    expressionUnordered <- dist(ngenes, ...)
    sort(expressionUnordered, decreasing = FALSE)
}



#' Generate ordered expression matrix from a provided distribution function
#'
#' @param ngenes Integer, number of genes to generate expression values for.
#' @param nsamples Integer, number of samples to generate expression values for.
#' @param dist Function, distribution function (rnorm, rnbinom, etc.), which
#'    will be passed the ngenes, and other keyword arguments, and should
#'    return a numeric vector of length ngenes.
#' @param ... Keyword arguments passed to distribution function
#'
#' @return Numeric matrix of ordered gene expression values, with rows
#'    representing genes, and columns samples
#' @export
#'
#' @examples
datagenOrderedMatrix <- function(ngenes, nsamples, dist, ...) {
    # Create an empty expression matrix
    expressionMatrix <- matrix(0, nrow = ngenes, ncol = nsamples)
    # Generate the expression for each sample of the ordered phenotype
    for (i in 1:nsamples) {
        expressionMatrix[, i] <- datagenOrderedVector(ngenes, dist, ...)
    }
    expressionMatrix
}


#' Generate unordered expression vector from a provided distribution function
#'
#' @param ngenes Integer, number of genes to generate expression values for.
#' @param dist Function, distribution function (rnorm, rnbinom, etc.), which
#'    will be passed the ngenes, and other keyword arguments, and should
#'    return a numeric vector of length ngenes.
#' @param ... Keyword arguments passed to distribution function
#'
#' @return Numeric vector of unordered gene expression values
#' @export
#'
#' @examples
datagenUnorderedVector <- function(ngenes, dist, ...) {
    dist(ngenes, ...)
}

#' Generate unordered expression matrix from a provided distribution function
#'
#' @param ngenes Integer, number of genes to generate expression values for.
#' @param nsamples Integer, number of samples to generate expression values for.
#' @param dist Function, distribution function (rnorm, rnbinom, etc.), which
#'    will be passed the ngenes, and other keyword arguments, and should
#'    return a numeric vector of length ngenes.
#' @param ... Keyword arguments passed to distribution function
#'
#' @return Numeric matrix of unordered gene expression values, with rows
#'    representing genes, and columns samples
#' @export
#'
#' @examples
datagenUnorderedMatrix <- function(ngenes, nsamples, dist, ...) {
    # Create an empty expression matrix
    expressionMatrix <- matrix(0, nrow = ngenes, ncol = nsamples)
    # Generate the expression for each sample of the ordered phenotype
    for (i in 1:nsamples) {
        expressionMatrix[, i] <- datagenUnorderedVector(ngenes, dist, ...)
    }
    expressionMatrix
}
