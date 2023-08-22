# Main Generating Functions -----------------------------------------------

#' Generates gene expression data with ordered sub network
#'
#' @param ngenes.ordered Integer, number of genes in the ordered sub network.
#' @param ngenes.total Integer, total number of genes desired.
#' @param nsamples.ordered Integer, number of samples in the ordered phenotype.
#' @param nsamples.total Integer, total number of samples
#' @param dist Function, distribution to use for generating data. It is passed
#'    ngenes (for number of values to generate), and any other keyword
#'    arguments. Must return a numeric vector.
#' @param reorder.genes Boolean, whether the order of genes within the ordered
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
#'    ordered.genes: Integer vector with the indices of the genes within
#'      the ordered sub network. These will be in ascending order, so the gene
#'      with the lowest expression values will be in the row represented by the
#'      first index, and the gene with the highest expression will be
#'      represented by the last index.
#'    ordered.samples: Integer vector containing the indices of the samples
#'      within the ordered phenotype.
#' @export
#'
#' @examples
datagen.generate <- function(ngenes.ordered, ngenes.total, nsamples.ordered,
                             nsamples.total, dist, reorder.genes = TRUE, ...) {
    expression.matrix <- datagen.unordered.matrix(
        ngenes.total,
        nsamples.total,
        dist, ...
    )
    ordered.expression <- datagen.ordered.matrix(
        ngenes.ordered, nsamples.ordered,
        dist, ...
    )
    ordered.genes <- sample(ngenes.total, ngenes.ordered, replace = FALSE)
    ordered.samples <- sample(nsamples.total, nsamples.ordered, replace = FALSE)
    if (!reorder.genes) {
        ordered.genes <- sort(ordered.genes)
    }
    for (i in 1:ngenes.ordered) {
        expression.matrix[ordered.genes[i], ordered.samples] <-
            ordered.expression[i, ]
    }
    list(
        expression = expression.matrix,
        ordered.genes = ordered.genes,
        ordered.samples = sort(ordered.samples)
    )
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
datagen.swap <- function(expression, nswaps) {
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
datagen.noise <- function(expression, dist, ...) {
    ngenes <- nrow(expression)
    nsamples <- ncol(expression)
    noise.matrix <- matrix(dist(ngenes * nsamples, ...), nrow = ngenes, ncol = nsamples)
    expression + noise.matrix
}

#' Dropout a proportion of expression data
#'
#' Sets dropout.rate proportion of data to 0.
#'
#' @param expression Numeric matrix, represents gene expression data with genes
#'    as rows, and samples as columns.
#' @param dropout.rate Double, representing the proportion of values to set to
#'    0, 0<=dropout.rate<=1.
#'
#' @return Numeric matrix, expression data with dropout.rate proportion of
#'    values set to 0.
#' @export
#'
#' @examples
datagen.dropout <- function(expression, dropout.rate) {
    ngenes <- nrow(expression)
    nsamples <- ncol(expression)
    nvalues.drop <- floor(dropout.rate * (ngenes * nsamples))
    drop.index <- sample(ngenes * nsamples, size = nvalues.drop, replace = FALSE)
    drop.vector <- logical(ngenes * nsamples)
    drop.vector[drop.index] <- TRUE
    expression[drop.vector] <- 0
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
datagen.ordered.vector <- function(ngenes, dist, ...) {
    expression.unordered <- dist(ngenes, ...)
    sort(expression.unordered, decreasing = FALSE)
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
datagen.ordered.matrix <- function(ngenes, nsamples, dist, ...) {
    # Create an empty expression matrix
    expression.matrix <- matrix(0, nrow = ngenes, ncol = nsamples)
    # Generate the expression for each sample of the ordered phenotype
    for (i in 1:nsamples) {
        expression.matrix[, i] <- datagen.ordered.vector(ngenes, dist, ...)
    }
    expression.matrix
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
datagen.unordered.vector <- function(ngenes, dist, ...) {
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
datagen.unordered.matrix <- function(ngenes, nsamples, dist, ...) {
    # Create an empty expression matrix
    expression.matrix <- matrix(0, nrow = ngenes, ncol = nsamples)
    # Generate the expression for each sample of the ordered phenotype
    for (i in 1:nsamples) {
        expression.matrix[, i] <- datagen.unordered.vector(ngenes, dist, ...)
    }
    expression.matrix
}
