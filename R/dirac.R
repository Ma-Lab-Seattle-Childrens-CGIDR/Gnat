# Based on: Eddy, James A., et al. “Identifying Tightly Regulated and Variably
# Expressed Networks by Differential Rank Conservation (DIRAC).” PLoS
# Computational Biology, vol. 6, no. 5, May 2010.
# go-gale-com.offcampus.lib.washington.edu,
# https://doi.org/10.1371/journal.pcbi.1000792.



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
#' dirac.rank_vector(c(4,2,1,3))
#' dirac.rank_vector(c(4,2,1,3),4)
#' @export
dirac.rank_vector <- function(expression,
                              expression.length=NULL,
                              matrix.index=NULL){
  # Calculate the expression.length if needed
  if(is.null(expression.length)){
    expression.length <- length(expression)
  }
  # Calculate the matrix.index if needed
  if(is.null(matrix.index)){
    matrix.index <- lower.tri(matrix(1, nrow = expression.length,
                                     ncol = expression.length))
  }
  # Create a square matrix from the gene expression data,
  # used for finding pairwise rank comparison
  expression.matrix <- matrix(expression,
                              nrow = expression.length,
                              ncol = expression.length,
                              byrow = TRUE)
  ((expression.matrix - t(expression.matrix))<0)[matrix.index]
}

#' Determine DIRAC ranking matrix
#'
#' `dirac.rank_matrix` creates a DIRAC ranking matrix from expression data.
#'
#' @param expression A gene expression matrix, rows are genes, columns are
#'    samples
#'
#' @return matrix whose columns represent the rank_vector for each sample
#'
#' @examples
#' dirac.rank_matrix(matrix(1:16, ncol=4, nrow=4))
#'
#' @export
dirac.rank_matrix <- function(expression){
  # Calculate the expression.length (number of genes)
  expression.length <- nrow(expression)
  # Calculate the matrix.index vector
  matrix.index <- as.vector(
    lower.tri(matrix(1, nrow = expression.length, ncol = expression.length)))
  # Apply the rank_vector function to find the rank_vector of each column
  apply(expression,2, dirac.rank_vector,
        expression.length=expression.length,
        matrix.index=matrix.index)
}



#' Generate rank template
#'
#' `dirac.rank_template` creates a DIRAC ranking template from expression data.
#'
#' Strict inequality is used for finding the rank template.
#'
#' @param expression A gene expression matrix, rows are genes, columns are
#'    samples
#'
#' @return Numeric vector, rank template
#'
#' @examples
#' dirac.rank_template(matrix(1:16,ncol=4, nrow=4))
#'
#' @export
dirac.rank_template <- function(expression){
  # From the expression data, create the rank matrix
  rank_matrix <- dirac.rank_matrix(expression = expression)
  # Take the row means to get the conditional probability
  means <- rowMeans(rank_matrix)
  # Determine which rows have probability below 0.5
  means>0.5
}

#' Find rank matching score
#'
#' `dirac.rank_matching_score` finds the rank matching score between a rank
#'    vector, and a rank template.
#'
#' @param rank_vector Boolean vector representing the rank vector, must be
#'    the same length as rank_template.
#' @param rank_template Boolean vector representing the rank template.
#'
#' @returns Float, representing the rank matching score
#'
#' @examples
#'  dirac.rank_matching_score(c(1,0,0,1,0), c(0,1,0,1,0))
#' @export
dirac.rank_matching_score <- function(rank_vector, rank_template){
  if(length(rank_vector)!=length(rank_template)){
    stop("rank_vector and rank_template must be same length")
  }
  mean((rank_vector == rank_template))
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
#' expression <-  matrix(1:16, ncol=4, nrow=4)
#' rank_matrix <-dirac.rank_matrix(expression)
#' rank_template <- dirac.rank_template(expression)
#' dirac.rank_matching_score.vector(rank_matrix, rank_template)
#' @export
dirac.rank_matching_score.vector <- function(rank_matrix, rank_template){
  apply(rank_matrix,2, dirac.rank_matching_score,
        rank_template=rank_template)
}

#' Find rank conservation index
#'
#' `dirac.rank_conservation_index` finds the rank rank conservation index in
#'    gene expression data
#'
#' @param rank_matrix Boolean matrix, columns represent each samples rank
#'      vector
#' @param rank_template Boolean vector, represents template to match against
#'
#' @returns Float representing the rank conservation index
#'
#' @examples
#' expression <-  matrix(1:16, ncol=4, nrow=4)
#' rank_matrix <-dirac.rank_matrix(expression)
#' rank_template <- dirac.rank_template(expression)
#' dirac.rank_conservation_index(rank_matrix, rank_template)
#' @export
dirac.rank_conservation_index <- function(rank_matrix, rank_template){
  mean(dirac.rank_matching_score.vector(rank_matrix = rank_matrix,
                                   rank_template = rank_template))
}
