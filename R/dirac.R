#' dirac.rank_vector: Determine DIRAC ranking vector
#'
#' Creates a DIRAC ranking vector from expression data
#'
#' @param expression A numeric vector of gene expression data
#' @param expression.legnth Optional, length of expression matrix, useful for
#' repeated calls to this function with same expression length
#' @param matrix.index Optional, a Boolean vector or matrix with TRUE for entries
#' of gene-by-gene square matrix to be included in the rank_vector. Mainly used
#' for repeated calls to this function so it doesn't have to be recalculated
#' each time
#'
#' @return Numeric rank vector
#'
#' @examples
#' # example code
#'  dirac.rank_vector(c(4,2,1,3), 4)

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

#' dirac.rank_matrix: Determine DIRAC ranking matrix
#'
#' Creates a DIRAC ranking matrix from expression data
#'
#' @param expression A gene expression matrix, rows are genes, columns are
#' samples
#'
#' @return matrix whose columns represent the rank_vector for each sample

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
