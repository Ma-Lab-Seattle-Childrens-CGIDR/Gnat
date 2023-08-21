# Main Function -----------------------------------------------------------
#' Title
#'
#' @param expression
#' @param gene_network_list
#' @param phenoytpe1
#' @param phenotype2
#' @param rank_fun
#' @param score_fun
#' @param bootstrap_iterations
#' @param parallel
#' @param cores
#' @param BPPARAM
#' @param replace
#' @param seed
#' @param as_frame
#'
#' @return
#' @export
#'
#' @examples
compare_phenotypes <- function(expression, gene_network_list, phenoytpe1,
                               phenotype2, rank_fun, score_fun,
                               bootstrap_iterations = 1000,
                               parallel = TRUE, cores = 2, BPPARAM=NULL,
                               replace=TRUE, seed=NULL, as_frame = TRUE){
  # Ensure list is named
  gene_network_list <- ensure_named(gene_network_list)
  # Run the bootstrapping
  res_list <- lapply(gene_network_list, compare_phenotypes_single,
                     expression=expression, phenotype1=phenotype1,
                     phenotype2=phenotype2,
                     bootstrap_iterations=bootstrap_iterations,
                     parallel = parallel, cores = cores, replace = replace,
                     seed=seed)
  if(!as_frame){
    names(res_list <- names(gene_network_list))
    resturn(res_list)
  }
  p1_score_list <- sapply(res_list, function(x) x$p1_score)
  p2_score_list <- sapply(res_list, function(x) x$p2_score)
  absolute_difference_list <- sapply(res_list,
                                     function(x) x$absolute_difference)
  p_values_list <- sappy(res_list, function(x) x$p_value)
  gene_networks <- names(gene_network_list)
  data.frame(
    gene_network = gene_networks,
    p1_score = p1_score_list,
    p2_score = p2_score_list,
    absolute_difference = absolute_difference_list,
    p_value = p_values_list
  )
}


# Helper Functions --------------------------------------------------------
#' Compare two phenotypes for a single gene network
#'
#' @param rank_matrix Numeric matrix With samples as the columns. See details
#'    for further information.
#' @param phenotype1,phenotype2 Index vector represent the indices of each
#'    phenotype within the `rank_matrix`.
#' @param score_fun Function which calculates a score form the
#'    rank matrix.
#'
#' @return `numeric` Absolute difference of the scores for the two phenotypes
#' @export
#'
#' @details
#' Rank matrix is normally generated from a passed rank_fun function, so the
#' columns should represent samples, but what the rows represent depends
#' on the specific rank function. For RACE, CRANE, and INFER the rows simply
#' represent genes, while for DIRAC they represent specific pairwise
#' comparison.
#'
#' The phenotype vectors can contain any values that can be used with the
#' `[` matrix function, and the `c` function.
#'
#' @examples
compare_phenotypes_single <- function(rank_matrix, phenotype1, phenotype2,
                                      score_fun){
  p1_rank_matrix <- rank_matrix[,phenotype1]
  p2_rank_matrix <- rank_matrix[,phenotype2]
  p1_score <- score_fun(p1_rank_matrix)
  p2_score <- score_fun(p2_rank_matrix)
  abs(p1_score-p2_score)
}

#' Compare two shuffled phenotypes
#'
#' @param iteration Value representing the current bootstrap iteration,
#'    used to make signature match the required signature for `bplapply`.
#' @param rank_matrix Numeric matrix with samples as the columns. See details
#'    for further information.
#' @param combined Vector representing the concatenated indices of
#'    phenotype 1 and 2.
#' @param p1_size,p2_size Integer representing the size of the two phenotypes
#' @param score_fun Function which calculates a score form the
#'    rank matrix.
#' @param replace Boolean representing Whether sampling should be done with
#'    replacement
#'
#' @return
#' @export
#'
#' @details
#' Rank matrix is normally generated from a passed rank_fun function, so the
#' columns should represent samples, but what the rows represent depends
#' on the specific rank function. For RACE, CRANE, and INFER the rows simply
#' represent genes, while for DIRAC they represent specific pairwise
#' comparison.
#'
#' The combined vector can contain any values that can be used with the
#' `[` matrix function, as well as the vector `c` and `sample` functions.
#'
#' The replace parameter dictates whether the sampling to create shuffled
#' phenotype index vectors should be done with replacement. If TRUE, sampling
#' with replacement if performed to create the two shuffled index vectors
#' of the appropriate sizes. If FALSE, random sampling without replacement
#' occurs once to select phenotype1 indices, and the remaining indices are used
#' for phenotype2.
#'
#' @examples
score_shuffled <- function(iteration, rank_matrix, combined, p1_size, p2_size,
                           score_fun, replace=TRUE){
  if(replace){
    p1_idx <- sample(combined, p1_size, replace = replace)
    p2_idx <- sample(combined, p2_size, repalce = replace)
  } else {
    p1_idx <- sample(combined, p1_size, replace = replace)
    p2_idx <- setdiff(combined, p1_idx)
  }
  compare_phenotypes_single(rank_matrix = rank_matrix, phenotype1 = p1_idx,
                     phenotype2 = p2_idx, score_fun = score_fun)

}

#' Title
#'
#' @param BPPARAM
#' @param parallel
#' @param cores
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
create_parallel_param <- function(BPPARAM, parallel, cores, seed) {
  # Create the BPPARAM param object if needed
  if(is.null(BPPARAM)){
    if(!parallel){
      BPPARAM <- SerialParam(RNGseed = NULL)
    } else {
      core_limit <- parallel::detectCores()
      cores <- if(cores<=core_limit) cores else core_limit
      os_type <- .Platform$OS.type
      if(os_type == "unix"){
        BPPARAM <- MulticoreParam(workers = cores, RNGseed = seed)
      } else if(os_type=="windows"){
        BPPARAM <- SnowParam(workers = cores, RNGseed = seed)
      } else {
        BPPARAM <- SnowParam(workers = cores, RNGseed = seed)
      }
    }
  }
  BPPARAM
}



#' Title
#'
#' @param gene_network
#' @param rank_fun
#' @param score_fun
#' @param expression
#' @param phenotype1
#' @param phenotype2
#' @param bootstrap_iterations
#' @param parallel
#' @param cores
#' @param BPPARAM
#' @param replace
#' @param seed
#'
#' @return
#' @importFrom BiocParallel bplapply MulticoreParam SnowParam
#' @export
#'
#' @examples
bootstrap_score <- function(gene_network, rank_fun, score_fun, expression,
                            phenotype1, phenotype2, bootstrap_iterations=1000,
                            parallel = TRUE, cores=2, BPPARAM=NULL,
                            replace=TRUE, seed=NULL){
  # Create the BPPARAM object
  BPPARAM <- create_parallel_param(BPPARAM =  BPPARAM, parallel = parallel,
                                   cores=cores, seed = seed)
  # Get the lengths of the phenotypes
  p1_size <- length(phenotype1)
  p2_size <- length(phenotype2)
  # Get the combined phenotype index vector
  combined_phenotypes <- c(phenotype1, phenotype2)
  # Compute the rank matrix using the provided method
  rank_matrix <- rank_fun(expression[gene_network,combined_phenotypes])
  # Compute the unshuffled score for the two phenotypes
  p1_score <- score_fun(rank_matrix[,phenotype1])
  p2_score <- score_fun(rank_matrix[,phenotype2])
  abs_diff <- abs(p1_score-p2_score)
  res <- unlist(bplapply(seq_len(bootstrap_iterations), score_shuffled,
                         rank_matrix = rank_matrix, p1_size=p1_size,
                         p2_size=p2_size, replace=replace, BPPARAM = BPPARAM))
  # Create the empirical CDF
  boot_cdf <- stats::ecdf(res)
  p_value <- 1-boot_cdf(abs_diff)
  # Return results list
  list(
    p1_score = p1_score,
    p2_score = p2_score,
    absolute_difference = abs_diff,
    p_value=p_value
  )

}
