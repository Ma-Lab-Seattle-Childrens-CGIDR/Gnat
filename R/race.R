# Main Function -----------------------------------------------------------

#' Compare phenotypes and gene networks using RACE
#'
#' @param gene_network_list Named list of integer vectors, each representing
#'    the indices of the genes in the network within the expression matrix.
#' @param expression Numeric matrix, represents gene expression with rows
#'    representing genes, and columns representing samples
#' @param phenotype1,phenotype2 Integer vectors representing the indices
#'    of the samples in the phenotype within the expression matrix
#' @param bootstrap_iterations Integer number of iterations for the
#'    bootstrapping
#' @param parallel Boolean, whether computations should be done in parallel.
#'    See details for more information.
#' @param cores Integer, number of cores to use for parallel operation.
#' @param replace Boolean, whether sampling should be done with replacement
#' @param seed Integer, seed for RNG (note, even for the same seed, the
#'    RNG is handled differently for parallel operation between different
#'    operating systems).
#' @param as.frame Boolean, whether the result should be returned as a
#'    dataframe or a named list
#'
#' For parallel operations, on unix operating systems a fork cluster will be
#' used, while on windows a psock cluster will be used. The RNG for both uses
#' the L'Ecuyer method. For serial operation, a Mersenne-Twister RNG is used.
#' This means that results will vary between serial and parallel operation.
#'
#' @return A named list or dataframe depending on `as.dataframe`.
#'    If as.dataframe is TRUE:
#'      Will return a dataframe with a row for each of the gene networks,
#'      and columns for gene_network, p1.mean_rank_correlation,
#'      p2.mean_rank_correlation, absolute_difference, and p.value
#'    If as.dataframe is FALSE:
#'      will return a named list, named according to the names of the
#'      `gene_network_list` argument. Each element will also be a named list,
#'      with p1.mean_rank_correlation, p2.mean_rank_correlation,
#'      absolute_difference, and p.value.
#' @export
#'
#' @examples
RACE <- function(expression, gene_network_list,
                 phenotype1, phenotype2, bootstrap_iterations=1000,
                 parallel=TRUE, cores=4, replace=TRUE, seed=NULL,
                 as.frame=TRUE){
  # Run race.compare_phenotypes for each of the gene networks in the list
  res_list <- lapply(gene_network_list, race.compare_phenotypes,
                     expression=expression,
                     phenotype1=phenotype1,
                     phenotype2=phenotype2,
                     bootstrap_iterations=bootstrap_iterations,
                     parallel=parallel, cores=cores, replace=replace,
                     seed=seed)
  if(!as.frame){
    names(res_list) <- names(gene_network_list)
    return(res_list)
  } else {
    p1.mean_rank_correlation.list <- sapply(
      res_list, function(x) x$p1.mean_rank_correlation)
    p2.mean_rank_correlation.list <- sapply(
      res_list, function(x) x$p2.mean_rank_correlation
    )
    abs_diff_list <- sapply(res_list, function(x) x$absolute_difference)
    p.value.list <- sapply(res_list, function(x) x$p.value)
    gene_network_names <- names(gene_network_list)
    res_frame <- data.frame(
      gene_network = gene_network_names,
      p1.mean_rank_correlation = p1.mean_rank_correlation.list,
      p2.mean_rank_correlation = p2.mean_rank_correlation.list,
      absolute_difference = abs_diff_list,
      p.value = p.value.list
    )
    return(res_frame)
  }
}


# Helper Functions --------------------------------------------------------

#' Find the mean Kendall correlation between all samples
#'
#' @param expression.filtered Numeric Matrix, rows are genes in the gene
#'    network, and columns are samples in the phenotype
#'
#' @return Double, the mean of the Kendall Rank Correlation
#' @export
#'
#' @examples
#' # Create example matrix
#' mat <- matrix(c(1,2,3,4,2,3,1,4,4,1,2,3,1,4,2,3), ncol=4)
#' # Print the result
#' print(race.mean_rank_corr(mat))
race.mean_rank_corr <- function(expression.filtered){
  correlation <- cor(expression.filtered, method="kendall")
  mean(correlation[lower.tri(correlation)])
}


#' Compare two phenotypes mean Kendall Correlation Coefficient
#'
#' `race.compare_phenotypes.single` computes the absolute difference between the
#' mean Kendall rank correlation coefficient of the two phenotypes
#'
#' @param expression.gene_network Numeric matrix, represents the gene
#'    expression, with rows representing genes within the network of interest,
#'    and columns representing the samples.
#' @param phenotype1,phenotype2 Integer vectors, representing the indices of the
#'    two phenotypes within the expression.gene_network matrix
#'
#' @return Double, the absolute difference between the mean Kendall Rank
#'    Correlation of the two phenotypes
#' @export
#'
#' @examples
race.compare_phenotypes.single <- function(expression.gene_network, phenotype1,
                                           phenotype2){
  # get the expression matrices for each phenotype
  p1.expression <- expression.gene_network[,phenotype1]
  p2.expression <- expression.gene_network[,phenotype2]
  # Get the rank correlation means for each phenotype
  p1.rank_corr_mean <- race.mean_rank_corr(p1.expression)
  p2.rank_corr_mean <- race.mean_rank_corr(p2.expression)
  # Return the absolute value of the differences between the mean rank
  #   correlations
  abs(p1.rank_corr_mean-p2.rank_corr_mean)
}


#' Compare shuffled Phenotypes
#'
#' `race.compare_phenotypes.shuffle` Computes the absolute difference in
#' mean Kendall Rank Correlation Coefficient between shuffled phenotypes
#'
#' @param i Bootstrap iteration
#' @param expression.gene_network Numeric matrix, represents the gene
#'    expression, with rows representing genes in the network of interest,
#'    and columns representing samples.
#' @param combined Integer vector, represents the indices of both phenotypes
#'    (concatenated) within the expression.gene_network matrix
#' @param p1.size,p2.size Integers representing the size of the two phenotypes
#' @param replace Boolean, whether sampling should be done with replacement
#'
#' @return Double, the absolute difference between the mean Kendall
#'    Rank Correlation Coefficient of the two shuffled phenotypes
#' @export
#'
#' @examples
race.compare_phenotypes.shuffle <- function(i, expression.gene_network,
                                            combined, p1.size, p2.size,
                                            replace=TRUE){
  if(replace){
    # Sample combined for p1.indices
    p1.idx <- sample(combined, p1.size, replace=replace)
    # Sample combined for the p2.indices
    p2.idx <- sample(combined, p2.size, replace = replace)
  } else {
    # Sample for p1 indices
    p1.idx <- sample(combined, p1.size, replace=replace)
    # Use the remaining indices for p2.idx
    p2.idx <- setdiff(combined, p1.idx)
  }
  # Compute the absolute difference in mean Kendall Correlation between the
  # shuffled phenotypes
  race.compare_phenotypes.single(
    expression.gene_network = expression.gene_network, phenotype1=p1.idx,
    phenotype2=p2.idx)
}


#' Compare phenotypes using the Absolute Difference in Mean Rank Correlation
#' Coefficient
#'
#' @param gene_network Integer vector representing the indices of the genes
#'    in the network of interest within the expression matrix
#' @param expression Numeric matrix, represents gene expression with rows
#'    representing genes, and columns representing samples
#' @param phenotype1,phenotype2 Integer vectors representing the indices
#'    of the samples in the phenotype within the expression matrix
#' @param bootstrap_iterations Integer number of iterations for the
#'    bootstrapping
#' @param parallel Boolean, whether computations should be done in parallel.
#'    See details for more information.
#' @param cores Integer, number of cores to use for parallel operation.
#' @param replace Boolean, whether sampling should be done with replacement
#' @param seed Integer, seed for RNG (note, even for the same seed, the
#'    RNG is handled differently for parallel operation between different
#'    operating systems).
#'
#' For parallel operations, on unix operating systems a fork cluster will be
#' used, while on windows a psock cluster will be used. The RNG for both uses
#' the L'Ecuyer method. For serial operation, a Mersenne-Twister RNG is used.
#' This means that results will vary between serial and parallel operation.
#'
#' @return Named list, with p1.mean_rank_correlation, p2.mean_rank_correlation,
#'    absolute_difference, and p.value.
#' @export
#'
#' @examples
race.compare_phenotypes <- function(gene_network, expression,
                                    phenotype1, phenotype2,
                                    bootstrap_iterations=1000, parallel=TRUE,
                                    cores=4, replace=TRUE, seed=NULL){
  # Get the size of each of the phenotypes
  p1.size <- length(phenotype1)
  p2.size <- length(phenotype2)
  # Get the combined indices of the two phenotypes
  combined <- c(phenotype1, phenotype2)
  # Find the mean kendall rank correlation coefficient for the two phenotypes
  p1.mean_rc <- race.mean_rank_corr(
    expression.filtered = expression[gene_network, phenotype1])
  p2.mean_rc <- race.mean_rank_corr(
    expression.filtered = expression[gene_network, phenotype2])
  # Find the absolute difference between the two phenotypes
  abs_diff_mean_rc <- abs(p1.mean_rc - p2.mean_rc)
  # Perform the bootstrapping for the null distribution
  if(parallel){
    cores=if(parallel::detectCores()>cores) cores else parallel::detectCores()
  }
  os_type <- .Platform$OS.type
  if(parallel){
    if(os_type=="unix"){
      # Set seed
      set.seed(seed, "L'Ecuyer")
      # Use multicore apply
      res <- unlist(parallel::mclapply(1:bootstrap_iterations,
                                       race.compare_phenotypes.shuffle,
                                       expression.gene_network =
                                         expression[gene_network,],
                                       combined=combined, p1.size=p1.size,
                                       p2.size=p2.size, replace=replace,
                                       mc.cores=cores))
    } else if(os_type=="windows"){
      # make the cluster
      cl <- parallel::makeCluster(cores)
      # Set the RNG stream seed
      parallel::clusterSetRNGStream(cl, seed)
      # Export the needed functions
      parallel::clusterExport(cl, list("race.mean_rank_corr",
                                       "race.compare_phenotypes.single",
                                       "race.compare_phenotypes.shuffle"))
      res <- tryCatch(expr={
        unlist(parallel::parLapply(cl, 1:bootstrap_iterations,
                                   race.compare_phenotypes.shuffle,
                                   expression.gene_network=
                                     expression[gene_network,],
                                   combined=combined, p1.size=p1.size,
                                   p2.size=p2.size, replace=replace))
      }, finally = {parallel::stopCluster(cl)})
    } else {
      stop("Unsupported OS for parallel operation")
    }
  } else {
    set.seed(seed, kind="Mersenne-Twister", normal.kind = "Inversion")
    res <- unlist(lapply(1:bootstrap_iterations,
                         race.compare_phenotypes.shuffle,
                         expression.gene_network =
                           expression[gene_network,],
                         combined=combined, p1.size=p1.size,
                         p2.size=p2.size, replace=replace))
  }
  # Now create the ecdf
  boot_cdf <- stats::ecdf(res)
  # Get the p-value for the absolute rank correlation difference
  p.value <- 1-boot_cdf(abs_diff_mean_rc)
  # Return named list with the results
  list(p1.mean_rank_correlation = p1.mean_rc,
       p2.mean_rank_correlation = p2.mean_rc,
       absolute_difference = abs_diff_mean_rc,
       p.value = p.value)
}




