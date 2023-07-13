# Main Function -----------------------------------------------------------

#' Compare phenotypes using a list of Gene Networks
#'
#' @param expression Numeric matrix representing gene expression, rows
#'    represent genes, and columns represent samples
#' @param gene_network_list Named list, with each element being a integer
#'    vector containing the indices of the genes within the network
#' @param phenotype1,phenotype2 Integer vectors of indices of the samples
#'    in each phenotype.
#' @param bootstrap_iterations Integer number of iterations to perform
#'    for the bootstrapping
#' @param parallel Boolean whether the calculation should be done in
#'    parallel
#' @param cores Integer number of cores to use for parallel computation
#' @param replace Boolean whether to use replacement when sampling during
#'    bootstrap operations.
#' @param seed Integer seed for RNG
#' @param as.frame Boolean whether return should be data.frame or named
#'    list
#'
#' @return data.frame or named list
#'    If as.frame is TRUE:
#'      Returns a dataframe where each row is a gene network, and the columns
#'      are gene_network, p1.mean_centroid_distance,
#'      p2.mean_centroid_distance, absolute_difference, and p.value.
#'    If as.frame is FALSE:
#'      Returns a named list, where the names are the gene_network,
#'      taken from the gene_network_list names. Each element is another
#'      named list, with p1.mean_centroid_distance,
#'      p2.mean_centroid_distance, absolute_difference, and p.value.
#' @export
#'
#' @examples
CRANE <- function(expression, gene_network_list, phenotype1, phenotype2,
                  bootstrap_iterations=1000, parallel=TRUE, cores=4,
                  replace=TRUE, seed=NULL, as.frame=TRUE){
  # Run the crane.compare_phenotypes
  res_list <- lapply(gene_network_list, crane.compare_phenotypes,
                     expression=expression, phenotype1=phenotype1,
                     phenotype2=phenotype2,
                     bootstrap_iterations = bootstrap_iterations,
                     parallel=parallel,cores=cores, replace=replace,
                     seed=seed)
  if(!as.frame){
    names(res_list) <- names(gene_network_list)
    return(res_list)
  } else {
    p1.mean_centroid_distance.list <-
      sapply(res_list, function(x) x$p1.mean_centroid_distance)
    p2.mean_centroid_distance.list <-
      sapply(res_list, function(x) x$p2.mean_centroid_distance)
    absolute_difference.list <- sapply(
      res_list, function(x) x$absolute_difference)
    p.values.list <- sapply(res_list, function(x) x$p.value)
    gene_network.names <- names(gene_network_list)
    res_frame <- data.frame(
      gene_network = gene_network.names,
      p1.mean_centroid_distance = p1.mean_centroid_distance.list,
      p2.mean_centroid_distance = p2.mean_centroid_distance.list,
      absolute_difference=absolute_difference.list,
      p.value=p.values.list
    )
    return(res_frame)
  }
}





# Helper Functions --------------------------------------------------------

#' Compute the rank matrix
#'
#' @param expression.filtered Numeric matrix of gene expression values, genes
#'    are rows, and samples are columns
#'
#' @return Integer matrix of ranks
#' @export
#'
#' @examples
crane.rank_matrix <- function(expression.filtered){
  apply(expression.filtered, MARGIN=2, rank, ties.method="first")
}

#' Compute the mean centroid distance
#'
#' @param rank_matrix Integer matrix with genes as rows, and samples as columns
#'
#' @return Double, mean distance of each samples rank vector from the centroid
#'    of all the samples rank vectors
#' @export
#'
#' @examples
crane.mean_centroid_distance <- function(rank_matrix){
  centroid <- apply(rank_matrix, MARGIN=1, mean, na.rm=TRUE)
  sd <- (rank_matrix-centroid)^2
  mean(sqrt(apply(sd, MARGIN=2, sum)))
}

#' Compare phenotypes using the mean rank centroid distance
#'
#' @param rank_matrix Integer matrix with genes as rows, and samples as columns
#' @param phenotype1,phenotype2 Integer vector of indices for samples in a
#'    phenotype
#'
#' @return Double, absolute difference between the mean centroid distance
#'    for the two phenotypes
#' @export
#'
#' @examples
crane.compare_phenotypes.single <- function(rank_matrix,
                                            phenotype1, phenotype2){
  # Get the rank matrix for each phenotype
  p1.rank_matrix <- rank_matrix[,phenotype1]
  p2.rank_matrix <- rank_matrix[,phenotype2]
  # Get the mean rank centroid difference
  p1.mean_centroid_distance <- crane.mean_centroid_distance(p1.rank_matrix)
  p2.mean_centroid_distance <- crane.mean_centroid_distance(p2.rank_matrix)
  # Return the absolute value of the difference in mean centroid distance
  # between the two phenotypes
  abs(p1.mean_centroid_distance-p2.mean_centroid_distance)
}

#' Compare shuffled phenotypes using mean rank centroid distance
#'
#' @param i Integer, bootstrap iteration
#' @param rank_matrix Integer matrix, ranks of gene expression with genes as
#'    rows, and samples as columns
#' @param combined Integer vector, indices of both phenotypes
#' @param p1.size,p2.size Integers, number of samples in each phenotype
#' @param replace Boolean whether replacement should be used when sampling
#'
#' @return Double, the absolute difference between the mean centroid
#'    distance for the shuffled phenotypes
#' @export
#'
#' @examples
crane.compare_phenotypes.shuffle <- function(i, rank_matrix, combined,
                                             p1.size, p2.size, replace=TRUE){
  if(replace){
    p1.idx <- sample(combined, p1.size, replace = replace)
    p2.idx <- sample(combined, p2.size, replace = replace)
  } else {
    p1.idx <- sample(combined, p1.size, replace=replace)
    p1.idx <- setdiff(combined, p1.idx)
  }
  crane.compare_phenotypes.single(rank_matrix = rank_matrix,
                                  phenotype1 = p1.idx,
                                  phenotype2 = p2.idx)
}

#' Compare phenotypes using a list of Gene Networks
#'
#' @param expression Numeric matrix representing gene expression, rows
#'    represent genes, and columns represent samples
#' @param gene_network Integer vector containing the indices of the genes
#'    within the network
#' @param phenotype1,phenotype2 Integer vectors of indices of the samples
#'    in each phenotype.
#' @param bootstrap_iterations Integer number of iterations to perform
#'    for the bootstrapping
#' @param parallel Boolean whether the calculation should be done in
#'    parallel
#' @param cores Integer number of cores to use for parallel computation
#' @param replace Boolean whether to use replacement when sampling during
#'    bootstrap operations.
#' @param seed Integer seed for RNG
#'
#' @return Named list, with p1.mean_centroid_distance,
#'    p2.mean_centroid_distance, absolute_difference, and p.value
#' @export
#'
#' @examples
crane.compare_phenotypes <- function(gene_network, expression, phenotype1,
                                     phenotype2, bootstrap_iterations=1000,
                                     parallel=TRUE, cores=4, replace=TRUE,
                                     seed=NULL){
  # Get the sizes of the phenotypes
  p1.size <- length(phenotype1)
  p2.size <- length(phenotype2)
  # Get the combined phenotype index vector
  combined <- c(phenotype1, phenotype2)
  # Compute the rank matrix
  rank_matrix <- crane.rank_matrix(expression[gene_network, combined])
  # Find the centroid difference for the unshuffled phenotypes
  p1.mean_centroid_distance <- crane.mean_centroid_distance(
    rank_matrix[,phenotype1])
  p2.mean_centroid_distance <- crane.mean_centroid_distance(
    rank_matrix[,phenotype2])
  abs_diff <- abs(p1.mean_centroid_distance-p2.mean_centroid_distance)
  # Perform bootstrapping to create the null distribution
  if(parallel){
    cores = if(parallel::detectCores()>cores) cores else parallel::detectCores()

  }
  os_type <- .Platform$OS.type
  if(parallel){
    if(os_type=="unix"){
      # Set Seed
      set.seed(seed, "L'Ecuyer")
      # Use multicore apply
      res <- unlist(parallel::mclapply(1:bootstrap_iterations,
                                       crane.compare_phenotypes.shuffle,
                                       rank_matrix=rank_matrix,
                                       combined=combined, p1.size=p1.size,
                                       p2.size=p2.size, replace=replace,
                                       mc.cores=cores))
    } else if (os_type=="windows"){
      # Create the cluster
      cl <- parallel::makeCluster(cores)
      # Set the RNG stream seed
      parallel::clusterSetRNGStream(cl, seed)
      # Export the needed functions
      parallel::clusterExport(cl, list("crane.mean_centroid_distance",
                                       "crane.compare_phenotypes.single",
                                       "crane.compare_phenotypes.shuffle"))
      # Run the bootstrap
      res <- tryCatch(expr={
        unlist(parallel::parLapply(cl, 1:bootstrap_iterations,
                                   crane.compare_phenotypes.shuffle,
                                   rank_matrix=rank_matrix,
                                   combined=combined, p1.size=p1.size,
                                   p2.size=p2.size,replace=replace))
      },
      finally = {parallel::stopCluster(cl)})
    } else {
        stop("Unsupported OS for parallel operation")
    }
  } else {
    set.seed(seed, kind="Mersenne-Twister", normal.kind = "Inversion")
    res <- unlist(lapply(1:bootstrap_iterations,
                         crane.compare_phenotypes.shuffle,
                         rank_matrix=rank_matrix,
                         combined=combined,
                         p1.size=p1.size,
                         p2.size=p2.size,
                         replace=replace))
  }
  # Now create the ECDF
  boot_cdf <- stats::ecdf(res)
  # Get the p-value for the absolute rank entropy difference
  p.value <- 1-boot_cdf(abs_diff)
  # Return results list
  list(p1.mean_centroid_distance=p1.mean_centroid_distance,
       p2.mean_centroid_distance=p2.mean_centroid_distance,
       absolute_difference=abs_diff,
       p.value=p.value)
}













