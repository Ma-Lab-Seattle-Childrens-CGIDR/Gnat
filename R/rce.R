# Main Function -----------------------------------------------------------

#' Calculate p-values for Rank Conservation Entropy
#'
#' For each of the gene networks described by the gene_index_list, calculates
#' the p-values for the difference in rank conservation entropy between
#' the two phenotypes.
#'
#' @param gene_expression Numeric matrix representing gene expression, rows are
#'    genes, and columns are samples
#' @param gene_index_list List of Integer vectors, each representing the indices
#'    of a gene network in the gene_expression matrix
#' @param phenotype1,phenotype2 Integer vectors representing the indices of the
#'    two phenotypes in the gene_expression matrix
#' @param iterations Integer, Number of templates to sample from each of the
#'    phenotypes to calculate the distribution of the Kendall Rank Correlation
#'    coefficient
#' @param parallel Boolean describing whether parallel operation is desired. On
#'    unix systems creates fork cluster, on windows creates a psock cluster
#' @param cores Integer describing number of cores to use for parallel
#'    operation.
#' @param replace Boolean, describing whether the phenotype samples should be
#'    taken with replacement
#' @param seed1,seed2 Seeds to use for sampling the respective phenotypes
#' @param as.frame Boolean, describing whether the results should be returned
#'    as a data.frame
#'
#' @return If as.frame is TRUE:
#'          returns a dataframe with columns of gene network names
#'          (`$network_names`), and p-values (`$pvals`)
#'         If as.frame is FALSE:
#'          returns a named list of p-values, with names corresponding to names
#'          of the gene_index_list
#' @export
#'
#' @examples
#' # Create example random expression data
#' expression <- matrix(runif(50*120, min=1, max=500000), ncol=50)
#' # Create list of integer vectors representing different gene networks
#' gene_index_list <- list(
#'   gene_network_A = c(34,35,36,40,47,101),
#'   gene_network_B = c(115,116,120),
#'   gene_network_C = c(1,5,17,19,20,25,31,33)
#' )
#' # Create integer vector representing indices of phenotype1 and phenotype2
#' phenotype1 <- sample(50, 20,replace=FALSE)
#' phenotype2 <- sample(50, 15,replace=FALSE)
#' # Run the RCE
#' results <- rce(expression, gene_index_list, phenotype1, phenotype2,
#'   iterations=100, parallel=TRUE, cores=2, replace=TRUE, seed1=42, seed2=314,
#'   as.frame=TRUE)
#' # View the results
#' print(results)
#' # Repeat using serial computation
#' results <- rce(expression, gene_index_list, phenotype1, phenotype2,
#'   iterations=100, parallel=FALSE, replace=TRUE, seed1=42, seed2=314,
#'   as.frame=TRUE)
#' # Show results, note likely not the same due to small number of iterations
#' # and the different ways that random number generation is handled between
#' # serial and parallel operation
#' print(results)
rce <- function(gene_expression, gene_index_list, phenotype1, phenotype2,
                iterations=1000, parallel=TRUE, cores=4, replace=TRUE, seed1,
                seed2, as.frame=TRUE){
  pvals <- unlist(lapply(gene_index_list, rce.compare_phenotypes,
                         phenotype1=phenotype1, phenotype2=phenotype2,
                         iterations=iterations, replace=replace, seed1=seed1,
                         seed2=seed2))
  if(!as.frame){
    names(pvals) <- names(gene_index_list)
    return(pvals)
  }
  network_names <- names(gene_index_list)
  data.frame(cbind(network_names, pvals))
}




# Helper Functions --------------------------------------------------------

#' Compute the rank correlation entropy between two vectors
#'
#' `rce.rank_correlation.vector` computes the Kendall Tau Rank Correlation
#' coefficient between two vectors
#'
#' @param x,y Numeric vectors between which the correlation coefficient is to be
#'    calculated.
#' @returns Float, the Kendall Rank Correlation Coefficient between `x` and `y`
rce.rank_correlation.vector <- function(x,y){
  stats::cor.test(x,y, method="kendall")$estimate[[1]]
}

#' Compute the rank correlation entropy for one template
#'
#' `rce.rank_correlation.matrix` computes the Kendall Tau Rank correlation
#' coefficient between each column of gene_expression and the template
#'
#' @param gene_expression Numeric matrix, rows are genes in the network,
#'    and columns are samples
#' @param template Numeric vector,
#' @returns Float, mean of the Kendall Tau Rank correlation coefficients
#'    between the template and all columns of gene_expression
rce.rank_correlation.matrix <- function(gene_expression, template){
  apply(gene_expression, 2, rce.rank_correlation.vector,
        y=template)
}

#' Compute the Kendall correlation coefficient between gene expression and
#' sampled template
#'
#' `rce.sample` Takes in gene expression matrix, and index for phenotype to be
#' used as template, then computes the Kendall Rank Correlation Coefficient
#' between all other samples and the template
#'
#' @param gene_expression Numeric matrix representing the gene expression, with
#'    genes in the network as rows, and samples as columns
#' @param template_index Integer, index of sample to use as template
#'
#' @returns Float, the mean of the Kendall Rank correlation coefficient
#'    between the sampled template and the other samples in the phenotype
rce.sample <- function(template_index, gene_expression){
  mean(rce.rank_correlation.matrix(gene_expression[,-template_index],
                              gene_expression[,template_index]))
}


#' Compute empirical distribution for the Kendall Correlation Coefficient
#'
#' `rce.distribution` calculates the empirical distribution of the Kendall
#' Correlation Coefficient.
#'
#' For each iteration a template sample is chosen, and the Kendall rank
#' correlation coefficient is calculated between that template and all other
#' samples. Then the mean of these values is taken.
#'
#' @param gene_expression Numeric matrix, rows represent genes, and columns
#'    represent samples
#' @param gene_index Integer vector, representing the indices in the gene
#'    network
#' @param phenotype Integer vector, representing the indices of the phenotype
#' @param iterations
#' @param parallel Boolean, whether parallel operation is desired
#' @param cores Integer, number of  cores to use for parallel operation,
#'    if greater than the number of available cores will be reduced to the
#'    number of available cores
#' @param replace Boolean, Whether replacement should be used when sampling
#'    phenotype for templates.
#' @param seed Integer, seed for random number generator
#'
#' @return Float vector, represents the empirical rank correlation entropy
#'    distribution with the phenotype.
rce.distribution.central <- function(gene_expression, gene_index, phenotype,
                             iterations=1000, parallel=TRUE, cores=4,
                             replace=TRUE, seed){
  os_type=.Platform$OS.type
  set.seed(seed)
  if(!replace){
    if(iterations>length(phenotype)){
      templates <- phenotype
    } else {
      templates <- sample(phenotype,iterations, replace=replace)
    }
  } else{
    templates <- sample(phenotype, iterations, replace=replace)
  }
  gene_expression.filtered <- gene_expression[gene_index,]
  if(!parallel){
    res <- unlist(
      lapply(templates, rce.sample,
                         gene_expression = gene_expression.filtered)
                  )
  } else if(os_type=="unix"){
    res <- unlist(
      parallel::mclapply(templates, rce.sample,
             gene_expression = gene_expression.filtered, mc.cores=cores)
    )
  } else if(os_type=="windows"){
    # Make cluster
    cl <- parallel::makeCluster(cores)
    # Export needed functions
    parallel::clusterExport(cl, list("rce.rank_correlation.vector",
                                     "rce.rank_correlation.matrix",
                                     "rce.sample"))
    # perform calculation
    res <- tryCatch(expr={unlist(
      parallel::parLapply(cl, templates, rce.sample,
                          gene_expression=gene_expression.filtered))},
      finally = {parallel::stopCluster(cl)})
  } else{
    stop("Unsupported OS for parallel operation")
  }
  return(res)
}


#' Compare phenotypes using Rank Correlation Entropy
#'
#' @param gene_expression Numeric matrix representing gene expression, rows
#'    represent genes, and columns represent samples.
#' @param gene_index Integer vector, represents the indices of genes within
#'    a gene network.
#' @param phenotype1,phenotype2 Integer vectors, representing the indices
#'    of the two phenotypes
#' @param boostrap_iterations Integer, number of iterations to perform for
#'    bootstrapping to create the distribution.
#' @param parallel Boolean, whether to use parallel computing. On unix systems
#'    will use a fork cluster, on windows will use a psock cluster.
#' @param cores Integer, number of cores to use for parallel computation.
#' @param replace Boolean, whether to sample with replacement when performing
#'    bootstrapping.
#' @param seed1 Seed for bootstrapping the phenotype1 distribution.
#' @param seed2 Seed for bootstrapping the phenotype2 distribution.
#'
#' @return Float representing the p-value,
#'    Null hypothesis: The two distributions of rank correlations are from the
#'      same distributions.
#'    Alternate hypothesis: The two distributions of rank correlations are from
#'      different distributions.
#' @export
#'
#' @examples
rce.compare_phenotypes <- function(gene_expression, gene_index, phenotype1,
                                   phenotype2, bootstrap_iterations=1000,
                                   parallel=TRUE,
                                   cores=4, replace=TRUE, seed1, seed2){
    p1.dist <- rce.distribution.central(gene_expression=gene_expression,
                                gene_index=gene_index, phenotype = phenotype1,
                                iterations = bootstrap_iterations,
                                replace=replace,
                                seed=seed1)
    p2.dist <- rce.distribution.central(gene_expression=gene_expression,
                                gene_index=gene_index, phenotype = phenotype2,
                                iterations=bootstrap_iterations,
                                replace=replace,
                                seed=seed2)
    stats::ks.test(p1.dist, p2.dist)$p.value
}
