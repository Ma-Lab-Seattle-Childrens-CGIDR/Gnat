# Citation ----------------------------------------------------------------

# Based on: Eddy, James A., et al. “Identifying Tightly Regulated and Variably
# Expressed Networks by Differential Rank Conservation (DIRAC).” PLoS
# Computational Biology, vol. 6, no. 5, May 2010.
# go-gale-com.offcampus.lib.washington.edu,
# https://doi.org/10.1371/journal.pcbi.1000792.


# Main Functions ----------------------------------------------------------

#' Compare classification rate for different gene networks
#'
#' `DIRAC.compare_network_classification` compares the classification rate
#'    between multiple gene networks for two phenotypes
#'
#' @param expression Numeric matrix, represents gene expression, rows are genes
#'    columns are samples
#' @param phenotype1,phenotype2 Numeric vectors, contain the indices for the
#'    samples corresponding to each phenotype.
#' @param gene_index List of numeric vectors, each containing the indices for
#'    one of the gene networks.
#' @param parallel Boolean indicating Whether the computation should be run in
#'   parallel.
#' @param cores Integer, indicates number of cores to use for parallel
#'    calculations. If greater than parallel::detectCores() return value,
#'    will instead be set to said return value.
#' @param as.frame Boolean, whether the return value should be a dataframe or
#'    a named list.
#'
#' @returns
#'    if as.frame is TRUE:
#'      data.frame with columns for gene_network name (determined by names of
#'      gene_index list), and classification rate.
#'    if as.frame is FALSE:
#'      Numeric vector containing the classification rate for each of the
#'      provided gene networks. Named based on names of gene_index argument.
#' @examples
#' # example code
#' @export
DIRAC.compare_network_classification <- function(expression, phenotype1,
                                                 phenotype2, gene_index,
                                                 parallel = TRUE, cores=4,
                                                 as.frame=FALSE){
  # Wrapper around dirac.classification_rate.compare
  dirac.classification_rate.compare(expression, phenotype1,
                                    phenotype2, gene_index,
                                    parallel, cores,
                                    as.frame)

}

#' Compare two phenotypes using the difference in rank conservation indices
#'
#' `DIRAC.compare_phenotypes` computes the absolute difference between the rank
#' conservation indices for two phenotypes based on a single gene network, and
#' the associated p-value using a bootstrapping approach.
#'
#' When using the parallel computation option, the bootstrapping is
#' parallelized, but testing each gene network is done serially.
#'
#' @param expression Numeric matrix representing gene expression, rows
#'    represent genes, and columns represent samples.
#' @param phenotype1,phenotype2 Numeric vectors representing the indices of
#'    two phenotypes
#' @param gene_index A list of numeric vectors, each representing the indices
#'    of a gene network.
#' @param bootstrap_iterations Number of bootstrap iterations to run for computing
#'    the null distribution for calculating the p-value
#' @param parallel Boolean determining if calculation should be run in parallel,
#'    on unix systems it uses the mclapply function, and on windows it creates
#'    a socket cluster and uses parLapply.
#' @param cores number of cores to use for parallel computation if desired.
#'    If it is greater than the number of cores available (as determined by
#'    parallel::detectCores), it is instead set to the number of cores
#'    available.
#' @param replace Whether the sampling should be done with replacement. If it is
#'    TRUE, then two samples of size p1.size and p2.size respectively are taken
#'    with replacement from combined (these two samples can overlap). If FALSE,
#'    then a sample of size p1.size is taken without replacement from combined,
#'    and the remainder of combined is taken as the other sample.
#' @param seed Integer, used to set the seed for the random number generator
#'    used for sampling.
#' @param as.frame Boolean, whether return value should be dataframe or named
#'    list, see return for more information.
#' @returns
#'    If as.frame is TRUE:
#'      a dataframe with columns for gene network name
#'      (determined by the names of the gene_index list), value of absolute
#'      rank conservation index, and p-value of the absolute rank conservation
#'      index.
#'    If as.frame is FALSE:
#'      A list named according to the names in the
#'      gene_index argument, with each entry being a named list, with `$value`
#'      equal to the absolute difference in rank conservation scores between the
#'      two phenotypes, and `$p.value` equal to the p-value for `$value` found with
#'      bootstrapping.
#' @examples
#' # example code
#' @export
DIRAC.compare_phenotypes <- function(expression, phenotype1, phenotype2,
                                     gene_index, bootstrap_iterations=1000,
                                     parallel=TRUE, cores=4, replace=TRUE,
                                     seed=NULL, as.frame=TRUE){
  # Run the dirac.commpare phenotype function for each of the gene networks
  res_list <- lapply(gene_index, dirac.compare_phenotype,
                     expression=expression,
                     phenotype1=phenotype1,
                     phenotype2=phenotype2,
                     bootstrap_iterations=bootstrap_iterations,
                     parallel=parallel, cores=cores, replace=replace,
                     seed=seed)
  if(!as.frame){
    names(res_list) <- names(gene_index)
    return(res_list)
  } else{
    # Get the values
    values.list <- sapply(res_list, function(x) x$value)
    p.values.list <- sapply(res_list, function(x) x$p.value)
    gene_networks.list <- names(gene_index)
    res_frame <- data.frame(
      gene_network=gene_networks.list,
      value=values.list,
      p.value=p.values.list
    )
    return(res_frame)
  }
}

#' Creates a phenotype classifier based on DIRAC
#'
#' `DIRAC.classifier` Creates a function which takes in a matrix of gene
#' expression values and returns a classification of phenotype1, or phenotype2
#' for each of the samples.
#'
#' @param expression Numeric matrix of gene expression values, rows are genes
#'    columns are samples
#' @param phenotype1,phenotype2 Numeric vectors with indices for each of the
#'    phenotypes
#' @param gene_index Numeric vector of indices for the genes within the gene
#'    network being used for classification
#' @param return String describing desired return value from created function:
#'    "name":  returns a character vector of the names of the predicted
#'      phenotypes
#'    "is phenotype1": returns a boolean vector indicating if the
#'      predicted phenotype is phenotype1
#'    "is phenotype2": returns a boolean vector indicating if the predicted
#'      phenotype is phenotype2
#'    "phenotype number": returns numeric vector, with 1 representing the
#'      predicted phenotype is phenotype 1, and 2 representing that the
#'      predicted phenotype is phenotype2
#'    Any other value will lead the function to just return the rank matching
#'      difference
#' @param names optional Character vector with the desired names for the phenotypes,
#'
#' @returns A function which takes as an argument a gene expression numeric
#'    matrix with genes as the rows and samples as the columns, and returns
#'    a vector predicting which samples are from which phenotype (see return
#'    argument).
#' @examples
#' # example code
#' @export
DIRAC.classifier <- function(expression, phenotype1, phenotype2, gene_index,
                             return="is phenotype1",
                             names=c("phenotype1", "phenotype2")){
  # get the rank matrices for the expression matrix
  p1.rank_matrix <- dirac.rank_matrix(expression[gene_index, phenotype1])
  p2.rank_matrix <- dirac.rank_matrix(expression[gene_index, phenotype2])
  # Get the rank templates for each of the phenotypes
  p1.rank_template <- dirac.rank_template(p1.rank_matrix)
  p2.rank_template <- dirac.rank_template(p2.rank_matrix)
  # Create an intermediate function which returns a vector of the rank
  #   matching score differences from an expression matrix
  rank_matching.diff <- function(expression){
    # Get the rank matrix
    rank_matrix <- dirac.rank_matrix(expression[gene_index,])
    # Get the p1 matching score
    p1.matching_score <- dirac.rank_matching_score.vector(rank_matrix,
                                                          p1.rank_template)
    p2.matching_score <- dirac.rank_matching_score.vector(rank_matrix,
                                                          p2.rank_template)
    # Return the difference between the two
    p1.matching_score-p2.matching_score
  }
  # Based on desired return type, define a function
  if(return == "name"){
    if(is.null(names)){
      names <- c("phenotype1", "phenotype2")
    }
    return(
      function(expression){
        # Get the rank matching score difference
        rmd <- rank_matching.diff(expression)
        ifelse(rmd>0, names[1], names[2])
      }
    )
  } else if(return=="is phenotype1"){
    return(
      function(expression){
        # Get the rank matching score difference
        rmd <- rank_matching.diff(expression)
        (rmd>0)
      }
    )
  } else if(return == "is phenotype2"){
    return(
      function(expression){
        # Get the rank matching score difference
        rmd <- rank_matching.diff(expression)
        (rmd<=0)
      }
    )
  } else if(return=="phenotype number"){
    return(
      function(expression){
        # Get the rank matching score difference
        rmd <- rank_matching.diff(expression)
        ifelse(rmd>0, 1, 2)
      }
    )
  } else {
    return(
      function(expression){
        # Get the rank matching score difference
        rmd <- rank_matching.diff(expression)
        rmd
      }
    )
  }
}


# Helper Functions --------------------------------------------------------




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
#' @param rank_matrix Rank matrix, from dirac.rank_matrix
#'
#' @return Numeric vector, rank template
#'
#' @examples
#' dirac.rank_template(matrix(1:16,ncol=4, nrow=4))
#'
#' @export
dirac.rank_template <- function(rank_matrix){
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
#' rank_template <- dirac.rank_template(rank_matrix)
#' dirac.rank_matching_score.vector(rank_matrix, rank_template)
#' @export
dirac.rank_matching_score.vector <- function(rank_matrix, rank_template){
  apply(rank_matrix, 2, dirac.rank_matching_score,
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
#' rank_template <- dirac.rank_template(rank_matrix)
#' dirac.rank_conservation_index(rank_matrix, rank_template)
#' @export
dirac.rank_conservation_index <- function(rank_matrix, rank_template){
  mean(dirac.rank_matching_score.vector(rank_matrix = rank_matrix,
                                   rank_template = rank_template))
}

#' Find the Rank Difference Score
#'
#' `dirac.rank_difference_score` finds the rank difference score for a network
#' between 2 phenotypes. The score against the rank template for the second
#' phenotype is subtracted from the score against the rank template for the
#' first phenotype.
#'
#' @param expression Numeric matrix with gene expression data, rows are
#'    genes, columns are samples.
#' @param phenotype1,phenotype2 Numeric vectors containing the indexes for the
#'    phenotypes.
#' @param gene_index Numeric vector of indices for the genes in the network
#'
#' @returns A list of two vectors, each containing the rank difference scores
#'    for the samples within the phenotypes.
#'
#' @examples
#' expression=matrix(1:16, ncol=4, nrow=4)
#' dirac.rank_difference_score(expression, phenotype1=c(1,4), phenotype2=
#'  c(2,3), gene_index=c(1,2,4))
#' @export
dirac.rank_difference_score <- function(expression, phenotype1, phenotype2,
                                        gene_index){
  # Find the rank matrix for each phenotype
  p1.rank_matrix <- dirac.rank_matrix(expression[gene_index, phenotype1])
  p2.rank_matrix <- dirac.rank_matrix(expression[gene_index, phenotype2])
  # Find the rank templates for the gene network in each of the phenotypes
  p1.rank_template <- dirac.rank_template(p1.rank_matrix)
  p2.rank_template <- dirac.rank_template(p2.rank_matrix)
  # Find the rank difference score for each phenotype
  p1.rank_difference <- dirac.rank_matching_score.vector(p1.rank_matrix,
                                                         p1.rank_template) -
    dirac.rank_matching_score.vector(p1.rank_matrix, p2.rank_template)
  p2.rank_difference <- dirac.rank_matching_score.vector(p2.rank_matrix,
                                                         p1.rank_template) -
    dirac.rank_matching_score.vector(p2.rank_matrix, p2.rank_template)
  list(p1.rank_difference, p2.rank_difference)
}


#' Find the classification rate for single network
#'
#' `dirac.classification_rate` Finds the classification rate between two
#'    phenotypes based on a gene network.
#'
#' Classification rate is defined as the average between sensitivity and
#'   specificity. classification_rate =
#'   Pr(rank_diff>0|phenotype=A)*Pr(phenotype=A) +
#'   Pr(rank_diff<=0|phenotype=B)*Pr(phenotype=B)
#'
#' Due to the way that the classification rate splits 0's, this value depends
#' on the ordering of the gene_index vector. The order of this vector should
#' therefor be kept consistent.
#'
#' @param expression Numeric matrix with gene expression data, rows are genes
#'     columns are samples.
#' @param phenotype1,phenotype2 Numeric vectors containing the indexes for the
#'    two phenotypes to compare.
#' @param gene_index Numeric vector of index for genes in network
#'
#' @returns Float representing the classification rate based on DIRAC
#' @examples
#' # example code
#' @export
dirac.classification_rate <- function(gene_index, expression, phenotype1,
                                      phenotype2){
  # Get the rank differences
  rank_diff <- dirac.rank_difference_score(expression, phenotype1 = phenotype1,
                                           phenotype2 = phenotype2,
                                           gene_index = gene_index)
  # Unpack the list
  p1.rank_diff <- rank_diff[[1]]
  p2.rank_diff <- rank_diff[[2]]
  # calculate Pr(diff_score>O|phenotype=A), number true positives
  tp <- mean(p1.rank_diff>0)
  # calculate Pr(diff_score<=O|phenotype=B), number true negative
  tn <- mean(p2.rank_diff<=0)
  # Calculate Pr(phenotype=A) and Pr(phenotype=B)
  p1.num <- length(phenotype1)
  p2.num <- length(phenotype2)
  total <- p1.num+p2.num
  p1.prob <- p1.num/total
  p2.prob <- p2.num/total
  # Calculate the classification rate
  (tp*p1.prob+tn*p2.prob)
}

#' Compare classification rate for different gene networks
#'
#' `dirac.classification_rate.compare` compares the classification rate between
#'    multiple gene networks for two phenotypes
#'
#' @param expression Numeric matrix, represents gene expression, rows are genes
#'    columns are samples
#' @param phenotype1,phenotype2 Numeric vectors, contain the indices for the
#'    samples corresponding to each phenotype.
#' @param gene_index List of numeric vectors, each containing the indices for
#'    one of the gene networks.
#' @param parallel Boolean indicating Whether the computation should be run in
#'   parallel. On unix makes use of the R core library parallel, on windows
#'   requires foreach, and doParallel.
#' @param cores Integer, indicates number of cores to use for parallel
#'    calculations. If greater than parallel::detectCores() return value,
#'    will instead be set to said return value.
#' @param as.frame Boolean, whether return value should be data.frame or named
#'    list, see return for more information.
#'
#' @returns
#'    if as.frame is TRUE:
#'      data.frame with gene_network names as one column, and classification
#'      rate as another.
#'    if as.frame is FALSE:
#'      Numeric vector containing the classification rate for each of the
#'      provided gene networks. Named based on names of gene_index argument.
#' @examples
#' # example code
#' @export
dirac.classification_rate.compare <- function(expression, phenotype1,
                                              phenotype2, gene_index,
                                              parallel = TRUE, cores=4,
                                              as.frame=FALSE){
  if(parallel){
    cores = if(parallel::detectCores()>cores) cores else parallel::detectCores()
  }
  os_type <- .Platform$OS.type
  if(parallel && os_type=="unix"){
    res <-  unlist(
      parallel::mclapply(gene_index, dirac.classification_rate,
                         expression=expression, phenotype1=phenotype1,
                         phenotype2=phenotype2,
                         mc.cores=cores)
      )
  } else if(parallel && os_type == "windows"){
    # make cluster
    cl <- parallel::makeCluster(cores)
    # Export needed functions to cluster
    parallel::clusterExport(cl, list("dirac.rank_matching_score",
                                     "dirac.rank_vector",
                                     "dirac.rank_template",
                                     "dirac.rank_matrix",
                                     "dirac.rank_matching_score.vector",
                                     "dirac.rank_difference_score",
                                     "dirac.classification_rate"))
    # Perform calculation
    res <- tryCatch(expr={unlist(parallel::parLapply(cl, gene_index, dirac.classification_rate,
                                      expression=expression,
                                      phenotype1=phenotype1,
                                      phenotype2=phenotype2))},
                    # Whether there is an error or not, stop the cluster
                    finally={parallel::stopCluster(cl)})
  } else if(parallel){
    stop("Unsupported OS for parallel operation")
  } else {
    res <- unlist(
      lapply(gene_index, dirac.classification_rate, expression=expression,
             phenotype1=phenotype1, phenotype2=phenotype2)
    )
  }
  if(!as.frame){
    # Rename the results based on the names of the gene_index list arg
    names(res) <- names(gene_index)
    # Return the renamed results
    return(res)
  }
  data.frame(gene_network=names(gene_index), classification_rate=res)
}


#' Compare classification rate for different gene networks
#'
#' `dirac.classification_rate.best` finds the network with the best
#'    classification rate between multiple gene networks for two phenotypes.
#'
#' @param expression Numeric matrix, represents gene expression, rows are genes
#'    columns are samples
#' @param phenotype1,phenotype2 Numeric vectors, contain the indices for the
#'    samples corresponding to each phenotype.
#' @param gene_index List of numeric vectors, each containing the indices for
#'    one of the gene networks.
#' @param parallel Boolean indicating Whether the computation should be run in
#'   parallel. On unix makes use of the R core library parallel, on windows
#'   requires foreach, and doParallel.
#' @param cores Integer, indicates number of cores to use for parallel
#'    calculations. If greater than parallel::detectCores() return value,
#'    will instead be set to said return value.
#'
#' @returns String, name of best gene network
#' @examples
#' # example code
#' @export
dirac.classification_rate.best <- function(expression, phenotype1, phenotype2,
                                           gene_index, parallel=TRUE, cores=4){
  names(which.max(dirac.classification_rate.compare(expression, phenotype1,
                                              phenotype2, gene_index, parallel,
                                              cores, as.frame=FALSE)))
}


#' Compare two phenotypes rankings using the conservation index
#'
#' `dirac.compare_phenotype.single`computes the absolute difference between the
#' rank conservation index for a gene network between two phenotypes.
#'
#' @param rank_matrix Boolean matrix, columns represent the rank vector for
#'    a given sample
#' @param phenotype1,phenotype2 Numeric vectors, indices of the samples in a
#'    particular phenotype
#' @returns Float, absolute difference in rank conservation index between the
#'    two phenotypes for a given gene network (implied by rank_matrix)
#' @examples
#' #example code
#' @export
dirac.compare_phenotype.single <- function(rank_matrix, phenotype1, phenotype2){
  # Get the rank matrices for each of the phenotypes
  p1.rank_matrix <- rank_matrix[, phenotype1]
  p2.rank_matrix <- rank_matrix[, phenotype2]
  # Get the rank templates for each of the phenotypes
  p1.rank_template <- dirac.rank_template(p1.rank_matrix)
  p2.rank_template <- dirac.rank_template(p2.rank_matrix)
  # Get the rank conservation indices
  p1.rank_conservation_index <- dirac.rank_conservation_index(p1.rank_matrix,
                                                              p1.rank_template)
  p2.rank_conservation_index <- dirac.rank_conservation_index(p2.rank_matrix,
                                                              p2.rank_template)
  # Return the absolute value of the difference between the two
  abs(p1.rank_conservation_index-p2.rank_conservation_index)
}


# Note: When replace is true, just samples combined twice,
#   when replace is false it samples it for p1, and the remainder becomes p2
# The i parameter is the bootstrap iteration, and is mostly for compatibility
#   with apply functions


#' Compare two shuffled phenotypes rankings using the conservation index
#'
#' `dirac.compare_phenotype.shuffle`computes the absolute difference between the
#' rank conservation index for a gene network between two shuffled phenotypes.
#'
#' @param i Integer, represents the bootstrap iteration currently being
#'    performed. Mostly ensures compatibility with apply function.
#' @param rank_matrix Boolean matrix, columns represent the rank vector for
#'    a given sample
#' @param combined Numeric Vector, represents the combined indices for
#'     both phenotypes.
#' @param p1.size,p2.size Integer, represents the size of phenotype1, and
#'    phenotype2 respectively
#' @param replace Whether the sampling should be done with replacement. If it is
#'    TRUE, then two samples of size p1.size and p2.size respectively are taken
#'    with replacement from combined (these two samples can overlap). If FALSE,
#'    then a sample of size p1.size is taken without replacement from combined,
#'    and the remainder of combined is taken as the other sample.
#' @returns Float, absolute difference in rank conservation index between the
#'    two phenotypes for a given gene network (implied by rank_matrix)
#' @examples
#' #example code
#' @export
dirac.compare_phenotype.shuffle <- function(i,rank_matrix, combined, p1.size,
                                            p2.size, replace=TRUE){
  if(replace){
    # Sample combined for p1 indices
    p1.idx <- sample(combined, p1.size, replace=replace)
    # Repeat for p2 indices
    p2.idx <- sample(combined, p2.size, replace=replace)
  } else {
    # sample for p1 indices
    p1.idx <- sample(combined, p1.size, replace=replace)
    # Get all the elements from combined not in p1, and put into p2.idx
    p2.idx <- setdiff(combined,p1.idx)
  }
  dirac.compare_phenotype.single(rank_matrix = rank_matrix, phenotype1 = p1.idx,
                                 phenotype2 = p2.idx)
}



#' Compare two phenotypes using the difference in rank conservation indices
#'
#' `dirac.compare_phenotype` computes the absolute difference between the rank
#' conservation indices for two phenotypes based on a single gene network, and
#' the associated p-value using a bootstrapping approach.
#'
#' @param expression Numeric matrix representing gene expression, rows
#'    represent genes, and columns represent samples.
#' @param phenotype1,phenotype2 Numeric vectors representing the indices of
#'    two phenotypes
#' @param gene_index Numeric vector representing the indices of the genes in
#'    the network.
#' @param bootstrap_iterations Number of bootstrap iterations to run for
#'   computing the null distribution for calculating the p-value
#' @param parallel Boolean determining if calculation should be run in parallel,
#'    on unix systems it uses the mclapply function, and on windows it creates
#'    a socket cluster and uses parLapply.
#' @param cores number of cores to use for parallel computation if desired.
#'    If it is greater than the number of cores available (as determined by
#'    parallel::detectCores), it is instead set to the number of cores
#'    available.
#' @param replace Whether the sampling should be done with replacement. If it is
#'    TRUE, then two samples of size p1.size and p2.size respectively are taken
#'    with replacement from combined (these two samples can overlap). If FALSE,
#'    then a sample of size p1.size is taken without replacement from combined,
#'    and the remainder of combined is taken as the other sample.
#' @param seed Integer, used to set the seed for the random number generator
#'    used for sampling.
#' @returns A named list, with value equal to the absolute difference in rank
#'    conservation scores between the two phenotypes, and p.value equal to the
#'    p-value found with bootstrapping.
#' @examples
#' # example code
#' @export
dirac.compare_phenotype <- function(gene_index, expression, phenotype1,
                                    phenotype2, bootstrap_iterations=1000,
                                    parallel=TRUE, cores=4, replace=TRUE,
                                    seed=NULL){
  # Get size of both of the phenotypes
  p1.size <- length(phenotype1)
  p2.size <- length(phenotype2)
  # Combine the phenotype vectors into one for resampling
  combined <- c(phenotype1, phenotype2)
  # Find the rank matrix
  rank_matrix <- dirac.rank_matrix(expression[gene_index,])
  # Find the difference value for the unshuffled phenotypes
  abs_rci_diff <-  dirac.compare_phenotype.single(expression[gene_index,],
                                                  phenotype1,
                                                  phenotype2)
  # Now perform the bootstrapping to find the null distribution
  if(parallel){
    cores = if(parallel::detectCores()>cores) cores else parallel::detectCores()
  }
  os_type <- .Platform$OS.type
  if(parallel && os_type=="unix"){
    # Set seed for
    set.seed(seed, "L'Ecuyer")
    res <- unlist(parallel::mclapply(1:bootstrap_iterations,
                                    dirac.compare_phenotype.shuffle,
                                    rank_matrix=rank_matrix, combined=combined,
                                    p1.size=p1.size, p2.size=p2.size,
                                    replace=replace, mc.cores = cores))
  } else if(parallel && os_type == "windows"){
    # make cluster
    cl <- parallel::makeCluster(cores)
    # Set the RNG stream seed
    parallel::clusterSetRNGStream(cl, seed)
    # Export Needed functions to cluster
    parallel::clusterExport(cl, list("dirac.rank_vector",
                                     "dirac.rank_matrix",
                                     "dirac.rank_matching_score",
                                     "dirac.rank_matching_score.vector",
                                     "dirac.rank_conservation_index",
                                     "dirac.rank_template",
                                     "dirac.compare_phenotype.single",
                                     "dirac.compare_phenotype.shuffle"))
    # Run the bootstrap
    res <-tryCatch(expr={unlist(parallel::parLapply(cl, 1:bootstrap_iterations,
                                      dirac.compare_phenotype.shuffle,
                                      rank_matrix=rank_matrix,
                                      combined=combined, p1.size=p1.size,
                                      p2.size=p2.size, replace=replace))},
                   finally={parallel::stopCluster(cl)})
  } else if(parallel){
    stop("Unsupported OS for parallel operation")
  } else {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
    res <- unlist(lapply(1:bootstrap_iterations,dirac.compare_phenotype.shuffle,
                         rank_matrix=rank_matrix, combined=combined,
                         p1.size=p1.size, p2.size=p1.size, replace=replace))
  }
  # Now, use the bootstrapped values to create an empirical cdf
  boot_cdf <- stats::ecdf(res)
  # Get the p-value for the value
  p.value <- 1-boot_cdf(abs_rci_diff)
  # Return named list with value being the difference, and p.value being
  #   the p.value calculated using the empirical cdf
  list(value=abs_rci_diff, p.value=p.value)
}











































