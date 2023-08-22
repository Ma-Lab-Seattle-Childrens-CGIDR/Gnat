# Main Function -----------------------------------------------------------


#' Compute the rank entropy and difference probability for a list of gene
#' networks
#'
#' @param expression Numeric matrix representing gene expression data,
#'    rows represent genes and columns represent samples
#' @param gene_network_list List of integer vectors, each vector is the
#'    indices for a genes network within the expression matrix. The names of
#'    this list are used for naming the gene networks.
#' @param phenotype1,phenotype2 Integer vectors representing the indices of the
#'    two phenotypes within the expression matrix.
#' @param bootstrap_iterations Integer, number of bootstrap_iterations to
#'    perform
#' @param parallel Boolean, whether the calculation should be performed in
#'    parallel. See details for more information.
#' @param cores Integer, number of cores to use for parallel operations
#' @param replace Boolean, whether to sample with replacement during the
#'    bootstrapping.
#' @param seed Integer, seed to use for the phenotype sampling.
#' @param as.frame Boolean, whether return should be a dataframe or not. If
#'    TRUE, will return a dataframe with gene_network, p1.rank_entropy,
#'    p2.rank_entropy, absolute_difference, and p.value as columns. If FALSE,
#'    will return this as a named list instead.
#'
#' @return dataframe or named list depending on as.frame.
#'    If as.frame is TRUE:
#'      returns a dataframe with gene_network, p1.rank_entropy,
#'      p2.rank_entropy, absolute_difference and p.value columns.
#'    If as.frame is FALSE:
#'      returns the data as a named list. Names are the same as names for
#'      gene_network_list argument, and each element is a named list with
#'      gene_network, p1.rank_entropy, p2.rank_entropy, absolute_difference, and
#'      p.value.
#'
#' For parallel operation, on unix the parallel mclapply is used, which employs
#' a fork cluster. On windows, a psock cluster is used.
#'
#' @export
#'
#' @examples
INFER <- function(expression, gene_network_list,
                  phenotype1, phenotype2,
                  bootstrap_iterations = 1000,
                  parallel = TRUE,
                  cores = 4, replace = TRUE, seed = NULL, as.frame = TRUE) {
    # Ensure that the gene_network_list is named
    gene_network_list <- ensure_named(gene_network_list, prefix = "gene_network_")
    # Run the infer.compare_phenotypes
    res_list <- lapply(gene_network_list, infer.compare_phenotypes,
        expression = expression,
        phenotype1 = phenotype1,
        phenotype2 = phenotype2,
        bootstrap_iterations = bootstrap_iterations,
        parallel = parallel, cores = cores, replace = replace,
        seed = seed
    )
    if (!as.frame) {
        names(res_list) <- names(gene_network_list)
        return(res_list)
    } else {
        p1.rank_entropy.list <- sapply(res_list, function(x) x$p1.rank_entropy)
        p2.rank_entropy.list <- sapply(res_list, function(x) x$p2.rank_entropy)
        values.list <- sapply(res_list, function(x) x$absolute_difference)
        p.values.list <- sapply(res_list, function(x) x$p.value)
        gene_networks.names <- names(gene_network_list)
        res_frame <- data.frame(
            gene_network = gene_networks.names,
            p1.rank_entropy = p1.rank_entropy.list,
            p2.rank_entropy = p2.rank_entropy.list,
            absolute_difference = values.list,
            p.value = p.values.list
        )
        return(res_frame)
    }
}




# Helper Functions --------------------------------------------------------



#' Determine Entropy of a Vector
#'
#' `entropy.vector` finds the information entropy of a vector
#'
#' @param vec Numeric vector (or other class which implements table and
#'    length methods), for which the entropy is calculated
#'
#' @return
#'    Float, representing the entropy of the vector
#' @export
#'
#' @examples
#' # Create an example vector
#' example_vector <- c(1, 2, 3, 2, 1, 2, 3, 4, 5, 6, 6, 7, 1)
#' # Print the entropy of the vector
#' print(entropy.vector(example_vector))
entropy.vector <- function(vec) {
    t <- table(vec) / length(vec)
    -sum(t * log2(t))
}


#' Compute entropy for a matrix
#'
#' @param mat Numeric matrix to calculate the entropy of
#' @param margin Integer, determines which axis the calculation is performed
#'    along, 1 for rows, 2 for columns
#'
#' @return Numeric Vector, each element representing the entropy of
#'    corresponding row or column
#' @export
#'
#' @examples
#' # Create an example matrix
#' example_matrix <- matrix(rnorm(4 * 4), ncol = 4)
#' # Print result vector
#' print(entropy.matrix(example_matrix, 2))
entropy.matrix <- function(mat, margin = 1) {
    apply(mat, MARGIN = margin, entropy.vector)
}

#' Compute the rank matrix
#'
#' @param expression.filtered Numeric matrix representing the gene expression,
#'    with genes as rows and samples as columns. Genes should only be those
#'    in the network of interest, and samples should only be those in phenotype
#'    of interest. Can function with the transpose of this by changing
#'    margin to 2, but other methods don't work with this alternate
#'    orientation.
#' @param margin Integer representing whether function rank should be applied
#'    along rows (1) or columns (2)
#'
#' @return Integer array representing the ranks of the gene expression within
#'    a sample.
#' @export
#'
#' @examples
infer.rank_matrix <- function(expression.filtered, margin = 2) {
    apply(expression.filtered, MARGIN = margin, rank, ties.method = "first")
}


#' Find the rank entropy for each gene in a network
#'
#' `infer.gene.entropy` Determines the rank entropy for each gene in a network
#'
#' @param expression Numeric matrix representing gene expression, with
#'    rows representing genes and columns representing samples.
#' @param gene_network Integer vector representing indices of genes in the
#'    network in the expression matrix.
#' @param phenotype Integer vector representing indices of samples in the
#'    phenotype in the expression matrix.
#'
#' @return Numeric vector representing the entropy for each gene in the network
#'    in the given phenotype.
#' @export
#'
#' @examples
infer.gene.entropy <- function(expression, gene_network,
                               phenotype) {
    entropy.matrix(infer.rank_matrix(expression[gene_network, phenotype]))
}


#' Find the entropy for gene network ranks
#'
#' `infer.gene_network.entropy` finds the rank entropy for a gene network
#'
#' @param expression Numeric matrix representing gene expression values,
#'    with rows representing genes, and columns representing samples
#' @param gene_network Integer vector representing the indices of the genes
#'    in the gene network within the expression matrix
#' @param phenotype Integer vector representing the indices of the samples
#'    in the phenotype within the expression matrix
#'
#' @return Float, representing the rank entropy of the gene_network in the
#'    given phenotype
#' @export
#'
#' @examples
infer.gene_network.entropy <- function(expression, gene_network, phenotype) {
    mean(infer.gene.entropy(expression, gene_network, phenotype))
}


#' Calculate the absolute difference in rank entropies between two phenotypes
#'
#' `infer.compare_phenotypes.single` computes the rank entropy for each of two
#' phenotypes, and returns the absolute value of the difference between them
#'
#' @param rank_matrix Integer matrix representing the gene expression ranks,
#'    rows representing genes, and columns representing samples
#' @param phenotype1,phenotype2 Integer vectors representing the indices of the
#'    two phenotypes within the rank_matrix
#'
#' @return Float, the absolute value of the difference between the rank
#'    entropies for each of the phenotypes.
#' @export
#'
#' @examples
infer.compare_phenotypes.single <- function(rank_matrix,
                                            phenotype1,
                                            phenotype2) {
    # Get the rank matrices for each phenotype
    p1.rank_matrix <- rank_matrix[, phenotype1]
    p2.rank_matrix <- rank_matrix[, phenotype2]
    # Get the rank entropies for each phenotype
    p1.rank_entropy <- mean(entropy.matrix(p1.rank_matrix, margin = 1))
    p2.rank_entropy <- mean(entropy.matrix(p2.rank_matrix, margin = 1))
    # Return the absolute value of the differences between the rank entropies
    abs(p1.rank_entropy - p2.rank_entropy)
}





#' Find the absolute difference between rank entropies for shuffled phenotypes
#'
#' @param i Integer, Bootstrap iteration
#' @param rank_matrix Integer Matrix, represents the within sample ranks of gene
#'    expression, rows representing genes, columns representing samples
#' @param combined Integer vector, represents the concatenation of the two
#'    phenotype index vectors
#' @param p1.size,p2.size Integers, the number of samples within each phenotype
#' @param replace Boolean, whether sampling should be done with replacement.
#'
#' If replace is TRUE, will take samples for each phenotype, so the two
#' phenotypes may overlap. If replace is FALSE, will take sample for phenotype1,
#' then use the remaining indices for phenotype2, so there won't be any overlap
#'
#' @return double, The absolute difference between the rank entropies of the
#'    shuffled phenotypes
#' @export
#'
#' @examples
infer.compare_phenotypes.shuffle <- function(i, rank_matrix, combined,
                                             p1.size, p2.size, replace = TRUE) {
    if (replace) {
        # Sample combined for p1 indices
        p1.idx <- sample(combined, p1.size, replace = replace)
        # Repeat for p2 indices
        p2.idx <- sample(combined, p2.size, replace = replace)
    } else {
        # Sample for p1 indices
        p1.idx <- sample(combined, p1.size, replace = replace)
        # Get remaining elements as those will be p2 indices
        p2.idx <- setdiff(combined, p1.idx)
    }
    # Compute the absolute difference in rank entropies between the shuffled
    #   phenotypes
    infer.compare_phenotypes.single(
        rank_matrix = rank_matrix,
        phenotype1 = p1.idx,
        phenotype2 = p2.idx
    )
}

#' Compute the rank entropy difference and p-value between two phenotypes for
#'    a gene network
#'
#' @param gene_network Integer vector, the indices of the gene network within
#'    the expression matrix
#' @param expression Numeric matrix, the gene expression data. Rows represent
#'    genes, and columns represent samples.
#' @param phenotype1,phenotype2 Integer vectors representing the indices of the
#'    phenotypes withint the expression matrix
#' @param bootstrap_iterations Integer, number of bootstrap iterations to
#'    perform
#' @param parallel Boolean, whether parallel operation is desired
#' @param cores Integer, number of cores to use for parallel operation
#' @param replace Boolean, whether to sample with replacement during
#'    bootstrapping
#' @param seed Integer, seed to use for random number generation and sampling
#'
#' @return Named list with p1.rank_entropy, p2.rank_entropy,
#'    absolute_difference, and p.value.
#' @export
#'
#' @examples
infer.compare_phenotypes <- function(gene_network, expression,
                                     phenotype1, phenotype2,
                                     bootstrap_iterations = 1000, parallel = TRUE,
                                     cores = 4, replace = TRUE, seed = NULL) {
    # Get the size of each of the phenotypes
    p1.size <- length(phenotype1)
    p2.size <- length(phenotype2)
    # Get the combined phenotype index vector
    combined <- c(phenotype1, phenotype2)
    # Compute the rank_matrix
    rank_matrix <- infer.rank_matrix(expression[gene_network, ])
    # Find the difference value for the unshuffled phenotype
    abs_re_diff <- infer.compare_phenotypes.single(
        rank_matrix = rank_matrix,
        phenotype1 = phenotype1,
        phenotype2 = phenotype2
    )
    p1.rank_entropy <- infer.gene_network.entropy(
        expression = expression,
        gene_network = gene_network,
        phenotype = phenotype1
    )
    p2.rank_entropy <- infer.gene_network.entropy(
        expression = expression,
        gene_network = gene_network,
        phenotype = phenotype2
    )
    # Perform bootstrapping for the null distribution
    if (parallel) {
        cores <- if (parallel::detectCores() > cores) cores else parallel::detectCores()
    }
    os_type <- .Platform$OS.type
    if (parallel) {
        if (os_type == "unix") {
            # Set seed
            set.seed(seed, "L'Ecuyer")
            # Use multicore apply
            res <- unlist(parallel::mclapply(1:bootstrap_iterations,
                infer.compare_phenotypes.shuffle,
                rank_matrix = rank_matrix,
                combined = combined, p1.size = p1.size,
                p2.size = p2.size, replace = replace,
                mc.cores = cores
            ))
        } else if (os_type == "windows") {
            # make the cluster
            cl <- parallel::makeCluster(cores)
            # Set the RNG stream seed
            parallel::clusterSetRNGStream(cl, seed)
            # Export needed functions
            parallel::clusterExport(cl, list(
                "entropy.vector",
                "entropy.matrix",
                "infer.compare_phenotypes.single",
                "infer.compare_phenotypes.shuffle"
            ))
            # Run the bootstrap
            res <- tryCatch(expr = {
                unlist(parallel::parLapply(cl, 1:bootstrap_iterations,
                    infer.compare_phenotypes.shuffle,
                    rank_matrix = rank_matrix,
                    combined = combined,
                    p1.size = p1.size,
                    p2.size = p2.size,
                    replace = replace
                ))
            }, finally = {
                parallel::stopCluster(cl)
            })
        } else {
            stop("Unsupported OS for parallel operation")
        }
    } else {
        set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
        res <- unlist(lapply(1:bootstrap_iterations,
            infer.compare_phenotypes.shuffle,
            rank_matrix = rank_matrix,
            combined = combined,
            p1.size = p1.size,
            p2.size = p2.size,
            replace = replace
        ))
    }
    # Now create the ECDF
    boot_cdf <- stats::ecdf(res)
    # Get the p-value for the absolute rank entropy difference
    p.value <- 1 - boot_cdf(abs_re_diff)
    # Return named list with value being the absolute difference, and
    #   p.value being the p-value calculated using the ecdf
    list(
        p1.rank_entropy = p1.rank_entropy, p2.rank_entropy = p2.rank_entropy,
        absolute_difference = abs_re_diff, p.value = p.value
    )
}
