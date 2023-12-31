% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/race.R
\name{race.compare_phenotypes}
\alias{race.compare_phenotypes}
\title{Compare phenotypes using the Absolute Difference in Mean Rank Correlation
Coefficient}
\usage{
race.compare_phenotypes(
  gene_network,
  expression,
  phenotype1,
  phenotype2,
  bootstrap_iterations = 1000,
  parallel = TRUE,
  cores = 4,
  replace = TRUE,
  seed = NULL
)
}
\arguments{
\item{gene_network}{Integer vector representing the indices of the genes
in the network of interest within the expression matrix}

\item{expression}{Numeric matrix, represents gene expression with rows
representing genes, and columns representing samples}

\item{phenotype1, phenotype2}{Integer vectors representing the indices
of the samples in the phenotype within the expression matrix}

\item{bootstrap_iterations}{Integer number of iterations for the
bootstrapping}

\item{parallel}{Boolean, whether computations should be done in parallel.
See details for more information.}

\item{cores}{Integer, number of cores to use for parallel operation.}

\item{replace}{Boolean, whether sampling should be done with replacement}

\item{seed}{Integer, seed for RNG (note, even for the same seed, the
RNG is handled differently for parallel operation between different
operating systems).

For parallel operations, on unix operating systems a fork cluster will be
used, while on windows a psock cluster will be used. The RNG for both uses
the L'Ecuyer method. For serial operation, a Mersenne-Twister RNG is used.
This means that results will vary between serial and parallel operation.}
}
\value{
Named list, with p1.mean_rank_correlation, p2.mean_rank_correlation,
absolute_difference, and p.value.
}
\description{
Compare phenotypes using the Absolute Difference in Mean Rank Correlation
Coefficient
}
