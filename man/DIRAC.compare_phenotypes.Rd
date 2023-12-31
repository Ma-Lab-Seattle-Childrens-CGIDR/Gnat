% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dirac.R
\name{DIRAC.compare_phenotypes}
\alias{DIRAC.compare_phenotypes}
\title{Compare two phenotypes using the difference in rank conservation indices}
\usage{
DIRAC.compare_phenotypes(
  expression,
  gene_network_list,
  phenotype1,
  phenotype2,
  bootstrap_iterations = 1000,
  parallel = TRUE,
  cores = 4,
  replace = TRUE,
  seed = NULL,
  as.frame = TRUE
)
}
\arguments{
\item{expression}{Numeric matrix representing gene expression, rows
represent genes, and columns represent samples.}

\item{gene_network_list}{A list of numeric vectors, each representing the
indices of a gene network.}

\item{phenotype1, phenotype2}{Numeric vectors representing the indices of
two phenotypes}

\item{bootstrap_iterations}{Number of bootstrap iterations to run for computing
the null distribution for calculating the p-value}

\item{parallel}{Boolean determining if calculation should be run in parallel,
on unix systems it uses the mclapply function, and on windows it creates
a socket cluster and uses parLapply.}

\item{cores}{number of cores to use for parallel computation if desired.
If it is greater than the number of cores available (as determined by
parallel::detectCores), it is instead set to the number of cores
available.}

\item{replace}{Whether the sampling should be done with replacement. If it is
TRUE, then two samples of size p1.size and p2.size respectively are taken
with replacement from combined (these two samples can overlap). If FALSE,
then a sample of size p1.size is taken without replacement from combined,
and the remainder of combined is taken as the other sample.}

\item{seed}{Integer, used to set the seed for the random number generator
used for sampling.}

\item{as.frame}{Boolean, whether return value should be dataframe or named
list, see return for more information.}
}
\value{
If as.frame is TRUE:
a dataframe with columns \code{gene_network} (determined by the names of the
gene_index list),  \code{p1.rank_conservation_index},
\code{p2.rank_conservation_index}, \code{absolute_difference}, and \code{p.value}.
If as.frame is FALSE:
A list named according to the names in the gene_index argument, with
each entry being a named list, with \verb{$p1.rank_conservation_index},
\verb{$p2.rank_conservation_index}, \verb{$absolute_difference}, and \verb{$p.value}
\verb{$p.value} equal to the p-value for \verb{$value} found with bootstrapping.
}
\description{
\code{DIRAC.compare_phenotypes} computes the absolute difference between the rank
conservation indices for two phenotypes based on a single gene network, and
the associated p-value using a bootstrapping approach.
}
\details{
When using the parallel computation option, the bootstrapping is
parallelized, but testing each gene network is done serially.
}
\examples{
# example code
}
