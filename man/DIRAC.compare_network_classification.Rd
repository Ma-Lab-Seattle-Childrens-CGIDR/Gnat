% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dirac.R
\name{DIRAC.compare_network_classification}
\alias{DIRAC.compare_network_classification}
\title{Compare classification rate for different gene networks}
\usage{
DIRAC.compare_network_classification(
  expression,
  phenotype1,
  phenotype2,
  gene_network_list,
  parallel = TRUE,
  cores = 4,
  as.frame = FALSE
)
}
\arguments{
\item{expression}{Numeric matrix, represents gene expression, rows are genes
columns are samples}

\item{phenotype1, phenotype2}{Numeric vectors, contain the indices for the
samples corresponding to each phenotype.}

\item{gene_network_list}{List of numeric vectors, each containing the indices
for one of the gene networks.}

\item{parallel}{Boolean indicating Whether the computation should be run in
parallel.}

\item{cores}{Integer, indicates number of cores to use for parallel
calculations. If greater than parallel::detectCores() return value,
will instead be set to said return value.}

\item{as.frame}{Boolean, whether the return value should be a dataframe or
a named list.}
}
\value{
if as.frame is TRUE:
data.frame with columns for gene_network name (determined by names of
gene_index list), and classification rate.
if as.frame is FALSE:
Numeric vector containing the classification rate for each of the
provided gene networks. Named based on names of gene_index argument.
}
\description{
\code{DIRAC.compare_network_classification} compares the classification rate
between multiple gene networks for two phenotypes
}
\examples{
# example code
}
