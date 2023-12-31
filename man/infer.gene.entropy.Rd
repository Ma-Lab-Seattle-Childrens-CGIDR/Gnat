% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infer.R
\name{infer.gene.entropy}
\alias{infer.gene.entropy}
\title{Find the rank entropy for each gene in a network}
\usage{
infer.gene.entropy(expression, gene_network, phenotype)
}
\arguments{
\item{expression}{Numeric matrix representing gene expression, with
rows representing genes and columns representing samples.}

\item{gene_network}{Integer vector representing indices of genes in the
network in the expression matrix.}

\item{phenotype}{Integer vector representing indices of samples in the
phenotype in the expression matrix.}
}
\value{
Numeric vector representing the entropy for each gene in the network
in the given phenotype.
}
\description{
\code{infer.gene.entropy} Determines the rank entropy for each gene in a network
}
