% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datagen.R
\name{datagen.swap}
\alias{datagen.swap}
\title{Perform nswaps on the expression data, and return the result}
\usage{
datagen.swap(expression, nswaps)
}
\arguments{
\item{expression}{Numeric matrix of gene expression data, rows represent
genes, and columns represent samples}

\item{nswaps}{Integer, number of swaps to perform}
}
\value{
Numeric matrix, expression data but with nswaps performed
}
\description{
For each swap, a random sample, and two random genes are chosen. The values
for gene expression for those two genes, in that sample are then swapped.
}
