% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infer.R
\name{entropy.matrix}
\alias{entropy.matrix}
\title{Compute entropy for a matrix}
\usage{
entropy.matrix(mat, margin = 1)
}
\arguments{
\item{mat}{Numeric matrix to calculate the entropy of}

\item{margin}{Integer, determines which axis the calculation is performed
along, 1 for rows, 2 for columns}
}
\value{
Numeric Vector, each element representing the entropy of
corresponding row or column
}
\description{
Compute entropy for a matrix
}
\examples{
  # Create an example matrix
  example_matrix <- matrix(rnorm(4*4), ncol=4)
  # Print result vector
  print(entropy.matrix(example_matrix, 2))
}
