% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infer.R
\name{entropy.vector}
\alias{entropy.vector}
\title{Determine Entropy of a Vector}
\usage{
entropy.vector(vec)
}
\arguments{
\item{vec}{Numeric vector (or other class which implements table and
length methods), for which the entropy is calculated}
}
\value{
Float, representing the entropy of the vector
}
\description{
\code{entropy.vector} finds the information entropy of a vector
}
\examples{
  # Create an example vector
  example_vector <- c(1,2,3,2,1,2,3,4,5,6,6,7,1)
  # Print the entropy of the vector
  print(entropy.vector(example_vector))
}
