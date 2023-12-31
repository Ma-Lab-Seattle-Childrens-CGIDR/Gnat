% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/race.R
\name{race.mean_rank_corr}
\alias{race.mean_rank_corr}
\title{Find the mean Kendall correlation between all samples}
\usage{
race.mean_rank_corr(expression.filtered)
}
\arguments{
\item{expression.filtered}{Numeric Matrix, rows are genes in the gene
network, and columns are samples in the phenotype}
}
\value{
Double, the mean of the Kendall Rank Correlation
}
\description{
Find the mean Kendall correlation between all samples
}
\examples{
# Create example matrix
mat <- matrix(c(1,2,3,4,2,3,1,4,4,1,2,3,1,4,2,3), ncol=4)
# Print the result
print(race.mean_rank_corr(mat))
}
