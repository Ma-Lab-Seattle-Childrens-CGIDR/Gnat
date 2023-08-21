
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Gnat (Gene Network Analysis Tools)

<!-- badges: start -->
<!-- badges: end -->

Gnat includes several methods for comparing how tightly gene networks
are controlled between two phenotypes by analyzing how the order of the
genes within a network of interest vary between conditions.

Currently implemented are:

- DIRAC (Differential Rank Conservation): Implemented based on Eddy et
  al.(2010). Creates a pairwise rank template between genes, and
  compares how that varies between phenotypes.
- INFER (Information Entropy of Ranks): Compares the average information
  entropy genes within the network between phenotypes.
- RACE (Rank Correlation Entropy): Compares the average Kendall Rank
  Correlation Coefficient between all samples within a phenotype.
- CRANE (Centroid Rank Entropy): Compares the average euclidean distance
  to the centroid of the rank vectors for all samples within a
  phenotype.

See below for further details.

## Installation

You can install the development version of Gnat from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Ma-Lab-Seattle-Childrens-CGIDR/Gnat")
```

## Example

``` r
library(Gnat)
# Create some example data, with one tightly regulated phenotype, and very
#   loosely regulated
# Create the tightly regulated gene network expression data
expression.p1 <- matrix(c(
  1, 2, 3, 4, 5, 6, 7, 8,
  1, 3, 2, 4, 5, 6, 8, 7,
  2, 1, 3, 4, 6, 5, 7, 8,
  1, 2, 3, 4, 5, 7, 6, 8,
  1, 3, 2, 4, 5, 6, 7, 8,
  3, 2, 1, 4, 5, 6, 7, 8,
  1, 2, 3, 4, 6, 5, 7, 8
), ncol = 7)
# Use random data to simulate a very loosely regulated network
# Note, the magnitudes of expression is very different between the two
#   phenotypes, but since only ranks are considered this won't affect the
#   outcome
expression.p2 <- matrix(runif(7 * 8, min = 1, max = 500), ncol = 7)
# Combine the two phenotypes into one matrix
expression <- cbind(expression.p1, expression.p2)
# Store the indices for both phenotypes
phenotype1 <- 1:7
phenotype2 <- 8:14
# Create a list to represent several gene networks
#   Note that overlaps between the networks are allowed
gene_network_list <- list(
  A = c(1, 3, 5, 6),
  B = c(2, 3, 7, 8),
  C = c(4, 5, 6)
)
# Run the analysis with different methods
DIRAC.results <- DIRAC.compare_phenotypes(
  expression=expression, phenotype1=phenotype1, phenotype2=phenotype2,
  gene_network_list=gene_network_list, bootstrap_iterations = 1000,
  parallel = FALSE, cores = 1,  replace = TRUE, seed = 42, as.frame = TRUE
)
INFER.results <- INFER(expression, gene_network_list, phenotype1, phenotype2,
  bootstrap_iterations = 1000, parallel = FALSE,
  cores = 1, replace = TRUE, seed = 42, as.frame = TRUE
)
CRANE.results <- CRANE(expression, gene_network_list, phenotype1, phenotype2,
  bootstrap_iterations = 1000, parallel = FALSE,
  cores = 1, replace = TRUE, seed = 42, as.frame = TRUE
)
RACE.results <- RACE(expression, gene_network_list, phenotype1, phenotype2,
  bootstrap_iterations = 1000, parallel = FALSE,
  cores = 1, replace = TRUE, seed = 42, as.frame = TRUE
)
print("DIRAC results")
#> [1] "DIRAC results"
print(DIRAC.results)
#>   gene_network p1.rank_conservation_index p2.rank_conservation_index
#> A            A                  0.9285714                  0.6904762
#> B            B                  0.9047619                  0.6428571
#> C            C                  0.9047619                  0.6666667
#>   absolute_difference p.value
#> A           0.2380952   0.001
#> B           0.2619048   0.018
#> C           0.2380952   0.031
print("INFER results")
#> [1] "INFER results"
print(INFER.results)
#>   gene_network p1.rank_entropy p2.rank_entropy absolute_difference p.value
#> A            A       0.7273967        1.753434           1.0260377   0.004
#> B            B       0.7884505        1.672554           0.8841031   0.018
#> C            C       0.5754137        1.378783           0.8033698   0.039
print("CRANE results")
#> [1] "CRANE results"
print(CRANE.results)
#>   gene_network p1.mean_centroid_distance p2.mean_centroid_distance
#> A            A                 0.7350120                  2.042387
#> B            B                 0.8111936                  2.091460
#> C            C                 0.5772300                  1.255934
#>   absolute_difference p.value
#> A           1.3073746   0.000
#> B           1.2802666   0.014
#> C           0.6787035   0.046
print("RACE results")
#> [1] "RACE results"
print(RACE.results)
#>   gene_network p1.mean_rank_correlation p2.mean_rank_correlation
#> A            A                0.7460317               0.04761905
#> B            B                0.7142857              -0.04761905
#> C            C                0.6825397              -0.01587302
#>   absolute_difference p.value
#> A           0.6984127   0.002
#> B           0.7619048   0.012
#> C           0.6984127   0.038
```

# Method Details

## DIRAC

For each sample, DIRAC creates a binary pairwise ordering vector,
$$Gene~1~ < Gene~2~, Gene~1~ < Gene~3~,...$$, with a 1 where TRUE, and 0
where FALSE. For population level DIRAC, a rank ordering template is
created, where if Gene<sub>i</sub> \< Gene<sub>j</sub> in the majority
of the samples, then the corresponding element in the rank template is
1, otherwise it is 0. For each sample, this can then be used to compute
a rank matching score, representing the fraction of values in the
samples ordering vector that match the rank templates. From this, a rank
conservation index is calculated, essentially the average of the rank
matching score for all samples within a phenotype. Additionally, DIRAC
can be used for classification, by comparing how well a sample’s
pairwise rank ordering vector matches the phenotypes rank template. See
Eddy et al., 2010 for further details.

## INFER

This method uses the information entropy of the rank for each gene to
compare how tightly regulated the network is between two phenotypes.
First a rank matrix is computed from the expression data. For each gene
within the network, the information entropy
$H(x) = -\sum_{x \in X} p(x) \log p(x)$ (where $X$ is all possible
values of the rank, i.e. 1,…, $n$, $n$ being the number of genes in the
network) of its ranks across samples is calculated, called the rank
entropy. The average of the rank entropy for each gene in the network is
taken to provide a score for the rank entropy of the network. This score
can then be compared between two phenotypes using a permutation test
approach, where the samples are assigned randomly to the two phenotypes
to create the null distribution, which can be used to calculate the
p-value of the difference in network rank entropy scores between the
phenotypes.

## RACE

This method employs the Kendall Rank Correlation Coefficient between all
of the samples within the phenotype. For each pair of samples, the rank
correlation between the genes within the network of interest is
calculated, and these values are averaged to yield the phenotypes rank
correlation score. The difference between two phenotypes can then be
assessed using a similar permutation testing approach to INFER, where
the samples are randomly assigned to two phenotypes and the difference
in phenotype ran correlation score is used to create a null
distribution. The distribution is then used to calculate the p-value for
the difference in phenotype rank correlations being different.

## CRANE

For each sample, a rank vector is calculated, and then the centroid of
all the rank vectors within a phenotype is found. The average Euclidean
distance between the samples and the centroid provides the phenotype’s
rank grouping score. Then, similar to RACE and INFER a permutation
testing approach can be used to calculate the statistical significance
of the difference in rank grouping scores between two phenotypes.

# References

\[1\] Eddy, J. A., Hood, L., Price, N. D., & Geman, D. (2010).
Identifying tightly regulated and variably expressed networks by
Differential Rank Conservation (DIRAC). PLoS Computational Biology,
6(5). <https://doi.org/10.1371/journal.pcbi.1000792>
