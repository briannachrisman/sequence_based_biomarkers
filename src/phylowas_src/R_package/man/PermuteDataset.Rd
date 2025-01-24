% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{PermuteDataset}
\alias{PermuteDataset}
\title{Computes a matrix illustrating how to permute the dataset.}
\usage{
PermuteDataset(ids, n_permute_iter, prob = c(0.5, 0.5))
}
\arguments{
\item{ids}{identifier of rows you want to keep together. IE if you have a family, and always want to flip their phenotypes together, those two rows should have the same ID.}

\item{n_permute_iter}{number of permutations you want to run.}
}
\value{
permute_switch_pheno a matrix with length(ids) rows and n_permute_iter columns, showing how to permute the dataset for each permutation.TrimTaxaVar
}
\description{
This function creates a matrix with length(ids) rows and n_permute_iter columns, showing how to permute the dataset
for each permutation (whether or not to flip the phenotype at row (row # for permutation #).
}
