% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{ComputePairedP}
\alias{ComputePairedP}
\title{Computes the p-value for a paired wilcoxon rank test.}
\usage{
ComputePairedP(person_groups, phenotypes, shuff_switch_sign, i_group)
}
\arguments{
\item{person_groups}{table of person (rows) vs abundance (columns) for groups of taxa or variants.}

\item{phenotypes}{list of phenotypes for each person (T/F).}

\item{shuff_switch_sign}{table of pair (row) vs iteration (col) of whether or not to switch the phenotypes of each pair for that permutation.}

\item{i_group}{which group (col) of person_groups you would like to compute the p-value for.}
}
\value{
p-value.
}
\description{
This function computes the p-value for a paired wilcoxon rank sum test, and computes p-values from the permuted data.
}
