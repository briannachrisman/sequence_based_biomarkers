% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{ComputeUnpairedP}
\alias{ComputeUnpairedP}
\title{Computes the p-value for an unpaired wilcoxon rank test.}
\usage{
ComputeUnpairedP(person_groups, phenotypes, shuff_phenos, i_group)
}
\arguments{
\item{person_groups}{table of person (rows) vs abundance (columns) for groups of taxa or variants.}

\item{phenotypes}{list of phenotypes for each person (T/F).}

\item{i_group}{which group (col) of person_groups you would like to compute the p-value for.}

\item{shuff_switch_sign}{table of person (row) vs iteration (col) of whether or not to switch the phenotypes of each person for that permutation.}
}
\value{
p-value.
}
\description{
This function computes the p-value for an unpaired wilcoxon rank sum test, and computes p-values from the permuted data.
}
