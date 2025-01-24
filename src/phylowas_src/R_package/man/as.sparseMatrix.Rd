% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{as.sparseMatrix}
\alias{as.sparseMatrix}
\title{Turns to sparse matrix.}
\usage{
as.sparseMatrix(dense_matrix, use_pointers = F)
}
\arguments{
\item{dense_matrix}{dense matrix.}

\item{use_pointers}{whether or not to use the actual values in the matrix or to just set to 1.}
}
\value{
sparse_matrix
}
\description{
Turns dense matrix to sparse matrix.
}
