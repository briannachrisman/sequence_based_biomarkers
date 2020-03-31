library(Matrix)
#' Compute N Primes.
#'
#' This function computes the first n primes.
#'
#' @param n number of primes.
#' @return sequence of primes.
#'
#' @export
GetPrimes = function(n) {
   n <- as.integer(n)
   if(n > 1e6) stop("n too large")
   primes <- rep(TRUE, n)
   primes[1] <- FALSE
   last.prime <- 2L
   for(i in last.prime:floor(sqrt(n)))
   {
      primes[seq.int(2L*last.prime, n, last.prime)] <- FALSE
      last.prime <- last.prime + min(which(primes[(last.prime+1):n]))
   }
   which(primes)
}

#' Recomputes the time.
#'
#' This function prints the time elapsed between now and the last time SetTime was called.
#'
#' @param t current time.
#' @return new time.
#'
#' @export
SetTime = function(t=NA) {
    if (!is.na(t)) {message(paste0("Time elapsed: ", round(Sys.time() - t, 2)))}
    return(Sys.time())
}

#' Turns to sparse matrix.
#'
#' Turns dense matrix to sparse matrix.
#'
#' @param dense_matrix dense matrix.
#' @param use_pointers whether or not to use the actual values in the matrix or to just set to 1.
#' @return sparse_matrix
#'
#' @export
as.sparseMatrix = function(dense_matrix, use_pointers=F) {
    ids = which(dense_matrix!=0, arr.ind = T)
    if (use_pointers) {
        pointers = sapply(1:nrow(ids), function(i) {dense_matrix[ids[i,'row'], ids[i,'col']]})
        sparse_matrix = sparseMatrix(i=ids[,'row'], j=ids[,'col'], x = pointers, dims=dim(dense_matrix))
    } else {
        sparse_matrix = sparseMatrix(i=ids[,'row'], j=ids[,'col'], dims=dim(dense_matrix))
    }
    return(sparse_matrix)
}


as.sparseMatrix = function(dense_matrix, use_pointers=F) {
    ids = which(dense_matrix!=0, arr.ind = T)
    if (use_pointers) {
        pointers = sapply(1:nrow(ids), function(i) {dense_matrix[ids[i,'row'], ids[i,'col']]})
        sparse_matrix = sparseMatrix(i=ids[,'row'], j=ids[,'col'], x = pointers, dims=dim(dense_matrix))
    } else {
        sparse_matrix = sparseMatrix(i=ids[,'row'], j=ids[,'col'], dims=dim(dense_matrix))
    }
    return(sparse_matrix)
}


as.sparseVector = function(denseVector, use_pointers=F) {
    ids = which(denseVector!=0, arr.ind = T)
    if (use_pointers) {
        pointers = denseVector[ids]
        sparse_vector = sparseVector(i=ids, x=pointers, length=length(denseVector))
    } else {
        sparse_vector = sparseVector(i=ids,length=length(denseVector))
    }
    return(sparse_vector)
}

duplicated.sparseMatrix = function (sparse_mat, margin, has_vals=F) {
  j = rep(1:ncol(sparse_mat), diff(sparse_mat@p))
  i = sparse_mat@i + 1
  if (has_vals) {x = sparse_mat@x}
  else {x = sapply(j, function(x){1})}
  if (margin == 1L) { # Duplicated rows.
    names(x) = j
    result = duplicated.default(split(x, i))
  } else if (margin == 2L) { # Duplicated columns.
    names(x) = i
    result = duplicated.default(split(x, j))
  }
  return(result)
}

ElementIdsFromDimensionIds = function(idxs, dims) {
    element_id = rep(0, nrow(idxs))
    for (i in seq(length(dims), 1, by=-1)) {
        if (i>1) {element_id = element_id + prod(dims[seq(i-1, 1, by=-1)]) %*% (idxs[,i]-1)}
        else {element_id = element_id + idxs[,1]}
    }
    return((element_id))
}



