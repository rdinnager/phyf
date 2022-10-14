## from: https://stackoverflow.com/questions/65635086/how-to-insert-row-in-a-sparse-matrix-efficiently
add_rows <- function(M, v, new_rows = NULL) {
    oldind <- seq(dim(M)[1])-1L   ## all rows (0-indexed)
    ## new row indices: count how many rows are inserted before a given row
    ## (maybe a clearer/better/more efficient way to do this?)
    newind <- oldind + as.integer(rowSums(outer(oldind,v,">=")))
    ## modify dimensions
    M@Dim <- M@Dim + c(length(v),0L)
    ## modify row indices (1L accounts for 0 vs 1-indexing)
    M@i <- newind[M@i+1L]
    
    if(!is.null(new_rows)) {
      M[added_rows(v), ] <- new_rows
    }
    
    M
}

added_rows <- function(v) {
  v + 0:(length(v)-1)
}

split_indexes <- function(x) {
  res <- split(x@i + 1, findInterval(seq_len(Matrix::nnzero(x)), x@p, left.open = TRUE))
  res
}

split_xs <- function(x) {
  res <- split(x@x + 1, findInterval(seq_len(Matrix::nnzero(x)), x@p, left.open = TRUE))
  res
}