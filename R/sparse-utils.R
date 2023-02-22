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

split_indexes <- function(x, empty = NULL) {
  
  if(!is.null(empty)) {
    empty2 <- Matrix::colSums(x) == 0
    if(any(empty2)) {
      x[1, empty2] <- 1e-20
    }
  }
  res <- split(x@i + 1, findInterval(seq_len(Matrix::nnzero(x)), x@p, left.open = TRUE))
  
  if(!is.null(empty) && any(empty2)) {
    res[empty2] <- list(empty)
  }
  res
}

split_xs <- function(x, empty = NULL) {
  ## deal with empty columns
  if(!is.null(empty)) {
    empty2 <- Matrix::colSums(x) == 0
    if(any(empty2)) {
      x[1, empty2] <- 1e-20
    }
  }
  res <- split(x@x, findInterval(seq_len(Matrix::nnzero(x)), x@p, left.open = TRUE))
  
  if(!is.null(empty) && any(empty2)) {
    res[empty2] <- list(empty)
  }
  res
}

force_dgCMatrix <- function(x) {
  if(inherits(x, "sparseVector")) {
    x <- Matrix::Matrix(x, nrow = 1)
  }
  if(!inherits(x, "dgCMatrix")) {
    return(as(x, "CsparseMatrix"))
  } else {
    return(x)
  }
}

dimnames_kron <- function(x, y) {
  
  rnames_x <- rep(rownames(x), each = nrow(y))
  cnames_x <- rep(colnames(x), each = ncol(y))
  rnames_y <- rep(rownames(y), nrow(x))
  cnames_y <- rep(colnames(y), ncol(x))
  rsep <- ifelse(is.null(rnames_x) || is.null(rnames_y), "", ":")
  csep <- ifelse(is.null(cnames_x) || is.null(cnames_y), "", ":")
  list(rownames = paste(rnames_x, rnames_y, sep = rsep),
       colnames = paste(cnames_x, cnames_y, sep = csep))
}

dimnames_row_kron <- function(x, y) {
  
  rnames_x <- rownames(x)
  cnames_x <- rep(colnames(x), each = ncol(y))
  rnames_y <- rownames(y)
  cnames_y <- rep(colnames(y), ncol(x))
  rsep <- ifelse(is.null(rnames_x) || is.null(rnames_y), "", ":")
  csep <- ifelse(is.null(cnames_x) || is.null(cnames_y), "", ":")
  list(rownames = paste(rnames_x, rnames_y, sep = rsep),
       colnames = paste(cnames_x, cnames_y, sep = csep))
}

#' @importClassesFrom Matrix "dgCMatrix"
#' @importFrom methods as new
empty_sparse <- function(nrow = 0, ncol = 0) {
  nrow <- as.integer(nrow)
  ncol <- as.integer(ncol)
  out <- new("dgCMatrix")
  out@Dim <- as.integer(c(nrow, ncol))
  out@p <- integer(ncol + 1L)
  out
}