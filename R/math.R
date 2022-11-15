#' Calculate a kronecker product when a `pfc` is the 
#' multiplicand or the multiplier.
#'
#' @param x The kronecker multiplicand.
#' @param y The kronecker multiplier.
#' @param ... Other arguments passed to or from other methods.
#'
#' @return A `pfc` object
#' @export
#'
#' @examples
#' pf_kronecker(rpfc(20), rpfc(20))
pf_kronecker <- function(x, y, ...) {
  UseMethod("pf_kronecker")
}

#' @method pf_kronecker pfc
#' @export
pf_kronecker.pfc <- function(x, y, ...) {
  UseMethod("pf_kronecker.pfc", y)
}

#' @method pf_kronecker.pfc default
#' @export
pf_kronecker.pfc.default <- function(x, y, ...) {
  stop_incompatible_op("pf_kronecker", x, y)
}

#' @method pf_kronecker Matrix
#' @export
pf_kronecker.Matrix <- function(x, y, ...) {
  UseMethod("pf_kronecker.Matrix", y)
}

#' @method pf_kronecker.Matrix default
#' @export
pf_kronecker.Matrix.default <- function(x, y, ...) {
  stop_incompatible_op("pf_kronecker", x, y)
}

#' @method pf_kronecker dgCMatrix
#' @export
pf_kronecker.dgCMatrix <- function(x, y, ...) {
  UseMethod("pf_kronecker.dgCMatrix", y)
}

#' @method pf_kronecker.dgCMatrix default
#' @export
pf_kronecker.dgCMatrix.default <- function(x, y, ...) {
  stop_incompatible_op("pf_kronecker", x, y)
}

#' @method pf_kronecker matrix
#' @export
pf_kronecker.matrix <- function(x, y, ...) {
  UseMethod("pf_kronecker.matrix", y)
}

#' @method pf_kronecker.matrix default
#' @export
pf_kronecker.matrix.default <- function(x, y, ...) {
  stop_incompatible_op("pf_kronecker", x, y)
}

#' @method pf_kronecker.pfc pfc
#' @export
pf_kronecker.pfc.pfc <- function(x, y, ...) {
  
  m1 <- pf_as_sparse(x)
  m2 <- pf_as_sparse(y)

  new_is_tip <- rep(field(x, "is_tip"), each = nrow(m2)) | rep(field(y, "is_tip"), nrow(m1))
  
  m <- force_dgCMatrix(kronecker(m1, m2, ...))
  new_names <- dimnames_kron(m1, m2)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, new_is_tip)
}

#' @method pf_kronecker.pfc dgCMatrix
#' @export
pf_kronecker.pfc.dgCMatrix <- function(x, y, ...) {
  
  m1 <- pf_as_sparse(x)
  
  new_is_tip <- rep(field(x, "is_tip"), each = nrow(y))
  
  m <- force_dgCMatrix(kronecker(m1, y, ...))
  new_names <- dimnames_kron(m1, y)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, new_is_tip)
}

#' @method pf_kronecker.pfc matrix
#' @export
pf_kronecker.pfc.matrix <- function(x, y, ...) {
  
  m1 <- pf_as_sparse(x)
  
  new_is_tip <- rep(field(x, "is_tip"), each = nrow(y))
  
  m <- force_dgCMatrix(kronecker(m1, y, ...))
  new_names <- dimnames_kron(m1, y)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, new_is_tip)
}

#' @method pf_kronecker.pfc Matrix
#' @export
pf_kronecker.pfc.Matrix <- function(x, y, ...) {
  
  m1 <- pf_as_sparse(x)
  
  new_is_tip <- rep(field(x, "is_tip"), each = nrow(y))
  
  m <- force_dgCMatrix(kronecker(m1, y, ...))
  new_names <- dimnames_kron(m1, y)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, new_is_tip)
}


#' @method pf_kronecker.dgCMatrix pfc
#' @export
pf_kronecker.dgCMatrix.pfc <- function(x, y, ...) {
  
  m2 <- pf_as_sparse(y)
  
  new_is_tip <- rep(field(y, "is_tip"), nrow(x))
  
  m <- force_dgCMatrix(kronecker(x, m2, ...))
  new_names <- dimnames_kron(x, m2)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, new_is_tip)
}

#' @method pf_kronecker.matrix pfc
#' @export
pf_kronecker.matrix.pfc <- function(x, y, ...) {
  
  m2 <- pf_as_sparse(y)
  
  new_is_tip <- rep(field(y, "is_tip"), nrow(x))
  
  m <- force_dgCMatrix(kronecker(x, m2, ...))
  new_names <- dimnames_kron(x, m2)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, new_is_tip)
}

#' @method pf_kronecker.Matrix pfc
#' @export
pf_kronecker.Matrix.pfc <- function(x, y, ...) {
  
  m2 <- pf_as_sparse(y)
  
  new_is_tip <- rep(field(y, "is_tip"), nrow(x))
  
  m <- force_dgCMatrix(kronecker(x, m2, ...))
  new_names <- dimnames_kron(x, m2)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, new_is_tip)
}


#' Calculate a rowwise kronecker product when the multiplicand or multiplier
#' is a `pfc`.
#'
#' @param x The rowwise kronecker multiplicand
#' @param y The rowwise kronecker multiplier
#' @param ... Other arguments passed to or from other methods.
#'
#' @return A `pfc` object.
#' @export
#'
#' @examples
#' pf_row_kron(rpfc(20), rpfc(20))
pf_row_kron <- function(x, y, ...) {
  UseMethod("pf_row_kron")
}

#' @method pf_row_kron pfc
#' @export
pf_row_kron.pfc <- function(x, y, ...) {
  UseMethod("pf_row_kron.pfc", y)
}

#' @method pf_row_kron.pfc default
#' @export
pf_row_kron.pfc.default <- function(x, y, ...) {
  stop_incompatible_op("pf_row_kron", x, y)
}

#' @method pf_row_kron Matrix
#' @export
pf_row_kron.Matrix <- function(x, y, ...) {
  UseMethod("pf_row_kron.Matrix", y)
}

#' @method pf_row_kron.Matrix default
#' @export
pf_row_kron.Matrix.default <- function(x, y, ...) {
  stop_incompatible_op("pf_row_kron", x, y)
}

#' @method pf_row_kron dgCMatrix
#' @export
pf_row_kron.dgCMatrix <- function(x, y, ...) {
  UseMethod("pf_row_kron.dgCMatrix", y)
}

#' @method pf_row_kron.dgCMatrix default
#' @export
pf_row_kron.dgCMatrix.default <- function(x, y, ...) {
  stop_incompatible_op("pf_row_kron", x, y)
}

#' @method pf_row_kron matrix
#' @export
pf_row_kron.matrix <- function(x, y, ...) {
  UseMethod("pf_row_kron.matrix", y)
}

#' @method pf_row_kron.matrix default
#' @export
pf_row_kron.matrix.default <- function(x, y, ...) {
  stop_incompatible_op("pf_row_kron", x, y)
}

#' @method pf_row_kron.pfc pfc
#' @export
pf_row_kron.pfc.pfc <- function(x, y, ...) {
  
  m1 <- pf_as_sparse(x)
  m2 <- pf_as_sparse(y)

  new_is_tip <- field(x, "is_tip") | field(y, "is_tip")
  new_internal <- attr(x, "internal") | attr(y, "internal")
  
  m <- force_dgCMatrix(inlabru::row_kron(m1, m2, ...))
  new_names <- dimnames_row_kron(m1, m2)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, is_tip = new_is_tip, internal = new_internal)
}

#' @method pf_row_kron.pfc dgCMatrix
#' @export
pf_row_kron.pfc.dgCMatrix <- function(x, y, ...) {
  
  m1 <- pf_as_sparse(x)
  
  new_is_tip <- field(x, "is_tip")
  new_internal <- attr(x, "internal")
  
  m <- force_dgCMatrix(inlabru::row_kron(m1, y, ...))
  new_names <- dimnames_row_kron(m1, y)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, is_tip = new_is_tip, internal = new_internal)
}

#' @method pf_row_kron.pfc matrix
#' @export
pf_row_kron.pfc.matrix <- function(x, y, ...) {
  
  m1 <- pf_as_sparse(x)
  
  new_is_tip <- field(x, "is_tip")
  new_internal <- field(x, "is_tip")
  
  m <- force_dgCMatrix(inlabru::row_kron(m1, y, ...))
  new_names <- dimnames_row_kron(m1, y)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, is_tip = new_is_tip, internal = new_internal)
}

#' @method pf_row_kron.pfc Matrix
#' @export
pf_row_kron.pfc.Matrix <- function(x, y, ...) {
  
  m1 <- pf_as_sparse(x)
  
  new_is_tip <- field(x, "is_tip")
  new_internal <- attr(x, "internal")
  
  m <- force_dgCMatrix(inlabru::row_kron(m1, y, ...))
  new_names <- dimnames_row_kron(m1, y)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, is_tip = new_is_tip, internal = new_internal)
}


#' @method pf_row_kron.dgCMatrix pfc
#' @export
pf_row_kron.dgCMatrix.pfc <- function(x, y, ...) {
  
  m2 <- pf_as_sparse(y)
  
  new_is_tip <- field(y, "is_tip")
  new_internal <- attr(y, "internal")
  
  m <- force_dgCMatrix(inlabru::row_kron(x, m2, ...))
  new_names <- dimnames_row_kron(x, m2)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, is_tip = new_is_tip, internal = new_internal)
}

#' @method pf_row_kron.matrix pfc
#' @export
pf_row_kron.matrix.pfc <- function(x, y, ...) {
  
  m2 <- pf_as_sparse(y)
  
  new_is_tip <- field(y, "is_tip")
  new_internal <- attr(y, "internal")
  
  m <- force_dgCMatrix(inlabru::row_kron(x, m2, ...))
  new_names <- dimnames_row_kron(x, m2)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, is_tip = new_is_tip, internal = new_internal)
}

#' @method pf_row_kron.Matrix pfc
#' @export
pf_row_kron.Matrix.pfc <- function(x, y, ...) {
  
  m2 <- pf_as_sparse(y)
  
  new_is_tip <- field(y, "is_tip")
  new_internal <- attr(y, "internal")
  
  m <- force_dgCMatrix(inlabru::row_kron(x, m2, ...))
  new_names <- dimnames_row_kron(x, m2)
  
  rownames(m) <- new_names$rownames
  colnames(m) <- new_names$colnames
  
  pf_as_pfc(m, is_tip = new_is_tip, internal = new_internal)
}


#' @method vec_arith pfc
#' @export
vec_arith.pfc <- function(op, x, y, ...) {
  UseMethod("vec_arith.pfc", y)
}

#' @method vec_arith.pfc default
#' @export
vec_arith.pfc.default <- function(op, x, y, ...) {
  stop_incompatible_op(op, x, y)
}

#' @method vec_arith.pfc pfc
#' @export
vec_arith.pfc.pfc <- function(op, x, y, ...) {
  
  if(!structure_equal(x, y)) {
    if(op == '*') {
      rlang::warn("pfcs have different structures, doing a rowwise kronecker product instead of elementwise.",
                  class = "warning_asterisk_kronecker")
      return(pf_row_kron(x, y))
    }
    rlang::warn("You are doing arithmetic on phylogenetic flows with a different structures. Are you sure this is what you want to do?")
  }
  
  if(!identical(field(x, "is_tip"), field(y, "is_tip"))) {
    rlang::abort("Both pfc objects must have the same elements as tips")
  }
  if(!identical(field(x, "pfn"), field(y, "pfn"))) {
    rlang::abort("Both pfc objects must have the same labels")
  }
  
  m1 <- pf_as_sparse(x)
  m2 <- pf_as_sparse(y)
  
  m <- vec_arith_sparse(op, m1, m2)
  pf_as_pfc(m, field(x, "is_tip"))
}

#' @method vec_arith.pfc numeric
#' @export
vec_arith.pfc.numeric <- function(op, x, y, ...) {
  
  m1 <- pf_as_sparse(x)
  
  if(length(y) > 1) {
    if(length(y) == nrow(m1)) {
      m <- vec_arith_sparse(op, m1, structural_sparse(m1) * y)
      return(pf_as_pfc(m, field(x, "is_tip")))  
    } else {
      rlang::abort("Arithmetic with a pfc and a vector only works if the length of the vector is 1 or equal to the number of flows in the `pfc`")
    }
  }
  
  m <- m1
  m@x <- vec_arith_base(op, m1@x, y)
  pf_as_pfc(m, field(x, "is_tip"))
  
}

#' @method vec_arith.pfc matrix
#' @export
vec_arith.pfc.matrix <- function(op, x, y, ...) {
  
  m1 <- pf_as_sparse(x)
  if(all(dim(y) == dim(m1))) {
    m <- vec_arith_sparse(op, m1, structural_sparse(m1) * y)
    return(pf_as_pfc(m, field(x, "is_tip")))  
  }
  if(!any(dim(y) == 1)) {
    rlang::abort("One dimension of the matrix must be one.")
  }
  if(dim(y)[1] == nrow(m1)) {
     m <- vec_arith_sparse(op, m1, structural_sparse(m1) * as.vector(y))
     return(pf_as_pfc(m, field(x, "is_tip")))  
  } else {
    if(dim(y)[2] == ncol(m1)) {
      m <- vec_arith_sparse(op, m1, Matrix::t(structural_sparse(Matrix::t(m1)) * as.vector(y)))
      return(pf_as_pfc(m, field(x, "is_tip"))) 
    } else {
      rlang::abort("Arithmetic with a pfc and a matrix only works if one of the dims is equal to one of the `pfc` dims")
    }
  }

  
  m <- m1
  m@x <- vec_arith_base(op, m1@x, y)
  pf_as_pfc(m, field(x, "is_tip"))
  
}

#' @method vec_arith.numeric pfc
#' @export
vec_arith.numeric.pfc <- function(op, x, y, ...) {
  
  m1 <- pf_as_sparse(y)
  
  if(length(x) > 1) {
    if(length(x) == nrow(m1)) {
      m <- vec_arith_sparse(structural_sparse(m1) * x, m1)
      return( pf_as_pfc(m, field(y, "is_tip")))  
    } else {
      rlang::abort("Arithmetic with a pfc and a vector only works if the length of the vector is 1 or eqal to the number of flows in the `pfc`")
    }
  }
  
  m <- m1
  m@x <- vec_arith_base(op, x, m1@x)
  pf_as_pfc(m, field(y, "is_tip"))
  
}

#' @method vec_arith.pfc dgCMatrix
#' @export
vec_arith.pfc.dgCMatrix <- function(op, x, y, ...) {
  
  if(!identical(field(x, "pfn"), rownames(y))) {
    rlang::abort("objects must have the same labels")
  }
  m1 <- pf_as_sparse(x)
  # if(!structure_equal(m1, y)) {
  #   rlang::warn("You are doing arithmetic on phylogenetic flows with a different structures. Are you sure this is what you want to do?")
  # }
  m <- vec_arith_sparse(op, m1, y)
  pf_as_pfc(m, field(x, "is_tip"))
  
}


vec_arith_sparse <- function(op, x, y) {
  op_fn <- getExportedValue("base", op)
  op_fn(x, y)
}

#' @export
vec_math.pfc <- function(.fn, .x, ...) {
  switch(.fn,
    sum = pf_apply_sparse_fn(.x, Matrix::colSums, 
                             is_tip = any(field(.x, "is_tip")), 
                             sparseResult = TRUE, na.rm = TRUE),
    mean = pf_apply_sparse_fn(.x, Matrix::colMeans,
                              is_tip = any(field(.x, "is_tip")), 
                             sparseResult = TRUE, na.rm = TRUE),
    pf_apply_sparse_fn(.x, getExportedValue("base", .fn), ...)
  )
}

vec_math_sparse <- function(fn, x) {
  math_fn <- getExportedValue("base", fn)
  math_fn(x)
}

#' Calculate phylogenetic variance covariance matrix from `pfc`
#'
#' @param x A `pfc` object
#' @param ... Other arguments for future extensions.
#'
#' @return A `sparseMatrix`
#' @export
#'
#' @examples
#' pf_vcv(rpfc(10))
pf_vcv <- function(x, ...) {
  m <- pf_as_sparse(sqrt(x))
  covar <- m %*% Matrix::t(m)
  covar
}
