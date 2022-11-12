#' Make a phylogenetic interpolation `pfc`
#' 
#' @param x A `pfc` or `pf` object 
#' @param ... Other arguments passed to or from other methods
#'
#' @export
make_interp <- function(x, ...) {
  UseMethod("make_interp")
}

#' @export
make_interp.pfc <- function(x, ...) {
  descs <- pf_desc(x)
  inds <- pf_desc(pf_indexes(x))
  
  ind_which <- pf_edge_apply(inds, function(x) (as.numeric(as.factor(x)) * 2) - 3)
  
  desc_sp <- pf_as_sparse(descs)
  ind_sp <- pf_as_sparse(ind_which)
  
  new_mat <- desc_sp * ind_sp
  
  zero_br <- Matrix::colSums(structural_sparse(new_mat)) != 0
  new_mat <- new_mat[ , zero_br]
  
  new_pfc <- pf_as_pfc(new_mat, is_tip = rep(FALSE, nrow(new_mat)),
                       internal = rep(TRUE, ncol(new_mat)))
  
  class(new_pfc) <- c("pfc_contrast", class(new_pfc))
  
  new_pfc
}

#' @export
make_interp.phylo <- function(x, ...) {
  make_interp(pf_as_pfc(x))
}