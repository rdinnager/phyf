new_pfc <- function(phy) {
  assertthat::assert_that(inherits(phy, "phylo"))

  n_nodes <- ape::Ntip(phy) + ape::Nnode(phy)
  phy2 <- add_internal_tips(phy)
  
  edge_ord <- rev(ape::postorder(phy2))
  node_ord <- c(ape::Ntip(phy2) + 1, phy2$edge[edge_ord, 2])
  
  rtp <- root2tip_binary(phy2)
  rtp <- rtp[node_ord, ]
  lens <- phy2$edge.length[edge_ord]
  
  rtp <- Matrix::t(rtp)[ , -1] %*% Matrix::Diagonal(length(lens), lens)
  rtp <- Matrix::drop0(rtp)
  
  zero_br <- Matrix::colSums(rtp) != 0
  rtp <- rtp[ , zero_br]
  
  new_lens <- lens[zero_br]
  np <- split_indexes(Matrix::t(rtp))
  
  el <- purrr::map(np,
                   ~ new_lens[.x])

  new_rcrd()
}

new_pfp <- function(pfp_list) {
  vctrs::vec_assert(pfp_list, list())
  new_vctr(pfp_list, class = "pfp")
}


# for compatibility with the S4 system
methods::setOldClass(c("pfp", "vctrs_vctr"))
methods::setOldClass(c("pfc", "vctrs_vctr"))
methods::setOldClass(c("pf", "vctrs_vctr"))

#' Create a new phylogenetic flow path object (`pfp`)
#'
#' @param x
#'  * For `pfp()`: A list of nodes (integer) the path passes through
#'  * For `is_percent()`: An object to test.
#'
#' @return
#' @export
#'
#' @examples
pfp <- function(x = list()) {
  x <- vec_cast(x, list())
  new_pfp(x)
}

#' @export
#' @rdname pfp
is_pfp <- function(x) {
  inherits(x, "pfp")
}

#' @export
format.pfp <- function(x, ...) {
  purrr::map_chr(x,
                 function(y) paste(y, collapse = "->"))
}
