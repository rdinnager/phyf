new_pfc <- function(phy) {
  assertthat::assert_that(inherits(phy, "phylo"))

  n_nodes <- ape::Ntip(phy) + ape::Nnode(phy)
  phy2 <- add_internal_tips(phy)
  
  edge_ord <- rev(ape::postorder(phy))
  node_ord <- c(ape::Ntip(phy) + 1, phy$edge[edge_ord, 2])
  
  if(!is.null(phy$edge.length)) {
    lens <- phy$edge.length[edge_ord]
  } else {
    lens <- rep(1, nrow(phy$edge))
  }
  
  if(!is.null(phy$root.edge)) {
    lens <- c(phy$root.edge, lens)
  } else {
    lens <- c(0, lens)
  }
  
  np <- ape::nodepath(phy)
  el <- purrr::map(np,
                   ~ lens[.x])

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
