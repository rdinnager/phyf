new_pfc <- function(phy) {
  assertthat::assert_that(inherits(phy, "phylo"))

  n_nodes <- ape::Ntip(phy) + ape::Nnode(phy)
  phy <- ape::makeNodeLabel(phy)
  phy2 <- add_internal_tips(phy)
  
  edge_ord <- rev(ape::postorder(phy2))
  node_ord <- c(ape::Ntip(phy2) + 1, phy2$edge[edge_ord, 2])
  
  rtp <- root2tip_binary(phy2)
  rtp <- rtp[node_ord, ]
  
  rownames(rtp) <- c(phy2$tip.label, phy2$node.label)[node_ord]
  colnames(rtp) <- phy2$tip.label
  
  lens <- phy2$edge.length[edge_ord]
  
  rtp <- Matrix::t(rtp)[ , -1] %*% Matrix::Diagonal(length(lens), lens)
  rtp <- Matrix::drop0(rtp)
  
  zero_br <- Matrix::colSums(rtp) != 0
  rtp <- rtp[ , zero_br]
  
  rtp <- rtp[-which(rownames(rtp) == "Node1"), ]
  
  new_lens <- lens[zero_br]
  np <- split_indexes(Matrix::t(rtp))
  
  el <- purrr::map(np,
                   ~ new_lens[.x])
  
  np <- new_pfp(np)
  
  edge_names <- colnames(rtp)
  tip_names <- rownames(rtp)

  new_rcrd(list(pfn = tip_names, pfp = np, pfl = el),
           edge_names = edge_names,
           sparse_rep = rtp,
           class = "pfc")
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
#'  * For `is_pfp()`: An object to test.
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

#' Create a new phylogenetic flow collection object (`pfc`)
#'
#' @param x
#'  * For `pfc()`: A phylogenetic tree in `ape::phylo` format
#'  * For `is_pfc()`: An object to test.
#'
#' @return
#' @export
#'
#' @examples
pfc <- function(x) {
  if(!inherits(x, "phylo")) {
    rlang::warn("x is not a phylo object. Attempting to coerce it.")
    x <- ape::as.phylo(x)
  }
  
  new_pfc(x)
}

#' @export
#' @rdname pfp
is_pfc <- function(x) {
  inherits(x, "pfc")
}

#' @export
format.pfc <- function(x, ...) {
  edge_names <- edge_names(x)
  ch <- setNames(purrr::map2_chr(vctrs::field(x, "pfp"), vctrs::field(x, "pfl"),
                 function(x, y) paste0("root", 
                                       paste0("-", stringr::str_pad(format(y, digits = 1, trim = TRUE),
                                                                    8, side = "both", pad = "-"),
                                              ">",
                                              edge_names[x],
                                              collapse = ""))), 
           vctrs::field(x, "pfn"))
  ch
}

#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.pfc <- function(x, ...) {
  out <- format(x, ...)
  pillar::new_pillar_shaft_simple(out, shorten = "mid", min_width = 15)
}

edge_names <- function(x) {
  attr(x, "edge_names")  
}

tip_names <- function(x) {
  attr(x, "tip_names")  
}