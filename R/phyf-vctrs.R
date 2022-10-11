#' Make a `pfc` object from a `phylo` object
#'
#' @param phy A phylogenetic tree in `ape::phylo`
#' format
#'
#' @return a `pfc`
#' @export
#'
#' @examples
#' as_pfc(ape::rtree(100))
as_pfc <- function(phy) {
  
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
  
  new_pfc(unname(tip_names),
          unname(np),
          unname(el),
          unname(edge_names),
          rtp)
  
}

new_pfc <- function(pfn = character(), 
                    pfpp = pfp(), 
                    pfl = list(),
                    edge_names = character(),
                    sparse_mat = NULL) {
  
  if(is.null(sparse_mat)) {
    if(length(pfn) > 0) {
      sparse_mat <- as_sparse(pfn, pfpp, pfl,
                              edge_names)
    } else {
      sparse_mat <- MatrixExtra::emptySparse(format = "C")
    }
  } 

  new_rcrd(list(pfn = pfn, pfp = pfpp, pfl = pfl),
           edge_names = edge_names,
           sparse_rep = sparse_mat,
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
pfc <- function(pfn = character(), 
                pfpp = pfp(), 
                pfl = list(),
                edge_names = character(),
                sparse_mat = NULL) {
  
  vec_assert(pfn, character())
  vec_assert(pfpp, pfp())
  vec_assert(pfl, list())
  vec_assert(edge_names, character())
  
  if(!is.null(sparse_mat)) {
    if(length(pfn) > 0) {
      assertthat::assert_that(inherits(sparse_mat,
                                     "CsparseMatrix"))
    } 
  } 
  
  new_pfc(pfn,
          pfpp,
          pfl,
          edge_names,
          sparse_mat)
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

#' @export
vec_restore.pfc <- function(x, to, ..., i = NULL) {
  new_pfc(field(x, "pfn"),
          field(x, "pfp"),
          field(x, "pfl"),
          edge_names(to))  
}

#' @export
vec_ptype2.pfc.pfc <- function(x, y, ...) new_pfc()

#' @export
vec_cast.pfc.pfc <- function(x, to, ...) x

#' @export
vec_ptype.pfc <- function(x, ...) {
  spm <- attr(x, "sparse_rep")
  if(any(dim(spm)) != 0) {
    spm <- spm[0, 0]
  }
  new_pfc(field(x, "pfn"), field(x, "pfp"),
          field(x, "pfl"), edge_names(x)[integer()],
          spm)
}

#' @export
vec_ptype_full.pfc <- function(x, ...) {
  paste0("pfc<e:", length(edge_names(x)), ">")
}

edge_names <- function(x) {
  attr(x, "edge_names")  
}

tip_names <- function(x) {
  attr(x, "tip_names")  
}

as_sparse <- function(pfn, pfp, pfl, edge_names) {
  rtp <- Matrix::t(build_rtp(pfp, length(edge_names),
                   lens = pfl))
  rownames(rtp) <- pfn
  colnames(rtp) <- edge_names
  
  rtp <- Matrix::drop0(rtp)
  
  rtp
}

#' @export
pf <- function(x = pfc(), pf_column = "phyf", ...) {
  vec_assert(x, pfc())
  out <- tibble::tibble(node_name = field(x, "pfn"),
                        "{pf_column}" := x)
  tibble::new_tibble(out,
                     pf_column = pf_column,
                     class = "pf")
}

#' @export
as_pf <- function(phy) {
  pf(as_pfc(phy))
}