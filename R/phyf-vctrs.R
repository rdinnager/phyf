#' Make a `pfc` object from a `phylo` object
#'
#' @param x A phylogenetic tree in `ape::phylo`
#' format
#' @param ... Arguments passed to or from other methods
#'
#' @return a `pfc`
#' @export
#'
#' @examples
#' pf_as_pfc(ape::rtree(100))
pf_as_pfc <- function(x, ...) {
  UseMethod("pf_as_pfc")
}

#' @importFrom methods setClass
#' @export
pf_as_pfc.phylo <- function(x, ...) {
  
  assertthat::assert_that(inherits(x, "phylo"))

  n_nodes <- ape::Ntip(x) + ape::Nnode(x)
  x <- ape::makeNodeLabel(x)
  x2 <- add_internal_tips(x)
  
  edge_ord <- rev(ape::postorder(x2))
  node_ord <- c(ape::Ntip(x2) + 1, x2$edge[edge_ord, 2])
  
  rtp <- root2tip_binary(x2)
  rtp <- rtp[node_ord, ]
  
  rownames(rtp) <- c(x2$tip.label, x2$node.label)[node_ord]
  colnames(rtp) <- x2$tip.label
  
  lens <- x2$edge.length[edge_ord]
  
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
  
  is_tip <- tip_names %in% x$tip.label
  
  internal <- edge_names %in% x$node.label
  
  new_pfc(tip_names,
          np,
          el,
          is_tip,
          edge_names,
          new_lens,
          internal,
          rtp)
  
}

#' @export
pf_as_pfc.dgCMatrix <- function(x, is_tip = NULL, internal = NULL, ...) {
  
  ## deal with possible empty cols
  empty <- Matrix::rowSums(x) == 0
  if(any(empty)) {
    x[empty, 1] <- 1e-20
  }
  
  np <- split_indexes(Matrix::t(x))
  el <- split_xs(Matrix::t(x))
  edge_names <- colnames(x)
  tip_names <- rownames(x)
  
  np[empty] <- list(integer(0))
  el[empty] <- list(numeric(0))
  
  new_pfc(tip_names,
          np,
          el,
          is_tip,
          edge_names,
          purrr::map_dbl(split_xs(x), mean),
          internal,
          Matrix::drop0(x, 1e-10))
  
}

new_pfc <- function(pfn = character(), 
                    pfpp = pfp(), 
                    pfl = list(),
                    is_tip = logical(),
                    edge_names = character(),
                    edge_lengths = double(),
                    internal = logical(),
                    sparse_mat = NULL) {
  
  if(is.null(sparse_mat)) {
    if(length(pfn) > 0) {
      sparse_mat <- as_sparse(pfn, pfpp, pfl,
                              edge_names)
    } else {
      sparse_mat <- MatrixExtra::emptySparse(format = "C")
    }
  } 

  new_rcrd(list(pfn = unname(pfn), pfp = unname(pfpp), pfl = unname(pfl),
                is_tip = unname(is_tip)),
           edge_names = edge_names,
           edge_lengths = edge_lengths,
           internal = internal,
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
#' @return A `pfp` object
#' @export
#'
#' @examples
#' pfp(ape::nodepath(ape::rtree(100)))
pfp <- function(x = list()) {
  x <- vec_cast(x, list())
  new_pfp(x)
}

#' @export
#' @rdname pfp
pf_is_pfp <- function(x) {
  inherits(x, "pfp")
}

#' @export
format.pfp <- function(x, ...) {
  purrr::map_chr(x,
                 function(y) paste(y, collapse = "->"))
}

#' Create a new phylogenetic flow collection object (`pfc`)
#'
#' @param pfn A character vector of names for each phylogenetic flow
#' @param pfpp A vector of phylogenetic flow paths of class `pfp` 
#' @param pfl A list of phylogenetic flow feature. Should be numeric.
#' @param is_tip A logical vector specifying whether the phylogenetic
#' flow reaches the phylogeny's tips 
#' @param edge_names A character vector of names for the phylogeny's
#' edges
#' @param edge_lengths A numeric vector with a representation of 
#' overall features on each edge 
#' @param internal A logical vector specifying whether each edge
#' is internal (not leading to a tip) 
#' @param sparse_mat A sparse matrix representation of the phylogenetic
#' flow collection. Can be left `NULL`, in which case it is constructed
#' from the other arguments.
#'
#' @return A `pfc` object
#' @export
pfc <- function(pfn = character(), 
                pfpp = pfp(), 
                pfl = list(),
                is_tip = logical(),
                edge_names = character(),
                edge_lengths = double(),
                internal = logical(),
                sparse_mat = NULL) {
  
  vec_assert(pfn, character())
  vec_assert(pfpp, pfp())
  vec_assert(pfl, list())
  vec_assert(edge_names, character())
  vec_assert(edge_lengths, double())
  vec_assert(is_tip, logical())
  vec_assert(internal, logical())
  
  if(!is.null(sparse_mat)) {
    if(length(pfn) > 0) {
      assertthat::assert_that(inherits(sparse_mat,
                                     "CsparseMatrix"))
    } 
  } 
  
  new_pfc(pfn,
          pfpp,
          pfl,
          is_tip,
          edge_names,
          edge_lengths,
          internal,
          sparse_mat)
}

#' @export
#' @param x An object to be tested
#' @rdname pfc
pf_is_pfc <- function(x) {
  inherits(x, "pfc")
}

#' @importFrom stats setNames
#' @export
format.pfc <- function(x, ...) {
  edge_names <- edge_names(x)
  ch <- setNames(purrr::map2_chr(vctrs::field(x, "pfp"), vctrs::field(x, "pfl"),
                 function(x, y) paste0("", 
                                       paste0("-", stringr::str_pad(format(y, digits = 3, trim = TRUE),
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
          field(x, "is_tip"),
          edge_names(to),
          edge_lengths(to),
          internal_edges(to))  
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
          field(x, "pfl"), field(x, "is_tip"), 
          edge_names(x)[integer()],
          edge_lengths(x)[integer()],
          internal_edges(x)[integer()],
          spm)
}

#' @export
vec_ptype_full.pfc <- function(x, ...) {
  paste0("pfc<e:", length(edge_names(x)), ">")
}

edge_names <- function(x) {
  attr(x, "edge_names")  
}

edge_lengths <- function(x) {
  attr(x, "edge_lengths")  
}

tip_names <- function(x) {
  field(x, "pfn")  
}

internal_edges <- function(x) {
  attr(x, "internal")
}


#' Extract or assign into the internal or terminal edges of a `pfc`
#' 
#' `pf_internal` extracts or assigns a `pfc` or `dgCMatrix` into the internal
#' edges of a `pfc`. 
#' 
#' `pf_terminal` extracts or assigns a `pfc` or `dgCMatrix` into the terminal
#' edges of a `pfc`. 
#'
#' @param x a `pfc` object
#' @param value a `pfc` or `dgCMatrix` with the
#' same dimensions as `pf_internal(x)` or `pf_terminal(x)`
#' @param ... Arguments passed to or from other methods.
#'
#' @return a `pfc` object with flows truncated to internal nodes.
#' @export
#'
#' @examples
#' pfc1 <- rpfc(100)
#' # Pagel's lambda transformation:
#' lambda <- 0.5
#' lambda_pfc <- pfc1
#' root2node_lens <- pf_flow_sum(lambda_pfc)
#' pf_internal(lambda_pfc) <- pf_internal(lambda_pfc) * lambda
#' pf_terminal(lambda_pfc) <- pf_terminal(lambda_pfc) * (root2node_lens / pf_flow_sum(lambda_pfc))
#' plot(lambda_pfc)
pf_internal <- function(x, ...) {
  UseMethod("pf_internal")
}

#' @export
pf_internal.pfc <- function(x, ...) {
  edges <- internal_edges(x)
  is_tip <- field(x, "is_tip")
  m <- pf_as_sparse(x)[ , edges]
  pf_as_pfc(m, is_tip = is_tip, internal = edges[edges])
}


#' @rdname pf_internal
#' @export
`pf_internal<-` <- function(x, value) {
  edges <- internal_edges(x)
  is_tip <- field(x, "is_tip")
  m <- pf_as_sparse(x)
  m[ , edges] <- pf_as_sparse(value)
  pf_as_pfc(m, is_tip = is_tip, internal = edges)
}

#' @rdname pf_internal
#' @export
pf_terminal <- function(x, ...) {
  UseMethod("pf_terminal")
}

#' @export
pf_terminal.pfc <- function(x, ...) {
  edges <- !internal_edges(x)
  is_tip <- field(x, "is_tip")
  m <- pf_as_sparse(x)[ , edges]
  pf_as_pfc(m, is_tip = is_tip, internal = edges[edges])
}


#' @rdname pf_internal
#' @export
`pf_terminal<-` <- function(x, value) {
  edges <- !internal_edges(x)
  is_tip <- field(x, "is_tip")
  m <- pf_as_sparse(x)
  m[ , edges] <- pf_as_sparse(value)
  pf_as_pfc(m, is_tip = is_tip, internal = edges)
}

as_sparse <- function(pfn, pfp, pfl, edge_names) {
  rtp <- Matrix::t(build_rtp(pfp, length(edge_names),
                   lens = pfl))
  rownames(rtp) <- pfn
  colnames(rtp) <- edge_names
  
  rtp <- Matrix::drop0(rtp)
  
  rtp
}

#' `pf` object constructor
#' 
#' @param x a `pfc` object
#' @param pf_column Name of the column to hold `x`. Default: `'phlo'` 
#' @param ... Reserved for future extensions. Not used. 
#'
#' @return a `pf` object
#'
#' @importFrom rlang :=
#' @export
#' @examples
#' pf(rpfc(100))
pf <- function(x = pfc(), pf_column = "phlo", ...) {
  vec_assert(x, pfc())
  out <- tibble::tibble(label = field(x, "pfn"),
                        is_tip = field(x, "is_tip"),
                        "{pf_column}" := x)
  tibble::new_tibble(out,
                     pf_column = pf_column,
                     label_column = "label",
                     is_tip_column = "is_tip", 
                     class = "pf")
}


#' Convert a `ape::phylo` object to a `pf` object
#'
#' @param phy A `phylo` object
#'
#' @return a `pf` object with branch lengths as features
#' @export
#'
#' @examples
#' pf_as_pf(ape::rtree(100))
pf_as_pf <- function(phy) {
  pf(pf_as_pfc(phy))
}


#' Convert a `pf` or `pfc` object to a `ape::phylo` object
#' 
#' This function attempts to convert a `pf` or `pfc` object
#' to a phylogenyin `phylo` format (from package `{ape}`).
#' It will use the feature as branch lengths, and if edges
#' have multiple feature values, they will be aggregated 
#' by averaging them. Note that this function can fail if the
#' `pf` or `pfc` does not have a tree-like structure. An example
#' of this would be an 'interaction' `pfc` (as generated by
#' `pf_interaction`)
#'
#' @param x a `pf` or `pfc` object to be converted
#' @param ... Other arguments passed to or from other methods 
#'
#' @return A `phylo` object
#' @export
#'
#' @examples
#' pf_as_phylo(rpfc(100))
pf_as_phylo <- function(x, ...) {
  UseMethod("pf_as_phylo")
}

#' @export
pf_as_phylo.pfc <- function(x, ...) {
  
  tips <- x[field(x, "is_tip")]
  nodes <- field(tips, "pfp")
  edge_mats <- purrr::map(nodes,
                          flow_to_edge_list)
  edge_mat <- unique(do.call(rbind, edge_mats))
  
  edge_nams <- c("root", edge_names(x))
  edges_present <- unique(as.vector(edge_mat))
  tip_nodes <- edges_present[edge_nams[edges_present] %in% field(tips, "pfn")]#  which(edge_nams %in% field(tips, "pfn")) + 1
  other_nodes <- edges_present[!edges_present %in% tip_nodes]

  tip_node_changer <- seq_along(tip_nodes)
  names(tip_node_changer) <- tip_nodes
  
  other_node_changer <- max(tip_node_changer) + seq_along(other_nodes)
  names(other_node_changer) <- other_nodes
  
  changer <- c(tip_node_changer, other_node_changer)
  
  edge_mat2 <- matrix(changer[as.character(edge_mat)], ncol = 2)
  
  tip_lab <- edge_nams[tip_nodes]
  node_lab <- edge_nams[other_nodes]

  new_phylo <- list()
  new_phylo$edge <- edge_mat2
  new_phylo$tip.label <- tip_lab
  new_phylo$node.label <- node_lab
  new_phylo$Nnode <- length(node_lab)
  new_phylo$edge.length <- edge_lengths(x)[edge_mat[ , 2] - 1]
  class(new_phylo) <- "phylo"
  
  new_phylo <- ape::collapse.singles(new_phylo)
  
  # ape::checkValidPhylo(new_phylo)
  # plot(new_phylo)
  # phytools::cotangleplot(test_tree, new_phylo, type = "phylogram")
  # plot(phytools::cophylo(test_tree, new_phylo, rotate = FALSE))
  
  new_phylo
  
}

#' @export
pf_as_phylo.pf <- function(x, ...) {
  pf_as_phylo(pf_phyloflow(x, ...))
}


#' Convert a `pf` or `pfc` object to a sparse matrix
#' representation.
#'
#' @param x A `pfc` or `pf` object
#' @param ... Arguments passed to or from other methods
#'
#' @return A `dgCMatrix` object (from the `{Matrix} package`). See
#' `Matrix::sparseMatrix()` for details.
#' @export
#'
#' @examples
#' pf_as_sparse(rpfc(100))
pf_as_sparse <- function(x, ...) {
  UseMethod("pf_as_sparse")
}

#' @export
pf_as_sparse.pfc <- function(x, ...) {
  attr(x, "sparse_rep")
}

#' @export
pf_as_sparse.pf <- function(x, ...) {
  pf_as_sparse(pf_phyloflow(x))
}

#' @export
pf_as_sparse.dgCMatrix <- function(x, ...) {
  x
}


#' Extracts the phylogenetic flow column of an `pf` object
#'
#' @param x A `pf` object
#'
#' @return A `pfc` object (phylogenetic flow collection)
#' @export
#'
#' @examples
#' pf_phyloflow(avonet)
pf_phyloflow <- function(x) {
  
  if(!inherits(x, "pf")) {
    rlang::abort("x must be a pf object")
  }
  
  pf_col <- attr(x, "pf_column")
  x[[pf_col]]
  
}

`pf_phyloflow<-` <- function(x, value) {
  
  if(!inherits(x, "pf")) {
    rlang::abort("x must be a pf object")
  }
  
  if(is.character(value)) {
    if(length(value) >= 1) {
      attr(x, "pf_column") <- value[1]
    }
    if(length(value) >= 2) {
      attr(x, "label") <- value[2]
    }
    if(length(value) >= 3) {
      attr(x, "is_tip") <- value[3]
    }
  } else {
    pf_col <- attr(x, "pf_column")
    x[[pf_col]] <- value
  }
  
  x
  
}

#' @importFrom utils head
#' @export
plot.pf <- function(x, n = 4, ...) {
  phy <- pf_as_phylo(x)
  labels <- c(attr(x, "label_column"))
  ignore <- c(labels,
              attr(x, "is_tip_column"),
              attr(x, "pf_column"))
  
  if((ncol(x) - 3 ) == 0) {
    p <- phytools::contMap(phy, main = "A phylogenetic tree with no features",
                           x = rep(0, ape::Ntip(phy)) %>% setNames(phy$tip.label),
                           method = "user", anc.states = rep(0, ape::Nnode(phy)) %>%
                                                               setNames(ape::Ntip(phy) +
                                                                          seq_len(ape::Nnode(phy)))
                           )
    plot(p)
    return(p)
  }
  
  pcols <- setdiff(names(x), ignore)
  
  if((ncol(x) - 3 ) == 1) {
    sel <- c(labels, pcols)
    df <- dplyr::tibble("{labels}" := c(phy$tip.label,
                                       phy$node.label)) %>%
      dplyr::left_join(x %>%
                         dplyr::select(dplyr::all_of(sel))) 
  }
  
  if((ncol(x) - 3 ) > 1) {
    sel <- c(labels, head(pcols, n))
    df <- dplyr::tibble("{labels}" := c(phy$tip.label,
                                       phy$node.label)) %>%
      dplyr::left_join(x %>%
                         dplyr::select(dplyr::all_of(sel)))
    withr::with_par(list(mfrow = c(2, 2)), 
                    {
                      "not implemented yet!"
                    })
  }
  
}

#' @export
plot.pfc <- function(x, scale_bar_len = 0.5, ...) {
  phy <- pf_as_phylo(x)
  plot(phy, ...) 
  ape::add.scale.bar(length = scale_bar_len)
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
      return(vec_arith_sparse(op, m1, structural_sparse(m1) * y))  
    } else {
      rlang::abort("Arithmetic with a pfc and a vector only works if the length of the vector is 1 or eqal to the number of flows in the `pfc`")
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
      return(vec_arith_sparse(structural_sparse(m1) * x, m1))  
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

#' @method vec_arith.numeric pfc
#' @export
vec_arith.numeric.pfc <- function(op, x, y, ...) {
  
  m1 <- pf_as_sparse(y)
  
  m <- m1
  m@x <- vec_arith_base(op, x, m1@x)
  pf_as_pfc(m, field(y, "is_tip"))
  
}

vec_arith_sparse <- function(op, x, y) {
  op_fn <- getExportedValue("base", op)
  op_fn(x, y)
}

structural_sparse <- function(m) {
  m@x <- 1
  m
}

# structure_equal <- function(pfc1, pfc2) {
#   m1 <- pf_as_sparse(pfc1)
#   m2 <- pf_as_sparse(pfc2)
#   m1@x <- 1
#   m2@x <- 1
#   identical(m1, m2)
# }

structure_equal <- function(pfc1, pfc2) {
  all(vec_equal(unlist(field(pfc1, "pfp")),
            unlist(field(pfc2, "pfp"))))
}

#' Calculate the sum of features for each flow
#' 
#' Sums the phylogenetic flow features for each flow
#' in a `pfc`, returning a vector of values.
#'
#' @param x A `pfc` object
#' @param ... Other arguments to pass to or from other methods
#'
#' @return A vector equal to the length of `x` with the 
#' sum of each flow's features
#' @export
#'
#' @examples
#' pf_flow_sum(rpfc(100))
pf_flow_sum <- function(x, ...) {
  UseMethod("pf_flow_sum")
}

#' @export
pf_flow_sum.default <- function(x, ...) {
  rlang::abort("flow_sums only works on pf and pfc objects")
}

#' @export
pf_flow_sum.pfc <- function(x, ...) {
  Matrix::rowSums(pf_as_sparse(x))
}

#' Calculate the sum of features for each flow
#' 
#' Cumulative sum the phylogenetic flow features for each flow
#' in a `pfc`, returning a `pfc`.
#'
#' @param x A `pfc` object
#' @param ... Other arguments to pass to or from other methods
#'
#' @return A `pfc` object with cumulative sums of features in
#' place of the original features
#' @export
#'
#' @examples
#' pf_flow_cumsum(rpfc(100))
pf_flow_cumsum <- function(x, ...) {
  UseMethod("pf_flow_cumsum")
}

#' @export
pf_flow_cumsum.default <- function(x, ...) {
  rlang::abort("flow_sums only works on pf and pfc objects")
}

#' @export
pf_flow_cumsum.pfc <- function(x, ...) {
  field(x, "pfl") <- purrr::map(field(x, "pfl"),
                                cumsum)
  refresh_features(x)
}

#' Return a `pfc` with ancestral features
#'  
#' @param x A `pfc` object
#' @param replace Value to replace edge feature with no
#' ancestral value (such as when the ancestor is the root)
#'
#' @return A new `pfc` with the same structure as `x`, but with
#' the features of each edge's ancestors instead
#' 
#' @export
#'
#' @examples
#' pf_anc(rpfc(100))
pf_anc <- function(x, replace = 0) {
  field(x, "pfl") <- purrr::map(field(x, "pfl"),
                                function(y) c(replace, 
                                              y[-length(y)]))
  refresh_features(x)
}

#' Generate a random tree and return it as a 
#' `pfc` object
#'
#' @param n Number of tips in the generated tree.
#' @param method A function to generate a tree. Default
#' is `ape::rcoal`. See `ape::rtree()` for more options.
#' @param ... Additional arguments to pass to `method` 
#'
#' @return a `pfc` object
#' @export
#'
#' @examples
#' plot(rpfc(100))
rpfc <- function(n, method = ape::rcoal, ...) {
  pf_as_pfc(method(n, ...))
}

refresh_features <- function(x) {
  p <- new_pfc(field(x, "pfn"),
               field(x, "pfp"),
               field(x, "pfl"),
               field(x, "is_tip"),
               edge_names(x),
               edge_lengths(x),
               internal_edges(x))
  attr(p, "edge_lengths") <- purrr::map_dbl(split_xs(pf_as_sparse(p)), 
                                            mean)
  p
}

#' Scale the phylogenetic flow features to a constant sum. 
#' 
#' Can be used to standardise branch lengths to reasonable values
#'
#' @param x a `pfc` object
#' @param scale_to The value to scale the sums to. Default is 1.
#'
#' @return a scaled `pfc` object
#' @export
#'
#' @examples
#' pf_scale_flow_sums(rpfc(100))
pf_scale_flow_sum <- function(x, scale_to = 1) {
  
  field(x, "pfl") <- purrr::map(field(x, "pfl"),
                                function(y) (y / sum(y)) * scale_to)
  refresh_features(x)
  
}
