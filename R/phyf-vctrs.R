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
  
  rtp <- Matrix::t(rtp)[ , -1] 
  #cnames <- colnames(rtp)
  rtp <- rtp %*% Matrix::Diagonal(length(lens), lens, names = colnames(rtp))
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
          internal,
          rtp)
  
}

#' @export
pf_as_pfc.dgCMatrix <- function(x, is_tip = NULL, internal = NULL, ...) {
  
  x <- force_dgCMatrix(x)
  
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
          internal,
          Matrix::drop0(x, 1e-10))
  
}

#' @export
pf_as_pfc.Matrix <- function(x, is_tip = NULL, internal = NULL, ...) {
  
  x <- force_dgCMatrix(x)
  
  pf_as_pfc(x, is_tip = is_tip, internal = internal, ...)
  
}

#' @export
pf_as_pfc.matrix <- function(x, is_tip = NULL, internal = NULL, ...) {
  
  x <- force_dgCMatrix(x)
  
  pf_as_pfc(x, is_tip = is_tip, internal = internal, ...)
  
}


new_pfc <- function(pfn = character(), 
                    pfpp = pfp(), 
                    pfl = list(),
                    is_tip = logical(),
                    edge_names = character(),
                    internal = logical(),
                    sparse_mat = NULL) {
  
  if(is.null(sparse_mat)) {
    if(length(pfn) > 0) {
      sparse_mat <- as_sparse(pfn, pfpp, pfl,
                              edge_names)
    } else {
      sparse_mat <- empty_sparse()
    }
  } 

  new_rcrd(list(pfn = unname(pfn), pfp = unname(pfpp), pfl = unname(pfl),
                is_tip = unname(is_tip)),
           edge_names = edge_names,
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
  if(rlang::is_installed("cli")) {
    purrr::map_chr(x,
                   function(y) paste(y, collapse = "->"))
  } else {
    purrr::map_chr(x,
                   function(y) paste(y, collapse = "->"))
  }
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
                internal = logical(),
                sparse_mat = NULL) {
  
  vec_assert(pfn, character())
  vec_assert(pfpp, pfp())
  vec_assert(pfl, list())
  vec_assert(edge_names, character())
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
format.pfc <- function(x, ansi = length(x) <= 100, colour = length(x) <= 25, ...) {
  
  edge_names <- edge_names(x)
  
  if(rlang::is_installed("cli") & ansi) {
    
    root_sym <- cli::format_inline(cli::col_red("{cli::symbol$circle_double}"))
    line_sym <- cli::format_inline(rep(cli::col_red("{cli::symbol$line}"), 2))
    arrow_sym <- cli::format_inline(line_sym, cli::col_red("{cli::symbol$arrow_right}"))
    
    if(colour) {
      col_range <- range(unlist(field(x, "pfl")))
      col_pal <- scales::col_bin(palette = "viridis",
                                 bins = 10,
                                 domain = col_range)
      
      colors <- purrr::map(field(x, "pfl"),
                           col_pal)
      all_cols <- unique(unlist(colors))
      names(all_cols) <- all_cols
    
      ansi_funs <- purrr::map(all_cols,
                              ~ cli::make_ansi_style(.x, bg = TRUE)) 
      
      ch <- purrr::pmap_chr(list(vctrs::field(x, "pfp"), vctrs::field(x, "pfl"), colors),
                   function(x, y, z) paste0(root_sym, 
                                            paste0(line_sym, 
                                                   purrr::map2_chr(z, format(y, digits = 2, trim = TRUE, width = I(5)),
                                                                   function(z, y) cli::col_white(ansi_funs[[z]](y))),
                                                   arrow_sym, " ",
                                                   if(length(x) > 1)  {
                                                     c(cli::col_blue(cli::style_bold(edge_names[x[-length(x)]])),
                                                       cli::bg_white(cli::style_bold(cli::style_italic(edge_names[x[length(x)]]))))
                                                   } else {
                                                       cli::bg_white(cli::style_bold(cli::style_italic(edge_names[x[length(x)]]))) 
                                                   }, 
                                                   " ",
                                                   collapse = "")))
    } else {
      ch <- purrr::pmap_chr(list(vctrs::field(x, "pfp"), vctrs::field(x, "pfl")),
                            function(x, y, z) paste0(root_sym, 
                                                     paste0(line_sym, 
                                                            cli::bg_white(format(y, digits = 2, trim = TRUE, width = I(5))),
                                                            arrow_sym, " ",
                                                            if(length(x) > 1)  {
                                                              c(cli::col_blue(cli::style_bold(edge_names[x[-length(x)]])),
                                                                cli::bg_white(cli::style_bold(cli::style_italic(edge_names[x[length(x)]]))))
                                                              } else {
                                                                cli::bg_white(cli::style_bold(cli::style_italic(edge_names[x[length(x)]]))) 
                                                              },
                                                            " ",
                                                            collapse = ""))) 
    }
    
  } else {
    root_sym <- "O"
    line_sym <- "=="
    arrow_sym <- "==>"
    ch <- purrr::map2_chr(vctrs::field(x, "pfp"), vctrs::field(x, "pfl"),
                   function(x, y) paste0(root_sym, 
                                         paste0(line_sym, format(y, digits = 2, trim = TRUE,
                                                                 width = I(5)),
                                                arrow_sym, " ",
                                                edge_names[x], " ",
                                                collapse = "")))
    
  }
  ch
}

#' @importFrom utils head
#' @export
obj_print_data.pfc <- function(x, max = 10, ansi = max <= 100, colour = max <= 25, ...) {
  ch <- format(head(x, max), ansi = ansi, colour = colour)
  ch <- paste0("[", seq_along(ch), "] ", ch)
  cat(ch, sep = "\n\n")
}

#' @export
obj_print_header.pfc <- function(x, max = 10, ...) {
  NextMethod("obj_print_header")
  cat("First ", min(max, length(x)), "phylogenetic flows: \n")
}

#' @method obj_print pfc
#' @importFrom vctrs obj_print
#' @export
obj_print.pfc <- function(x, max = 10, ansi = max <= 100, colour = max <= 25, ...) {
  obj_print_header(x, ...)
  obj_print_data(x, max = max, ansi = ansi, colour = colour, ...)
}

#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.pfc <- function(x, ...) {
  out <- format(x, ansi = length(x) <= 100, colour = length(x) <= 25, ...)
  pillar::new_pillar_shaft_simple(out, shorten = "mid", min_width = 15)
}


#' @export
vec_restore.pfc <- function(x, to, ..., i = NULL) {
  
  new_sm <- as_sparse(field(x, "pfn"),
                      field(x, "pfp"),
                      field(x, "pfl"),
                      edge_names(to))
  
  new_pfc(field(x, "pfn"),
          field(x, "pfp"),
          field(x, "pfl"),
          field(x, "is_tip"),
          edge_names(to),
          internal_edges(to),
          new_sm)

}

#' @method vec_ptype2 pfc
#' @export
vec_ptype2.pfc <- function(x, y, ...) {
  UseMethod("vec_ptype2.pfc", y)
}

#' @method vec_ptype2.pfc default
#' @export
vec_ptype2.pfc.default <- function(x, y, ...) {
  vec_ptype(x, y)
}

#' @method vec_ptype2.pfc pfc
#' @export
vec_ptype2.pfc.pfc <- function(x, y, ...) {
  
  edge_names1 <- edge_names(x)
  edge_names2 <- edge_names(y)
  
  if(length(edge_names1) != length(edge_names2)) {
    rlang::abort("Cannot combine pfc objects with different number of edges.")
  }
  if(any(edge_names1 != edge_names2)) {
    rlang::abort("Cannot combine pfc object with different edge names.")
  }
  
  new_pfc(field(x, "pfn"), field(x, "pfp"),
          field(x, "pfl"), field(x, "is_tip"),
          edge_names(x),
          internal_edges(x),
          pf_as_sparse(x))
}

#' @export
vec_cast.pfc.pfc <- function(x, to, ...) x

#' @export
# vec_ptype.pfc <- function(x, ...) {
#   spm <- attr(x, "sparse_rep")
#   if(any(dim(spm)) != 0) {
#     spm <- spm[0, 0]
#   }
#   new_pfc(field(x, "pfn"), field(x, "pfp"),
#           field(x, "pfl"), field(x, "is_tip"), 
#           edge_names(x)[integer()],
#           edge_lengths(x)[integer()],
#           internal_edges(x)[integer()],
#           spm)
# }

#' @export
vec_ptype_full.pfc <- function(x, ...) {
  paste0("pfc<e:", length(edge_names(x)), ">")
}

edge_names <- function(x) {
  attr(x, "edge_names")  
}

edge_lengths <- function(x, fun, empty = NULL) {
  purrr::map_dbl(split_xs(pf_as_sparse(x), empty = empty), fun)
}

tip_names <- function(x) {
  field(x, "pfn")  
}

internal_edges <- function(x) {
  attr(x, "internal")
}

#' Get only the tip elements of a `pfc`
#'
#' @param x A `pfc` object 
#'
#' @return A `pfc` object with only tip elements
#' @export
#'
#' @examples
#' pf_tips(rpfc(100))
pf_tips <- function(x) {
  x[field(x, "is_tip")]
}

#' Get only the (internal) node elements of a `pfc`
#'
#' @param x A `pfc` object 
#'
#' @return A `pfc` object with only internal node elements
#' @export
#'
#' @examples
#' pf_nodes(rpfc(100))
pf_nodes <- function(x) {
  x[!field(x, "is_tip")]
}

#' Return a vector of labels for a `pfc` object
#'
#' @param x A `pfc` object
#'
#' @return A character vector of labels
#' @export
#'
#' @examples
#' pf_labels(rpfc(100))
pf_labels <- function(x) {
  field(x, "pfn")
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
  
  if(length(purrr::compact(pfp)) == 0) {
    rtp <- empty_sparse(ncol = length(edge_names),
                        nrow = length(pfn))
    rownames(rtp) <- pfn
    colnames(rtp) <- edge_names
    return(rtp)
  }
  
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
  #vec_assert(x, pfc())
  assertthat::assert_that(inherits(x, "pfc"))
  out <- tibble::tibble(label = field(x, "pfn"),
                        is_tip = field(x, "is_tip"),
                        "{pf_column}" := x)
  tibble::new_tibble(out,
                     pf_column = pf_column,
                     label_column = "label",
                     is_tip_column = "is_tip", 
                     class = "pf")
}

#' Convert an object to a `pf` object
#'
#' @param x An object to convert
#' @param ... Other arguments pass to or from other methods.
#'
#' @return a `pf` object with branch lengths as features
#' @export
#'
#' @examples
#' pf_as_pf(ape::rtree(100))
pf_as_pf <- function(x, ...) {
  UseMethod("pf_as_pf")
}

#' @export
pf_as_pf.phylo <- function(x, ...) {
  pf(pf_as_pfc(x, ...))
}

#' @export
pf_as_pf.pfc <- function(x, ...) {
  pf(x, ...)
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
  new_phylo$edge.length <- purrr::map_dbl(split_xs(pf_as_sparse(x)[ , edge_mat[ , 2] - 1]), mean)
  
  class(new_phylo) <- "phylo"
  
  new_phylo <- ape::collapse.singles(new_phylo, root.edge = TRUE)
  
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

#' @export
plot.pfc <- function(x, scale_bar_len = 0.5, ...) {
  phy <- pf_as_phylo(x)
  plot(phy, ...) 
  ape::add.scale.bar(length = scale_bar_len)
}


structural_sparse <- function(m) {
  m@x[] <- 1
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
  if(pf_nedges(pfc1) != pf_nedges(pfc2)) {
    return(FALSE)
  }
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

#' Calculate the cumulative sum of features for each flow
#' 
#' Cumulative sum the phylogenetic flow features for each flow
#' in a `pfc`, returning a `pfc`.
#'
#' @param x A `pfc` object
#' @param direction Which direction to take the cumulative sum along? 
#' "from_root", the default, goes from the root to the terminal nodes. 
#' "to_root" goes in the opposite direction, from the terminal nodes to
#' the root.
#' @param ... Other arguments to pass to or from other methods
#'
#' @return A `pfc` object with cumulative sums of features in
#' place of the original features
#' @export
#'
#' @examples
#' pf_flow_cumsum(rpfc(100))
pf_flow_cumsum <- function(x, direction = c("from_root", "to_root"), ...) {
  UseMethod("pf_flow_cumsum")
}

#' @export
pf_flow_cumsum.default <- function(x, direction = c("from_root", "to_root"), ...) {
  rlang::abort("flow_sums only works on pf and pfc objects")
}

#' @export
pf_flow_cumsum.pfc <- function(x, direction = c("from_root", "to_root"), ...) {
  direction <- match.arg(direction)
  if(direction == "from_root") {
    field(x, "pfl") <- purrr::map(field(x, "pfl"),
                                  cumsum)
  } else {
    field(x, "pfl") <- purrr::map(field(x, "pfl"),
                                  function(x) rev(cumsum(rev(x))))
  }
  refresh_features(x)
}


#' Return a `pfc` with ancestral features
#'  
#' @param x A `pfc` object
#' @param replace Value to replace edge feature with no
#' ancestral value (such as when the ancestor is the root)
#'
#' @return A new `pfc` with the same structure as `x`, but with
#' the features of each edge's ancestor instead
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

#' Return a `pfc` with descendent's features
#' 
#' Note that in the context of a phylogenetic flow, each element of the flow
#' only has a single descendent, which is the next edge on the sequence between
#' the root and the terminal node of the flow in question.
#'  
#' @param x A `pfc` object
#' @param replace Value to replace edge feature with no
#' descendent values (such as when the edges are the tip edges).
#'
#' @return A new `pfc` with the same structure as `x`, but with
#' the features of each edge's descendant instead
#' 
#' @export
#'
#' @examples
#' pf_desc(rpfc(100))
pf_desc <- function(x, replace = 0) {
  field(x, "pfl") <- purrr::map(field(x, "pfl"),
                                function(y) c(y[-1],
                                              replace))
  refresh_features(x)
}

#' Replace feature with edge index in a `pfc`
#' 
#' This is similar to base R's `col()` function, it returns a `pfc` with features
#' that are just the index of the edge. Useful for complex indexing procedures 
#'
#' @param x A `pfc` object
#'
#' @return A new `pfc` with the same structure as `x`, but with
#' the features of each edge equal to the edge index (the column number in the 
#' sparse matrix representation)
#' @export
#'
#' @examples
#' pf_indexes(rpfc(100))
pf_indexes <- function(x) {
  field(x, "pfl") <- field(x, "pfp")
  refresh_features(x)
}

#' Replace feature with edge position in flow in a `pfc`
#' 
#' @param x A `pfc` object
#'
#' @return A new `pfc` with the same structure as `x`, but with
#' the features of each edge equal to the position along the flow as an integer
#' index
#' @export
#'
#' @examples
#' pf_position(rpfc(100))
pf_position <- function(x) {
  field(x, "pfl") <- purrr::map(field(x, "pfp"),
                                seq_along)
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
               internal_edges(x))

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
#' pf_scale_flow_sum(rpfc(100))
pf_scale_flow_sum <- function(x, scale_to = 1) {
  
  field(x, "pfl") <- purrr::map(field(x, "pfl"),
                                function(y) (y / sum(y)) * scale_to)
  refresh_features(x)
  
}

#' Standardise the phylogenetic flow features to an implied typical
#' variance of 1. 
#' 
#' Can be used to make random effects based on `pfc`s comparable.
#'
#' @param x a `pfc` object
#'
#' @return a scaled `pfc` object
#' @export
#'
#' @examples
#' pf_standardise(rpfc(100))
pf_standardise <- function(x) {
  
  fact <- mean(pf_flow_sum(x))
  
  field(x, "pfl") <- purrr::map(field(x, "pfl"),
                                function(y) y / fact)
  refresh_features(x)
  
}

#' Return the number of edges in a `pfc`
#'
#' @param x a `pfc` object
#'
#' @return The number of edges
#' @export
#'
#' @examples
#' pf_nedges(rpfc(100))
pf_nedges <- function(x) {
  length(edge_names(x))
}

#' @export
pf_as_pfc.factor <- function(x, ...) {
  
  levs <- levels(x) 
  m <- Matrix::.sparseDiagonal(length(levs), shape = "g")
  colnames(m) <- rownames(m) <- levs
  m_expand <- m[as.character(x), ] 
  pf_as_pfc(m_expand, is_tip = rep(1, length(x)))
  
}

#' @export
pf_as_pfc.character <- function(x, ...) {
  
  x <- as.factor(x)
  pf_as_pfc(x)
  
}

#' Replace features with ones
#'
#' @param x A `pfc` object
#'
#' @return A `pfc` object with all feature values replaced with ones
#' @export
#'
#' @examples
#' pf_ones(rpfc(100))
pf_ones <- function(x) {
  m <- structural_sparse(pf_as_sparse(x))
  pf_as_pfc(m, field(x, "is_tip"))
}

#' Replace features with zeros
#'
#' @param x A `pfc` object
#'
#' @return A `pfc` object with all feature values replaced with zeros
#' @export
#'
#' @examples
#' pf_zeros(rpfc(100))
pf_zeros <- function(x) {
  m <- pf_as_sparse(x)
  m@x[] <- 0
  pf_as_pfc(m, field(x, "is_tip"))
}

#' Apply a function along edges within a `pfc` object
#'
#' Features of the returned `pfc` will be the output of the function applied
#' to the original features. Therefore the functionm provided must return a
#' vector the same length of the input.
#'
#' @param x A `pfc` object
#' @param fun A function to apply to the edge vectors, which will be passed as 
#' the first argument. Must return a vector the same length and type as the input.
#' @param ... Any additional arguments to pass to `fun`
#'
#' @return A `pfc` object with the same structure as `x` but with features 
#' replaced with the output of `fun`
#' @export
#'
#' @examples
#' ## average features of descendents of each edge as a pfc:
#' pf_edge_apply(pf_desc(pf_indexes(rpfc(100))), function(x) rep(mean(x), length(x)))
pf_edge_apply <- function(x, fun, ...) {
  x_sp <- pf_as_sparse(x)
  xs <- split_xs(x_sp)
  xs <- purrr::map(xs,
                   fun, 
                   ...)
  xs <- unlist(xs)
  if(length(xs) != length(x_sp@x)) {
    rlang::abort("Length of output does not match length of input!")
  }
  x_sp@x <- xs
  m <- pf_as_pfc(x_sp, is_tip = field(x, "is_tip"),
                 internal = attr(x, "internal"))
  m
}

pf_apply_sparse_fn <- function(.x, .fn,
                               is_tip = field(.x, "is_tip"),
                               internal = attr(.x, "internal"),                            
                               ...) {
  m <- pf_as_sparse(.x)
  m2 <- force_dgCMatrix(.fn(m, ...))
  if(is.null(colnames(m2))) {
    colnames(m2) <- colnames(m)[seq_len(ncol(m2))]
  }
  if(is.null(rownames(m2))) {
    rownames(m2) <- rownames(m)[seq_len(nrow(m2))]
  }
  pf_as_pfc(m2, is_tip = is_tip, internal = internal)
}

#' Extract edge names from `pfc` object 
#'
#' @param x `pfc` object
#'
#' @return A character vector of edge names
#' @export
#'
#' @examples
#' pf_edge_names(rpfc(100))
pf_edge_names <- function(x) {
  edge_names(x)
}


#' Extract mean edge features from `pfc` object 
#'
#' @param x `pfc` object
#'
#' @return A numeric vector of mean edge features
#' @export
#'
#' @examples
#' pf_mean_edge_features(rpfc(100))
pf_mean_edge_features <- function(x) {
  edge_lengths(x, mean, empty = 0)
}

#' @export
`[.pfc` <- function(x, i, j = NULL) {
  if(is.character(i)) {
    i <- vec_match(i, field(x, "pfn"))
    return(NextMethod("["))
  } else {
    if(is.null(j)) {
      return(NextMethod("["))
    }
  }
}

#' Extract `pfp` object from `pfc`, the paths of
#' each flow from root to terminal node.
#'
#' @param x A `pfc` object
#'
#' @return A `pfp` object
#' @export
#'
#' @examples
#' pf_path(rpfc(100))
pf_path <- function(x) {
  field(x, "pfp")
}

#' Return a logical vector which is `TRUE` for the elements
#' of a `pfc` whci reprsent tips of a phylogeny
#' @param x A `pfc` object.
#' @param ... Other arguments passed to or from other methods
#'
#' @export
#' @examples
#' pf_is_tips(rpfc(100)) 
pf_is_tips <- function(x, ...) {
  
  field(x, "is_tip")
  
}

drop_zero_edges <- function(x, ...) {
  m <- pf_as_sparse(x)
  internal <- attr(x, "internal")
  zero_br <- Matrix::colSums(m) != 0
  m <- m[ , zero_br]
  pf_as_pfc(m, internal = internal[zero_br],
            is_tip = field(x, "is_tip"))
}

drop_root_edges <- function(x, ...) {
  mrca <- pf_mrca(x)
  m <- pf_as_sparse(x)
  internal <- attr(x, "internal")
  m <- m[ , -1:-mrca]
  pf_as_pfc(m, internal = internal[-1:-mrca],
            is_tip = field(x, "is_tip"))
}
  