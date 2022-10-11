add_internal_tips <- function(phy) {
  nodes <- ape::Ntip(phy) + seq_len(ape::Nnode(phy))
  if(!is.null(phy$node.label)) {
    node_labels <- phy$node.label
  } else {
    node_labels <- as.character(node_labels)
  }
  new_phy <- phangorn::add.tips(phy, node_labels,
                                nodes, 0)
  new_phy
}

root2tip_binary <- function(phy) {

  np <- ape::nodepath(phy)
  rtp <- build_rtp(np, nrow(phy$edge) + 1)

  rtp

}

build_rtp <- function(paths, n_nodes,
                      lens = NULL,
                      sparse = TRUE) {
  
  if(is.null(lens)) {
    lens <- 1
  } else {
    lens <- unlist(lens)
  }

  if(!sparse) {

    rtp_mat <- lapply(paths, function(x) tabulate(x, n_nodes))
    do.call(rbind, rtp_mat)

  } else {

    js <- rep(seq_along(paths),
              lengths(paths))

    nj <- length(paths)

    ig_out <- unlist(paths)

    nn <- seq_len(n_nodes)

    Matrix::sparseMatrix(
      j = js,
      i = fastmatch::fmatch(ig_out, nn, nomatch = 0),
      x = lens,
      dims = c(n_nodes, nj),
    )
  }

}

split_indexes <- function(x) {
  res <- split(x@i + 1, findInterval(seq_len(Matrix::nnzero(x)), x@p, left.open = TRUE))
  res
}
