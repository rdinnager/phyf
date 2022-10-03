add_internal_tips <- function(phy) {
  node_labels <- seq_len(ape::Ntip(phy) + ape::Nnode(phy))
  new_phy <- phangorn::add.tips(phy, as.character(node_labels),
                                node_labels, 0)
  new_phy
}
