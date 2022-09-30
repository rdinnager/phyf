add_internal_tips <- function(phy) {
  node_labels <- length(phy$tip.label) + seq_len(phy$Nnode)
  new_phy <- phangorn::add.tips(phy, as.character(node_labels),
                                node_labels, 0)
}
