new_pfc <- function(phy) {
  assertthat::assert_that(inherits(phy, "phylo"))

  phy <- add_internal_tips(phy)


  new_rcrd()
}
