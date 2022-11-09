#' Return the edge leading up to the most recent common ancestor of a set
#' of phylogenetic flows
#'
#' @param x A `pfc` object
#' @param name Should the edge name be returned? If `FALSE`, the default,
#' the edge number is returned
#' @param ... Other arguments passed to or from other methods
#'
#' @return A string with the edge name is `name = TRUE`, or an integer 
#' index otherwise
#' @export
#'
#' @examples
#' require(dplyr)
#' avonet %>%
#'   filter(label %in% c("Platalea_minor", "Pelecanus_occidentalis")) %>%
#'   pull(phlo) %>%
#'   pf_mrca()
pf_mrca <- function(x, name = FALSE, ...) {
  if(length(x) < 2) {
    rlang::abort("pfc object x must have a length of at least 2 for pf_mrca")
  }
  if(length(x) == 2) {
    shared_edges <- do.call(intersect, field(x, "pfp")) 
  } else {
    shared_edges <- Reduce(intersect, field(x, "pfp"))
  }
  max(shared_edges)
}

#' Return a logical vector determining if a flow's terminal node is
#' a descendant of an edge.
#'
#' @param x A `pfc` object
#' @param edge An edge number or name
#'
#' @return A logical vector
#' @export
#'
#' @examples
#' pf_is_desc(rpfc(100), "Node12")
pf_is_desc <- function(x, edge) {
  
  if(is.character(edge)) {
    nams <- edge_names(x)
    descs <- purrr::map_lgl(field(x, "pfp"),
                          ~ edge %in% nams[.x])
  } else {
    descs <- purrr::map_lgl(field(x, "pfp"),
                            ~ edge %in% .x)
  }
  
  descs
  
}