#' Filter descendents of a MRCA
#' 
#' Filter a `pf` or `pfc` object by taking all
#' descendents of the most recent common ancestor (mrca)
#' of `descendents`
#'
#' @param x A `pf` or `pfc` object
#' @param descendents An index into a `pfc`, using any way
#' that you can index a `pfc`. `filter_with_mrca` will find
#' the most recent common ancestor of everything in
#' `descendents` and then it will return a filtered `pf`
#' or `pfc` with all descendents of the mrca. If of length one,
#' it will instead be assumed to be the node name or number representing
#' the mrca itself.
#' @param drop_null_edges If `TRUE` (the default), drop edges with
#' no descendants after the filtering.
#' @param drop_root_edges If `TRUE` (the default), drop all edges
#' leading up to the mrca (e.g. when combined together these would make
#' a root edge leading to the desired clade)
#' @param ... Other arguments for future extensions.
#'
#' @return A `pf` or `pfc` object
#' @export
#'
#' @examples
#' avonet %>% 
#'     pf_filter_with_mrca(label %in% c("Platalea_minor", "Pelecanus_occidentalis"))
pf_filter_with_mrca <- function(x, descendents, drop_null_edges = TRUE, 
                                drop_root_edges = TRUE, ...) {
  UseMethod("pf_filter_with_mrca")
}

#' @export
pf_filter_with_mrca.pf <- function(x, descendents, drop_null_edges = TRUE,
                                   drop_root_edges = TRUE, ...) {
  pf_col <- attr(x, "pf_column")[1]
  desc <- rlang::eval_tidy(rlang::enquo(descendents), x)
  
  if(length(desc) == 1) {
    logi <- pf_is_desc(desc)
  } else {
    if(pf_is_pfc(desc)) {
      mrca <- pf_mrca(desc)
      
    } else {
      mrca <- pf_mrca(x[[pf_col]][which(desc)])
    } 
    
    logi <- pf_is_desc(x[[pf_col]], 
                         mrca)
  }
  
  
  
  res <- x %>%
    dplyr::filter(logi)
  
  if(drop_null_edges) {
    res[[pf_col]] <- drop_zero_edges(res[[pf_col]])  
  }
  
  if(drop_root_edges) {
    res[[pf_col]] <- drop_root_edges(res[[pf_col]])
  }
  
  return(res)

}

#' @export
pf_filter_with_mrca.pfc <- function(x, descendents, drop_null_edges = TRUE,
                                    drop_root_edges = TRUE, ...) {
  
  if(length(descendents) == 1) {
    mrca <- descendents
  } else {
    mrca <- pf_mrca(x[descendents])
  }
  res <- x[pf_is_desc(x, mrca)]
  if(drop_null_edges) {
    res <- drop_zero_edges(res)  
  }
  if(drop_root_edges) {
    res <- drop_root_edges(res)
  }
  res
}