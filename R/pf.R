#' Filter descendents of a MRCA
#' 
#' Filter a `pf` or `pfc` object by taking all
#' descendents of the most recent common ancestor (mrca)
#' of `descendents`
#'
#' @param x A `pf` or `pfc` object
#' @param descendents An index into a `pfc`, using any way
#' that you can index a `pfc`. `filter_with_mrca` will find
#' the most recent common ancestor of everyhting in
#' `descendents` and then it will return a filtered `pf`
#' or `pfc` with all descendents of the mrca.
#' @param drop_null_edges If `TRUE` (the default), drop edges with
#' no descendants after the filtering.
#' @param ... Other arguments for future
#' extensions.
#'
#' @return A `pf` or `pfc` object
#' @export
#'
#' @examples
#' avonet %>% 
#'     pf_filter_with_mrca(label %in% c("Platalea_minor", "Pelecanus_occidentalis"))
pf_filter_with_mrca <- function(x, descendents, drop_null_edges = TRUE, ...) {
  UseMethod("pf_filter_with_mrca")
}

#' @export
pf_filter_with_mrca.pf <- function(x, descendents, drop_null_edges = TRUE, ...) {
  pf_col <- attr(x, "pf_column")[1]
  desc <- rlang::eval_tidy(rlang::enquo(descendents), x)
  if(pf_is_pfc(desc)) {
    logi <- pf_is_desc(x[[pf_col]], 
                               pf_mrca(desc))
    
  } else {
    logi <- pf_is_desc(x[[pf_col]], 
                       pf_mrca(x[[pf_col]][which(desc)]))
  }
  
  res <- x %>%
    dplyr::filter(logi)
  
  if(drop_null_edges) {
    res[[pf_col]] <- drop_zero_edges(res[[pf_col]])  
  }
  
  return(res)

}

#' @export
pf_filter_with_mrca.pfc <- function(x, descendents, drop_null_edges = TRUE, ...) {
  res <- x[pf_is_desc(x, pf_mrca(x[descendents]))]
  if(drop_null_edges) {
    res <- drop_zero_edges(res)  
  }
  res
}