#' Calculate a 'second order' `pfc`, which represents summed
#' branch length from each node to the terminal node along each flow.
#' When used in modelling traits, a second order `pfc` implies that
#' rates of evolution of the trait are themselves evolving according to
#' Brownian motion.
#'
#' @param x A `pfc` object
#' @param ... Additional arguments passed to or other methods.
#'
#' @return A `pfc` object
#' @export
#'
#' @examples
#' pf_second_order(rpfc(100))
pf_second_order <- function(x, ...) {
  x2 <- pf_flow_cumsum(x, "to_root")
  x2
}