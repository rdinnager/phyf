flow_to_edge_list <- function(x) {
  x <- c(0, x) + 1
  cbind(x[-length(x)], x[-1])
}