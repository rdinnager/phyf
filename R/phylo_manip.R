#' Segments the edges of a phylogeny by splitting them at particular positions
#'
#' @param x A `pfc` object
#' @param edges A vector of edge names
#' @param positions A vector of position to cut the edges (must be same length as `edges`)
#'
#' @return A `pfc` object with edges segmented
#' @export
#'
#' @examples
#' pf_edge_segmentize(rpfc(100), "t1", 0.01)
pf_edge_segmentize  <- function(x, edges, positions) {
  enames <- pf_edge_names(x)
  edge_nums <- fastmatch::fmatch(edges, enames) 
  split_df <- dplyr::tibble(edge = edge_nums, position = positions) %>%
    dplyr::group_by(edge) %>%
    dplyr::summarise(pos = list(unname(position)))
  pft <- pf_as_pftibble(pf_tips(x)) %>%
    dplyr::left_join(split_df) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pos = list(c(pos, val))) %>%
    dplyr::mutate(pos = list(c(pos[1], diff(pos))))
  new_pft <- pft %>%
    dplyr::select(-val) %>%
    tidyr::unnest(pos) %>%
    dplyr::mutate(edge_name = enames[edge]) %>%
    dplyr::group_by(path, edge) %>%
    dplyr::mutate(edge_name = paste0(edge_name, paste0("_", 1:dplyr::n()))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(edge_fact = reorder(as.factor(edge_name), edge)) %>%
    dplyr::mutate(edge_num = as.numeric(edge_fact))
  
  new_pf <- new_pft %>%
    dplyr::select(path, edge_num, pos) %>%
    dplyr::group_by(path) %>%
    dplyr::group_nest()
  new_pf_data <- purrr::list_transpose(new_pf$data)
  new_enames <- new_pft %>%
    dplyr::group_by(edge_num) %>%
    dplyr::summarise(edge_name = edge_name[1]) %>%
    dplyr::arrange(edge_num)
  
  last <- purrr::map_int(new_pf_data$edge_num,
                         ~ tail(.x, 1))
  new_tips <- new_enames$edge_name[last]
  internal <- !new_enames$edge_num %in% unique(last)
  
  pfc(pfn = new_tips,
      pfpp = pfp(new_pf_data$edge_num),
      pfl = new_pf_data$pos,
      is_tip = rep(TRUE, length(new_pf$path)),
      edge_names = new_enames$edge_name,
      internal = internal)
}

#' Get edges and positions along edges where a set of epoch times intersect
#'
#' @param x An object of class pfc
#' @param times A vector of epoch times to slice the pfc along
#'
#' @return A tibble with edge labels in the first column and position in the second column
#' @export
#'
#' @examples
#' pf_epoch_info(rpfc(100), c(1, 2, 3))
pf_epoch_info <- function(x, times) {
  max_len <- pf_max_length(x)
  if(any(times > max_len | times < 0)) {
    rlang::warn("Some times are outside the range of zero to the maximum length of the phyf ({max_len}), these elements will be NULL in the return value")
  }
  cpf <- pf_flow_cumsum(pf_tips(x))
  edge_cs <- pf_mean_edge_features(cpf)
  names(edge_cs) <- pf_edge_names(cpf)
  ends <- pf_ends(x)
  edges <- ivs::iv(edge_cs[ends$start], edge_cs[ends$end])
  edges_bw <- ivs::iv_locate_between(times, edges)
  
  edges_df <- ivs::iv_align(times, edges, locations = edges_bw) %>%
    dplyr::mutate(edge = edges_bw$haystack) %>%
    dplyr::mutate(position = needles - ivs::iv_start(haystack)) %>%
    dplyr::select(edge, position) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(edge = ends$end[edge])
  
  edges_df
  
}

pf_max_length <- function(x) {
  max(pf_flow_sum(x))
}