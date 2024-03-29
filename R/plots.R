#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' Make an automatic `ggplot2` plot for a `pf` object
#'
#' @param object A `pf` object to plot
#' @param columns Columns to plot along with the phylogeny.
#' Can use bare column names or any other `tidyselect` syntax 
#' @param layout `ggtree::ggtree()` layout to use.
#' @param suppress_tiplabels If `TRUE`, don't draw tip labels.
#' @param edge_traits A logical indicating whether the continuous avriable
#' refers to edge traits, where the edge is determines by the terminal node.
#' (default is FALSE, which means the variable refers to a node trait).
#' @param continuous continuous transition for selected aesthethic 
#' ('size' or 'color'('colour')). It should be one of 'color' (or 'colour'), 
#' 'size', 'all' and 'none', default is 'colour'.
#' @param ... Other arguments passed to or from other methods.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' autoplot(pf(rpfc(100)) %>% dplyr::mutate(trait = rnorm(dplyr::n())), trait,
#' layout = "rectangular")
autoplot.pf <- function(object, columns = NULL, layout = "circular",
                        suppress_tiplabels = FALSE,
                        suppress_tippoints = FALSE,
                        edge_traits = FALSE,
                        continuous = "colour",
                        size = 1.4,
                        outline_size = 1.4 * size,
                        ...) {
  
  sel <- dplyr::select(object, {{ columns }})
  
  ob <- phyf::pf_phyloflow(object)
  tree <- phyf::pf_as_phylo(ob)
  
  tree_df <- ggtree::fortify(tree)
  
  if(length(sel) == 0) {
    p1 <- ggtree::ggtree(tree) 
    if(length(tree_df$label[tree_df$isTip]) < 200 && !suppress_tiplabels) {
      p1 <- p1 + ggtree::geom_tiplab()
    }
    return(p1)
  }
  
  tree_df <- tree_df %>%
      dplyr::left_join(sel %>%
                         dplyr::mutate(label = phyf::pf_labels(ob)),
                       by = "label")
  
  if(ncol(sel) == 1) {
    
    # && (nrow(tree_df) == length(ob) + 1 || nrow(tree_df) == length(ob))
    
    if(all(is.na(unlist(sel)[!pf_is_tips(ob)]))) {
      p1 <- ggtree::ggtree(tree_df, layout = layout, ladderize = FALSE, size = size,
                           ...) + 
        ggplot2::scale_color_viridis_c() +
        ggplot2::theme(legend.position = c(.05, .85))
      if(!suppress_tippoints) {
        p1 <- p1 +
          ggtree::geom_tippoint(colour = "black",
                                size = 1.4 * outline_size) +
          ggtree::geom_tippoint(ggplot2::aes(color = {{ columns }}),
                                size = 1.4 * size)
      }
      return(p1)
    }
    
    if(edge_traits) {
      continuous <- "none"
    }
    
    p1 <- ggtree::ggtree(tree_df, layout = layout, ladderize = FALSE, size = outline_size) + 
      ggtree::geom_tree(ggplot2::aes(color = {{ columns }}), 
                        continuous = continuous, size = size, ...) +  
      ggplot2::scale_color_viridis_c() +
      ggplot2::theme(legend.position = c(.05, .85)) 
    
    if((length(tree_df$label[tree_df$isTip]) > 200 || suppress_tiplabels)) {
      if(!suppress_tippoints) {
        p1 <- p1 + ggtree::geom_tippoint(colour = "black",
                                         size = 1.4 * outline_size) +
          ggtree::geom_tippoint(ggplot2::aes(color = {{ columns }}),
                                size = 1.4 * size)
      }
    } else {
      p1 <- p1 + ggtree::geom_tiplab(ggplot2::aes(color = {{ columns }}), hjust = -.1)
    }
    
  }
  
  p1
  
}


#' Make a plot for a `pf` object
#'
#' @param x A `pf` obect to plot
#' @param columns Bare column names of variables to plot with tree.
#' @param layout Either 'phylogram' or 'fan'
#' @param ... Other arguments passed to `phytools::contMap()`
#'
#' @return A `phytools::contMap()` object.
#' @export 
plot.pf <- function(x, columns = NULL, 
                    layout = "fan", 
                    suppress_tiplabels = FALSE,
                    ...) {
  
  object <- x
  
  sel <- dplyr::select(object, {{ columns }})
  
  ob <- phyf::pf_phyloflow(object)
  tree <- phyf::pf_as_phylo(ob)
  
  if(length(sel) == 0) {
    if(length(tree$tip.label) > 200) {
      p1 <- plot(tree, show.tip.label = FALSE, ...)
    } else {
      p1 <- plot(tree, ...)
    }
    return(invisible(p1))
  }
  
  total_nodes <- tree$Nnode + length(tree$tip.label)
  
  if(ncol(sel) == 1) {
    
    # && total_nodes == length(ob) + 1
    
    if(length(tree$tip.label) > 200 || suppress_tiplabels){
      ftype <- c("off", "reg")  
    } else {
      ftype <- NULL
    }
    
    x <- sel %>%
      dplyr::pull({{ columns }})
    
    names(x) <- pf_labels(ob)
    
    tips <- x[field(ob, "is_tip")]
    ancs <- x[!field(ob, "is_tip")]
    
    if(all(is.na(ancs))) {
      method = "fastAnc"
    } else {
      method = "user"
      ancs <- ancs[!is.na(ancs)]
    }
    
    pal <- scales::viridis_pal()
    
    p1 <- phytools::contMap(tree, tips, 
                            anc.states = ancs,
                            method = method,
                            plot = FALSE,
                            ...)
    p1 <- phytools::setMap(p1, colors = pal(50))
    plot(p1, ftype = ftype, fsize = 1.2,
         type = layout,
         sig = 1,
         leg.txt = rlang::as_label(rlang::enquo(columns)))
    
  }
  
  
}

#' @importFrom ggplot2 fortify
#' @export
ggplot2::fortify

#' @export
fortify.pfc <- function(model, data, ...) {
  
  tree <- pf_as_phylo(model)
  ggtree::fortify(tree, ...)
  
}

#' @export
fortify.pf <- function(model, data, ...) {
  
  tree <- pf_as_phylo(model)
  pf_col <- attr(model, "pf_column")[1]
  dat <- model[ , -which(colnames(model) == pf_col)]
  if(!"label" %in% colnames(model)) {
     dat <- dat %>%
      dplyr::mutate(label = pf_labels(model[[pf_col]]))
  }
  tree_df <- ggtree::fortify(tree, ...) %>%
    dplyr::left_join(dat, by = "label")
  tree_df
  
}

