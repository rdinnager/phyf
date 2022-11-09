#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' Make an automatic `ggplot2` plot for a `pf` object
#'
#' @param object 
#' @param element 
#' @param ... 
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' autoplot(pf(rpfc(100)) %>% dplyr::mutate(trait = rnorm(dplyr::n())), trait,
#' layout = "rectangular")
autoplot.pf <- function(object, columns = dplyr::everything(), layout = "circular",
                        suppress_tiplabels = FALSE, ...) {
  
  sel <- dplyr::select(object, {{ columns }})
  
  ob <- phyf::pf_phyloflow(object)
  tree <- phyf::pf_as_phylo(ob)
  
  tree_df <- ggtree::fortify(tree)
  
  if(ncol(sel) == 1 && (nrow(tree_df) == length(ob) + 1 || nrow(tree_df) == length(ob))) {
    tree_df <- tree_df %>%
      dplyr::left_join(sel %>%
                         dplyr::mutate(label = phyf::pf_labels(ob)),
                       by = "label")
    
    p1 <- ggtree::ggtree(tree_df, layout = layout, ladderize = FALSE, size = 2.2) + 
      ggtree::geom_tree(ggplot2::aes(color = {{ columns }}), continuous = 'colour', size = 1.4) +  
      ggplot2::scale_color_viridis_c() +
      ggplot2::theme(legend.position = c(.05, .85)) 
    
    if(length(tree_df$label[tree_df$isTip]) > 200 || suppress_tiplabels) {
      p1 <- p1 + ggtree::geom_tippoint(colour = "black",
                                       size = 2.8) +
        ggtree::geom_tippoint(ggplot2::aes(color = {{ columns }}),
                              size = 2)
    } else {
      p1 <- p1 + ggtree::geom_tiplab(ggplot2::aes(color = {{ columns }}), hjust = -.1)
    }
    
  }
  
  p1
  
}


#' Make a plot for a `pf` object
#'
#' @param object A `pf` obect to plot
#' @param columns Bare column names of variables to plot with tree.
#' @param layout Either 'phylogram' or 'fan'
#' @param direction plotting direction for type="phylogram".
#' @param ... Other arguments passed to `phytools::contMap()`
#'
#' @return A `phytools::contMap()` object.
#' @export 
#'
#' @examples
#' plot(pf(rpfc(100)) %>% dplyr::mutate(trait = rnorm(dplyr::n())), trait,
#' layout = "rectangular")
plot.pf <- function(object, columns = dplyr::everything(), 
                        layout = "fan", ...) {
  
  sel <- dplyr::select(object, {{ columns }})
  
  ob <- phyf::pf_phyloflow(object)
  tree <- phyf::pf_as_phylo(ob)
  
  total_nodes <- tree$Nnode + length(tree$tip.label)
  
  if(ncol(sel) == 1 && total_nodes == length(ob) + 1) {
    
    if(length(tree$tip.label) > 200){
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
    
    p1 <- phytools::contMap(tree, tips, 
                            anc.states = ancs,
                            method = method,
                            plot = FALSE,
                            ...)
    p1 <- phytools::setMap(p1, colors = viridis::viridis(50))
    plot(p1, ftype = ftype, fsize = 1.2,
         type = layout,
         sig = 1,
         leg.txt = rlang::as_label(rlang::enquo(columns)))
    
  }
  
  p1
  
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
