discrete_colours <- function(n = 1) {
  ghibli::ghibli_palettes[c("TotoroMedium")]
}

totoro_palette <- function() {
  #col2rgb(ghibli::ghibli_palettes$TotoroMedium[c(7, 2, 1, 4, 5)])
}

continuous_colors <- function(x) {
  # colourvalues::color_values(x, totoro_palette(), 
  #                            na_colour = "grey")
}

#ghibli::ghibli_palette("TotoroMedium", n = 20, type = "continuous")