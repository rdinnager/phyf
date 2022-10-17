#' AVONET Bird Trait Data with Phylogeny
#'
#' The AVONET Bird Trait Database joined to a `pf` object (for `{phyf}`)
#'
#' @format ## `avonet`
#' A 'pf' data frame (subclasses `tibble`) with 13,338 rows and 38 columns:
#' \describe{
#'   \item{label}{Node labels including species name for the tip labels}
#'   \item{phlo}{The phylogenetic flow column which stores the phylogenetic information}
#'   \item{Species3, Family3, Order3}{Taxonomic names -- Names of the species, family and order, respectively}
#'   \item{Total.individuals}{Number of individuals used to measure the data}
#'   \item{Beak.Length_Culmen:Species.Status}{Various traits of the bird species, see Source section to get more detailed information}
#' }
#' @source <https://figshare.com/articles/dataset/AVONET_morphological_ecological_and_geographical_data_for_all_birds_Tobias_et_al_2021_Ecology_Letters_/16586228>
#' @references Tobias JA, Sheard C, Pigot AL, Devenish AJ, Yang J, Sayol F, Neate-Clegg MH, Alioravainen N, Weeks TL, Barber RA, Walkden PA. AVONET: morphological, ecological and geographical data for all birds. Ecology Letters. 2022 Mar 1;25(3):581-97.
"avonet"

#' Vertebrate Base Metabolic Rates with Phylogeny
#'
#' Data on vertebrate base Metabolic rates joined to a `pf` object (for `{phyf}`)
#'
#' @format ## `vert_bmr`
#' A 'pf' data frame (subclasses `tibble`) with 1,712 rows and 8 columns:
#' \describe{
#'   \item{label}{Node labels including species name for the tip labels}
#'   \item{phlo}{The phylogenetic flow column which stores the phylogenetic information}
#'   \item{lnBMR}{Natural log of the base metabolic rate}
#'   \item{lnMass}{Natural log of body mass}
#'   \item{lnMass2}{Squared natural log of body mass}
#'   \item{lnGS}{Natural log of genome size}
#'   \item{endo}{Is the species endothermic? 1 for yes, 0 for no}
#' }
#' @source <https://datadryad.org/stash/dataset/doi:10.5061/dryad.3c6d2>
#' @references Uyeda JC, Pennell MW, Miller ET, Maia R, McClain CR. The evolution of energetic scaling across the vertebrate tree of life. The American Naturalist. 2017 Aug 1;190(2):185-99.
"vert_bmr"