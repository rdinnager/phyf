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

#' Primate Diet Diversity and Threat Status Data with Phylogeny
#'
#' Data on primates diets and their threat status joined to a `pf` object (for `{phyf}`)
#'
#' @format ## `primate_diet`
#' A 'pf' data frame (subclasses `tibble`) with 504 rows and 55 columns:
#' \describe{
#'   \item{label}{Node labels including species name for the tip labels}
#'   \item{phlo}{The phylogenetic flow column which stores the phylogenetic information}
#'   \item{Threat status}{IUCN threat category}
#'   \item{Threat status source}{Threat status data source -- see `primate_refs` dataset}
#'   \item{Body mass (g)}{Primate species' mean body mass in grams}
#'   \item{Body mass source)}{Body mass data source -- see `primate_refs` dataset}
#'   \item{Range size (km2)}{Primate species' range size in kilometers}
#'   \item{Range size source)}{Range size data source -- see `primate_refs` dataset}
#'   \item{Diet disparity (PSV)}{Primate species' disparity of diet items as 
#'   measured by Phylogenetic Species Variability (PSV) metric.}
#'   \item{Diet breadth}{Primate species' breadth of diet items as measured by 
#'   the number of different diet items (richness).}
#'   \item{Diet diversity (DDI)}{Primate species' diversity of diet items as measured by 
#'   the DDI metric.}
#'   \item{Trophic guild}{Primate species' trophic guild. Possible values: 
#'   "Omnivore", "Frugivore", "Gummivore", "Insectivore", or "Folivore-frugivore"}
#'   \item{diet_item:item (40 columns)}{Each of the 40 columns starting with "diet_item:"
#'   represents a different type of item in primates' diets. These are binary integer
#'   columns with a 1 if the species feeds on that diet item or 0 if it does not.}
#' }
#' @source <https://zslpublications.onlinelibrary.wiley.com/doi/full/10.1111/acv.12823>
#' @references Machado, F. F., Jardim, L., Dinnage, R., Brito, D., & Cardillo, M. (2022). Diet disparity and diversity predict extinction risk in primates. Animal Conservation.
"primate_diet"

#' Primate Diet Items Hierarchy
#' 
#' A `data.frame` representing an hierarchical categorization of primate diet items
#'
#' @format ## `primate_hierarchy`
#' A data frame with 40 rows and 6 columns:
#' \describe{
#'   \item{FOOD ITEM}{Food item name. This matches the 'item' 'diet_item:item' 
#'   columns in  the `primate_diet` data set}
#'   \item{level 2}{A category categorizing the food items immediately above
#'   the items themselves}
#'   \item{level 3:level 6 (4 columns)}{Additional categories arranged in an 
#'   heirarchy}
#' }
#' @source <https://zslpublications.onlinelibrary.wiley.com/doi/full/10.1111/acv.12823>
#' @references Machado, F. F., Jardim, L., Dinnage, R., Brito, D., & Cardillo, M. (2022). Diet disparity and diversity predict extinction risk in primates. Animal Conservation.
"primate_diet_hierarchy"

#' Primate Diet Data References
#' 
#' A `data.frame` with references for the Primate Diet Data. See the `primate_diet`
#' data set.
#'
#' @format ## `primate_diet_refs`
#' A 'pf' data frame (subclasses `tibble`) with 4 rows and 3 columns:
#' \describe{
#'   \item{FOOD ITEM}{Food item name. This matches the 'item' 'diet_item:item' 
#'   columns in  the `primate_diet` data set}
#'   \item{level 2}{A category categorizing the food items immediately above
#'   the items themselves}
#'   \item{level 3:level 6 (4 columns)}{Additional categories arranged in an 
#'   heirarchy}
#' }
#' @source <https://zslpublications.onlinelibrary.wiley.com/doi/full/10.1111/acv.12823>
#' @references Machado, F. F., Jardim, L., Dinnage, R., Brito, D., & Cardillo, M. (2022). Diet disparity and diversity predict extinction risk in primates. Animal Conservation.
"primate_diet_hierarchy"


#' Terrestrial Mammal Bioegeography Data with Phylogeny
#'
#' Data on terrestrial mammals biogeographic distributions across the world's ecoregions 
#' joined to a `pf` object (for `{phyf}`)
#'
#' @format ## `mammal_biogeo`
#' A 'pf' data frame (subclasses `tibble`) with 10,648 rows and 844 columns:
#' \describe{
#'   \item{label}{Node labels including species name for the tip labels}
#'   \item{phlo}{The phylogenetic flow column which stores the phylogenetic information}
#'   \item{IUCN_binomial}{Species name used by the IUCN, matches with IUCN range polygons}
#'   \item{body_mass_median}{Mammal species' body mass}
#'   \item{litter_clutch_size}{Mammal species' average clutch size 
#'   (# of offspring in a litter)}
#'   \item{activity)}{Mammal species' primary time of activity}
#'   \item{hab_breadth}{...}
#'   \item{volant}{...}
#'   \item{diet_5cat}{...}
#'   \item{range_size_km2}{Mammal species' range size in kilometers squared}
#'   \item{threat}{Mammal species' IUCN threat category}
#'   \item{ecoregion:eco_id}{Each of the 832 columns starting with "ecoregion:"
#'   represents the proportion of the mammal species' range that fall in the ecoregion 
#'   with id equal to "eco_id".}
#' }
#' @source <https://ecoregions.appspot.com>, <https://www.iucnredlist.org/resources/spatial-data-download>
#' @references None yet.
"mammal_biogeo"

#' Ecoregion IDs
#' 
#' A `data.frame` with ecoregion ids and their name. This can be matched to the
#' ecoregions referred to in the dataset `mammal_biogeo`. 
#'
#' @format ## `ecoregion_ids`
#' A 'pf' data frame (subclasses `tibble`) with 847 rows and 2 columns:
#' \describe{
#'   \item{ECO_ID}{ID numbers for 847 ecoregions across the world}
#'   \item{ECO_NAME}{The name for the corresponding ecoregion ID}
#' }
#' @source <https://ecoregions.appspot.com>
"ecoregion_ids"