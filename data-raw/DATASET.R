## create AVONET dataset
library(phyf)
library(dplyr)
library(sf)
library(ape)
library(fasterize)
library(raster)

avonet_tree <- ape::read.nexus("extdata/HackettStage1_0001_1000_MCCTreeTargetHeights.nex")
avonet_dat <- readr::read_csv("extdata/AVONET3_BirdTree.csv")

avonet <- pf_as_pf(avonet_tree)
avonet <- avonet %>%
  left_join(avonet_dat %>%
              mutate(label = gsub(" ", "_", Species3)))

usethis::use_data(avonet, overwrite = TRUE)

#### Uyeda et al data

vert_tree <- ape::read.tree("extdata/vertTree.tre")
vert_dat <- readr::read_csv("extdata/vertData.csv")

vert_bmr <- pf_as_pf(vert_tree)
vert_bmr <- vert_bmr %>%
  left_join(vert_dat %>%
              rename(label = `...1`))

usethis::use_data(vert_bmr, overwrite = TRUE)


#### Mammal Biogeography

mammal_tree <- ape::read.nexus("extdata/terrestrial_mammal_tree_matching_IUCN.nexus")
mammal_maps <- sf::read_sf("extdata/MAMMALS")
ecoregions <- sf::st_read("extdata/Ecoregions")

ecoregions <- ecoregions %>% 
  dplyr::select(ECO_NAME, ECO_ID) %>%
  st_make_valid()

ecoregions_rast <- fasterize(ecoregions, raster::raster(ecoregions, resolution = 0.1),
                             field = "ECO_ID")

mammal_maps <- mammal_maps %>%
  dplyr::select(binomial, presence, origin, seasonal) %>%
  filter(origin %in% c(1, 2)) %>%
  filter(presence %in% c(1, 4, 5)) %>%
  dplyr::select(binomial) %>%
  st_make_valid()

not_valid <- !st_is_valid(mammal_maps)
sum(not_valid)
plot(mammal_maps[not_valid, "geometry"])

mammal_maps <- mammal_maps[!not_valid, ]

mammal_maps <- mammal_maps %>%
  group_by(binomial) %>%
  summarise()

nrow(mammal_maps)

readr::write_rds(mammal_maps, "exdata/mammal_ranges_processed.rds")

ecoregion_ex <- raster::extract(ecoregions_rast, mammal_maps)

mammal_biogeo <- pf_as_pf(mammal_tree)
mammal_biogeo <- mammal_biogeo %>%
  left_join(vert_dat %>%
              rename(label = `...1`))

usethis::use_data(vert_bmr, overwrite = TRUE)


