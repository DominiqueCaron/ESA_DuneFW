library(dplyr)
library(rgbif)
library(purrr)
library(missForest)
library(brms)
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# prepare traits ----------------------------------------------------------
tetrapod_traits <- read.csv("data/raw/TetrapodTraits.csv", row.names = 1)
tetrapod_traits$Order <- firstup(tetrapod_traits$Order)
tetrapod_traits$Family <- firstup(tetrapod_traits$Family)

Dune_sp <- c(
  "Vulpes macrotis",
  "Microdipodops megacephalus",
  "Lepus californicus",
  "Terrapene ornata",
  "Buteo jamaicensis",
  "Glaucidium brasilianum",
  "Aquila chrysaetos",
  "Strix hadorami")
 
Dune_taxo <- 
  map_df(Dune_sp, name_backbone)

Dune_traits <- filter(tetrapod_traits, Species %in% Dune_sp) %>%
  rbind(filter(tetrapod_traits, Family %in% Dune_taxo$family))

traits <- Dune_traits %>%
  select_at(vars(Order:Family, Body_mass_g:Other.Unknown, Adult_svl_cm:Generation_length_d)) %>%
  mutate_at(vars(Order:Family, Trophic_level:Diel_activity, Artificial_habitat_use:Other.Unknown), as.factor) %>%
  droplevels() %>%
  missForest()

Dune_traits <- cbind(
  select(Dune_taxo, Species = species, Order = order, Family = family, Genus = genus),
  select(traits$ximp, -Order, -Family, -Adult_svl_cm, -Longevity_d, -Generation_length_d, -Maturity_d, -Artificial_habitat_use, -Other.Unknown)
)[c(1:length(Dune_sp)),] %>%
  mutate(Habitat_breadth_IUCN = log(Habitat_breadth_IUCN),
         logBM = log(Body_mass_g),
         logLongevity = log(Max_longevity_d),
         logClutchSize = log(Litter_clutch_size),
         Herbivore = ifelse(Trophic_level == "Herbivore", 1,0),
         Omnivore = ifelse(Trophic_level == "Omnivore", 1,0),
         Carnivore = ifelse(Trophic_level == "Carnivore", 1, 0)) %>%
  select(-Body_mass_g, -Trophic_level, -Max_longevity_d, -Litter_clutch_size)

write.csv(Dune_traits, "data/clean/DuneTraits.csv")

# Get predictors ----------------------------------------------------------
source("functions.R")
Dune_traits <- read.csv("data/clean/DuneTraits.csv", row.names = 1)
DuneFW <- get_predictors(Dune_traits$Species, Dune_traits)

var2scale <- c("Habitat_breadth.predator", "BM.predator", "Longevity.predator", 
               "ClutchSize.predator", "Habitat_breadth.prey", "BM.prey",
               "Longevity.prey", "ClutchSize.prey", "Habitat.match", "BM.match")

# Europe stats and models
EuropeModel <- readRDS("~/OneDrive/Chapt-NicheSpace/models/EuroModel_brms.rds")
Europe_stats <- read.csv("data/clean/europe_stats.csv", row.names = 1)

# scale predictors
DuneFWscaled <- DuneFW
DuneFWscaled[,var2scale] <- sweep(DuneFWscaled[,var2scale], MARGIN = 2, Europe_stats$mean)
DuneFWscaled[,var2scale] <- sweep(DuneFWscaled[,var2scale], MARGIN = 2, FUN = "/", 2*Europe_stats$sd)


# Make predictions --------------------------------------------------------
DuneFWpredictions <- make_predictions(EuropeModel, newdata = DuneFWscaled, 
                                     allow_new_levels = TRUE, ndraws = 100, extrapolation = T)

write.csv(DuneFWpredictions, "data/output/DuneFWpredictions.csv")
