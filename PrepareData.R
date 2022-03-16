library(dplyr)
library(rgbif)
library(purrr)
library(missForest)
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
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
  "Strix hadorami",
  "Antrozous pallidus")
 
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
#####
Dune_traits <- read.csv("data/clean/DuneTraits.csv", row.names = 1)
DuneFW <- expand.grid(Dune_traits$Species, Dune_traits$Species) 
colnames(DuneFW) <- c("Predator", "Prey")

DuneFW <- left_join(DuneFW, Dune_traits, by = c("Prey" = "Species")) %>%
  left_join(Dune_traits, by = c("Predator" = "Species"))

DuneFW <- traits2predictors(DuneFW)
DuneFW[, columns_to_scale] <- sweep(DuneFW[, columns_to_scale], MARGIN = 2, Predictor_means)
DuneFW[, columns_to_scale] <- sweep(DuneFW[, columns_to_scale], MARGIN = 2, 2*Predictor_sd, FUN = "/")

predictors <- select(DuneFW, 
                     -Predator, -Prey, -Order.predator, -Order.prey,
                     -Herbivore.predator, -Herbivore.prey)

predictors <- cbind(rep(1, nrow(DuneFW)), predictors) %>% as_data()
predator_order <- as.numeric(factor(DuneFW$Order.pred, levels = unique(FuncTraits$Order)))
linear_predictor <- rowSums(predictors * t(coef[,predator_order]))

p <- ilogit(linear_predictor)

predictions <- calculate(p, values = GLMM, nsim = 100)
predictions <- colMeans(predictions[[1]])
predictions <- data.frame(Predator = as.factor(DuneFW$Predator), Prey = as.factor(DuneFW$Prey), predictions = predictions)

ggplot(predictions) +
  geom_tile(aes(Predator, Prey, fill = predictions)) +
  scale_fill_distiller(palette = "YlGn") +
  theme(axis.text.x = element_text(angle = 90))

