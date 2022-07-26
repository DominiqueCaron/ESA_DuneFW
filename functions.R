# Functions used throughout the different scripts

#' Jaccard similarity
#' 
#' @param pred_df dataframe. Traits of the predator.
#' @param prey_df dataframe. Traits of the prey.
#' @return a vector indicating the similarity between the traits of the predator and the prey.
Jaccard <- function(pred_df, prey_df){
  jaccard_s <- c()
  a <- rowSums(pred_df == 1 & prey_df == 1)
  b <- rowSums(pred_df == 1 & prey_df == 0)
  c <- rowSums(pred_df == 0 & prey_df == 1)
  jaccard_s <- a/(a+b+c)
  jaccard_s[is.na(jaccard_s)] <- 0.5
  return(jaccard_s)
}

#' Transform traits into predictors
#' 
#' @param FW dataframe. Traits of the the predators and the prey.
#' @param prey_suffix character. Suffix after all column names for the traits of the prey.
#' @param predator_suffix character. Suffix after all column names for the traits of the predator.
#' @return a dataframe with the predictors of species interactions.
traits2predictors <- function(FW, prey_suffix = ".x", predator_suffix = ".y"){
  FW <- FW %>%
    rename_with(~ gsub(paste0("\\", predator_suffix, "$"), ".predator", .x)) %>%
    rename_with(~ gsub(paste0("\\", prey_suffix, "$"), ".prey", .x))
  foraging_predictor <- FW %>%
    select(Herbivore.predator, Omnivore.predator, Carnivore.predator,
           Habitat_breadth.predator = Habitat_breadth_IUCN.predator, 
           Order.predator, BM.predator = logBM.predator, Longevity.predator = logLongevity.predator,
           ClutchSize.predator = logClutchSize.predator)
  
  vulnerability_predictor <- FW %>%
    select(Herbivore.prey, Omnivore.prey, Carnivore.prey,
           Habitat_breadth.prey = Habitat_breadth_IUCN.prey, 
           Order.prey, BM.prey = logBM.prey, Longevity.prey = logLongevity.prey, 
           ClutchSize.prey = logClutchSize.prey)
  
  match_predictor <- data.frame(
    ActivityTime.match = FW$Diel_activity.prey == FW$Diel_activity.predator,
    Habitat.match = Jaccard(
      select_at(FW, vars(Forest.predator:Introduced.vegetation.predator)), 
      select_at(FW, vars(Forest.prey:Introduced.vegetation.prey))),
    BM.match = (FW$logBM.predator - FW$logBM.prey)^2
  )
  predictors <- data.frame(Predator = FW$Predator,
                           Prey = FW$Prey) %>%
    cbind(foraging_predictor) %>%
    cbind(vulnerability_predictor) %>%
    cbind(match_predictor) %>%
    mutate(ActivityTime.match = as.integer(ActivityTime.match))
  return(predictors)
}


scale2 <- function(x){(x - mean(x)) / (2*sd(x))}

get_predictors <- function(Species_List, FuncTraits){
  FW <- expand.grid(Species_List, Species_List)
  colnames(FW) <- c("Predator", "Prey")
  FW <- left_join(FW, FuncTraits, by = c("Prey" = "Species")) %>%
    left_join(FuncTraits, by = c("Predator" = "Species"))
  FW <- traits2predictors(FW)
  return(FW)
}

make_predictions <- function(Model, newdata, ndraws = 100, allow_new_levels = TRUE, extrapolation = F){
  predictions <- predict(Model, newdata = newdata, allow_new_levels = allow_new_levels, ndraws = ndraws, summary = F)
  rownames(predictions) <- paste0("draws", c(1:ndraws))
  predictions <- as.data.frame(t(predictions))
  predictions$Estimate <- apply(predictions, MARGIN = 1, mean)
  predictions$Est.Error <- apply(predictions[,-ncol(predictions)], MARGIN = 1, sd)
  predictions <- select(newdata, Predator, Prey) %>%
    bind_cols(predictions)
  if (!extrapolation){
    predictions$training <- ifelse(c(1:nrow(newdata) %in% as.numeric(rownames(Model$data))), 1, 0)
  }
  return(predictions)
}

