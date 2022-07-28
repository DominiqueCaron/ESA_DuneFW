library(dplyr)
library(ggplot2)
library(igraph)
library(NetIndices)

DuneFW <- read.csv("data/output/DuneFWpredictions.csv", row.names = 1)

n_draws <- 100 # number of fw in global web
col_names <- c("size", "links", "connectance","meanTL", "maxTL","pl_exp", paste0("motif", c(1:13)))
duneFW_stats <- matrix(NA, nrow = n_draws, ncol = 19, dimnames = list(c(1:n_draws), col_names)) %>% as.data.frame()

for (i in c(1:n_draws)){
  fw <- DuneFW[,c("Prey", "Predator", paste0("draws", i))]
  fw <- fw[fw[,3] == 1,]
  g <- graph_from_edgelist(as.matrix(fw[,c(1,2)]))
  size <- vcount(g)
  links <- ecount(g)
  connectance <- links / (size^2)
  trophiclevels <- TrophInd(as_adjacency_matrix(g, sparse = F))$TL
  meanTL <- mean(trophiclevels)
  maxTL <- max(trophiclevels)
  pl_exp <- fit_power_law(degree(g), implementation = "R.mle")@coef
  motif_profile <- motifs(g)[c(5, 8, 12, 7, 3, 9, 14, 6, 10, 13, 16, 15, 11)]
  motif_profile <- motif_profile / sum(motif_profile)
  duneFW_stats[i, ] <- c(size, links, connectance, meanTL, maxTL, pl_exp, motif_profile)
}

# let's compare Dune FW to earth's FW
globalWeb_stats <- read.csv("data/output/GlobalWebStructure.csv", row.names = 1) %>%
  drop_na()
globalWeb_stats$logSize <- log10(globalWeb_stats$size)
globalWeb_stats$logLinks <- log10(globalWeb_stats$links)

# links-species relationship?
m <- lm(logLinks ~ logSize, data = globalWeb_stats)

mpi <- cbind(globalWeb_stats, predict(m, interval = "prediction"))

ggplot() +
  geom_point(aes(x = logSize, y = logLinks), data = globalWeb_stats, alpha = 0.3) +
  geom_smooth(aes(x = logSize, y = logLinks), data = globalWeb_stats, method = "lm") +
  geom_ribbon(aes(x = logSize, y = logLinks, ymin = lwr, ymax = upr), data = mpi,
              fill = "blue", alpha = 0.2) +
  geom_point(aes(x = log10(size), y = log10(links)), data = duneFW_stats, alpha = 0.7, size = 2,color = "gold")

# connectance relationship?
ggplot() +
  geom_point(aes(x = logSize, y = connectance), data = globalWeb_stats, alpha = 0.3) +
  geom_smooth(aes(x = logSize, y = connectance), data = globalWeb_stats, method = "loess", fill = "orange") +
  geom_point(aes(x = log10(size), y = connectance), data = duneFW_stats, alpha = 0.7, size = 2,color = "gold")

# trophic level
globalWeb_stats <- filter(globalWeb_stats, meanTL < 10)
ggplot() +
  geom_point(aes(x = logSize, y = meanTL), data = globalWeb_stats, alpha = 0.3) +
  geom_smooth(aes(x = logSize, y = meanTL), data = globalWeb_stats, method = "loess", fill = "orange") +
  geom_point(aes(x = log10(size), y = meanTL), data = duneFW_stats, alpha = 0.7, size = 2,color = "gold")

# power - law exponent
ggplot() +
  geom_point(aes(x = logSize, y = pl_exp), data = globalWeb_stats, alpha = 0.3) +
  geom_smooth(aes(x = logSize, y = pl_exp), data = globalWeb_stats, method = "loess", fill = "orange") +
  geom_point(aes(x = log10(size), y = pl_exp), data = duneFW_stats, alpha = 0.7, size = 2,color = "gold")

# motifs frequency
dune_motifs <- dplyr::select(duneFW_stats, starts_with("motif")) %>%
  pivot_longer(everything(), names_to = "motif", values_to = "frequency") %>%
  mutate(fw = "dune")

global_motifs <- dplyr::select(globalWeb_stats, starts_with("motif")) %>%
  pivot_longer(everything(), names_to = "motif", values_to = "frequency") %>%
  mutate(fw = "earth")

ggplot(aes(x = motif, y = frequency, colour = fw), data = rbind(dune_motifs, global_motifs)) +
  geom_boxplot()
