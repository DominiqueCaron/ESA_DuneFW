library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)
library(NetIndices)

earth_color <- "#f3ecdf"
dune_color <- "#ec8c14"
background_color <- "#537479"

DuneFW <- read.csv("data/output/DuneFWpredictions.csv", row.names = 1)

n_draws <- 1000 # number of fw in global web
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
globalWeb_stats$logConnectance <- log10(globalWeb_stats$connectance)

# Get quantiles for dune stats
quantiles <- c(0.01,0.05,0.25,0.5,0.75,0.95,0.99)
duneFW_summary <- apply(duneFW_stats, 2, quantile, quantiles) %>% t() %>% as.data.frame()

duneFW_summary$metric <- rownames(duneFW_summary)

# links-species relationship?
m <- lm(logLinks ~ logSize, data = globalWeb_stats)
mpi <- cbind(globalWeb_stats, predict(m, interval = "prediction"))


ggplot() +
  geom_point(aes(x = logSize, y = logLinks), data = globalWeb_stats, alpha = 0.3, color = earth_color) +
  geom_smooth(aes(x = logSize, y = logLinks), data = globalWeb_stats, method = "lm", se = F, color = earth_color) +
  geom_pointrange(aes(x = log10(8), y = log10(get("50%")), ymin = log10(get("5%")), ymax = log10(get("95%"))), data = duneFW_summary["links",],color = dune_color, size = 1.25, fatten = 2) +
  geom_ribbon(aes(x = logSize, y = logLinks, ymin = lwr, ymax = upr), data = mpi,
              fill = "transparent", color = earth_color, linetype = "dashed") +
  theme_classic() +
  labs(y = "", x = "") +
  theme(plot.background = element_rect(fill = background_color, colour = background_color),
        panel.background = element_rect(fill = background_color, colour = background_color),
        axis.text = element_text(size = 15, colour = earth_color),
        axis.line = element_line(colour = earth_color),
        axis.ticks = element_line(colour = earth_color))

ggsave("figures/links.png", width =6, height = 6)

# connectance relationship?
m <- lm(logConnectance ~ logSize, data = globalWeb_stats)
mpi <- cbind(globalWeb_stats, predict(m, interval = "prediction"))

ggplot() +
  geom_point(aes(x = logSize, y = logConnectance), data = globalWeb_stats, alpha = 0.3, color = earth_color) +
  geom_smooth(aes(x = logSize, y = logConnectance), data = globalWeb_stats, method = "lm", se = F, color = earth_color) +
  geom_pointrange(aes(x = log10(8), y = log10(get("50%")), ymin = log10(get("5%")), ymax = log10(get("95%"))), data = duneFW_summary["connectance",],color = dune_color, size = 0.75) +
  geom_ribbon(aes(x = logSize, y = logConnectance, ymin = lwr, ymax = upr), data = mpi,
              fill = "transparent", color = earth_color, linetype = "dashed") +
  theme_classic() +
  labs(y = "log10(connectance)", x = "log10(# of species)")

ggsave("figures/connectance.png", width =6, height = 6)

# mean trophic level
globalWeb_stats_tl <- filter(globalWeb_stats, meanTL < 10)
m <- lm(meanTL ~ logSize, data = globalWeb_stats_tl)
mpi <- cbind(globalWeb_stats_tl, predict(m, interval = "prediction"))

ggplot() +
  geom_point(aes(x = logSize, y = meanTL), data = globalWeb_stats_tl, alpha = 0.3, color = earth_color) +
  geom_smooth(aes(x = logSize, y = meanTL), data = globalWeb_stats_tl, method = "lm", se = F, color = earth_color) +
  geom_pointrange(aes(x = log10(8), y = get("50%"), ymin = get("5%"), ymax = get("95%")), data = duneFW_summary["meanTL",],color = dune_color, size = 0.75) +
  geom_ribbon(aes(x = logSize, y = meanTL, ymin = lwr, ymax = upr), data = mpi,
              fill = "transparent", color = earth_color, linetype = "dashed") +
  lims(y=c(0,6)) +
  theme_classic() +
  labs(y = "Mean trophic level", x = "log10(# of species)")

ggsave("figures/meanTL.png", width =6, height = 6)

# max trophic level
m <- lm(maxTL ~ logSize, data = globalWeb_stats_tl)
mpi <- cbind(globalWeb_stats_tl, predict(m, interval = "prediction"))

ggplot() +
  geom_point(aes(x = logSize, y = maxTL), data = globalWeb_stats_tl, alpha = 0.3, color = earth_color) +
  geom_smooth(aes(x = logSize, y = maxTL), data = globalWeb_stats_tl, method = "lm", se = F, color = earth_color) +
  geom_pointrange(aes(x = log10(8), y = get("50%"), ymin = get("5%"), ymax = get("95%")), data = duneFW_summary["maxTL",],color = dune_color, size = 0.75) +
  geom_ribbon(aes(x = logSize, y = maxTL, ymin = lwr, ymax = upr), data = mpi,
              fill = "transparent", color = earth_color, linetype = "dashed") +
  lims(y=c(0,8)) +
  theme_classic() +
  labs(y = "Max trophic level", x = "log10(# of species)")

ggsave("figures/meanTL.png", width =6, height = 6)

# power - law exponent
ggplot() +
  geom_point(aes(x = logSize, y = pl_exp), data = globalWeb_stats, alpha = 0.3, color = earth_color) +
  geom_smooth(aes(x = logSize, y = pl_exp), data = globalWeb_stats, method = "loess", fill = earth_color, color = earth_color) +
  geom_pointrange(aes(x = log10(8), y = get("50%"), ymin = get("5%"), ymax = get("95%")), data = duneFW_summary["pl_exp",],color = dune_color, size = 0.75) +
  theme_classic() +
  labs(y = "", x = "") +
  theme(plot.background = element_rect(fill = background_color, colour = background_color),
        panel.background = element_rect(fill = background_color, colour = background_color),
        axis.text = element_text(size = 15, colour = earth_color),
        axis.line = element_line(colour = earth_color),
        axis.ticks = element_line(colour = earth_color))

ggsave("figures/pl_exponent.png", width = 6, height = 6)

# motifs frequency
dune_motifs <- dplyr::select(duneFW_stats, starts_with("motif")) %>%
  pivot_longer(everything(), names_to = "motif", values_to = "frequency") %>%
  mutate(fw = "dune", motif = factor(motif, levels = paste0("motif", c(1:13))))

global_motifs <- dplyr::select(globalWeb_stats, starts_with("motif")) %>%
  pivot_longer(everything(), names_to = "motif", values_to = "frequency") %>%
  mutate(fw = "earth", motif = factor(motif, levels = paste0("motif", c(1:13))))

plot_data = rbind(global_motifs, dune_motifs) %>%
  mutate(fw = factor(fw, levels = c("earth", "dune")))

ggplot(aes(x = motif, y = frequency, fill = fw), color = "black", data = plot_data) +
  geom_boxplot() +
  scale_fill_manual(values = c(earth_color, dune_color)) +
  theme_classic() +
  labs(x = "", y = "") +
  theme(plot.background = element_rect(fill = background_color, colour = background_color),
        panel.background = element_rect(fill = background_color, colour = background_color),
        axis.text = element_text(size = 15, colour = earth_color),
        axis.line = element_line(colour = earth_color),
        axis.ticks = element_line(colour = earth_color),
        legend.background = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.7,0.7),
        legend.key.size = unit(30, "pt"))

ggsave("figures/motifs.png", width = 10, height = 4)
