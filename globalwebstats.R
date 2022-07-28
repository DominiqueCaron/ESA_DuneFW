library(dplyr)
library(tidyr)
library(fread)
library(igraph)
library(NetIndices)

to_edge_list <- function(fw){
  out <- c()
  for (j in c(1:ncol(fw))){
    for (i in c(1:nrow(fw))){
      if (fw[i,j] > 0){
        out <- rbind(out,
                     c(rownames(fw)[i], colnames(fw)[j]))
      }
    }
  }
  return(out)
}

n_web <- 359 # number of fw in global web
col_names <- c("ref", "size", "links", "connectance","meanTL", "maxTL","pl_exp", paste0("motif", c(1:13)))
globalweb_stats <- matrix(NA, nrow = n_web, ncol = 20, dimnames = list(c(1:n_web), col_names)) %>% as.data.frame()

for (i in c(1:n_web)){
  path2file <- paste0("data/raw/GlobalWeb/FoodsWeb_35079406/WEB",i,".csv")
  if( file.exists(path2file)){
    globalweb_stats[i,"ref"] <- as.character(fread(path2file, nrows = 1, header = F)[,1])
    row2drop <- which(is.na(fread(path2file, skip = 1)[,2]))
    fw <- fread(path2file, check.names = F, skip = 1) %>%
      drop_na() %>% as.data.frame()
    fw <- fw[!duplicated(fw[,1]),]
    fw <- fw[,!duplicated(colnames(fw))]
    if (nrow(fw) >1 & i != 333){
      fw <- data.frame(fw, row.names = 1, check.names = F)
      fw[fw == "<0.01" | fw == "<0.001"] <- "0"
      fw[fw == "1a" |fw == "1b" |fw == "1c" |fw == "2a" |fw == "2b" |fw == "2c"] <- "1"
      fw <- fw %>% mutate_all(as.numeric)
      fw <- to_edge_list(fw)
      g <- graph_from_edgelist(fw)
      size <- vcount(g)
      links <- ecount(g)
      connectance <- links / (size^2)
      trophiclevels <- TrophInd(as_adjacency_matrix(g, sparse = F))$TL
      meanTL <- mean(trophiclevels)
      maxTL <- max(trophiclevels)
      pl_exp <- fit_power_law(degree(g), implementation = "R.mle")@coef
      motif_profile <- motifs(g)[c(5, 8, 12, 7, 3, 9, 14, 6, 10, 13, 16, 15, 11)]
      motif_profile <- motif_profile / sum(motif_profile)
      globalweb_stats[i, c(2:20)] <- c(size, links, connectance, meanTL, maxTL, pl_exp, motif_profile)
    }
  }
}

write.csv(globalweb_stats, "data/output/GlobalWebStructure.csv")
