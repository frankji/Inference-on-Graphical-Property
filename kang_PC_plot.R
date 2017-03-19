# Load Library #
library('WGCNA')
library('igraph')
library('Hmisc')
library('clusterSim')
library('clime')
options(digits=16)
# Load data #
load('~/Project/hongyu/Project/brain/kang.rda')
load('~/Project/hongyu/Project/brain/data_all.rda')
load('~/Project/hongyu/Project/brain/gene_symbol.rda')
load('~/Project/hongyu/Project/brain/brain_regions.rda')
load('~/Project/hongyu/Project/brain/time_p.rda')
load('~/Project/hongyu/Project/brain/pat_mat.rda')

kang <- datET$datETcollapsed
#kang <- apply(kang, 1, function(x) (x - mean(x))/sd(x))
# Get corresponding time period for each column
ctime <- unlist(mapply(function(x, y) y[colnames(x)], data_all, time_p))
# Get corresponding regions for each column
cregion <- unlist(mapply(function(x, y) rep(y, ncol(x)), data_all, brain_regions))
# Function to get network information
getnets <- function(datai, pathw, i, j, P = TRUE){
  group <- rnorm(length(pathw[, i] == 1 | pathw[, j] == 1))
  X <- datai[pathw[, i] == 1 | pathw[, j] == 1, ]
  X_r <- rcorr(t(X))
  X_cor <- X_r$r
  if(P){
    X_cor <- X_r$P
    X_cor[X_cor > 0.05/(0.5*nrow(X_cor)*nrow(X_cor))] <- 0
  }
  diag(X_cor) <- 0
  g <- graph_from_adjacency_matrix(X_cor, mode = 'undirected', weighted = TRUE)
  group[pathw[, i] == 1] <- 1
  group[pathw[, j] == 1] <- 3
  group[pathw[, i] == 1 & pathw[, j] == 1] <- 2
  group <- group[pathw[, i] == 1 | pathw[, j] == 1]
  g$group <- group
  g$cor <- X_r$r
  g$X <- X
  return(g)
}

# Notch 170
which(colnames(pat_mat) == '157118')
# EGFR 35
which(colnames(pat_mat) == '177929')

# Get data given time and region
dat <- list()
for(i in 3:10){
  # J: spatial index
  j <- c('V1C', 'ITC', 'IPC', 'A1C', 'STC')
  i_time <- ctime %in% i
  j_region <- cregion %in% j
  Xi <- kang[, i_time & j_region]
  Xi <- t(apply(Xi, 1, function(x) (x - mean(x))/sd(x)))
  net <- getnets(Xi, pat_mat, 170, 35)
  dat[[i-2]] <- net$X
}

Ets <- read.csv('~/Project/hongyu/Pathway_Connectivity/Ets_VEts_V1C_ITC_IPC_A1C_STC_cross_walk_edges.csv', sep = ';')

library(networkD3)
library(randomNames) 
library(igraph)
plotG <- function(x, group){
  vnames <- unique(c(x[, 1], x[, 2]))
  vindex <- (1:length(vnames)) - 1
  names(vindex) <- vnames
  vertices <- data.frame(name = vnames,
                         group = group)
  g_edges <- data.frame(x[, c(1, 2)], stringsAsFactors = FALSE)
  g_edges$value <- x[, 3]
  colnames(g_edges) <- c("source", "target", "value")
  g_edges$source <- vindex[g_edges$source]
  g_edges$target <- vindex[g_edges$target]
  forceNetwork(Links = g_edges, Nodes = vertices, Source = "source",
               Target = "target", Value = "value", NodeID = "name",
               Group = "group", opacity = 1, zoom = TRUE, fontSize = 24)
}

group <- net$group

par(mfrow=c(3,3))
for(i in 3:10){
Et <- as.matrix(Ets[Ets[, i] == 1, c(1,2)])
tmp <- matrix(0, nrow = 429, ncol = 429)
tmp[Et] <- 1
rownames(tmp) <- 1:429
colnames(tmp) <- 1:429
net <- graph_from_adjacency_matrix(tmp, mode = 'undirected', weighted = TRUE)
V(net)$color <- c("gray50", "tomato", "gold")[group]
V(net)$size <- 3
V(net)$label <- NA
E(net)$width <- 1
E(net)$edge.color <- "black"
plot(net, layout=layout_with_fr, margin = c(0,0,0,0))
}


