# Load Library #b
library('WGCNA')
library('igraph')
library('Hmisc')
library('clusterSim')
library('clime')

for(i in dir('result/')){
  load(file = paste('result/', i, sep = ''))
  assign(gsub('\\.rda', '', i), Ets)
}

# Load data #
load('~/Project/hongyu/Project/brain/kang.rda')
load('~/Project/hongyu/Project/brain/data_all.rda')
load('~/Project/hongyu/Project/brain/gene_symbol.rda')
load('~/Project/hongyu/Project/brain/brain_regions.rda')
load('~/Project/hongyu/Project/brain/time_p.rda')
load('~/Project/hongyu/Project/brain/pat_mat.rda')

library(randomNames) 
library(igraph)
kang <- datET$datETcollapsed

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

net <- getnets(datai = kang, pathw = pat_mat, i = 170, j = 35)


for(k in 1:length(dir('mirror_reflect/'))){
  if(length(grep('rda', dir('mirror_reflect/')[k])) == 0) next
  load(paste('mirror_reflect/', dir('mirror_reflect/')[k], sep = ''))
  fnamek <- gsub('\\.rda', '', dir('mirror_reflect/')[k])
  fnamek <- gsub('_+', '_', fnamek)
  fnamek <- gsub('_$', '', fnamek)
  fnamesk <- strsplit(fnamek, split = '_')[[1]]
  nfnamesk <- length(fnamesk)
  path1 <- as.numeric(fnamesk[nfnamesk - 1])
  path2 <- as.numeric(fnamesk[nfnamesk])
  myregion <- paste(fnamesk[2:(nfnamesk-4)], collapse = ',')
  mytask <- paste('', fnamesk[nfnamesk-3], fnamesk[nfnamesk-2], sep = ' ')
  jpeg(paste('mirror_reflect/', fnamek, '.jpg', sep = ''), width = 1600, height = 1200)
#  Ets <- get(gsub('\\.rda', '', dir('mirror_reflect/')[k]))
  Ets <- cbind(do.call(rbind, Ets), rep(3:10, lapply(Ets, nrow)))
  net <- getnets(datai = kang, pathw = pat_mat, i = path1, j = path2)
  X <- net$X
  group <- net$group
  par(mfrow=c(2,4))
  par(oma=c(0,0,4,0))
  for(i in 3:10){
    Et <- as.matrix(Ets[Ets[, 3] == i, c(1,2)])
    tmp <- matrix(0, nrow = nrow(X), ncol = nrow(X))
    tmp[Et] <- 1
    rownames(tmp) <- 1:nrow(X)
    colnames(tmp) <- 1:nrow(X)
    net <- graph_from_adjacency_matrix(tmp, mode = 'undirected', weighted = TRUE)
    V(net)$color <- c("gray50", "tomato", "gold")[group]
    if(path1 == path2 & path1 == 170){
      V(net)$color <- rep('gray50', length(group))
    }
    if(path1 == path2 & path1 == 35){
      V(net)$color <- rep('gold', length(group))
    }
    V(net)$frame.color <- V(net)$color
    V(net)$size <- 5
    V(net)$label <- NA
    E(net)$width <- 3
    E(net)$color <- "deepskyblue4"
    plot(net, layout=layout_with_fr, margin = rep(-0.08, 4))
  }
  if(path1 == path2 & path1 == 35){
    mytask <- paste(' EGFR', mytask, sep = '')
  }
  if(path1 == path2 & path1 == 170){
    mytask <- paste(' Notch', mytask, sep = '')
  }
  if(path1 != path2){
    mytask <- paste(' EGFR/Notch', mytask, sep = '')
  }
  mtext(paste(myregion, mytask, sep = ''), outer = TRUE, cex = 3)
  dev.off()
}