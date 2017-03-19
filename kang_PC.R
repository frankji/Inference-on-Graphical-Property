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

#Exploration before calculation. To see if the covariance matrix is close to each other
#`1covdist <- matrix(0, nrow = 8, ncol = 8)
#`1for(i in 1:8){
#`1  for(j in 1:8){
#`1    cvi <- dat[[i]] %*% t(dat[[i]])
#`1    cvj <- dat[[j]] %*% t(dat[[j]])
#`1    covdist[i,j] <- sum((cvi %*% chol2inv(cvj) - diag(nrow(cvi)))^2)
#`1  }
#`1}




add <- function(x) Reduce("+", x)
# Tune the parameter of h
kh <- function(h, zi, z) 0.75*(1 - ((zi - z)/h)^2)/h*(abs((zi - z)/h) <= 1)
hs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2)

n <- as.numeric(lapply(dat, ncol))
ind <- seq(1, sum(n))
ts <- (rep(3:10, n) - 3) / 7
inds <- sapply(((3:10)-3)/7, function(x) ind[ts == x])
X <- do.call(cbind, dat)
groups <- matrix(0, nrow = nrow(X), ncol = 2)
groups[net$group == 1, 1] <- 1
groups[net$group == 3, 2] <- 1
groups[net$group == 2, ] <- 2
#errs <- list()
#for(u in 1:11){
#  print(u)
#  errors <- matrix(0, nrow = 20, ncol = 8)
#  for(v in 3:10){
#    z <- (v - 3)/7
#    ksi <- kh(hs[u], ts, z)
#    for(m in 1:20){
#      indm <- do.call(c, lapply(inds, function(x) sample(x, length(x)/2)))
#      sigmai <- (ksi[indm] * X[, indm]) %*% t(X[, indm]) / sum(ksi[indm])
#      thetai <- chol2inv(sigmai)
#      sigmai2 <- (ksi[-indm] * X[, -indm]) %*% t(X[, -indm]) / sum(ksi[-indm])
#      errors[m, v-2] <- sum((sigmai2 %*% thetai - diag(nrow(sigmai)))^2)
#    }
#  }
#  errs[[u]] <- errors 
#}
#
#boxplot(lapply(errs, log), names = hs)
#
#save(errs, file = '~/Project/hongyu/Pathway_Connectivity/tmp/errs.rda')

# Set h as 0.5
h <- 0.5

sig <- function(z, X, ts, h){
  ks <- kh(h, ts, z)
  (t(replicate(nrow(X), ks)) * X) %*% t(X) / sum(ks)
}

## Estimation of precision matrix
sigmas <- list()
precs <- list()
for(t in 3:10){
  z <- (t - 3)/7 
  sigmai <- sig(z, X, ts, h)
  thetai <- chol2inv(sigmai)
  sigmas[[t-2]] <- sigmai 
  #  km(thetai[1, ], sigmai, thetai) / as.numeric(t(thetai[1, ]) %*% sigmai[1, ])
  # thetadi <- thetai - (t(thetai) %*% (sigmai %*% thetai - diag(nrow(thetai)))) / replicate(nrow(X), apply(thetai * sigmai, 1, sum))
  precs[[t-2]] <- thetai 
}

precs <- lapply(precs, function(x) {
  x[(x - t(x) < 0)] <- t(x)[(x - t(x) < 0)]
  x})
## New estimator TE
n <- as.numeric(lapply(dat, ncol))

epsilon <- rnorm(length(ts))

Te <- function(x, h, precs, ts, epsilon, z, E){
  t <- 7 * z + 3
  k <- kh(h, ts, z)
  ke <- k * epsilon
  tmat <- sqrt(length(ts)*h)*abs(t(precs[[t-2]]) %*% ((x * t(replicate(nrow(x), ke))) %*% t(x) %*% precs[[t-2]] - sum(ke)*diag(nrow(x)))) / length(ts)
  adj <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  adj[E] <- 1 
  max(tmat[adj > 0])
}


# temp data to build igraph to initialize edge set
Ets <- list()

for(t in 3:10){
  tmp <- sigmas[[1]]
  rownames(tmp) <- 1:nrow(X)
  colnames(tmp) <- 1:nrow(X)
  diag(tmp) <- 0
  E_igraph <- graph_from_adjacency_matrix(tmp, mode = 'undirected', weighted = TRUE)
  Ec <- get.edgelist(E_igraph)
  Ec <- apply(Ec, 2, as.numeric)
  E0 <- Ec
  Et_ind <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  Et <- matrix(NA, nrow = 0, ncol =2)
  Vt <- 1:nrow(X)
  
  iter <- 500
  a <- 0.95 * iter
  #ks <- kh(h, ts, z)
  Rtnorm <- 1
  # Start iteration, stop when |Rt| = 0
  while(Rtnorm != 0){
    teb <- c()
    for(i in 1:iter){
      epsilon <- rnorm(length(ts))
      teb_temp <- c()
      for(zi in 3:10){
        z <- (zi - 3)/7
        teb_temp[zi-2] <- Te(X, h, precs, ts, epsilon, z, Ec)
      }
      teb[i] <- max(teb_temp)
    }
    
    tebq <- sort(teb)[a]
    # Get the edge above the quantile
    Rt <- which(sqrt(length(ts)*h)*abs(precs[[t - 2]])*sum(kh(h, ts, (t - 3) / 7))/length(ts) > tebq & Et_ind == 0, arr.ind = TRUE)
    Rt <- Rt[Rt[, 1] < Rt[, 2], ]
    Et <- unique(rbind(Et, Rt))
    Et_ind[Et] <- 1
    Et_ind[which(Et_ind %*% Et_ind != 0, arr.ind = TRUE)] <- 1
    Ec <- which(Et_ind == 0, arr.ind = TRUE)
    Ec <- Ec[Ec[, 2] > Ec[, 1], ]
    Rtnorm <- ifelse(class(Rt) == 'numeric', 1, nrow(Rt))
  }
  
  Ets[[t-2]] <- Et
  
}
Ets_V1C_ITC_IPC_A1C_STC <- Ets
#save(Ets_V1C_ITC_IPC_A1C_STC, file = '~/Project/hongyu/Pathway_Connectivity/Ets_V1C_ITC_IPC_A1C_STC.rda')
tmp <- matrix(1, nrow = 429, ncol = 429)
diag(tmp) <- 0
E_all <- get.edgelist(graph_from_adjacency_matrix(tmp, weighted = TRUE, mode = 'undirected'))
Ets_V1C_ITC_IPC_A1C_STC_edges <- lapply(Ets_V1C_ITC_IPC_A1C_STC_cross_walk, function(x, tmp) {
  tmp[x] <- 0
  as.numeric(tmp[E_all] == 0)
  }, tmp = tmp)
Ets_V1C_ITC_IPC_A1C_STC_cross_walk_edges <- do.call(cbind, Ets_V1C_ITC_IPC_A1C_STC_edges)
Ets_V1C_ITC_IPC_A1C_STC_cross_walk_edges <- cbind(E_all, Ets_V1C_ITC_IPC_A1C_STC_cross_walk_edges)
write.csv2(Ets_V1C_ITC_IPC_A1C_STC_cross_walk_edges, file = '~/Project/hongyu/Pathway_Connectivity/Ets_VEts_V1C_ITC_IPC_A1C_STC_cross_walk_edges.csv', row.names = FALSE, col.names = FALSE)
#Connected giant for individual group
Ets <- list()
precs1 <- lapply(precs, function(x) x[net$group != 3, net$group != 3])
for(t in 3:10){
  tmp <- sigmas[[1]][net$group != 3, net$group != 3]
  rownames(tmp) <- 1:nrow(tmp)
  colnames(tmp) <- 1:nrow(tmp)
  diag(tmp) <- 0
  E_igraph <- graph_from_adjacency_matrix(tmp, mode = 'undirected', weighted = TRUE)
  Ec <- get.edgelist(E_igraph)
  Ec <- apply(Ec, 2, as.numeric)
  E0 <- Ec
  Et_ind <- matrix(0, nrow = nrow(tmp), ncol = nrow(tmp))
  Et <- matrix(NA, nrow = 0, ncol =2)
  Vt <- 1:nrow(tmp)
  
  iter <- 500
  a <- 0.95 * iter
  #ks <- kh(h, ts, z)
  Rtnorm <- 1
  # Start iteration, stop when |Rt| = 0
  while(Rtnorm != 0){
    teb <- c()
    for(i in 1:iter){
      epsilon <- rnorm(length(ts))
      teb_temp <- c()
      for(zi in 3:10){
        z <- (zi - 3)/7
        teb_temp[zi-2] <- Te(X[net$group != 3, ], h, precs1, ts, epsilon, z, Ec)
      }
      teb[i] <- max(teb_temp)
    }
    
    tebq <- sort(teb)[a]
    # Get the edge above the quantile
    Rt <- which(sqrt(length(ts)*h)*abs(precs1[[t - 2]])*sum(kh(h, ts, (t - 3) / 7))/length(ts) > tebq & Et_ind == 0, arr.ind = TRUE)
    Rt <- Rt[Rt[, 1] < Rt[, 2], ]
    Et <- unique(rbind(Et, Rt))
    Et_ind[Et] <- 1
    Et_ind[which(Et_ind %*% Et_ind != 0, arr.ind = TRUE)] <- 1
    Ec <- which(Et_ind == 0, arr.ind = TRUE)
    Ec <- Ec[Ec[, 2] > Ec[, 1], ]
    Rtnorm <- ifelse(class(Rt) == 'numeric', 1, nrow(Rt))
  }
  
  Ets[[t-2]] <- Et
  
}


### Cross walk ###

Ets <- list()

for(t in 3:10){
  
  tmp <- sigmas[[1]]
  rownames(tmp) <- 1:nrow(X)
  colnames(tmp) <- 1:nrow(X)
  diag(tmp) <- 0
  E_igraph <- graph_from_adjacency_matrix(tmp, mode = 'undirected', weighted = TRUE)
  Ec <- get.edgelist(E_igraph)
  Ec <- apply(Ec, 2, as.numeric)
  Ec <- Ec[apply(groups[Ec[, 1], ] * groups[Ec[, 2], ], 1, sum) != 1, ] 
  E0 <- Ec
  Et_ind <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  Et <- matrix(NA, nrow = 0, ncol =2)
  Vt <- 1:nrow(X)
  
  iter <- 500
  a <- 0.95 * iter
  #ks <- kh(h, ts, z)
  Rtnorm <- 1
  # Start iteration, stop when |Rt| = 0
  while(Rtnorm != 0){
    teb <- c()
    for(i in 1:iter){
      epsilon <- rnorm(length(ts))
      teb_temp <- c()
      for(zi in 3:10){
        z <- (zi - 3)/7
        teb_temp[zi-2] <- Te(X, h, precs, ts, epsilon, z, Ec)
      }
      teb[i] <- max(teb_temp)
    }
    
    tebq <- sort(teb)[a]
    Ec_ind <- matrix(1, nrow = nrow(X), ncol = nrow(X))
    Ec_ind[Ec] <- 0
    # Get the edge above the quantile
    Rt <- which(sqrt(length(ts)*h)*abs(precs[[t - 2]])*sum(kh(h, ts, (t - 3) / 7))/length(ts) > tebq & Ec_ind == 0, arr.ind = TRUE)
    Rt <- Rt[Rt[, 1] < Rt[, 2], ]
    Et <- unique(rbind(Et, Rt))
    Et_ind[Et] <- 1
    #    Et_ind[which(Et_ind %*% Et_ind != 0, arr.ind = TRUE)] <- 1
    Ec <- which(Et_ind == 0, arr.ind = TRUE)
    Ec <- Ec[Ec[, 2] > Ec[, 1], ]
    Ec <- Ec[apply(groups[Ec[, 1], ] * groups[Ec[, 2], ], 1, sum) != 1, ]
    Rtnorm <- ifelse(class(Rt) == 'integer', 1, nrow(Rt))
  }
  
  Ets[[t-2]] <- Et
}
Ets_V1C_ITC_IPC_A1C_STC_cross_walk <- Ets
save(Ets_V1C_ITC_IPC_A1C_STC_cross_walk, file = '~/Project/hongyu/Pathway_Connectivity/Ets_V1C_ITC_IPC_A1C_STC_cross_walk.rda')
