# Load Library #
library('glasso')
library('WGCNA')
library('igraph')
library('Hmisc')
library('clusterSim')
library('clime')
options(digits=16)

args <- commandArgs(trailingOnly = TRUE)

if(length(args)!=2){
  stop('Rscript Kang_PC_clean <cross_talk/connected_giant> <regions> <path1> <path2>')
}

print(args)

# Load data #
load('~/Project/hongyu/Project/brain/kang.rda')
load('~/Project/hongyu/Project/brain/data_all.rda')
load('~/Project/hongyu/Project/brain/gene_symbol.rda')
load('~/Project/hongyu/Project/brain/brain_regions.rda')
load('~/Project/hongyu/Project/brain/time_p.rda')
load('~/Project/hongyu/Project/brain/pat_mat.rda')

kang <- datET$datETcollapsed

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

get_pathway_data <- function(x, timespan, regionspan, ctime, cregion, pathway, path1, path2){
  dat <- list()
  for(i in 1:length(timespan)){
    i_time <- ctime %in% timespan[i]
    j_region <- cregion %in% regionspan
    xi <- x[, i_time & j_region]
    xi <- t(apply(xi, 1, function(x) (x - mean(x))/sd(x)))
    net <- getnets(xi, pathway, path1, path2)
    dat[[i]] <- net$X
  }
  return(dat)
}


################################################################################
# Tune the parameter of h
# Kernel Function
kh <- function(h, zi, z) {
  output <- 0.75*(1 - ((zi - z)/h)^2)/h*(abs((zi - z)/h) <= 1)
  output <- output * ((zi != z) + 1)
  return(output)
}
#hs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2)

#ind <- seq(1, sum(n))
#inds <- sapply(((3:10)-3)/7, function(x) ind[ts == x])


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
################################################################################


# Function to calculate sigma
sig_est <- function(t, X, ts, h, mirror = FALSE){
  sigmas <- list()
  precs <- list()
  for(i in 1:length(t)){
    ks <- kh(h, ts, t[i])
    sigmas[[i]] <- (t(replicate(nrow(X), ks)) * X) %*% t(X) / sum(ks)
    precs[[i]] <- glasso(sigmas[[i]], 0.01)$wi
    precs[[i]] <- precs[[i]] - (t(precs[[i]]) %*% (sigmas[[i]] %*% precs[[i]] - diag(nrow(precs[[i]])))) / replicate(nrow(X), apply(precs[[i]] * sigmas[[i]], 1, sum))
  }
  return(list(sigmas = sigmas, precs = precs))
}

# Estimator Te
Te <- function(x, h, precs, ts, epsilon, zs, i_z, E){
  # Whole dataset, p by n matrix in which p is the number of gene and n is sample size
  z <- zs[i_z]
  k <- kh(h, ts, z)
  ke <- k * epsilon
  tmat <- sqrt(length(ts)*h)*abs(t(precs[[i_z]]) %*% ((x * t(replicate(nrow(x), ke))) %*% t(x) %*% precs[[i_z]] - sum(ke)*diag(nrow(x)))) / length(ts)
  adj <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  adj[E] <- 1 
  max(tmat[adj > 0])
}

# Quantile Tc
Tc <- function(x, h, precs, ts, zspan, Ec, iter = 500, a = 0.95){
  teb <- c()
  for(i in 1:iter){
    epsilon <- rnorm(length(ts))
    teb_temp <- c()
    for(j in 1:length(zspan)){
      z <- zspan[j]
      teb_temp[j] <- Te(x, h, precs, ts, epsilon, zspan, j, Ec)
    }
    teb[i] <- max(teb_temp)
  }
  return(sort(teb)[as.integer(a*iter)])
}

# Functions to define critical edges

cross_talk <- function(Et_ind, groups){
  Ec <- which(Et_ind == 0, arr.ind = TRUE)
  Ec <- Ec[Ec[, 2] > Ec[, 1], ]
  Ec <- Ec[apply(groups[Ec[, 1], ] * groups[Ec[, 2], ], 1, sum) != 1, ]
  results <- list(Ec = Ec, Et_ind = Et_ind)
  return(results)
}

connected_giant <- function(Et_ind, groups = NULL){
  Et_ind[which(Et_ind %*% Et_ind != 0, arr.ind = TRUE)] <- 1
  Ec <- which(Et_ind == 0, arr.ind = TRUE)
  Ec <- Ec[Ec[, 2] > Ec[, 1], ]
  results <- list(Ec = Ec, Et_ind = Et_ind)
  return(results)
}

# Function to infer the graph by critical edges definition

pathggm <- function(x, h, precs, ts, E0_ind, cric_func, zspan, groups){
  # Et_ind contains edges considered as connected
  Ets <- list()
  for(i in 1:length(zspan)){
    z <- zspan[i]
    Rtnorm <- 1
    crics0 <- cric_func(E0_ind, groups)
    Ec <- crics0$Ec
    Ec_ind <- matrix(1, nrow = nrow(x), ncol = nrow(x))
    Ec_ind[Ec] <- 0
    Et_ind <- crics0$Et_ind
    Et <- matrix(NA, nrow = 0, ncol = 2)
    while(Rtnorm != 0){
      tebq <- Tc(x, h, precs, ts, zspan, Ec)
      Tb <- sqrt(length(ts)*h)*abs(precs[[i]])*sum(kh(h, ts, z))/length(ts)
      Rt <- which(Tb > tebq & Ec_ind == 0, arr.ind = TRUE)
      Rt <- matrix(Rt[Rt[, 1] < Rt[, 2], ], ncol = 2)
      Et <- rbind(Et, Rt)
      Et_ind[Et] <- 1
      crics <- cric_func(Et_ind, groups)
      Et_ind <- crics$Et_ind
      Ec <- crics$Ec
      Ec_ind <- matrix(1, nrow = nrow(x), ncol = nrow(x))
      Ec_ind[Ec] <- 0
      Rtnorm <- nrow(Rt)
    }
    Ets[[i]] <- Et
  }
  return(Ets)
}

h <- 0.5 


theregion <- strsplit(args[2], split = ',')[[1]]
dat <- get_pathway_data(kang, 3:10, theregion, ctime, cregion, pat_mat, args[3], args[4])
n <- as.numeric(lapply(dat, ncol))
ts <- (rep(3:10, n) - 3) / 7
X <- do.call(cbind, dat)
groups <- matrix(0, nrow = nrow(X), ncol = 2)
net <- getnets(datai = kang, pathw = pat_mat, i = args[3], j = args[4])
groups[net$group == 1, 1] <- 1
groups[net$group == 3, 2] <- 1
groups[net$group == 2, ] <- 2

sig_prec <- sig_est((3:10-3)/7, X, ts, h)
sigmas <- sig_prec$sigmas
precs <- sig_prec$precs

zspan <- (3:10 - 3)/7
E0_ind <- matrix(0, nrow = nrow(X), ncol = nrow(X))
Ets <- pathggm(X, h, precs, ts, E0_ind, get(args[1]), zspan, groups)
fname <- paste(c('Ets', theregion, args[1], args[3], args[4], '.rda'), collapse = '_')
fname <- gsub('_\\.rda', '\\.rda', fname)
save(Ets, file = fname)