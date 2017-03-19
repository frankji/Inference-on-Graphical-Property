# Function to compute MI for given mu and Sigma
mi <- function(mu, Sigma, partition){
  V <- chol2inv(Sigma)
  Sprime <- matrix(0, nrow=nrow(Sigma), ncol=ncol(Sigma))
  Vprime <- matrix(0, nrow=nrow(V), ncol=ncol(V))
  V12 <- V[partition, !partition]
  V21 <- V[!partition, partition]
  Sprime[partition, partition] <- Sigma[partition, partition]
  Sprime[!partition, !partition] <- Sigma[!partition, !partition]
  Vprime[partition, !partition] <- V12
  Vprime[!partition, partition] <- V21
  #  MI <- 0.5*log(det(Sprime)) - 0.5*log(det(Sigma))-0.5*sum(diag(Vprime%*%Sigma))
  MI <- 0.5*log(det(Sprime))-0.5*sum(diag(Vprime%*%Sigma))
  return(MI)
}

# Function to compute probability of MI by permutation.
pmi <- function(mu, Sigma, partition, n = 500){
  mi0 <- mi(mu, Sigma, partition)
  mis <- rep(0, n)
  for(i in 1:n){
    mis[i] <- mi(mu, Sigma, sample(partition))
  }
  return(sum(mis >= mi0, na.rm = TRUE)/n)
}

load('~/Project/hongyu/Project/brain/kang.rda')
load('~/Project/hongyu/Project/brain/data_all.rda')
load('~/Project/hongyu/Project/brain/gene_symbol.rda')
load('~/Project/hongyu/Project/brain/brain_regions.rda')
load('~/Project/hongyu/Project/brain/time_p.rda')
load('~/Project/hongyu/Project/brain/pat_mat.rda')

kang <- datET$datETcollapsed
getnets <- function(datai, pathw, i, j, P = TRUE){
  group <- rnorm(length(pathw[, i] == 1 | pathw[, j] == 1))
  X <- datai[pathw[, i] == 1 | pathw[, j] == 1, ]
  group[pathw[, i] == 1] <- 1
  group[pathw[, j] == 1] <- 3
  group[pathw[, i] == 1 & pathw[, j] == 1] <- 2
  group <- group[pathw[, i] == 1 | pathw[, j] == 1]
  g <- list()
  g$group <- group
  g$X <- X
  return(g)
}

# Get corresponding time period for each column
ctime <- unlist(mapply(function(x, y) y[colnames(x)], data_all, time_p))
# Get corresponding regions for each column
cregion <- unlist(mapply(function(x, y) rep(y, ncol(x)), data_all, brain_regions))
theregion <- c('M1C', 'S1C', 'VFC', 'MFC', 'DFC', 'OFC')

react <- read.csv('~/Project/hongyu/Project/brain/NCBI2Reactome_All_Levels.txt', as.is = TRUE, header = FALSE, sep = '\t')
react <- react[react$V6 == 'Homo sapiens', ]
react_relation <- read.csv('~/Project/hongyu/Project/brain/ReactomePathwaysRelation.txt', as.is = TRUE, header = FALSE, sep = '\t')
react_hsa <- c(react_relation[, 1][grep('HSA' ,react_relation[, 1])], react_relation[, 2][grep('HSA' ,react_relation[, 2])])
react_hsa <- unique(gsub('[^0-9]+', '', react_hsa))
react_relation <- apply(react_relation, 2, gsub, pattern = '[^0-9]+', replacement = '', perl = TRUE)
path_names <- unique(as.character(react_relation))
path_child <- unique(react_relation[, 2][!(react_relation[, 2] %in% react_relation[, 1])])
path_child <- path_child[path_child %in% react_hsa]
colnames(react) <- c('EntrezID', 'Pathway', 'Webinfo', 'Description', 'Evidentce', 'Species')
react$ID <- gsub('R-HSA-', '', react[, 2])
react_map <- unique(react[, c('Description', 'ID')])
rownames(react_map) <- react_map$ID

notch <- '157118'
EGFR <- '177929'
X <- getnets(kang[, cregion %in% theregion], pat_mat, notch, EGFR)
paths <- as.numeric(rownames(kang) %in% rownames(X$X))
paths <- pat_mat[, t(pat_mat) %*% paths / apply(pat_mat, 2, sum) == 1]
tmp <- t(paths) %*% paths
diag(tmp) <- 0
tmp <- tmp > 0
tmp <- tmp[rownames(tmp) %in% path_child, rownames(tmp) %in% path_child]
pat_shared <- tmp
pathsi <- as.numeric(rownames(kang) %in% rownames(X$X)[X$group!=3])
pathsi <- pat_mat[, t(pat_mat) %*% pathsi / apply(pat_mat, 2, sum) == 1]
tmp <- t(pathsi) %*% pathsi
diag(tmp) <- 0
tmp <- tmp > 0
tmp <- tmp[rownames(tmp) %in% path_child, rownames(tmp) %in% path_child]

react_selected <- (react_map[rownames(pat_shared), ])
paths_label <- rep('EGFR', nrow(pat_shared))
names(paths_label) <- rownames(pat_shared)
paths_label[rownames(tmp)] <- 'Notch'
react_selected$Pathway <- paths_label


results <- list()
for(t in 3:10){
  result <- matrix(-1, nrow=nrow(pat_shared), ncol=ncol(pat_shared))
  for(i in 2:nrow(pat_shared)){
    for(j in 1:(i-1)){
      if(pat_shared[i,j] == 0){
        pati <- rownames(pat_shared)[i]
        patj <- rownames(pat_shared)[j]
        X <- getnets(kang[, ctime %in% t & cregion %in% theregion], pat_mat, pati, patj)
        x <- X$X
        partition <- X$group == 1
        p <- nrow(x)
        N <- ncol(x)
        Sigma <- cov(t(x))
        result[i, j] <- 1 - pmi(apply(x, 1, mean), Sigma = Sigma, partition = partition)
      }
    }
  }
  results[[t-2]] <- result
}


output <- lapply(results, function(x) {x[which(x<=0)] <- 0
return(x)})
par(mfrow=c(2,4))
for(i in 1:8){
  net <- graph_from_adjacency_matrix(output[[i]], mode = 'undirected', weighted = TRUE)
  E(net)$width <- (E(net)$weight) ^ 3 * 10
  V(net)$color <- ifelse(paths_label == 'Notch', 'brown1', 'deepskyblue3')
  plot(net, layout=layout_with_fr, margin = rep(-0.08, 4))
}

