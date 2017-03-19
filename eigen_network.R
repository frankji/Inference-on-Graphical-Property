library('WGCNA')
library('topGO')
library('biomaRt')
library('ReactomePA')
load('~/Project/hongyu/Project/brain/kang.rda')
load('~/Project/hongyu/Project/brain/data_all.rda')
load('~/Project/hongyu/Project/brain/gene_symbol.rda')
load('~/Project/hongyu/Project/brain/brain_regions.rda')
load('~/Project/hongyu/Project/brain/time_p.rda')
load('~/Project/hongyu/Project/brain/pat_mat.rda')

# Get corresponding time period for each column
ctime <- unlist(mapply(function(x, y) y[colnames(x)], data_all, time_p))
# Get corresponding regions for each column
cregion <- unlist(mapply(function(x, y) rep(y, ncol(x)), data_all, brain_regions))


myregions <- list()
myregions[[1]] <- c('M1C','S1C','VFC','DFC','OFC')
myregions[[2]] <- c('MD', 'CBC')
myregions[[3]] <- c('STR', 'HTP', 'AMY')
myregions[[4]] <- c('V1C','ITC','IPC','A1C','STC')
i <- 1
kang <- datET$datETcollapsed
mydata <- kang[, ctime <= 10 & cregion %in% myregions[[i]]]
# Get the coexpression network, power 6
myadj <- adjacency(datExpr = t(mydata))
# Get the Topological Overlap Map
mytop <- TOMsimilarity(myadj, TOMDenom = 'min')
# Get the dissimilarity
mydis <- 1 - mytop
# Average linkage hierarchical clustering
myhier <- hclust(as.dist(mydis), method = 'average')
# Get modules by dynamic cutting of the branches
mymodules<- cutreeDynamic(myhier, method="tree", minClusterSize = 30, minSplitHeight = 0.15)
# Merge modules closed to each other
unmergedColors <- labels2colors(mymodules)
mymerged<- mergeCloseModules(t(mydata), mymodules, cutHeight = 0.15)
mymerged <- as.numeric(as.factor(mymerged$colors))
### Get GO first
# Load package for annotation functions
library('AnnotationFuncs')
# Load microarray index information database
library('hgu95av2.db')
# Load Gene Ontology database
library('GO.db')
# Load Genome wide annotation for human
library('org.Hs.eg.db')
# Get gene ontology data from gene symbol
gene_GOs<-translate(rownames(kang), remove.missing=FALSE, from=org.Hs.egSYMBOL2EG, to=org.Hs.egGO)
### Perform Enrichment analysis on genes in modules###
toenrich <- rownames(kang)[mymerged == unique(mymerged)[1]]
# Get access to Emsembl database
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl_dat <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'hgnc_symbol', 
                     mart = ensembl, values = rownames(kang), bmHeader = TRUE)
gene_ensembl <- ensembl_dat$`Gene ID`
names(gene_ensembl) <- ensembl_dat$`HGNC symbol`
# Get access to Entrez ID database
entrez <- as.list(org.Hs.egALIAS2EG)
gene_entrez <- unlist(lapply(entrez[rownames(kang)], function(x) ifelse(length(x) > 0, x[1], NA)))
gene_entrez <- gene_entrez[!is.na(gene_entrez)]
entrez_symbol <- names(gene_entrez)
names(entrez_symbol) <- gene_entrez
# Run GO enrichment analysis using topG0
tmp1 <- rep(0, nrow(kang))
names(tmp1) <- rownames(kang)
tmp1[toenrich] <- 1
zeroorone <- function(x) x > 0
sampleGOdata <- new('topGOdata',
                    description = 'Module gene GO enrichment', ontology = 'BP',
                    allGenes = tmp1, geneSel = zeroorone,
                    nodeSize = 10,
                    annotationFun = annFUN.org,
                    mapping = 'org.Hs.eg.db',
                    ID = 'symbol')
GO_EA <- runTest(sampleGOdata, algorithm = 'classic', statistic = 'fisher')
GO_EA@score <- p.adjust(GO_EA@score, method = 'none')
GO_result <- GenTable(sampleGOdata, Fis = GO_EA, topNodes = 500)

myeigens <-list()
mysvds <- list()
i <- 1
myregion <- cregion %in% myregions[[i]]
for(t in 3:10){
  tmp <- kang[, myregion & ctime == t]
  mysvd <- list()
  myeigen <- matrix(NA, nrow = ncol(tmp), ncol = length(unique(mymerged)))
  for(i in 1:length(unique(mymerged))){
    mysvd[[i]] <- svd(tmp[mymerged==i, ])
    myeigen[, i] <- mysvd[[i]]$v[, 1]
  }
  mysvds[[t-2]] <- mysvd
  myeigens[[t-2]] <- myeigen
}

kh <- function(h, zi, z) as.numeric(z == zi)

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


sig_est <- function(t, X, ts, h, mirror = FALSE){
  sigmas <- list()
  precs <- list()
  for(i in 1:length(t)){
    ks <- kh(h, ts, t[i])
    sigmas[[i]] <- (t(replicate(nrow(X), ks)) * X) %*% t(X) / sum(ks)
    precs[[i]] <- chol2inv(sigmas[[i]])
  }
  return(list(sigmas = sigmas, precs = precs))
}

connected_giant <- function(Et_ind, groups = NULL){
  Et_ind[which(Et_ind %*% Et_ind != 0, arr.ind = TRUE)] <- 1
  Ec <- which(Et_ind == 0, arr.ind = TRUE)
  Ec <- Ec[Ec[, 2] > Ec[, 1], ]
  results <- list(Ec = Ec, Et_ind = Et_ind)
  return(results)
}

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
dat <- lapply(myeigens, t)
n <- as.numeric(lapply(dat, ncol))
ts <- (rep(3:10, n) - 3) / 7
X <- do.call(cbind, dat)
sig_prec <- sig_est((3:10-3)/7, X, ts, h)
sigmas <- sig_prec$sigmas
precs <- sig_prec$precs
zspan <- (3:10 - 3)/7
E0_ind <- matrix(0, nrow = nrow(X), ncol = nrow(X))

Ets <- pathggm(X, h, precs, ts, E0_ind, connected_giant, zspan, NULL)
