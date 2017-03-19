library('WGCNA')
library('igraph')
library('Hmisc')
library('clusterSim')
library('clime')
library('tsne')
library('plotly')
options(digits=16)

args <- commandArgs(trailingOnly = TRUE)

if(length(args)!=2){
  stop('Rscript Kang_PC_clean <cross_talk/connected_giant> <regions> <path1> <path2>')
}

print(args)
args <- c('connected_giant', 'V1C,ITC,IPC,A1C,STC', 170, 35)


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

theregion <- strsplit(args[2], split = ',')[[1]]
dat <- get_pathway_data(kang, 3:10, theregion, ctime, cregion, pat_mat, as.numeric(args[3]), as.numeric(args[4]))
n <- as.numeric(lapply(dat, ncol))
ts <- (rep(3:10, n) - 3) / 7
X <- do.call(cbind, dat)
groups <- matrix(0, nrow = nrow(X), ncol = 2)
net <- getnets(datai = kang, pathw = pat_mat, i = as.numeric(args[3]), j = as.numeric(args[4]))
groups[net$group == 1, 1] <- 1
groups[net$group == 3, 2] <- 1
groups[net$group == 2, ] <- 2
regions <- gsub('.+_(.+)_.+', '\\1', colnames(X))

plot_tsne <- function(mydata, axis.name = c('tsne1', 'tsne2'), col = 'brown1', mode = 'groups', group = 'group', shape = 16, fill = NA, gradient = NULL){
  shapes <- c('\u25CF', '\u25B2', '\u25C4', '\u25A0', '\u25AA', "\u25BC")
  temp <- shapes[1:length(unique(shape))]
  names(temp) <- unique(shape)
  shape_value <- temp
  toplot <- data.frame(tsne1 = mydata[, 1], tsne2 = mydata[, 2])
  colnames(toplot) <- axis.name
  if(mode == 'groups'){
    toplot[, group] <- col
    myplot <- ggplot(data = toplot)
    if(class(fill)!='numeric'){
      myplot <- myplot + geom_point(aes(x = tsne1, y = tsne2, shape = shape), size = 4, fill = 'black')
    }
    if(class(fill)=='numeric'){
      myplot <- myplot + geom_point(aes(x = tsne1, y = tsne2, shape = shape, col = fill), size = 4)
      myplot <- myplot + scale_color_gradient('Time', low = 'deepskyblue', high = 'brown1')
    }
    myplot <- myplot + theme_classic() + theme(text = element_text(size = 15))
    myplot <- myplot + scale_shape_manual('Regions', values = shape_value)
    return(myplot)
  }
}

# Explorations of outliers

#theregion <- c('V1C', 'ITC', 'IPC', 'A1C', 'STC')
#theregion <- c('MD', 'CBC')
#theregion <- c('M1C', 'S1C', 'VFC', 'MFC', 'DFC', 'OFC')
theregion <- c('STR', 'HTP', 'AMY')
dat <- get_pathway_data(kang, 3:10, theregion, ctime, cregion, pat_mat, as.numeric(args[3]), as.numeric(args[4]))
n <- as.numeric(lapply(dat, ncol))
ts <- (rep(3:10, n) - 3) / 7
X <- do.call(cbind, dat)
groups <- matrix(0, nrow = nrow(X), ncol = 2)
net <- getnets(datai = kang, pathw = pat_mat, i = as.numeric(args[3]), j = as.numeric(args[4]))
groups[net$group == 1, 1] <- 1
groups[net$group == 3, 2] <- 1
groups[net$group == 2, ] <- 2
regions <- gsub('.+_(.+)_.+', '\\1', colnames(X))
tsne_all <- tsne(t(X), k = 2)
plots <- plot_tsne(tsne_all, shape = regions, fill = ts)
plot_name <- paste(c('tsne_all', theregion), collapse = '_')
plot_name <- paste(plot_name, '.jpg', sep = '')
ggsave(plot = plots, filename = plot_name, width = 8, height = 6)

# Explorations of global expression pattern
kh <- function(h, zi, z) {
  as.numeric(zi == z)
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
theregion <- c('V1C', 'ITC', 'IPC', 'A1C', 'STC', 'MD', 'CBC', 'M1C', 'S1C', 'VFC', 'MFC', 'DFC', 'OFC', 'STR', 'HTP', 'AMY')
h <- 0.05
t <- (3 - 3) / 7
dat <- get_pathway_data(kang, 3:10, theregion, ctime, cregion, pat_mat, as.numeric(args[3]), as.numeric(args[4]))
n <- as.numeric(lapply(dat, ncol))
ts <- (rep(3:10, n) - 3) / 7
X <- do.call(cbind, dat)
groups <- matrix(0, nrow = nrow(X), ncol = 2)
net <- getnets(datai = kang, pathw = pat_mat, i = as.numeric(args[3]), j = as.numeric(args[4]))
groups[net$group == 1, 1] <- 1
groups[net$group == 3, 2] <- 1
groups[net$group == 2, ] <- 2
regions <- gsub('.+_(.+)_.+', '\\1', colnames(X))
t <- (3 - 3) / 7
sig_prec <- sig_est(t, X, ts, h)
heatmap(sig_prec$sigmas[[1]])
sigmas <- lapply(sig_prec$sigmas, function(x) as.matrix(1 - abs(x)))
sigmas <- lapply(sigmas, as.dist)
tsne_all <- tsne(sigmas[[1]], k = 2)
plot(tsne_all)