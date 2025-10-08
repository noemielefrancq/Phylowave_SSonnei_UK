########################################################################################################################################
## Extract detected lineages from tree: output: Lineages_detected_20251007.Rdata
########################################################################################################################################
## Packages used
library(stringr)
library(ape)
library(phangorn)
library(BactDating)
library(phytools)
library(coda)
library(thd)
library(vcfR)
library(lubridate)
library(ggplot2)
library(ggtree)
library(extrafont)
loadfonts(device="all")
library(MetBrewer)

setwd('~/Dropbox/Projects/2025_Phylowave_SSonnei/Phylowave_SSonnei/')
load('2_analysis_index/1_index_computations/Initial_index_computation_and_parameters_20251007.Rdata')

########################################################################################################################################
## Useful functions
########################################################################################################################################
## Other functions
mean.and.ci <-function(v){ return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))}
axisPhylo_NL = function (side = 1, root.time = NULL, backward = TRUE, at_axis = NULL, lab_axis = NULL, ...){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  type <- lastPP$type
  if (type == "unrooted")
    stop("axisPhylo() not available for unrooted plots; try add.scale.bar()")
  if (type == "radial")
    stop("axisPhylo() not meaningful for this type of plot")
  if (is.null(root.time))
    root.time <- lastPP$root.time
  if (type %in% c("phylogram", "cladogram")) {
    xscale <- if (lastPP$direction %in% c("rightwards", "leftwards"))
      range(lastPP$xx)
    else range(lastPP$yy)
    tmp <- lastPP$direction %in% c("leftwards", "downwards")
    tscale <- c(0, xscale[2] - xscale[1])
    if (xor(backward, tmp))
      tscale <- tscale[2:1]
    if (!is.null(root.time)) {
      tscale <- tscale + root.time
      if (backward)
        tscale <- tscale - xscale[2]
    }
    beta <- diff(xscale)/diff(tscale)
    alpha <- xscale[1] - beta * tscale[1]
    if(is.null(at_axis) == T){
      x <- beta * lab + alpha
      lab <- pretty(tscale)
    }
    # if(is.null(at_axis) != F){
    x <- at_axis
    lab <- lab_axis
    # }
    axis(side = side, at = x, labels = lab, ...)
  }
  else {
    n <- lastPP$Ntip
    xx <- lastPP$xx[1:n]
    yy <- lastPP$yy[1:n]
    r0 <- max(sqrt(xx^2 + yy^2))
    alpha <- sort(setNames(rect2polar(xx, yy)$angle, 1:n))
    angles <- c(diff(alpha), 2 * pi - alpha[n] + alpha[1L])
    j <- which.max(angles)
    i <- if (j == 1L)
      n
    else j - 1L
    firstandlast <- as.integer(names(angles[c(i, j)]))
    theta0 <- mean(atan2(yy[firstandlast], xx[firstandlast]))
    x0 <- r0 * cos(theta0)
    y0 <- r0 * sin(theta0)
    inc <- diff(pretty(c(0, r0))[1:2])
    srt <- 360 * theta0/(2 * pi)
    coef <- -1
    if (abs(srt) > 90) {
      srt <- srt + 180
      coef <- 1
    }
    len <- 0.025 * r0
    r <- r0
    while (r > 1e-08) {
      x <- r * cos(theta0)
      y <- r * sin(theta0)
      if (len/r < 1) {
        ra <- sqrt(len^2 + r^2)
        thetaa <- theta0 + coef * asin(len/r)
        xa <- ra * cos(thetaa)
        ya <- ra * sin(thetaa)
        segments(xa, ya, x, y)
        text(xa, ya, r0 - r, srt = srt, adj = c(0.5,
                                                1.1), ...)
      }
      r <- r - inc
    }
    segments(x, y, x0, y0)
  }
}
read.chains.from.table = function(table){
  if(typeof(table) != 'double') {
    print('Changing type of table to matrix')
    table = as.matrix(table)
  }
  Chains = list()
  col_variables = sapply(colnames(table), function(x)str_split(x, pattern = "[.]")[[1]][1])
  variable_names = unique(col_variables)
  nchains = nrow(table)
  for(i in 1:length(variable_names)){
    a = match(col_variables, variable_names[i])
    if(length(which(is.na(a) == F)) == 1) {
      Chains[[i]] = table[,match(variable_names[i], col_variables)]
    }
    else{
      a = match(col_variables, variable_names[i])
      tmp = colnames(table)[which(is.na(a) == F)]
      ndims = 0
      dims = NULL
      empty = F
      while(empty == F & ndims < 10){
        l = sapply(tmp, function(x)str_split(x, pattern = "[.]")[[1]][1+ndims+1])
        if(all(is.na(l))) empty = T
        if(all(is.na(l)) == F) {
          ndims = ndims + 1
          dims = c(dims, max(as.numeric(l)))
        }
      }
      if(ndims >5) print('Error: this function only supports arrays of <=5 dimensions')
      Chains[[i]] = array(NA, dim = c(nchains, dims))
      a = which(is.na(match(col_variables, variable_names[i]))==F)
      
      if(ndims == 1){
        Chains[[i]] = table[,a]
        colnames(Chains[[i]]) = NULL
      }
      
      if(ndims == 2){
        j = 1
        for(d2 in 1:dims[2]){
          for(d1 in 1:dims[1]){
            Chains[[i]][,d1,d2] = table[,a[j]]
            j = j+1
          }
        }
      }
      
      if(ndims == 3){
        j = 1
        for(d3 in 1:dims[3]){
          for(d2 in 1:dims[2]){
            for(d1 in 1:dims[1]){
              Chains[[i]][,d1,d2,d3] = table[,a[j]]
              j = j+1
            }
          }
        }
      }
      
      if(ndims == 4){
        j = 1
        for(d4 in 1:dims[4]){
          for(d3 in 1:dims[3]){
            for(d2 in 1:dims[2]){
              for(d1 in 1:dims[1]){
                Chains[[i]][,d1,d2,d3,d4] = table[,a[j]]
                j = j+1
              }
            }
          }
        }
      }
      
      if(ndims == 5){
        j = 1
        for(d5 in 1:dims[5]){
          for(d4 in 1:dims[4]){
            for(d3 in 1:dims[3]){
              for(d2 in 1:dims[2]){
                for(d1 in 1:dims[1]){
                  Chains[[i]][,d1,d2,d3,d4] = table[,a[j]]
                  j = j+1
                }
              }
            }
          }
        }
      }
    }
    names(Chains)[i] = variable_names[i]
  }
  return(Chains)
}
########################################################################################################################################

######################################################################################################################################
## Read best splits
######################################################################################################################################
k_smooth = -1
min_group_size = 30

potential_splits = readRDS(paste0('2_analysis_index/2_find_index_groups/potential_splits_weighting_202510_', k_smooth, '_min_group_size_', min_group_size, '.rds'))$best_nodes_names[[9]]
split = merge.groups(timed_tree = tree, 
                          metadata = dataset_with_nodes, 
                          initial_splits = potential_splits,
                          group_count_threshold = 5, group_freq_threshold = 0)

min_group_size = 30
potential_splits_STN = readRDS(paste0('2_analysis_index/2_find_index_groups/potential_splits_weighting_STN_202510_', k_smooth, '_min_group_size_', min_group_size, '.rds'))$best_nodes_names[[3]]
split_STN = merge.groups(timed_tree = tree_STN, 
                     metadata = dataset_with_nodes_STN, 
                     initial_splits = potential_splits_STN,
                     group_count_threshold = 5, group_freq_threshold = 0)
potential_splits_nonSTN = readRDS(paste0('2_analysis_index/2_find_index_groups/potential_splits_weighting_nonSTN_202510_', k_smooth, '_min_group_size_', min_group_size, '.rds'))$best_nodes_names[[8]]
split_nonSTN = merge.groups(timed_tree = tree_nonSTN, 
                         metadata = dataset_with_nodes_nonSTN, 
                         initial_splits = potential_splits_nonSTN,
                         group_count_threshold = 5, group_freq_threshold = 0)
min_group_size = 30
potential_splits_pMSM = readRDS(paste0('2_analysis_index/2_find_index_groups/potential_splits_weighting_pMSM_202510_', k_smooth, '_min_group_size_', min_group_size, '.rds'))$best_nodes_names[[8]]
split_pMSM = merge.groups(timed_tree = tree_pMSM, 
                         metadata = dataset_with_nodes_pMSM, 
                         initial_splits = potential_splits_pMSM,
                         group_count_threshold = 5, group_freq_threshold = 0)
potential_splits_nonpMSM = readRDS(paste0('2_analysis_index/2_find_index_groups/potential_splits_weighting_nonpMSM_202510_', k_smooth, '_min_group_size_', min_group_size, '.rds'))$best_nodes_names[[7]]
split_nonpMSM = merge.groups(timed_tree = tree_nonpMSM, 
                            metadata = dataset_with_nodes_nonpMSM, 
                            initial_splits = potential_splits_nonpMSM,
                            group_count_threshold = 5, group_freq_threshold = 0)
potential_splits_MSM_pandemic = readRDS(paste0('2_analysis_index/2_find_index_groups/potential_splits_weighting_MSM_pandemic_202510_', k_smooth, '_min_group_size_', min_group_size, '.rds'))$best_nodes_names[[5]]
split_MSM_pandemic = merge.groups(timed_tree = tree_MSM_pandemic, 
                                  metadata = dataset_with_nodes_MSM_pandemic, 
                                  initial_splits = potential_splits_MSM_pandemic,
                                  group_count_threshold = 5, group_freq_threshold = 0)
######################################################################################################################################

#######################################################################################################################################
## Assign colour to each group - ALL
#######################################################################################################################################
dataset_with_nodes$groups = as.factor(split$groups)
name_groups = levels(dataset_with_nodes$groups)
n_groups <- length(name_groups)

## Reorder labels by time of emergence
name_groups = levels(dataset_with_nodes$groups)
time_groups = NULL
for(i in 1:length(name_groups)){
  time_groups = c(time_groups, min(dataset_with_nodes$time[which(dataset_with_nodes$groups == name_groups[i] &
                                                                                 dataset_with_nodes$is.node == 'no')]))
}
levels(dataset_with_nodes$groups) = match(name_groups, order(time_groups, decreasing = T))
dataset_with_nodes$groups = as.numeric(as.character(dataset_with_nodes$groups))
dataset_with_nodes$groups = as.factor(dataset_with_nodes$groups)
## Update split
split$tip_and_nodes_groups = match(split$tip_and_nodes_groups, order(time_groups, decreasing = T))
names(split$tip_and_nodes_groups) = 1:length(split$tip_and_nodes_groups)
split$groups = as.factor(split$groups)
levels(split$groups) = match(name_groups, order(time_groups, decreasing = T))
split$groups = as.numeric(as.character(split$groups))

## Choose color palette
colors_groups = met.brewer(name="Cross", n=n_groups, type="continuous")

## Color each group
dataset_with_nodes$group_color = dataset_with_nodes$groups
levels(dataset_with_nodes$group_color) = colors_groups
dataset_with_nodes$group_color = as.character(dataset_with_nodes$group_color)
#######################################################################################################################################

#######################################################################################################################################
## Assign colour to each group - STN
#######################################################################################################################################
dataset_with_nodes_STN$groups = as.factor(split_STN$groups)
name_groups_STN = levels(dataset_with_nodes_STN$groups)
n_groups_STN <- length(name_groups_STN)

## Reorder labels by time of emergence
name_groups_STN = levels(dataset_with_nodes_STN$groups)
time_groups_STN = NULL
for(i in 1:length(name_groups_STN)){
  time_groups_STN = c(time_groups_STN, min(dataset_with_nodes_STN$time[which(dataset_with_nodes_STN$groups == name_groups_STN[i] &
                                                                             dataset_with_nodes_STN$is.node == 'no')]))
}
levels(dataset_with_nodes_STN$groups) = match(name_groups_STN, order(time_groups_STN, decreasing = T))
dataset_with_nodes_STN$groups = as.numeric(as.character(dataset_with_nodes_STN$groups))
dataset_with_nodes_STN$groups = as.factor(dataset_with_nodes_STN$groups)
## Update split
split_STN$tip_and_nodes_groups = match(split_STN$tip_and_nodes_groups, order(time_groups_STN, decreasing = T))
names(split_STN$tip_and_nodes_groups) = 1:length(split_STN$tip_and_nodes_groups)
split_STN$groups = as.factor(split_STN$groups)
levels(split_STN$groups) = match(name_groups_STN, order(time_groups_STN, decreasing = T))
split_STN$groups = as.numeric(as.character(split_STN$groups))

## Choose color palette
colors_groups_STN = met.brewer(name="Cross", n=n_groups_STN, type="continuous")

## Color each group
dataset_with_nodes_STN$group_color = dataset_with_nodes_STN$groups
levels(dataset_with_nodes_STN$group_color) = colors_groups_STN
dataset_with_nodes_STN$group_color = as.character(dataset_with_nodes_STN$group_color)
#######################################################################################################################################

#######################################################################################################################################
## Assign colour to each group - non STN
#######################################################################################################################################
dataset_with_nodes_nonSTN$groups = as.factor(split_nonSTN$groups)
name_groups_nonSTN = levels(dataset_with_nodes_nonSTN$groups)
n_groups_nonSTN <- length(name_groups_nonSTN)

## Reorder labels by time of emergence
name_groups_nonSTN = levels(dataset_with_nodes_nonSTN$groups)
time_groups_nonSTN = NULL
for(i in 1:length(name_groups_nonSTN)){
  time_groups_nonSTN = c(time_groups_nonSTN, min(dataset_with_nodes_nonSTN$time[which(dataset_with_nodes_nonSTN$groups == name_groups_nonSTN[i] &
                                                                               dataset_with_nodes_nonSTN$is.node == 'no')]))
}
levels(dataset_with_nodes_nonSTN$groups) = match(name_groups_nonSTN, order(time_groups_nonSTN, decreasing = T))
dataset_with_nodes_nonSTN$groups = as.numeric(as.character(dataset_with_nodes_nonSTN$groups))
dataset_with_nodes_nonSTN$groups = as.factor(dataset_with_nodes_nonSTN$groups)
## Update split
split_nonSTN$tip_and_nodes_groups = match(split_nonSTN$tip_and_nodes_groups, order(time_groups_nonSTN, decreasing = T))
names(split_nonSTN$tip_and_nodes_groups) = 1:length(split_nonSTN$tip_and_nodes_groups)
split_nonSTN$groups = as.factor(split_nonSTN$groups)
levels(split_nonSTN$groups) = match(name_groups_nonSTN, order(time_groups_nonSTN, decreasing = T))
split_nonSTN$groups = as.numeric(as.character(split_nonSTN$groups))

## Choose color palette
colors_groups_nonSTN = met.brewer(name="Cross", n=n_groups_nonSTN, type="continuous")

## Color each group
dataset_with_nodes_nonSTN$group_color = dataset_with_nodes_nonSTN$groups
levels(dataset_with_nodes_nonSTN$group_color) = colors_groups_nonSTN
dataset_with_nodes_nonSTN$group_color = as.character(dataset_with_nodes_nonSTN$group_color)
#######################################################################################################################################

#######################################################################################################################################
## Assign colour to each group - pMSM
#######################################################################################################################################
dataset_with_nodes_pMSM$groups = as.factor(split_pMSM$groups)
name_groups_pMSM = levels(dataset_with_nodes_pMSM$groups)
n_groups_pMSM <- length(name_groups_pMSM)

## Reorder labels by time of emergence
name_groups_pMSM = levels(dataset_with_nodes_pMSM$groups)
time_groups_pMSM = NULL
for(i in 1:length(name_groups_pMSM)){
  time_groups_pMSM = c(time_groups_pMSM, min(dataset_with_nodes_pMSM$time[which(dataset_with_nodes_pMSM$groups == name_groups_pMSM[i] &
                                                                               dataset_with_nodes_pMSM$is.node == 'no')]))
}
levels(dataset_with_nodes_pMSM$groups) = match(name_groups_pMSM, order(time_groups_pMSM, decreasing = T))
dataset_with_nodes_pMSM$groups = as.numeric(as.character(dataset_with_nodes_pMSM$groups))
dataset_with_nodes_pMSM$groups = as.factor(dataset_with_nodes_pMSM$groups)
## Update split
split_pMSM$tip_and_nodes_groups = match(split_pMSM$tip_and_nodes_groups, order(time_groups_pMSM, decreasing = T))
names(split_pMSM$tip_and_nodes_groups) = 1:length(split_pMSM$tip_and_nodes_groups)
split_pMSM$groups = as.factor(split_pMSM$groups)
levels(split_pMSM$groups) = match(name_groups_pMSM, order(time_groups_pMSM, decreasing = T))
split_pMSM$groups = as.numeric(as.character(split_pMSM$groups))

## Choose color palette
colors_groups_pMSM = met.brewer(name="Cross", n=n_groups_pMSM, type="continuous")

## Color each group
dataset_with_nodes_pMSM$group_color = dataset_with_nodes_pMSM$groups
levels(dataset_with_nodes_pMSM$group_color) = colors_groups_pMSM
dataset_with_nodes_pMSM$group_color = as.character(dataset_with_nodes_pMSM$group_color)
#######################################################################################################################################

#######################################################################################################################################
## Assign colour to each group - non pMSM
#######################################################################################################################################
dataset_with_nodes_nonpMSM$groups = as.factor(split_nonpMSM$groups)
name_groups_nonpMSM = levels(dataset_with_nodes_nonpMSM$groups)
n_groups_nonpMSM <- length(name_groups_nonpMSM)

## Reorder labels by time of emergence
name_groups_nonpMSM = levels(dataset_with_nodes_nonpMSM$groups)
time_groups_nonpMSM = NULL
for(i in 1:length(name_groups_nonpMSM)){
  time_groups_nonpMSM = c(time_groups_nonpMSM, min(dataset_with_nodes_nonpMSM$time[which(dataset_with_nodes_nonpMSM$groups == name_groups_nonpMSM[i] &
                                                                                        dataset_with_nodes_nonpMSM$is.node == 'no')]))
}
levels(dataset_with_nodes_nonpMSM$groups) = match(name_groups_nonpMSM, order(time_groups_nonpMSM, decreasing = T))
dataset_with_nodes_nonpMSM$groups = as.numeric(as.character(dataset_with_nodes_nonpMSM$groups))
dataset_with_nodes_nonpMSM$groups = as.factor(dataset_with_nodes_nonpMSM$groups)
## Update split
split_nonpMSM$tip_and_nodes_groups = match(split_nonpMSM$tip_and_nodes_groups, order(time_groups_nonpMSM, decreasing = T))
names(split_nonpMSM$tip_and_nodes_groups) = 1:length(split_nonpMSM$tip_and_nodes_groups)
split_nonpMSM$groups = as.factor(split_nonpMSM$groups)
levels(split_nonpMSM$groups) = match(name_groups_nonpMSM, order(time_groups_nonpMSM, decreasing = T))
split_nonpMSM$groups = as.numeric(as.character(split_nonpMSM$groups))

## Choose color palette
colors_groups_nonpMSM = met.brewer(name="Cross", n=n_groups_nonpMSM, type="continuous")

## Color each group
dataset_with_nodes_nonpMSM$group_color = dataset_with_nodes_nonpMSM$groups
levels(dataset_with_nodes_nonpMSM$group_color) = colors_groups_nonpMSM
dataset_with_nodes_nonpMSM$group_color = as.character(dataset_with_nodes_nonpMSM$group_color)
#######################################################################################################################################

#######################################################################################################################################
## Assign colour to each group - MSM pandemic
#######################################################################################################################################
dataset_with_nodes_MSM_pandemic$groups = as.factor(split_MSM_pandemic$groups)
name_groups_MSM_pandemic = levels(dataset_with_nodes_MSM_pandemic$groups)
n_groups_MSM_pandemic <- length(name_groups_MSM_pandemic)

## Reorder labels by time of emergence
name_groups_MSM_pandemic = levels(dataset_with_nodes_MSM_pandemic$groups)
time_groups_MSM_pandemic = NULL
for(i in 1:length(name_groups_MSM_pandemic)){
  time_groups_MSM_pandemic = c(time_groups_MSM_pandemic, min(dataset_with_nodes_MSM_pandemic$time[which(dataset_with_nodes_MSM_pandemic$groups == name_groups_MSM_pandemic[i] &
                                                                                        dataset_with_nodes_MSM_pandemic$is.node == 'no')]))
}
levels(dataset_with_nodes_MSM_pandemic$groups) = match(name_groups_MSM_pandemic, order(time_groups_MSM_pandemic, decreasing = T))
dataset_with_nodes_MSM_pandemic$groups = as.numeric(as.character(dataset_with_nodes_MSM_pandemic$groups))
dataset_with_nodes_MSM_pandemic$groups = as.factor(dataset_with_nodes_MSM_pandemic$groups)
## Update split
split_MSM_pandemic$tip_and_nodes_groups = match(split_MSM_pandemic$tip_and_nodes_groups, order(time_groups_MSM_pandemic, decreasing = T))
names(split_MSM_pandemic$tip_and_nodes_groups) = 1:length(split_MSM_pandemic$tip_and_nodes_groups)
split_MSM_pandemic$groups = as.factor(split_MSM_pandemic$groups)
levels(split_MSM_pandemic$groups) = match(name_groups_MSM_pandemic, order(time_groups_MSM_pandemic, decreasing = T))
split_MSM_pandemic$groups = as.numeric(as.character(split_MSM_pandemic$groups))

## Choose color palette
colors_groups_MSM_pandemic = rev(met.brewer(name="Java", n=n_groups_MSM_pandemic, type="continuous"))

## Color each group
dataset_with_nodes_MSM_pandemic$group_color = dataset_with_nodes_MSM_pandemic$groups
levels(dataset_with_nodes_MSM_pandemic$group_color) = colors_groups_MSM_pandemic
dataset_with_nodes_MSM_pandemic$group_color = as.character(dataset_with_nodes_MSM_pandemic$group_color)
#######################################################################################################################################


########################################################################################################################################
## Reconstruction STN and non-STN groups on full phylogeny
########################################################################################################################################
dataset_with_nodes_reconstructed = dataset_with_nodes
dataset_with_nodes_reconstructed$group_color_STN = dataset_with_nodes_reconstructed$groups_STN = dataset_with_nodes_reconstructed$index_STN = NA
dataset_with_nodes_reconstructed$group_color_nonSTN = dataset_with_nodes_reconstructed$groups_nonSTN = dataset_with_nodes_reconstructed$index_nonSTN = NA
dataset_with_nodes_reconstructed$group_color_pMSM = dataset_with_nodes_reconstructed$groups_pMSM = dataset_with_nodes_reconstructed$index_pMSM = NA
dataset_with_nodes_reconstructed$group_color_nonpMSM = dataset_with_nodes_reconstructed$groups_nonpMSM = dataset_with_nodes_reconstructed$index_nonpMSM = NA

idx = match(dataset_with_nodes$name_seq2, dataset_with_nodes_STN$name_seq2)
dataset_with_nodes_reconstructed$index_STN[which(!is.na(idx))] = dataset_with_nodes_STN$index[idx[which(!is.na(idx))]]
dataset_with_nodes_reconstructed$groups_STN[which(!is.na(idx))] = dataset_with_nodes_STN$groups[idx[which(!is.na(idx))]]
dataset_with_nodes_reconstructed$group_color_STN[which(!is.na(idx))] = dataset_with_nodes_STN$group_color[idx[which(!is.na(idx))]]

idx = match(dataset_with_nodes$name_seq2, dataset_with_nodes_nonSTN$name_seq2)
dataset_with_nodes_reconstructed$index_nonSTN[which(!is.na(idx))] = dataset_with_nodes_nonSTN$index[idx[which(!is.na(idx))]]
dataset_with_nodes_reconstructed$groups_nonSTN[which(!is.na(idx))] = dataset_with_nodes_nonSTN$groups[idx[which(!is.na(idx))]]
dataset_with_nodes_reconstructed$group_color_nonSTN[which(!is.na(idx))] = dataset_with_nodes_nonSTN$group_color[idx[which(!is.na(idx))]]

idx = match(dataset_with_nodes$name_seq2, dataset_with_nodes_pMSM$name_seq2)
dataset_with_nodes_reconstructed$index_pMSM[which(!is.na(idx))] = dataset_with_nodes_pMSM$index[idx[which(!is.na(idx))]]
dataset_with_nodes_reconstructed$groups_pMSM[which(!is.na(idx))] = dataset_with_nodes_pMSM$groups[idx[which(!is.na(idx))]]
dataset_with_nodes_reconstructed$group_color_pMSM[which(!is.na(idx))] = dataset_with_nodes_pMSM$group_color[idx[which(!is.na(idx))]]

idx = match(dataset_with_nodes$name_seq2, dataset_with_nodes_nonpMSM$name_seq2)
dataset_with_nodes_reconstructed$index_nonpMSM[which(!is.na(idx))] = dataset_with_nodes_nonpMSM$index[idx[which(!is.na(idx))]]
dataset_with_nodes_reconstructed$groups_nonpMSM[which(!is.na(idx))] = dataset_with_nodes_nonpMSM$groups[idx[which(!is.na(idx))]]
dataset_with_nodes_reconstructed$group_color_nonpMSM[which(!is.na(idx))] = dataset_with_nodes_nonpMSM$group_color[idx[which(!is.na(idx))]]

## nonSTN
snp_data = dataset_with_nodes_reconstructed$groups_nonSTN[which(dataset_with_nodes_reconstructed$is.node == 'no')]
snp_data = as.numeric(as.factor(snp_data))
tree_tmp = tree
tree_tmp$node.label = NULL
rec = ace(snp_data, tree_tmp, type = 'discrete')
rec_nodes = apply(rec$lik.anc, MAR = 1, function(x)which.max(x))
rec_all = c(as.numeric(as.character(snp_data)), rec_nodes) ## List of all states: first all tips, then all nodes
dataset_with_nodes_reconstructed$groups_nonSTN = rec_all

## STN
snp_data = dataset_with_nodes_reconstructed$groups_STN[which(dataset_with_nodes_reconstructed$is.node == 'no')]
snp_data = as.numeric(as.factor(snp_data))
tree_tmp = tree
tree_tmp$node.label = NULL
rec = ace(snp_data, tree_tmp, type = 'discrete')
rec_nodes = apply(rec$lik.anc, MAR = 1, function(x)which.max(x))
rec_all = c(as.numeric(as.character(snp_data)), rec_nodes) ## List of all states: first all tips, then all nodes
dataset_with_nodes_reconstructed$groups_STN = rec_all

## nonpMSM
snp_data = dataset_with_nodes_reconstructed$groups_nonpMSM[which(dataset_with_nodes_reconstructed$is.node == 'no')]
snp_data = as.numeric(as.factor(snp_data))
tree_tmp = tree
tree_tmp$node.label = NULL
rec = ace(snp_data, tree_tmp, type = 'discrete')
rec_nodes = apply(rec$lik.anc, MAR = 1, function(x)which.max(x))
rec_all = c(as.numeric(as.character(snp_data)), rec_nodes) ## List of all states: first all tips, then all nodes
dataset_with_nodes_reconstructed$groups_nonpMSM = rec_all

## pMSM
snp_data = dataset_with_nodes_reconstructed$groups_pMSM[which(dataset_with_nodes_reconstructed$is.node == 'no')]
snp_data = as.numeric(as.factor(snp_data))
tree_tmp = tree
tree_tmp$node.label = NULL
rec = ace(snp_data, tree_tmp, type = 'discrete')
rec_nodes = apply(rec$lik.anc, MAR = 1, function(x)which.max(x))
rec_all = c(as.numeric(as.character(snp_data)), rec_nodes) ## List of all states: first all tips, then all nodes
dataset_with_nodes_reconstructed$groups_pMSM = rec_all
########################################################################################################################################

#######################################################################################################################################
## Save the detection
#######################################################################################################################################
save.image('2_analysis_index/2_find_index_groups/Lineages_detected_20251007.Rdata')
#######################################################################################################################################






