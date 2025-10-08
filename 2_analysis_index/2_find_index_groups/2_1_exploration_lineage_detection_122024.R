########################################################################################################################################
## Detect lineages, input: Initial_index_computation_and_parameters_20251007.Rdata
########################################################################################################################################
## Packages used
library(ape)
library(phangorn)
library(phytools)
library(coda)
library(thd)
library(snow)
library(doParallel)

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
## Parameters
######################################################################################################################################
## Parameter in common to all algos
time_window_initial = 2030;
time_window_increment = 100;
p_value_smooth = 0.05
weight_by_time = 1
k_smooth = -1
plot_screening = F
min_descendants_per_tested_node = 30
min_group_size = 30
weighting_transformation = c('inv_sqrt')

parallelize_code = T
number_cores = 3

max_stepwise_deviance_explained_threshold = 0
max_groups_found = 20
stepwise_AIC_threshold = 0

keep_track = T
######################################################################################################################################

######################################################################################################################################
## Detection lineages - full phylogeny
######################################################################################################################################
start_time = Sys.time()
potential_splits_weighting = find.groups.by.index.dynamics(timed_tree = tree,
                                                           metadata = dataset_with_nodes,
                                                           node_support = NULL,
                                                           threshold_node_support = NA, 
                                                           time_window_initial = time_window_initial,
                                                           time_window_increment = time_window_increment, 
                                                           min_descendants_per_tested_node = min_descendants_per_tested_node,
                                                           min_group_size = min_group_size,
                                                           p_value_smooth = p_value_smooth, 
                                                           stepwise_deviance_explained_threshold = max_stepwise_deviance_explained_threshold, 
                                                           stepwise_AIC_threshold = stepwise_AIC_threshold,
                                                           weight_by_time = weight_by_time,
                                                           weighting_transformation = weighting_transformation,
                                                           k_smooth = k_smooth,
                                                           plot_screening = plot_screening,
                                                           parallelize_code = parallelize_code,
                                                           number_cores = number_cores, 
                                                           max_groups_found = max_groups_found, 
                                                           keep_track = keep_track)
end_time = Sys.time()
print(end_time - start_time)
saveRDS(potential_splits_weighting, paste0('2_analysis_index/2_find_index_groups/potential_splits_weighting_202510_', k_smooth, '_min_group_size_', min_group_size, '.rds'))
######################################################################################################################################



######################################################################################################################################
### STN and non-STN, separatly ####
######################################################################################################################################

######################################################################################################################################
## Parameters
######################################################################################################################################
## Parameter in common to all algos
time_window_initial = 2030;
time_window_increment = 100;
p_value_smooth = 0.05
weight_by_time = 1
k_smooth = -1
plot_screening = F
min_descendants_per_tested_node = 30
min_group_size = 30
weighting_transformation = c('inv_sqrt')

parallelize_code = T
number_cores = 3

max_stepwise_deviance_explained_threshold = 0
max_groups_found = 20
stepwise_AIC_threshold = 0

keep_track = T
######################################################################################################################################

######################################################################################################################################
## Detection lineages - STN
######################################################################################################################################
start_time = Sys.time()
potential_splits_weighting = find.groups.by.index.dynamics(timed_tree = tree_STN,
                                                               metadata = dataset_with_nodes_STN,
                                                               node_support = NULL,
                                                               threshold_node_support = 0.5,
                                                               time_window_initial,
                                                               time_window_increment,
                                                               min_descendants_per_tested_node,
                                                               min_group_size,
                                                               p_value_smooth,
                                                               stepwise_deviance_explained_threshold = max_stepwise_deviance_explained_threshold,
                                                               stepwise_AIC_threshold,
                                                               weight_by_time,
                                                               weighting_transformation,
                                                               k_smooth,
                                                               plot_screening,
                                                               parallelize_code,
                                                               number_cores = number_cores, 
                                                               max_groups_found = max_groups_found, 
                                                               keep_track = keep_track)
end_time = Sys.time()
print(end_time - start_time)
saveRDS(potential_splits_weighting, paste0('2_analysis_index/2_find_index_groups/potential_splits_weighting_STN_202510_', k_smooth, '_min_group_size_', min_group_size, '.rds'))
######################################################################################################################################

######################################################################################################################################
## Detection lineages - non STN
######################################################################################################################################
start_time = Sys.time()
potential_splits_weighting = find.groups.by.index.dynamics(timed_tree = tree_nonSTN,
                                                               metadata = dataset_with_nodes_nonSTN,
                                                               node_support = NULL,
                                                               threshold_node_support = 0.5,
                                                               time_window_initial,
                                                               time_window_increment,
                                                               min_descendants_per_tested_node,
                                                               min_group_size,
                                                               p_value_smooth,
                                                               stepwise_deviance_explained_threshold = max_stepwise_deviance_explained_threshold,
                                                               stepwise_AIC_threshold,
                                                               weight_by_time,
                                                               weighting_transformation,
                                                               k_smooth,
                                                               plot_screening,
                                                               parallelize_code,
                                                               number_cores = number_cores, 
                                                               max_groups_found = max_groups_found, 
                                                               keep_track = keep_track)
end_time = Sys.time()
print(end_time - start_time)
saveRDS(potential_splits_weighting, paste0('2_analysis_index/2_find_index_groups/potential_splits_weighting_nonSTN_202510_', k_smooth, '_min_group_size_', min_group_size, '.rds'))
######################################################################################################################################




######################################################################################################################################
### pMSM and non-pMSM ####
######################################################################################################################################

######################################################################################################################################
## Parameters
######################################################################################################################################
## Parameter in common to all algos
time_window_initial = 2030;
time_window_increment = 100;
p_value_smooth = 0.05
weight_by_time = 1
k_smooth = -1
plot_screening = F
min_descendants_per_tested_node = 30
min_group_size = 30
weighting_transformation = c('inv_sqrt')

parallelize_code = T
number_cores = 3

max_stepwise_deviance_explained_threshold = 0
max_groups_found = 20
stepwise_AIC_threshold = 0

keep_track = T
######################################################################################################################################

######################################################################################################################################
## Detection lineages - pMSM
######################################################################################################################################
start_time = Sys.time()
potential_splits_weighting = find.groups.by.index.dynamics(timed_tree = tree_pMSM,
                                                               metadata = dataset_with_nodes_pMSM,
                                                               node_support = NULL,
                                                               threshold_node_support = 0.5,
                                                               time_window_initial,
                                                               time_window_increment,
                                                               min_descendants_per_tested_node,
                                                               min_group_size,
                                                               p_value_smooth,
                                                               stepwise_deviance_explained_threshold = max_stepwise_deviance_explained_threshold,
                                                               stepwise_AIC_threshold,
                                                               weight_by_time,
                                                               weighting_transformation,
                                                               k_smooth,
                                                               plot_screening,
                                                               parallelize_code,
                                                               number_cores = number_cores, 
                                                               max_groups_found = max_groups_found, 
                                                               keep_track = keep_track)
end_time = Sys.time()
print(end_time - start_time)
saveRDS(potential_splits_weighting, paste0('2_analysis_index/2_find_index_groups/potential_splits_weighting_pMSM_202510_', k_smooth, '_min_group_size_', min_group_size, '.rds'))
######################################################################################################################################

######################################################################################################################################
## Detection lineages - non pMSM
######################################################################################################################################
start_time = Sys.time()
potential_splits_weighting = find.groups.by.index.dynamics(timed_tree = tree_nonpMSM,
                                                               metadata = dataset_with_nodes_nonpMSM,
                                                               node_support = NULL,
                                                               threshold_node_support = 0.5,
                                                               time_window_initial,
                                                               time_window_increment,
                                                               min_descendants_per_tested_node,
                                                               min_group_size,
                                                               p_value_smooth,
                                                               stepwise_deviance_explained_threshold = max_stepwise_deviance_explained_threshold,
                                                               stepwise_AIC_threshold,
                                                               weight_by_time,
                                                               weighting_transformation,
                                                               k_smooth,
                                                               plot_screening,
                                                               parallelize_code,
                                                               number_cores = number_cores, 
                                                               max_groups_found = max_groups_found, 
                                                               keep_track = keep_track)
end_time = Sys.time()
print(end_time - start_time)
saveRDS(potential_splits_weighting, paste0('2_analysis_index/2_find_index_groups/potential_splits_weighting_nonpMSM_202510_', k_smooth, '_min_group_size_', min_group_size, '.rds'))
######################################################################################################################################

######################################################################################################################################
## Detection lineages - MSM pandemic
######################################################################################################################################
start_time = Sys.time()
potential_splits_weighting = find.groups.by.index.dynamics(timed_tree = tree_MSM_pandemic,
                                                               metadata = dataset_with_nodes_MSM_pandemic,
                                                               node_support = NULL,
                                                               threshold_node_support = 0.5,
                                                               time_window_initial,
                                                               time_window_increment,
                                                               min_descendants_per_tested_node,
                                                               min_group_size,
                                                               p_value_smooth,
                                                               stepwise_deviance_explained_threshold = max_stepwise_deviance_explained_threshold,
                                                               stepwise_AIC_threshold,
                                                               weight_by_time,
                                                               weighting_transformation,
                                                               k_smooth,
                                                               plot_screening,
                                                               parallelize_code,
                                                               number_cores = number_cores, 
                                                               max_groups_found = max_groups_found, 
                                                               keep_track = keep_track, log_y = F)
end_time = Sys.time()
print(end_time - start_time)
saveRDS(potential_splits_weighting, paste0('2_analysis_index/2_find_index_groups/potential_splits_weighting_MSM_pandemic_202510_', k_smooth, '_min_group_size_', min_group_size, '.rds'))
######################################################################################################################################eqw rt
 
