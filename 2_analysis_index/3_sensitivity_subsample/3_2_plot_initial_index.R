########################################################################################################################################
## Plot index: input: Initial_index_computation_and_parameters_20260313.Rdata
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

load('2_analysis_index/3_sensitivity_subsample/Initial_index_computation_and_parameters_subsampled_20260313.Rdata')

########################################################################################################################################
## Plot tree & index below, with colors from full genotypes
########################################################################################################################################
pdf('2_analysis_index/3_sensitivity_subsample/Index_SSonnei_UK_subsampled_20260313.pdf', width = 6/2.54, height = 6/2.54)
par(mfcol = c(2,1), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)

min_year = 1980
max_year = 2020
max_index = 0.25

## World
tree = tree
root_height = root_height
dataset_with_nodes = dataset_with_nodes
## Plot tree
tree_to_plot = ladderize(tree, right = F)
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$Genotype_color, cex = 0.15)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 10)-root_height,
             lab_axis = seq(min_year, max_year, 10), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$Genotype_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15, #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$Genotype_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, #0.15
       pch = 16)
axis(1, at = seq(min_year, max_year, 10), labels = seq(min_year, max_year, 10), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. sonnei", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

idx = which(!duplicated(dataset_with_nodes$Genotype))
legend('topleft', legend = dataset_with_nodes$Genotype[idx], col = dataset_with_nodes$Genotype_color[idx], pch = 16,
       cex = 0.15, ncol = 3, bty = 'n')

dev.off()
########################################################################################################################################

########################################################################################################################################
## Plot same dynamics for pMSM and non-pMSM, with colors from full genotypes
########################################################################################################################################
pdf('2_analysis_index/3_sensitivity_subsample/Index_SSonnei_UK_subsampled_20260313_non_pMSM.pdf', width = 10/2.54, height = 5/2.54)
par(mfcol = c(2,2), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)

min_year = 1980
max_year = 2020
max_index = 0.6

## pMSM
tree = tree_pMSM
dataset_with_nodes = dataset_with_nodes_pMSM
root_height = min(dataset_with_nodes$time)
## Plot tree
tree_to_plot = tree
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$Genotype_color, cex = 0.15)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 2)-root_height,
             lab_axis = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$Genotype_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15, #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$Genotype_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, #0.15
       pch = 16)
axis(1, at = seq(min_year, max_year, 2), labels = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. sonnei - pMSM", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

# idx = which(!duplicated(dataset_with_nodes$Genotype))
# legend('topleft', legend = dataset_with_nodes$Genotype[idx], col = dataset_with_nodes$Genotype_color[idx], pch = 16)

## Non pMSM
tree = tree_nonpMSM
dataset_with_nodes = dataset_with_nodes_nonpMSM
root_height = min(dataset_with_nodes$time)
## Plot tree
tree_to_plot = tree
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$Genotype_color, cex = 0.15)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 2)-root_height,
             lab_axis = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$Genotype_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15, #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$Genotype_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, #0.15
       pch = 16)
axis(1, at = seq(min_year, max_year, 2), labels = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. sonnei - non pMSM", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

dev.off()
########################################################################################################################################

########################################################################################################################################
## Plot same dynamics for travel pMSM, with colors from full genotypes
########################################################################################################################################
pdf('2_analysis_index/3_sensitivity_subsample/Index_SSonnei_UK_subsampled_20260313_travel.pdf', width = 10/2.54, height = 5/2.54)
par(mfcol = c(2,2), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)

min_year = 1980
max_year = 2020
max_index = 0.6

## travel pMSM
tree = tree_travelpMSM
dataset_with_nodes = dataset_with_nodes_travelpMSM
root_height = min(dataset_with_nodes$time)
## Plot tree
tree_to_plot = tree
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey',
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$Genotype_color, cex = 0.15)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 2)-root_height,
             lab_axis = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')],
     dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'yes')],
     col = adjustcolor(dataset_with_nodes$Genotype_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15, #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index),
     # main = paste0('World'),
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')],
       dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'no')],
       col = adjustcolor(dataset_with_nodes$Genotype_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, #0.15
       pch = 16)
axis(1, at = seq(min_year, max_year, 2), labels = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. sonnei - travel pMSM", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

dev.off()
########################################################################################################################################


#### Genotype simplified 
load('2_analysis_index/3_sensitivity_subsample/Initial_index_computation_and_parameters_subsampled_20260313.Rdata')

########################################################################################################################################
## Plot tree & index below, with colors from simplified genotypes
########################################################################################################################################
pdf('2_analysis_index/3_sensitivity_subsample/Index_SSonnei_UK_subsampled_20260313_simplified.pdf', width = 6/2.54, height = 6/2.54)
par(mfcol = c(2,1), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)

min_year = 1980
max_year = 2020
max_index = 0.25

## World
tree = tree
root_height = root_height
dataset_with_nodes = dataset_with_nodes
## Plot tree
tree_to_plot = ladderize(tree, right = F)
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tmp = dataset_with_nodes$Genotype_simplified_color
tmp[which(tmp != 'grey')] = NA
tiplabels(pch = 16, col = tmp, cex = 0.15)
tmp = dataset_with_nodes$Genotype_simplified_color
tmp[which(tmp == 'grey')] = NA
tiplabels(pch = 16, col = tmp, cex = 0.15)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 10)-root_height,
             lab_axis = seq(min_year, max_year, 10), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$Genotype_simplified_color == 'grey')], 
     dataset_with_nodes$index[which(dataset_with_nodes$Genotype_simplified_color == 'grey')], 
     col = adjustcolor(dataset_with_nodes$Genotype_simplified_color[which(dataset_with_nodes$Genotype_simplified_color == 'grey')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15, #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$Genotype_simplified_color != 'grey')], 
       dataset_with_nodes$index[which(dataset_with_nodes$Genotype_simplified_color != 'grey')], 
       col = adjustcolor(dataset_with_nodes$Genotype_simplified_color[which(dataset_with_nodes$Genotype_simplified_color != 'grey')], alpha.f = 1),
       cex = 0.15, #0.15
       pch = 16)
axis(1, at = seq(min_year, max_year, 10), labels = seq(min_year, max_year, 10), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. sonnei", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

idx = which(!duplicated(dataset_with_nodes$Genotype_simplified))
legend('topleft', legend = dataset_with_nodes$Genotype_simplified[idx], col = dataset_with_nodes$Genotype_simplified_color[idx], pch = 16,
       cex = 0.15, ncol = 3, bty = 'n')

dev.off()
########################################################################################################################################

########################################################################################################################################
## Plot same dynamics for pMSM and non-pMSM
########################################################################################################################################
pdf('2_analysis_index/3_sensitivity_subsample/Index_SSonnei_UK_subsampled_20260313_non_pMSM_simplified.pdf', width = 10/2.54, height = 5/2.54)
par(mfcol = c(2,2), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)

min_year = 1980
max_year = 2020
max_index = 0.6

## pMSM
tree = tree_pMSM
dataset_with_nodes = dataset_with_nodes_pMSM
root_height = min(dataset_with_nodes$time)
## Plot tree
tree_to_plot = tree
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$Genotype_simplified_color, cex = 0.15)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 2)-root_height,
             lab_axis = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$Genotype_simplified_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15, #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$Genotype_simplified_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, #0.15
       pch = 16)
axis(1, at = seq(min_year, max_year, 2), labels = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. sonnei - pMSM", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

# idx = which(!duplicated(dataset_with_nodes$Genotype))
# legend('topleft', legend = dataset_with_nodes$Genotype[idx], col = dataset_with_nodes$Genotype_simplified_color[idx], pch = 16)

## Non pMSM
tree = tree_nonpMSM
dataset_with_nodes = dataset_with_nodes_nonpMSM
root_height = min(dataset_with_nodes$time)
## Plot tree
tree_to_plot = tree
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$Genotype_simplified_color, cex = 0.15)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 2)-root_height,
             lab_axis = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$Genotype_simplified_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15, #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$Genotype_simplified_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, #0.15
       pch = 16)
axis(1, at = seq(min_year, max_year, 2), labels = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. sonnei - non pMSM", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

dev.off()
########################################################################################################################################

########################################################################################################################################
##  Plot same dynamics for MSM pandemic, with colors from full genotypes
########################################################################################################################################
pdf('2_analysis_index/3_sensitivity_subsample/Index_SSonnei_UK_subsampled_20260313_MSM_pandemic.pdf', width = 6/2.54, height = 6/2.54)
par(mfcol = c(2,1), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)

min_year = 2014
max_year = 2022
max_index = 0.6

## World
tree = tree_MSM_pandemic
root_height_plot = root_height_MSM_pandemic
dataset_with_nodes = dataset_with_nodes_MSM_pandemic
## Plot tree
tree_to_plot = ladderize(tree, right = F)
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height_plot)
tiplabels(pch = 16, col = dataset_with_nodes$Genotype_color, cex = 0.15)
axisPhylo_NL(side = 1, root.time = root_height_plot, backward = F,
             at_axis = seq(min_year, max_year, 2)-root_height_plot,
             lab_axis = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$Genotype_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15, #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$Genotype_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, #0.15
       pch = 16)
axis(1, at = seq(min_year, max_year, 2), labels = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. sonnei", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

idx = which(!duplicated(dataset_with_nodes$Genotype))
legend('topleft', legend = dataset_with_nodes$Genotype[idx], col = dataset_with_nodes$Genotype_color[idx], pch = 16,
       cex = 0.15, ncol = 3, bty = 'n')

dev.off()
########################################################################################################################################


########################################################################################################################################
#### Coloured by groups detected in the full analysis
########################################################################################################################################

###################################################
## Map groups to datasets
###################################################
load('2_analysis_index/3_sensitivity_subsample/Initial_index_computation_and_parameters_subsampled_20260313.Rdata')
data = read.csv('dataset_with_nodes_20251218_all.csv')
dataset_with_nodes$groups = data$groups[match(dataset_with_nodes$name_seq, data$name_seq)]

## Reconstruct groups on nodes
snp_data = dataset_with_nodes$groups[which(dataset_with_nodes$is.node == 'no')]
snp_data = as.numeric(as.factor(snp_data))
tree_tmp = tree
tree_tmp$node.label = NULL
rec = ace(snp_data, tree_tmp, type = 'discrete')
rec_nodes = apply(rec$lik.anc, MAR = 1, function(x)which.max(x))
rec_nodes_score = apply(rec$lik.anc, MAR = 1, function(x)max(x))
rec_nodes[which(rec_nodes_score < 0.75)] = NA
rec_all = c(as.numeric(as.character(snp_data)), rec_nodes) ## List of all states: first all tips, then all nodes
dataset_with_nodes$groups = rec_all

colours = data$group_color[match(sort(unique(data$groups), decreasing = F), data$groups)]
dataset_with_nodes$group_color = dataset_with_nodes$groups
dataset_with_nodes$group_color = factor(dataset_with_nodes$group_color)
levels(dataset_with_nodes$group_color) = colours
dataset_with_nodes$group_color = as.character(dataset_with_nodes$group_color)

## Reconstruct for pMSM and non-pMSM
dataset_with_nodes_reconstructed = dataset_with_nodes
dataset_with_nodes_reconstructed$group_color_pMSM = dataset_with_nodes_reconstructed$groups_pMSM = dataset_with_nodes_reconstructed$index_pMSM = NA
dataset_with_nodes_reconstructed$group_color_nonpMSM = dataset_with_nodes_reconstructed$groups_nonpMSM = dataset_with_nodes_reconstructed$index_nonpMSM = NA

idx = match(dataset_with_nodes$name_seq2, dataset_with_nodes_pMSM$name_seq2)
dataset_with_nodes_reconstructed$index_pMSM[which(!is.na(idx))] = dataset_with_nodes_pMSM$index[idx[which(!is.na(idx))]]
dataset_with_nodes_reconstructed$groups_pMSM[which(!is.na(idx))] = dataset_with_nodes_reconstructed$groups[which(!is.na(idx))]
dataset_with_nodes_reconstructed$group_color_pMSM[which(!is.na(idx))] = dataset_with_nodes_reconstructed$group_color[which(!is.na(idx))]

idx = match(dataset_with_nodes$name_seq2, dataset_with_nodes_nonpMSM$name_seq2)
dataset_with_nodes_reconstructed$index_nonpMSM[which(!is.na(idx))] = dataset_with_nodes_nonpMSM$index[idx[which(!is.na(idx))]]
dataset_with_nodes_reconstructed$groups_nonpMSM[which(!is.na(idx))] = dataset_with_nodes_reconstructed$groups[which(!is.na(idx))]
dataset_with_nodes_reconstructed$group_color_nonpMSM[which(!is.na(idx))] = dataset_with_nodes_reconstructed$group_color[which(!is.na(idx))]
dataset_with_nodes_reconstructed_subsampled = dataset_with_nodes_reconstructed
###################################################

###################################################
## Plot full tree
###################################################
pdf('2_analysis_index/3_sensitivity_subsample/Index_SSonnei_with_groups_subsampled_20260313.pdf', width = 7/2.54, height = 5/2.54, fonts = 'Arial', pointsize = 12)
par(mfcol = c(2,1), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.5, cex.sub=0.3)

dataset_with_nodes = dataset_with_nodes
root_height = root_height

min_year = 1990
max_year = 2020
max_index = 0.25

## Plot tree
tree_to_plot = ladderize(tree, right = F)
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'no')], cex = 0.15) # 0.15
# nodelabels(pch = 16, col = dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'yes')], cex = 0.5)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 5)-root_height,
             lab_axis = seq(min_year, max_year, 5), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15,  #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, pch = 16)
axis(1, at = seq(min_year, max_year, 5), labels = seq(min_year, max_year, 5), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. Sonnei", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

idx = which(!duplicated(dataset_with_nodes$group_color))
legend('topleft', legend = dataset_with_nodes$groups[idx], col = dataset_with_nodes$group_color[idx], pch = 16, bty = 'n', 
       cex = 0.25)
dev.off()
###################################################

###################################################
## Plot with groups pMSM and non-pMSM, group from full phylogeny
###################################################
pdf('2_analysis_index/3_sensitivity_subsample/Index_SSonnei_with_groups_subsampled_20260313_non_pMSM_groupsALL_2000.pdf', width = 8/2.54, height = 5/2.54, fonts = 'Arial', pointsize = 12)
par(mfcol = c(2,2), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.5, cex.sub=0.3)

dataset_with_nodes = dataset_with_nodes_reconstructed
root_height = min(dataset_with_nodes_pMSM$time)

min_year = 2000
max_year = 2020
max_index = 0.5

## Plot tree
tree_to_plot = ladderize(tree_pMSM, right = F)
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$group_color[match(tree_pMSM$tip.label, dataset_with_nodes$name_seq)], cex = 0.15) # 0.15
# nodelabels(pch = 16, col = dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'yes')], cex = 0.5)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 5)-root_height,
             lab_axis = seq(min_year, max_year, 5), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index_pMSM[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15,  #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index_pMSM[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, pch = 16)
axis(1, at = seq(min_year, max_year, 5), labels = seq(min_year, max_year, 5), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. Sonnei pMSM", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)


dataset_with_nodes = dataset_with_nodes_reconstructed
root_height = min(dataset_with_nodes_nonpMSM$time)

min_year = 2000
max_year = 2020
max_index = 0.5

tree = ladderize(tree, right = F)
tree_pMSM2 = drop.tip(tree, tip = names_seqs[which(is.na(match(names_seqs, names_seqs_pMSM)))])
tree_nonpMSM2 = drop.tip(tree, tip = names_seqs[which(is.na(match(names_seqs, names_seqs_nonpMSM)))])

## Plot tree
tree_to_plot = tree_nonpMSM2
# tree_to_plot = ape::rotate(tree_to_plot, node = tree_to_plot$Nnode+1+5)
tree_to_plot = phytools::rotateNodes(tree_to_plot, node = tree_to_plot$Nnode+1+5)
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$group_color[match(tree_to_plot$tip.label, dataset_with_nodes$name_seq)], cex = 0.15) # 0.15
# nodelabels(pch = 16, col = dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'yes')], cex = 0.5)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 5)-root_height,
             lab_axis = seq(min_year, max_year, 5), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index_nonpMSM[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15,  #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index_nonpMSM[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, pch = 16)
axis(1, at = seq(min_year, max_year, 5), labels = seq(min_year, max_year, 5), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. Sonnei non pMSM", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

dev.off()
###################################################

###################################################
## Plot index full tree vs subsampled
###################################################
pdf('2_analysis_index/3_sensitivity_subsample/Index_SSonnei_with_groups_subsampled_comparison_20260313.pdf', width = 10/2.54, height = 3.5/2.54, fonts = 'Arial', pointsize = 12)
par(mfcol = c(1,3), oma = c(0,0,0,0.5), mar = c(1,1,1,1), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=1, cex.lab=1, cex.main=0.5, cex.sub=0.3)
load('2_analysis_index/2_find_index_groups/Lineages_detected_20260301.Rdata')

idx = match(dataset_with_nodes_reconstructed_subsampled$name_seq[which(dataset_with_nodes_reconstructed_subsampled$is.node == 'no')], 
            dataset_with_nodes_reconstructed$name_seq)
plot(dataset_with_nodes_reconstructed_subsampled$index[which(dataset_with_nodes_reconstructed_subsampled$is.node == 'no')], 
     dataset_with_nodes_reconstructed$index[idx], 
     xlim = c(0, 0.25), ylim = c(0, 0.25), bty = 'n', yaxt = 'n', pch = 16,cex = 0.5,
     col = dataset_with_nodes_reconstructed$group_color[idx], 
     xlab = 'Index in subsampled dataset', ylab  = 'Index in full dataset', main = 'Full')
axis(2, las = 2)

idx = match(dataset_with_nodes_reconstructed_subsampled$name_seq[which(dataset_with_nodes_reconstructed_subsampled$is.node == 'no' & dataset_with_nodes_reconstructed_subsampled$pmsmdemdef == 'pMSM')], 
            dataset_with_nodes_reconstructed$name_seq)
plot(dataset_with_nodes_reconstructed_subsampled$index_pMSM[which(dataset_with_nodes_reconstructed_subsampled$is.node == 'no' & dataset_with_nodes_reconstructed_subsampled$pmsmdemdef == 'pMSM')], 
     dataset_with_nodes_reconstructed$index_pMSM[idx], 
     xlim = c(0, 0.5), ylim = c(0, 0.5), bty = 'n', yaxt = 'n', pch = 16, cex = 0.5,
     col = dataset_with_nodes_reconstructed$group_color[idx], 
     xlab = 'Index in subsampled dataset', ylab  = 'Index in full dataset', main = 'pMSM')
axis(2, las = 2)

idx = match(dataset_with_nodes_reconstructed_subsampled$name_seq[which(dataset_with_nodes_reconstructed_subsampled$is.node == 'no' & dataset_with_nodes_reconstructed_subsampled$pmsmdemdef == 'non-pMSM')], 
            dataset_with_nodes_reconstructed$name_seq)
plot(dataset_with_nodes_reconstructed_subsampled$index_nonpMSM[which(dataset_with_nodes_reconstructed_subsampled$is.node == 'no' & dataset_with_nodes_reconstructed_subsampled$pmsmdemdef == 'non-pMSM')], 
     dataset_with_nodes_reconstructed$index_nonpMSM[idx], 
     xlim = c(0, 0.5), ylim = c(0, 0.5), bty = 'n', yaxt = 'n', pch = 16,cex = 0.5,
     col = dataset_with_nodes_reconstructed$group_color[idx], 
     xlab = 'Index in subsampled dataset', ylab  = 'Index in full dataset', main = 'non-pMSM')
axis(2, las = 2)
dev.off()
###################################################





