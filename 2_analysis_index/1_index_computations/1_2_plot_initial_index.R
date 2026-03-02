########################################################################################################################################
## Plot index: input: Initial_index_computation_and_parameters_20260301.Rdata
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

load('2_analysis_index/1_index_computations/Initial_index_computation_and_parameters_20260301.Rdata')

########################################################################################################################################
## Plot tree & index below, with colors from full genotypes
########################################################################################################################################
pdf('2_analysis_index/1_index_computations/Index_SSonnei_UK_20260301.pdf', width = 6/2.54, height = 6/2.54)
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
pdf('2_analysis_index/1_index_computations/Index_SSonnei_UK_20260301_non_pMSM.pdf', width = 10/2.54, height = 5/2.54)
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
pdf('2_analysis_index/1_index_computations/Index_SSonnei_UK_20260301_travel.pdf', width = 10/2.54, height = 5/2.54)
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

########################################################################################################################################
##  Plot same dynamics for MSM pandemic
########################################################################################################################################
pdf('2_analysis_index/1_index_computations/Index_SSonnei_UK_20260301_MSM_pandemic_simplified.pdf', width = 6/2.54, height = 6/2.54)
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
tiplabels(pch = 16, col = dataset_with_nodes$Genotype_simplified_color, cex = 0.15)
axisPhylo_NL(side = 1, root.time = root_height_plot, backward = F,
             at_axis = seq(min_year, max_year, 2)-root_height_plot,
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
title(main="S. sonnei", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

idx = which(!duplicated(dataset_with_nodes$Genotype))
legend('topleft', legend = dataset_with_nodes$Genotype[idx], col = dataset_with_nodes$Genotype_simplified_color[idx], pch = 16,
       cex = 0.15, ncol = 3, bty = 'n')

dev.off()
########################################################################################################################################


#### Genotype simplified 
load('2_analysis_index/1_index_computations/Initial_index_computation_and_parameters_20260301.Rdata')

########################################################################################################################################
## Plot tree & index below, with colors from simplified genotypes
########################################################################################################################################
pdf('2_analysis_index/1_index_computations/Index_SSonnei_UK_20260301_simplified.pdf', width = 6/2.54, height = 6/2.54)
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
pdf('2_analysis_index/1_index_computations/Index_SSonnei_UK_20260301_non_pMSM_simplified.pdf', width = 10/2.54, height = 5/2.54)
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
pdf('2_analysis_index/1_index_computations/Index_SSonnei_UK_20260301_MSM_pandemic.pdf', width = 6/2.54, height = 6/2.54)
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


