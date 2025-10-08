########################################################################################################################################
## Plot detected lineages
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

setwd('~/Dropbox/Projects/2025_Phylowave_SSonnei/Phylowave_SSonnei/')
load('2_analysis_index/2_find_index_groups/Lineages_detected_20251007.Rdata')

########################################################################################################################################
## Plot with groups - full phylogeny
########################################################################################################################################
pdf('2_analysis_index/2_find_index_groups/Index_SSonnei_with_groups_20251007.pdf', width = 7/2.54, height = 5/2.54, fonts = 'Arial', pointsize = 12)
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
########################################################################################################################################

## Write dataset, full phylogeny
write.tree(tree, file = 'Tree_SSonnei_full_20251007.tree')
write.csv(dataset_with_nodes, file = 'dataset_with_nodes_20251007.csv')

########################################################################################################################################
## Plot with groups STN and non-STN, groups detected independently
########################################################################################################################################
pdf('2_analysis_index/2_find_index_groups/Index_SSonnei_with_groups_20251007_non_STN_0.5_2000.pdf', width = 12/2.54, height = 5/2.54, fonts = 'Arial', pointsize = 12)
par(mfcol = c(2,2), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.5, cex.sub=0.3)

dataset_with_nodes = dataset_with_nodes_STN
root_height = min(dataset_with_nodes_STN$time)

min_year = 2000
max_year = 2020
max_index = 0.5

## Plot tree
tree_to_plot = tree_STN
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'no')], cex = 0.15) # 0.15
# nodelabels(pch = 16, col = dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'yes')], cex = 0.5)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 2)-root_height,
             lab_axis = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
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
axis(1, at = seq(min_year, max_year, 2), labels = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. Sonnei STN", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

idx = which(!duplicated(dataset_with_nodes$group_color))
legend('topleft', legend = dataset_with_nodes$groups[idx], col = dataset_with_nodes$group_color[idx], pch = 16, bty = 'n', 
       cex = 0.25)


dataset_with_nodes = dataset_with_nodes_nonSTN
root_height = min(dataset_with_nodes_nonSTN$time)

min_year = 2000
max_year = 2020
max_index = 0.5

## Plot tree
tree_to_plot = tree_nonSTN
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'no')], cex = 0.15) # 0.15
# nodelabels(pch = 16, col = dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'yes')], cex = 0.5)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 2)-root_height,
             lab_axis = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
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
axis(1, at = seq(min_year, max_year, 2), labels = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. Sonnei non STN", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

idx = which(!duplicated(dataset_with_nodes$group_color))
legend('topleft', legend = dataset_with_nodes$groups[idx], col = dataset_with_nodes$group_color[idx], pch = 16, bty = 'n', 
       cex = 0.25)

dev.off()
########################################################################################################################################

########################################################################################################################################
## Plot with groups STN and non-STN, group from full phylogeny
########################################################################################################################################
pdf('2_analysis_index/2_find_index_groups/Index_SSonnei_with_groups_20251007_non_STN_groupsALL_2000.pdf', width = 8/2.54, height = 5/2.54, fonts = 'Arial', pointsize = 12)
par(mfcol = c(2,2), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.5, cex.sub=0.3)

dataset_with_nodes = dataset_with_nodes_reconstructed
root_height = min(dataset_with_nodes_STN$time)

min_year = 2000
max_year = 2020
max_index = 0.5

## Plot tree
tree_to_plot = ladderize(tree_STN, right = F)
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$group_color[match(tree_STN$tip.label, dataset_with_nodes$name_seq)], cex = 0.15) # 0.15
# nodelabels(pch = 16, col = dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'yes')], cex = 0.5)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 5)-root_height,
             lab_axis = seq(min_year, max_year, 5), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index_STN[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15,  #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index_STN[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, pch = 16)
axis(1, at = seq(min_year, max_year, 5), labels = seq(min_year, max_year, 5), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. Sonnei STN", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)


dataset_with_nodes = dataset_with_nodes_reconstructed
root_height = min(dataset_with_nodes_nonSTN$time)

min_year = 2000
max_year = 2020
max_index = 0.5

tree = ladderize(tree, right = F)
tree_STN2 = drop.tip(tree, tip = names_seqs[which(is.na(match(names_seqs, names_seqs_STN)))])
tree_nonSTN2 = drop.tip(tree, tip = names_seqs[which(is.na(match(names_seqs, names_seqs_nonSTN)))])

## Plot tree
tree_to_plot = tree_nonSTN2
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
     dataset_with_nodes$index_nonSTN[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.15,  #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index_nonSTN[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.15, pch = 16)
axis(1, at = seq(min_year, max_year, 5), labels = seq(min_year, max_year, 5), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. Sonnei non STN", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

dev.off()
########################################################################################################################################

########################################################################################################################################
## Plot with groups pMSM and non-pMSM, groups detected independently
########################################################################################################################################
pdf('2_analysis_index/2_find_index_groups/Index_SSonnei_with_groups_20251007_non_pMSM_0.5_2000.pdf', width = 12/2.54, height = 5/2.54, fonts = 'Arial', pointsize = 12)
par(mfcol = c(2,2), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.5, cex.sub=0.3)

dataset_with_nodes = dataset_with_nodes_pMSM
root_height = min(dataset_with_nodes_pMSM$time)

min_year = 2000
max_year = 2020
max_index = 0.5

## Plot tree
tree_to_plot = tree_pMSM
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'no')], cex = 0.15) # 0.15
# nodelabels(pch = 16, col = dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'yes')], cex = 0.5)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 2)-root_height,
             lab_axis = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
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
axis(1, at = seq(min_year, max_year, 2), labels = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. Sonnei pMSM", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

idx = which(!duplicated(dataset_with_nodes$group_color))
legend('topleft', legend = dataset_with_nodes$groups[idx], col = dataset_with_nodes$group_color[idx], pch = 16, bty = 'n', 
       cex = 0.25)


dataset_with_nodes = dataset_with_nodes_nonpMSM
root_height = min(dataset_with_nodes_nonpMSM$time)

min_year = 2000
max_year = 2020
max_index = 0.5

## Plot tree
tree_to_plot = tree_nonpMSM
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'no')], cex = 0.15) # 0.15
# nodelabels(pch = 16, col = dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'yes')], cex = 0.5)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 2)-root_height,
             lab_axis = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
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
axis(1, at = seq(min_year, max_year, 2), labels = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. Sonnei non pMSM", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

idx = which(!duplicated(dataset_with_nodes$group_color))
legend('topleft', legend = dataset_with_nodes$groups[idx], col = dataset_with_nodes$group_color[idx], pch = 16, bty = 'n', 
       cex = 0.25)

dev.off()
########################################################################################################################################

########################################################################################################################################
## Plot with groups pMSM and non-pMSM, group from full phylogeny
########################################################################################################################################
pdf('2_analysis_index/2_find_index_groups/Index_SSonnei_with_groups_20251007_non_pMSM_groupsALL_2000.pdf', width = 8/2.54, height = 5/2.54, fonts = 'Arial', pointsize = 12)
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
########################################################################################################################################

## Write dataset, full phylogeny
write.csv(dataset_with_nodes_reconstructed, file = 'dataset_with_nodes_20251007_all.csv')

########################################################################################################################################
## Plot with groups - MSM pandemic
########################################################################################################################################
pdf('2_analysis_index/2_find_index_groups/Index_SSonnei_with_groups_MSM_pandemic_20251007.pdf', width = 4/2.54, height = 5/2.54, fonts = 'Arial', pointsize = 12)
par(mfcol = c(2,1), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.5, cex.sub=0.3)

dataset_with_nodes = dataset_with_nodes_MSM_pandemic
root_height = root_height_MSM_pandemic

min_year = 2014
max_year = 2022
max_index = 0.75

## Plot tree
tree_to_plot = ladderize(tree_MSM_pandemic, right = F)
plot(tree_to_plot, show.tip.label = FALSE, edge.width = 0.15, edge.color = 'grey', 
     x.lim = c(min_year, max_year)-root_height)
tiplabels(pch = 16, col = dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'no')], cex = 0.2) # 0.15
# nodelabels(pch = 16, col = dataset_with_nodes$genotype_color[which(dataset_with_nodes$is.node == 'yes')], cex = 0.5)
axisPhylo_NL(side = 1, root.time = root_height, backward = F,
             at_axis = seq(min_year, max_year, 2)-root_height,
             lab_axis = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
## Plot thd
plot(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'yes')], 
     dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'yes')], 
     col = adjustcolor(dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'yes')], alpha.f = 1),
     bty = 'n', xlim = c(min_year, max_year), cex = 0.2,  #0.15
     pch = 16, bty = 'n', ylim = c(0, max_index), 
     # main = paste0('World'), 
     ylab = '', xlab = '', yaxt = 'n', xaxt = 'n')
points(dataset_with_nodes$time[which(dataset_with_nodes$is.node == 'no')], 
       dataset_with_nodes$index[which(dataset_with_nodes$is.node == 'no')], 
       col = adjustcolor(dataset_with_nodes$group_color[which(dataset_with_nodes$is.node == 'no')], alpha.f = 1),
       cex = 0.2, pch = 16)
axis(1, at = seq(min_year, max_year, 2), labels = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.02, mgp = c(0,-0.4,0))
axis(2, las = 2, tck=-0.01, lwd = 0.5)
title(main="S. Sonnei", line=-0.5, outer = F)
title(ylab="Index", line=0.5, outer = F)
title(xlab="Time (years)", line=0, outer = F)

idx = which(!duplicated(dataset_with_nodes$group_color))
legend('topleft', legend = dataset_with_nodes$groups[idx], col = dataset_with_nodes$group_color[idx], pch = 16, bty = 'n', 
       cex = 0.25)

dev.off()
########################################################################################################################################

## Write dataset, full phylogeny
write.csv(dataset_with_nodes_MSM_pandemic, file = 'dataset_with_nodes_MSM_pandemic_20251007.csv')





########################################################################################################################################
## Plot model fit
## FULL
########################################################################################################################################
load('2_analysis_index/2_find_index_groups/Lineages_detected_20251007.Rdata')
source('~/Dropbox/Projects/202312_Index_paper_methods_Sc2_H3N2_BP_TB/Phylowave_Learning-fitness-dynamics-pathogens-in-phylogenies/2_Functions/2_1_Index_computation_20240909.R')
source('~/Dropbox/Projects/202312_Index_paper_methods_Sc2_H3N2_BP_TB/Phylowave_Learning-fitness-dynamics-pathogens-in-phylogenies/2_Functions/2_2_Lineage_detection_20240909.R')
library(cowplot)

par(mfcol = c(2,1), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)

min_year = 1990
max_year = 2020
max_index = 0.25

potential_splits_model = compute_fit_for_node_set(timed_tree = tree, metadata = dataset_with_nodes, 
                                                  node_set = potential_splits, weight_by_time = NULL, k_smooth = -1, plot = F, log_y = T) 
plot_fit_world = function(){
  par(oma = c(0,0,0,0), mar =  c(2,2,1,0), mgp = c(2,0.5,0))
  plot(dataset_with_nodes$time, 
       dataset_with_nodes$index,
       col = adjustcolor(dataset_with_nodes$group_color, alpha.f = 0.2),
       bty = 'n', xlim = c(min_year, max_year), cex = 0.5, 
       pch = 16, ylim = c(0, max_index), 
       ylab = 'Index', xlab = 'Time (years)', yaxt = 'n', xaxt = 'n')
  points(dataset_with_nodes$time, 
         exp(predict(potential_splits_model$mod)), 
         col = adjustcolor(dataset_with_nodes$group_color, alpha.f = 1),
         bty = 'n', xlim = c(min_year, max_year), cex = 0.5, 
         pch = 16)
  axis(1, at = seq(min_year, max_year, 5), labels = seq(min_year, max_year, 5), lwd = 0.5, tck=-0.015)
  axis(2, las = 2, tck=-0.015, lwd = 0.5)
}
plot_fit_obs_pred = function(){
  par(oma = c(0,0,0,0.5), mar = c(2,2,1,0), mgp = c(2,0.5,0))
  idx = sample(1:nrow(dataset_with_nodes))
  plot(dataset_with_nodes$index[idx], 
       exp(predict(potential_splits_model$mod))[idx],
       col = adjustcolor(dataset_with_nodes$group_color[idx], alpha.f = 1),
       bty = 'n', xlim = c(0, max_index), cex = 0.5, xaxt = 'n',
       pch = 16, bty = 'n', ylim = c(0, max_index), 
       ylab = 'Predicted', xlab = 'Observed', yaxt = 'n')
  abline(a = 0, b = 1, lty = 2)
  axis(1, at = seq(0, max_index, 0.05), labels = seq(0, max_index, 0.05), lwd = 0.5, tck=-0.015)
  axis(2, at = seq(0, max_index, 0.05), labels = seq(0, max_index, 0.05), las = 2, tck=-0.015, lwd = 0.5)
}

p = plot_grid(ggdraw(plot_fit_world), ggdraw(plot_fit_obs_pred),
              rel_widths = c(2, 1), labels = c('A', 'B'), ncol = 2)
ggsave(filename = '2_analysis_index/2_find_index_groups/Fit_Index_SSonnei_with_groups_20251007.pdf',
       plot = p, device = 'pdf', scale = 1, width = 15, height = 5.5, units = 'cm')
########################################################################################################################################

########################################################################################################################################
## Plot model fit
## MSM pandemic
########################################################################################################################################
load('2_analysis_index/2_find_index_groups/Lineages_detected_20251007.Rdata')
library(cowplot)

par(mfcol = c(2,1), oma = c(0,0,0,0.5), mar = c(1,1,0,0), mgp = c(0,0.1,0), family = 'Arial',
    cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)

min_year = 2014
max_year = 2022
max_index = 0.6

potential_splits_model = compute_fit_for_node_set(timed_tree = tree_MSM_pandemic, metadata = dataset_with_nodes_MSM_pandemic, 
                                                  node_set = potential_splits_MSM_pandemic, weight_by_time = NULL, k_smooth = -1, plot = F, log_y = T) 
plot_fit_world = function(){
  par(oma = c(0,0,0,0), mar =  c(2,2,1,0), mgp = c(2,0.5,0))
  plot(dataset_with_nodes_MSM_pandemic$time, 
       dataset_with_nodes_MSM_pandemic$index,
       col = adjustcolor(dataset_with_nodes_MSM_pandemic$group_color, alpha.f = 0.3),
       bty = 'n', xlim = c(min_year, max_year), cex = 0.5, 
       pch = 16, ylim = c(0, max_index), 
       ylab = 'Index', xlab = 'Time (years)', yaxt = 'n', xaxt = 'n')
  points(dataset_with_nodes_MSM_pandemic$time, 
         exp(predict(potential_splits_model$mod)), 
         col = adjustcolor(dataset_with_nodes_MSM_pandemic$group_color, alpha.f = 1),
         bty = 'n', xlim = c(min_year, max_year), cex = 0.5, 
         pch = 16)
  axis(1, at = seq(min_year, max_year, 2), labels = seq(min_year, max_year, 2), lwd = 0.5, tck=-0.015)
  axis(2, las = 2, tck=-0.015, lwd = 0.5)
}
plot_fit_obs_pred = function(){
  par(oma = c(0,0,0,0.5), mar = c(2,2,1,0), mgp = c(2,0.5,0))
  idx = sample(1:nrow(dataset_with_nodes_MSM_pandemic))
  plot(dataset_with_nodes_MSM_pandemic$index[idx], 
       exp(predict(potential_splits_model$mod))[idx],
       col = adjustcolor(dataset_with_nodes_MSM_pandemic$group_color[idx], alpha.f = 1),
       bty = 'n', xlim = c(0, max_index), cex = 0.5, xaxt = 'n',
       pch = 16, bty = 'n', ylim = c(0, max_index), 
       ylab = 'Predicted', xlab = 'Observed', yaxt = 'n')
  abline(a = 0, b = 1, lty = 2)
  axis(1, at = seq(0, max_index, 0.1), labels = seq(0, max_index, 0.1), lwd = 0.5, tck=-0.015)
  axis(2, at = seq(0, max_index, 0.1), labels = seq(0, max_index, 0.1), las = 2, tck=-0.015, lwd = 0.5)
}

p = plot_grid(ggdraw(plot_fit_world), ggdraw(plot_fit_obs_pred),
              rel_widths = c(2, 1), labels = c('A', 'B'), ncol = 2)
ggsave(filename = '2_analysis_index/2_find_index_groups/Fit_Index_SSonnei_with_groups_MSM_pandemic_20251007.pdf',
       plot = p, device = 'pdf', scale = 1, width = 15, height = 5.5, units = 'cm')
########################################################################################################################################






