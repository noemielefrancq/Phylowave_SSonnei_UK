########################################################################################################################################
## Comparison detected phylowave groups vs. genotypes
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
library(cowplot)
library(extrafont)
library(ComplexHeatmap)
library(circlize)
loadfonts(device="all")

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
"%||%" <- devtools:::`%||%`
########################################################################################################################################

########################################################################################################################################
## Load data
########################################################################################################################################
setwd('~/Dropbox/Projects/2025_Phylowave_SSonnei/Phylowave_SSonnei/')
load('2_analysis_index/2_find_index_groups/Lineages_detected_20251007.Rdata')
########################################################################################################################################

########################################################################################################################################
## Find colors genotypes
########################################################################################################################################
eq_Gen_Genotypes = matrix(NA, ncol = 2, nrow = length(unique(dataset_with_nodes$Genotype)))
eq_Gen_Genotypes[,1] = unique(dataset_with_nodes$Genotype)
eq_Gen_Genotypes[,2] = dataset_with_nodes$Gen[match(unique(dataset_with_nodes$Genotype), dataset_with_nodes$Genotype)]

dataset_with_nodes$Genotype = eq_Gen_Genotypes[match(dataset_with_nodes$Genotype, eq_Gen_Genotypes[,1]),2]
dataset_with_nodes$Genotype = factor(dataset_with_nodes$Genotype)
dataset_with_nodes$Genotype = as.character(dataset_with_nodes$Genotype)
colors_groups_full_Genotype = MetBrewer::met.brewer(name="Cross", n=length(levels(as.factor(dataset_with_nodes$Gen))), type="continuous")
names(colors_groups_full_Genotype) = levels(as.factor(dataset_with_nodes$Gen))

dataset_with_nodes_STN$Genotype = eq_Gen_Genotypes[match(dataset_with_nodes_STN$Genotype, eq_Gen_Genotypes[,1]),2]
dataset_with_nodes_STN$Genotype = factor(dataset_with_nodes_STN$Genotype)
dataset_with_nodes_STN$Genotype = as.character(dataset_with_nodes_STN$Genotype)
colors_groups_STN_Genotype = MetBrewer::met.brewer(name="Cross", n=length(levels(as.factor(dataset_with_nodes_STN$Genotype))), type="continuous")
names(colors_groups_STN_Genotype) = levels(as.factor(dataset_with_nodes_STN$Genotype))
# names(colors_groups_STN_Genotype) = dataset_with_nodes$Gen[match(names(colors_groups_STN_Genotype), dataset_with_nodes$Genotype)]

dataset_with_nodes_nonSTN$Genotype = eq_Gen_Genotypes[match(dataset_with_nodes_nonSTN$Genotype, eq_Gen_Genotypes[,1]),2]
dataset_with_nodes_nonSTN$Genotype = factor(dataset_with_nodes_nonSTN$Genotype)
dataset_with_nodes_nonSTN$Genotype = as.character(dataset_with_nodes_nonSTN$Genotype)
colors_groups_nonSTN_Genotype = MetBrewer::met.brewer(name="Cross", n=length(levels(as.factor(dataset_with_nodes_nonSTN$Genotype))), type="continuous")
names(colors_groups_nonSTN_Genotype) = levels(as.factor(dataset_with_nodes_nonSTN$Genotype))
# names(colors_groups_nonSTN_Genotype) = dataset_with_nodes$Gen[match(names(colors_groups_nonSTN_Genotype), dataset_with_nodes$Genotype)]

dataset_with_nodes_MSM_pandemic$Genotype = eq_Gen_Genotypes[match(dataset_with_nodes_MSM_pandemic$Genotype, eq_Gen_Genotypes[,1]),2]
dataset_with_nodes_MSM_pandemic$Genotype = factor(dataset_with_nodes_MSM_pandemic$Genotype)
dataset_with_nodes_MSM_pandemic$Genotype = as.character(dataset_with_nodes_MSM_pandemic$Genotype)
colors_groups_MSM_pandemic_Genotype = MetBrewer::met.brewer(name="Cross", n=length(levels(as.factor(dataset_with_nodes_MSM_pandemic$Genotype))), type="continuous")
names(colors_groups_MSM_pandemic_Genotype) = levels(as.factor(dataset_with_nodes_MSM_pandemic$Genotype))
# names(colors_groups_MSM_pandemic_Genotype) = dataset_with_nodes$Gen[match(names(colors_groups_MSM_pandemic_Genotype), dataset_with_nodes$Genotype)]
########################################################################################################################################


########################################################################################################################################
## Plot comparison - full data
########################################################################################################################################
dataset_with_nodes$Genotype_simplified = as.character(dataset_with_nodes$Genotype_simplified)
dataset_with_nodes$Genotype_simplified[which(dataset_with_nodes$Genotype_simplified == 'NA' | dataset_with_nodes$Genotype_simplified == 'ZOthers')] = NA
dataset_with_nodes$Genotype_simplified = factor(dataset_with_nodes$Genotype_simplified)
Genotypes = matrix(as.character(dataset_with_nodes$Genotype_simplified), ncol = 1)
colnames(Genotypes) = 'groups'
rownames(Genotypes) = dataset_with_nodes$name_seq
plot_tree_Genotypes <- ggtree(tree, mrsd=lubridate::date_decimal(max(dataset_with_nodes$time)), 
                              size = 0.1,
                              aes(color = as.character(dataset_with_nodes$Genotype_simplified))) + 
  scale_color_manual(values = colors_clades_simplified, na.value = 'grey90')+
  theme_tree2()
plot_tree_Genotypes = gheatmap(plot_tree_Genotypes, Genotypes, offset=0.1, width=0.10, 
                               colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = colors_clades_simplified, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

groups = matrix(dataset_with_nodes$groups, ncol = 1)
colnames(groups) = 'groups'
rownames(groups) = dataset_with_nodes$name_seq
cols = as.character(colors_groups)
names(cols) = as.character(1:max(as.numeric(name_groups)))
plot_tree_groups <- ggtree(tree, mrsd=lubridate::date_decimal(max(dataset_with_nodes$time)), 
                           size = 0.1,
                           aes(color = as.character(split$groups))) + 
  scale_color_manual(values = cols, na.value = 'grey90')+
  theme_tree2()
plot_tree_groups = gheatmap(plot_tree_groups, groups, offset=0.1, width=0.10, 
                            colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = cols, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

comparison = plot_grid(plot_tree_groups, plot_tree_Genotypes, 
                       rel_widths = c(1, 1), labels = '', ncol = 2)

ggsave(filename = '2_analysis_index/2_find_index_groups/Comparison_groups_genotypes_full_simplified_20251007.pdf', plot = comparison, device = 'pdf', scale = 1,
       width = 10, height = 10, units = 'cm')
########################################################################################################################################

########################################################################################################################################
## Plot comparison - STN data
########################################################################################################################################
dataset_with_nodes_STN$Genotype_simplified = as.character(dataset_with_nodes_STN$Genotype_simplified)
dataset_with_nodes_STN$Genotype_simplified[which(dataset_with_nodes_STN$Genotype_simplified == 'NA' | dataset_with_nodes_STN$Genotype_simplified == 'ZOthers')] = NA
dataset_with_nodes_STN$Genotype_simplified = factor(dataset_with_nodes_STN$Genotype_simplified)
Genotypes_STN = matrix(as.character(dataset_with_nodes_STN$Genotype_simplified), ncol = 1)
colnames(Genotypes_STN) = 'groups'
rownames(Genotypes_STN) = dataset_with_nodes_STN$name_seq
plot_tree_Genotypes_STN <- ggtree(tree_STN, mrsd=lubridate::date_decimal(max(dataset_with_nodes_STN$time)), 
                                  size = 0.15,
                                  aes(color = as.character(dataset_with_nodes_STN$Genotype_simplified))) + 
  scale_color_manual(values = colors_clades_simplified, na.value = 'grey90')+
  theme_tree2()
plot_tree_Genotypes_STN = gheatmap(plot_tree_Genotypes_STN, Genotypes_STN, offset=0.1, width=0.10, 
                                   colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = colors_clades_simplified, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

groups_STN = matrix(as.character(dataset_with_nodes_reconstructed$groups[match(dataset_with_nodes_STN$name_seq2, dataset_with_nodes_reconstructed$name_seq2)]), ncol = 1)
colnames(groups_STN) = 'groups'
rownames(groups_STN) = dataset_with_nodes_STN$name_seq
cols = as.character(colors_groups)
names(cols) = as.character(1:max(as.numeric(name_groups)))
plot_tree_groups_STN <- ggtree(tree_STN, mrsd=lubridate::date_decimal(max(dataset_with_nodes_STN$time)), 
                               size = 0.15,
                               aes(color = as.character(dataset_with_nodes_reconstructed$groups[match(dataset_with_nodes_STN$name_seq2, dataset_with_nodes_reconstructed$name_seq2)]))) + 
  scale_color_manual(values = cols, na.value = 'grey90')+
  theme_tree2()
plot_tree_groups_STN = gheatmap(plot_tree_groups_STN, groups_STN, offset=0.1, width=0.10, 
                                colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = cols, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

comparison = plot_grid(plot_tree_groups_STN, plot_tree_Genotypes_STN, 
                       rel_widths = c(1, 1), labels = '', ncol = 2)

ggsave(filename = '2_analysis_index/2_find_index_groups/Comparison_groups_genotypes_STN_simplified_20251007.pdf', plot = comparison, device = 'pdf', scale = 1,
       width = 7, height = 7, units = 'cm')
########################################################################################################################################

########################################################################################################################################
## Plot comparison - non STN data
########################################################################################################################################
dataset_with_nodes_nonSTN$Genotype_simplified = as.character(dataset_with_nodes_nonSTN$Genotype_simplified)
dataset_with_nodes_nonSTN$Genotype_simplified[which(dataset_with_nodes_nonSTN$Genotype_simplified == 'NA' | dataset_with_nodes_nonSTN$Genotype_simplified == 'ZOthers')] = NA
dataset_with_nodes_nonSTN$Genotype_simplified = factor(dataset_with_nodes_nonSTN$Genotype_simplified)
Genotypes_nonSTN = matrix(as.character(dataset_with_nodes_nonSTN$Genotype_simplified), ncol = 1)
colnames(Genotypes_nonSTN) = 'groups'
rownames(Genotypes_nonSTN) = dataset_with_nodes_nonSTN$name_seq
plot_tree_Genotypes_nonSTN <- ggtree(tree_nonSTN, mrsd=lubridate::date_decimal(max(dataset_with_nodes_nonSTN$time)), 
                                     size = 0.15,
                                     aes(color = as.character(dataset_with_nodes_nonSTN$Genotype_simplified))) + 
  scale_color_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  # geom_text(aes(label=node), size = 3)+ ## to check for rotation
  theme_tree2()
plot_tree_Genotypes_nonSTN = rotate(plot_tree_Genotypes_nonSTN, 1656) ## to put yellow clade up
plot_tree_Genotypes_nonSTN = rotate(plot_tree_Genotypes_nonSTN, 1657)
plot_tree_Genotypes_nonSTN = rotate(plot_tree_Genotypes_nonSTN, 2435)
plot_tree_Genotypes_nonSTN = gheatmap(plot_tree_Genotypes_nonSTN, Genotypes_nonSTN, offset=0.1, width=0.10, 
                                      colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

groups_nonSTN = matrix(as.character(dataset_with_nodes_reconstructed$groups[match(dataset_with_nodes_nonSTN$name_seq2, dataset_with_nodes_reconstructed$name_seq2)]), ncol = 1)
colnames(groups_nonSTN) = 'groups'
rownames(groups_nonSTN) = dataset_with_nodes_nonSTN$name_seq
cols = as.character(colors_groups)
names(cols) = as.character(1:max(as.numeric(name_groups)))
plot_tree_groups_nonSTN <- ggtree(tree_nonSTN, mrsd=lubridate::date_decimal(max(dataset_with_nodes_nonSTN$time)), 
                                  size = 0.15,
                                  aes(color = as.character(dataset_with_nodes_reconstructed$groups[match(dataset_with_nodes_nonSTN$name_seq2, dataset_with_nodes_reconstructed$name_seq2)]))) + 
  scale_color_manual(values = cols, na.value = 'grey90')+
  # geom_text(aes(label=node), size = 3)+ ## to check for rotation
  theme_tree2()
# plot_tree_groups_nonSTN = rotate(plot_tree_groups_nonSTN, 1656) ## to put yellow clade up
# plot_tree_groups_nonSTN = rotate(plot_tree_groups_nonSTN, 1657)
# plot_tree_groups_nonSTN = rotate(plot_tree_groups_nonSTN, 2435)
plot_tree_groups_nonSTN = gheatmap(plot_tree_groups_nonSTN, groups_nonSTN, offset=0.1, width=0.10, 
                                   colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = cols, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

comparison = plot_grid(plot_tree_groups_nonSTN, plot_tree_Genotypes_nonSTN, 
                       rel_widths = c(1, 1), labels = '', ncol = 2)

ggsave(filename = '2_analysis_index/2_find_index_groups/Comparison_groups_genotypes_nonSTN_simplified_20251007.pdf', plot = comparison, device = 'pdf', scale = 1,
       width = 7, height = 7, units = 'cm')
########################################################################################################################################

########################################################################################################################################
## Plot comparison - pMSM data
########################################################################################################################################
dataset_with_nodes_pMSM$Genotype_simplified = as.character(dataset_with_nodes_pMSM$Genotype_simplified)
dataset_with_nodes_pMSM$Genotype_simplified[which(dataset_with_nodes_pMSM$Genotype_simplified == 'NA' | dataset_with_nodes_pMSM$Genotype_simplified == 'ZOthers')] = NA
dataset_with_nodes_pMSM$Genotype_simplified = factor(dataset_with_nodes_pMSM$Genotype_simplified)
Genotypes_pMSM = matrix(as.character(dataset_with_nodes_pMSM$Genotype_simplified), ncol = 1)
colnames(Genotypes_pMSM) = 'groups'
rownames(Genotypes_pMSM) = dataset_with_nodes_pMSM$name_seq
plot_tree_Genotypes_pMSM <- ggtree(tree_pMSM, mrsd=lubridate::date_decimal(max(dataset_with_nodes_pMSM$time)), 
                                  size = 0.15,
                                  aes(color = as.character(dataset_with_nodes_pMSM$Genotype_simplified))) + 
  scale_color_manual(values = colors_clades_simplified, na.value = 'grey90')+
  theme_tree2()
plot_tree_Genotypes_pMSM = gheatmap(plot_tree_Genotypes_pMSM, Genotypes_pMSM, offset=0.1, width=0.10, 
                                   colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = colors_clades_simplified, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

groups_pMSM = matrix(as.character(dataset_with_nodes_reconstructed$groups[match(dataset_with_nodes_pMSM$name_seq2, dataset_with_nodes_reconstructed$name_seq2)]), ncol = 1)
colnames(groups_pMSM) = 'groups'
rownames(groups_pMSM) = dataset_with_nodes_pMSM$name_seq
cols = as.character(colors_groups)
names(cols) = as.character(1:max(as.numeric(name_groups)))
plot_tree_groups_pMSM <- ggtree(tree_pMSM, mrsd=lubridate::date_decimal(max(dataset_with_nodes_pMSM$time)), 
                               size = 0.15,
                               aes(color = as.character(dataset_with_nodes_reconstructed$groups[match(dataset_with_nodes_pMSM$name_seq2, dataset_with_nodes_reconstructed$name_seq2)]))) + 
  scale_color_manual(values = cols, na.value = 'grey90')+
  theme_tree2()
plot_tree_groups_pMSM = gheatmap(plot_tree_groups_pMSM, groups_pMSM, offset=0.1, width=0.10, 
                                colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = cols, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

comparison = plot_grid(plot_tree_groups_pMSM, plot_tree_Genotypes_pMSM, 
                       rel_widths = c(1, 1), labels = '', ncol = 2)

ggsave(filename = '2_analysis_index/2_find_index_groups/Comparison_groups_genotypes_pMSM_simplified_20251007.pdf', plot = comparison, device = 'pdf', scale = 1,
       width = 7, height = 7, units = 'cm')
########################################################################################################################################

########################################################################################################################################
## Plot comparison - non pMSM data
########################################################################################################################################
dataset_with_nodes_nonpMSM$Genotype_simplified = as.character(dataset_with_nodes_nonpMSM$Genotype_simplified)
dataset_with_nodes_nonpMSM$Genotype_simplified[which(dataset_with_nodes_nonpMSM$Genotype_simplified == 'NA' | dataset_with_nodes_nonpMSM$Genotype_simplified == 'ZOthers')] = NA
dataset_with_nodes_nonpMSM$Genotype_simplified = factor(dataset_with_nodes_nonpMSM$Genotype_simplified)
Genotypes_nonpMSM = matrix(as.character(dataset_with_nodes_nonpMSM$Genotype_simplified), ncol = 1)
colnames(Genotypes_nonpMSM) = 'groups'
rownames(Genotypes_nonpMSM) = dataset_with_nodes_nonpMSM$name_seq
plot_tree_Genotypes_nonpMSM <- ggtree(tree_nonpMSM, mrsd=lubridate::date_decimal(max(dataset_with_nodes_nonpMSM$time)), 
                                     size = 0.15,
                                     aes(color = as.character(dataset_with_nodes_nonpMSM$Genotype_simplified))) + 
  scale_color_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  # geom_text(aes(label=node), size = 3)+ ## to check for rotation
  theme_tree2()
plot_tree_Genotypes_nonpMSM = rotate(plot_tree_Genotypes_nonpMSM, 1656) ## to put yellow clade up
plot_tree_Genotypes_nonpMSM = rotate(plot_tree_Genotypes_nonpMSM, 1657)
plot_tree_Genotypes_nonpMSM = rotate(plot_tree_Genotypes_nonpMSM, 2435)
plot_tree_Genotypes_nonpMSM = gheatmap(plot_tree_Genotypes_nonpMSM, Genotypes_nonpMSM, offset=0.1, width=0.10, 
                                      colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

groups_nonpMSM = matrix(as.character(dataset_with_nodes_reconstructed$groups[match(dataset_with_nodes_nonpMSM$name_seq2, dataset_with_nodes_reconstructed$name_seq2)]), ncol = 1)
colnames(groups_nonpMSM) = 'groups'
rownames(groups_nonpMSM) = dataset_with_nodes_nonpMSM$name_seq
cols = as.character(colors_groups)
names(cols) = as.character(1:max(as.numeric(name_groups)))
plot_tree_groups_nonpMSM <- ggtree(tree_nonpMSM, mrsd=lubridate::date_decimal(max(dataset_with_nodes_nonpMSM$time)), 
                                  size = 0.15,
                                  aes(color = as.character(dataset_with_nodes_reconstructed$groups[match(dataset_with_nodes_nonpMSM$name_seq2, dataset_with_nodes_reconstructed$name_seq2)]))) + 
  scale_color_manual(values = cols, na.value = 'grey90')+
  # geom_text(aes(label=node), size = 3)+ ## to check for rotation
  theme_tree2()
# plot_tree_groups_nonpMSM = rotate(plot_tree_groups_nonpMSM, 1656) ## to put yellow clade up
# plot_tree_groups_nonpMSM = rotate(plot_tree_groups_nonpMSM, 1657)
# plot_tree_groups_nonpMSM = rotate(plot_tree_groups_nonpMSM, 2435)
plot_tree_groups_nonpMSM = gheatmap(plot_tree_groups_nonpMSM, groups_nonpMSM, offset=0.1, width=0.10, 
                                   colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = cols, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

comparison = plot_grid(plot_tree_groups_nonpMSM, plot_tree_Genotypes_nonpMSM, 
                       rel_widths = c(1, 1), labels = '', ncol = 2)

ggsave(filename = '2_analysis_index/2_find_index_groups/Comparison_groups_genotypes_nonpMSM_simplified_20251007.pdf', plot = comparison, device = 'pdf', scale = 1,
       width = 7, height = 7, units = 'cm')
########################################################################################################################################


## HEATMAPS FULL DATASET
#######################################################################################################################################
## Compute groups vs ESBLR
#######################################################################################################################################
dataset_with_nodes$ESBLR = as.character(dataset_with_nodes$ESBLR)
dataset_with_nodes$ESBLR = factor(dataset_with_nodes$ESBLR)

grouping1 = levels(as.factor(dataset_with_nodes$ESBLR))
grouping2 = levels(as.factor(dataset_with_nodes$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes$ESBLR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,  #TRUE,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)

correspondance_groups_N_EBSL = correspondance_groups_N
#######################################################################################################################################

#######################################################################################################################################
## Compute groups vs CiprofloxacinR
#######################################################################################################################################
dataset_with_nodes$CiprofloxacinR = as.character(dataset_with_nodes$CiprofloxacinR)
dataset_with_nodes$CiprofloxacinR = factor(dataset_with_nodes$CiprofloxacinR)

grouping1 = levels(as.factor(dataset_with_nodes$CiprofloxacinR))
grouping2 = levels(as.factor(dataset_with_nodes$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes$CiprofloxacinR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

correspondance_groups_N_CIPRO = correspondance_groups_N
#######################################################################################################################################

#######################################################################################################################################
## Compute groups vs Az
#######################################################################################################################################
dataset_with_nodes$AzithromycinR = as.character(dataset_with_nodes$AzithromycinR)
dataset_with_nodes$AzithromycinR = factor(dataset_with_nodes$AzithromycinR)

grouping1 = levels(as.factor(dataset_with_nodes$AzithromycinR))
grouping2 = levels(as.factor(dataset_with_nodes$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes$AzithromycinR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  # t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

correspondance_groups_N_Az = correspondance_groups_N
#######################################################################################################################################

#######################################################################################################################################
## Combine all AMR one Matrix and plot
#######################################################################################################################################
correspondance_groups_N_AMR = cbind(correspondance_groups_N_Az[,1], correspondance_groups_N_CIPRO[,1])
correspondance_groups_N_AMR = cbind(correspondance_groups_N_AMR, correspondance_groups_N_EBSL[,1])

correspondance_groups_N_AMR = correspondance_groups_N_AMR[nrow(correspondance_groups_N_AMR):1,]
plot_AMR_compaarison_numbers = Heatmap((correspondance_groups_N_AMR),  col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('white', RColorBrewer::brewer.pal(name = 'Reds', n = 9)[c(2,3,6,9)])),
                                       heatmap_legend_param = list(at = c(0, 0.5, 1), 
                                                                   title_gp = gpar(fontsize = 3),
                                                                   labels_gp = gpar(fontsize = 2),
                                                                   legend_height = unit(0.25, "cm"),
                                                                   legend_width = unit(0.025, "cm")), 
                                       show_row_dend = F, show_column_dend = F, 
                                       cluster_rows = F, cluster_columns = F,
                                       show_heatmap_legend = T, name = 'Proportion', column_names_rot = 0, 
                                       column_names_gp = grid::gpar(fontsize = 3),
                                       row_names_gp = grid::gpar(fontsize = 3),
                                       cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                         if(!is.na(correspondance_groups_N_AMR[i, j])){
                                           grid.text(round(correspondance_groups_N_AMR[i, j], digits = 2), x, y, gp = gpar(fontsize = 4))
                                         }
                                       })
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_AMR_20251007_with_numbers.pdf', width = 3/2.54, height = 2.8/2.54)
par(oma = c(2,2,2,2), mar = c(0,0,0,0), mgp = c(2,0.5,0), mfrow = c(4,4),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
plot_AMR_compaarison_numbers
dev.off()

colnames(correspondance_groups_N_AMR) = c('Az', 'CIPRO', 'ESBL')
write.csv(correspondance_groups_N_AMR, '2_analysis_index/2_find_index_groups/Heatmap_SSonnei_AMR_20251007_with_numbers_data.csv')
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs STN
#######################################################################################################################################
dataset_with_nodes$STN = as.character(dataset_with_nodes$STN)
dataset_with_nodes$STN = factor(dataset_with_nodes$STN)

grouping1 = levels(as.factor(dataset_with_nodes$STN))
grouping2 = levels(as.factor(dataset_with_nodes$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes$STN == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  # t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
correspondance_groups_N[,3] = correspondance_groups_N[,3]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Heatmap_SSonnei_STN_20251007_row_data.csv')

library(fossil)
library(circlize)
library(ComplexHeatmap)
correspondance_groups_N = correspondance_groups_N[nrow(correspondance_groups_N):1,]
plot_STN_comparison_numbers = Heatmap((correspondance_groups_N), col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('white', RColorBrewer::brewer.pal(name = 'Reds', n = 9)[c(2,3,6,9)])),
                                       heatmap_legend_param = list(at = c(0, 0.5, 1), 
                                                                   title_gp = gpar(fontsize = 3),
                                                                   labels_gp = gpar(fontsize = 2),
                                                                   legend_height = unit(0.25, "cm"),
                                                                   legend_width = unit(0.025, "cm")), 
                                       show_row_dend = F, show_column_dend = F, 
                                       cluster_rows = F, cluster_columns = F,
                                       show_heatmap_legend = T, name = 'Proportion', column_names_rot = 0, 
                                       column_names_gp = grid::gpar(fontsize = 3),
                                       row_names_gp = grid::gpar(fontsize = 3),
                                       cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                         if(!is.na(correspondance_groups_N[i, j])){
                                           grid.text(round(correspondance_groups_N[i, j], digits = 2), x, y, gp = gpar(fontsize = 4))
                                         }
                                       })

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_STN_20251007_row_with_numbers.pdf', width = 3/2.54, height = 2.8/2.54)
par(oma = c(2,2,2,2), mar = c(0,0,0,0), mgp = c(2,0.5,0), mfrow = c(4,4),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
plot_STN_comparison_numbers
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs pMSM
#######################################################################################################################################
dataset_with_nodes$pmsmdemdef = as.character(dataset_with_nodes$pmsmdemdef)
dataset_with_nodes$pmsmdemdef = factor(dataset_with_nodes$pmsmdemdef)

grouping1 = levels(as.factor(dataset_with_nodes$pmsmdemdef))
grouping2 = levels(as.factor(dataset_with_nodes$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes$pmsmdemdef == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  # t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
correspondance_groups_N[,3] = correspondance_groups_N[,3]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Heatmap_SSonnei_pMSM_20251007_row_with_numbers_data.csv')

library(fossil)
library(circlize)
library(ComplexHeatmap)
correspondance_groups_N = correspondance_groups_N[nrow(correspondance_groups_N):1,]
plot_pMSM_comparison_numbers = Heatmap((correspondance_groups_N), col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('white', RColorBrewer::brewer.pal(name = 'Reds', n = 9)[c(2,3,6,9)])),
                                        heatmap_legend_param = list(at = c(0, 0.5, 1), 
                                                                    title_gp = gpar(fontsize = 3),
                                                                    labels_gp = gpar(fontsize = 2),
                                                                    legend_height = unit(0.25, "cm"),
                                                                    legend_width = unit(0.025, "cm")), 
                                        show_row_dend = F, show_column_dend = F, 
                                        cluster_rows = F, cluster_columns = F,
                                        show_heatmap_legend = T, name = 'Proportion', column_names_rot = 0, 
                                        column_names_gp = grid::gpar(fontsize = 3),
                                        row_names_gp = grid::gpar(fontsize = 3),
                                        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                          if(!is.na(correspondance_groups_N[i, j])){
                                            grid.text(round(correspondance_groups_N[i, j], digits = 2), x, y, gp = gpar(fontsize = 4))
                                          }
                                        })

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_pMSM_20251007_row_with_numbers.pdf', width = 3/2.54, height = 2.8/2.54)
par(oma = c(2,2,2,2), mar = c(0,0,0,0), mgp = c(2,0.5,0), mfrow = c(4,4),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
plot_pMSM_comparison_numbers
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs Genotype_simplified
#######################################################################################################################################
dataset_with_nodes$Genotype_simplified = as.character(dataset_with_nodes$Genotype_simplified)
dataset_with_nodes$Genotype_simplified = factor(dataset_with_nodes$Genotype_simplified)

grouping1 = levels(as.factor(dataset_with_nodes$Genotype_simplified))
grouping2 = levels(as.factor(dataset_with_nodes$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes$Genotype_simplified == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
for(i in 1:ncol(correspondance_groups_N)){
  correspondance_groups_N[,i] = correspondance_groups_N[,i]/tmp
}
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = TRUE,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups#[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N#[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_genotypes_simplified_row_20251007.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5)
dev.off()

colnames(correspondance_groups_N) = paste0('_', colnames(correspondance_groups_N))
write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Heatmap_SSonnei_genotypes_simplified_row_20251007_data.csv')
#######################################################################################################################################

#######################################################################################################################################
## Table, dominant simplified genotype per group
#######################################################################################################################################
dataset_with_nodes$Genotype_simplified = as.character(dataset_with_nodes$Genotype_simplified)
dataset_with_nodes$Genotype_simplified = factor(dataset_with_nodes$Genotype_simplified)

grouping1 = levels(as.factor(dataset_with_nodes$Genotype_simplified))
grouping2 = levels(as.factor(dataset_with_nodes$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes$Genotype_simplified == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
for(i in 1:ncol(correspondance_groups_N)){
  correspondance_groups_N[,i] = correspondance_groups_N[,i]/tmp
}
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = TRUE,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups#[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N#[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)

t = paste0('   ', colnames(correspondance_groups_N)[apply(correspondance_groups_N, MARGIN = 1, function(x)which.max(x))])
names(t) = 1:length(t)

write.csv(t, file = '2_analysis_index/2_find_index_groups/Table_most_prevalent_simplified_genotype_group_20251007.csv')
#######################################################################################################################################

#######################################################################################################################################
## Table, dominant true genotype per group
#######################################################################################################################################
dataset_with_nodes$Genotype = as.character(dataset_with_nodes$Genotype)
dataset_with_nodes$Genotype = factor(dataset_with_nodes$Genotype)

grouping1 = levels(as.factor(dataset_with_nodes$Genotype))
grouping2 = levels(as.factor(dataset_with_nodes$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes$Genotype == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
for(i in 1:ncol(correspondance_groups_N)){
  correspondance_groups_N[,i] = correspondance_groups_N[,i]/tmp
}
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = TRUE,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups#[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N#[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)

t_name = paste0('   ', colnames(correspondance_groups_N)[apply(correspondance_groups_N, MARGIN = 1, function(x)which.max(x))])
t_name = matrix(t_name, nrow = length(t_name))
rownames(t_name) = 1:length(t_name)

t_percent = apply(correspondance_groups_N, MARGIN = 1, function(x)max(x))
t_percent = matrix(t_percent, nrow = length(t_name))
rownames(t_percent) = 1:length(t_percent)

plot_main_genotype = Heatmap((t_percent),  col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('white', RColorBrewer::brewer.pal(name = 'Reds', n = 9)[c(2,3,6,9)])),
                                       heatmap_legend_param = list(at = c(0, 0.5, 1), 
                                                                   title_gp = gpar(fontsize = 3),
                                                                   labels_gp = gpar(fontsize = 2),
                                                                   legend_height = unit(0.25, "cm"),
                                                                   legend_width = unit(0.025, "cm")), 
                                       show_row_dend = F, show_column_dend = F, 
                                       cluster_rows = F, cluster_columns = F,
                                       show_heatmap_legend = T, name = 'Proportion', column_names_rot = 0, 
                                       column_names_gp = grid::gpar(fontsize = 3),
                                       row_names_gp = grid::gpar(fontsize = 3),
                                       cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                         if(!is.na(t_name[i, j])){
                                           grid.text(t_name[i, j], x, y, gp = gpar(fontsize = 10))
                                         }
                                       })
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_most_prevalent_genotype_group_20251007.pdf', width = 3/2.54, height = 2.8/2.54)
par(oma = c(2,2,2,2), mar = c(0,0,0,0), mgp = c(2,0.5,0), mfrow = c(4,4),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
plot_main_genotype
dev.off()

write.csv(plot_main_genotype, file = '2_analysis_index/2_find_index_groups/Table_most_prevalent_genotype_group_20251007.csv')
#######################################################################################################################################


#######################################################################################################################################
## ARI genotype simplified
#######################################################################################################################################
library(fossil)
idx = which(!is.na(dataset_with_nodes$Genotype_simplified))
group1 = as.numeric(dataset_with_nodes$Genotype_simplified[idx])
group2 = as.numeric(dataset_with_nodes$groups[idx])
rand.treeannotator_partitions = fossil::adj.rand.index(group1, group2)
rand.treeannotator_partitions
# 0.6496448
#######################################################################################################################################





































## Supp, to delete


########################################################################################################################################
## Plot comparison - full data
########################################################################################################################################
Genotypes = matrix(as.character(dataset_with_nodes$Genotype), ncol = 1)
colnames(Genotypes) = 'groups'
rownames(Genotypes) = dataset_with_nodes$name_seq
plot_tree_Genotypes <- ggtree(tree, mrsd=lubridate::date_decimal(max(dataset_with_nodes$time)), 
                           size = 0.1,
                            aes(color = as.character(dataset_with_nodes$Genotype))) + 
  scale_color_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  theme_tree2()
plot_tree_Genotypes = gheatmap(plot_tree_Genotypes, Genotypes, offset=0.1, width=0.10, 
                    colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

groups = matrix(dataset_with_nodes$groups, ncol = 1)
colnames(groups) = 'groups'
rownames(groups) = dataset_with_nodes$name_seq
cols = as.character(colors_groups)
names(cols) = as.character(1:max(as.numeric(name_groups)))
plot_tree_groups <- ggtree(tree, mrsd=lubridate::date_decimal(max(dataset_with_nodes$time)), 
                           size = 0.1,
                           aes(color = as.character(split$groups))) + 
  scale_color_manual(values = cols, na.value = 'grey90')+
  theme_tree2()
plot_tree_groups = gheatmap(plot_tree_groups, groups, offset=0.1, width=0.10, 
                                      colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = cols, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

comparison = plot_grid(plot_tree_groups, plot_tree_Genotypes, 
          rel_widths = c(1, 1), labels = '', ncol = 2)

ggsave(filename = '2_analysis_index/2_find_index_groups/Comparison_groups_genotypes_full_20250522.pdf', plot = comparison, device = 'pdf', scale = 1,
       width = 10, height = 10, units = 'cm')
########################################################################################################################################

########################################################################################################################################
## Plot comparison - full data, with legend 
########################################################################################################################################
Genotypes = matrix(as.character(dataset_with_nodes$Genotype), ncol = 1)
colnames(Genotypes) = 'groups'
rownames(Genotypes) = dataset_with_nodes$name_seq
plot_tree_Genotypes <- ggtree(tree, mrsd=lubridate::date_decimal(max(dataset_with_nodes$time)), 
                              size = 0.1,
                              aes(color = as.character(dataset_with_nodes$Genotype))) + 
  scale_color_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  theme_tree2()
plot_tree_Genotypes = gheatmap(plot_tree_Genotypes, Genotypes, offset=0.1, width=0.10, 
                               colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0, 0.3))

groups = matrix(dataset_with_nodes$groups, ncol = 1)
colnames(groups) = 'groups'
rownames(groups) = dataset_with_nodes$name_seq
cols = as.character(colors_groups)
names(cols) = as.character(1:max(as.numeric(name_groups)))
plot_tree_groups <- ggtree(tree, mrsd=lubridate::date_decimal(max(dataset_with_nodes$time)), 
                           size = 0.1,
                           aes(color = as.character(split$groups))) + 
  scale_color_manual(values = cols, na.value = 'grey90')+
  theme_tree2()
plot_tree_groups = gheatmap(plot_tree_groups, groups, offset=0.1, width=0.10, 
                            colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = cols, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

comparison = plot_grid(plot_tree_groups, plot_tree_Genotypes, 
                       rel_widths = c(1, 1), labels = '', ncol = 2)

ggsave(filename = '2_analysis_index/2_find_index_groups/Comparison_groups_genotypes_full_20250522_with_legend.pdf', plot = comparison, device = 'pdf', scale = 1,
       width = 10, height = 10, units = 'cm')
########################################################################################################################################

########################################################################################################################################
## Plot comparison - STN data
########################################################################################################################################
Genotypes_STN = matrix(as.character(dataset_with_nodes_STN$Genotype), ncol = 1)
colnames(Genotypes_STN) = 'groups'
rownames(Genotypes_STN) = dataset_with_nodes_STN$name_seq
plot_tree_Genotypes_STN <- ggtree(tree_STN, mrsd=lubridate::date_decimal(max(dataset_with_nodes_STN$time)), 
                              size = 0.15,
                              aes(color = as.character(dataset_with_nodes_STN$Genotype))) + 
  scale_color_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  theme_tree2()
plot_tree_Genotypes_STN = gheatmap(plot_tree_Genotypes_STN, Genotypes_STN, offset=0.1, width=0.10, 
                               colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

groups_STN = matrix(dataset_with_nodes_STN$groups, ncol = 1)
colnames(groups_STN) = 'groups'
rownames(groups_STN) = dataset_with_nodes_STN$name_seq
cols = as.character(colors_groups_STN)
names(cols) = as.character(1:max(as.numeric(name_groups_STN)))
plot_tree_groups_STN <- ggtree(tree_STN, mrsd=lubridate::date_decimal(max(dataset_with_nodes_STN$time)), 
                           size = 0.15,
                           aes(color = as.character(split_STN$groups))) + 
  scale_color_manual(values = cols, na.value = 'grey90')+
  theme_tree2()
plot_tree_groups_STN = gheatmap(plot_tree_groups_STN, groups_STN, offset=0.1, width=0.10, 
                            colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = cols, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

comparison = plot_grid(plot_tree_groups_STN, plot_tree_Genotypes_STN, 
                       rel_widths = c(1, 1), labels = '', ncol = 2)

ggsave(filename = '2_analysis_index/2_find_index_groups/Comparison_groups_genotypes_STN_20250405.pdf', plot = comparison, device = 'pdf', scale = 1,
       width = 7, height = 7, units = 'cm')
########################################################################################################################################

########################################################################################################################################
## Plot comparison - non STN data
########################################################################################################################################
Genotypes_nonSTN = matrix(as.character(dataset_with_nodes_nonSTN$Genotype), ncol = 1)
colnames(Genotypes_nonSTN) = 'groups'
rownames(Genotypes_nonSTN) = dataset_with_nodes_nonSTN$name_seq
plot_tree_Genotypes_nonSTN <- ggtree(tree_nonSTN, mrsd=lubridate::date_decimal(max(dataset_with_nodes_nonSTN$time)), 
                                  size = 0.15,
                                  aes(color = as.character(dataset_with_nodes_nonSTN$Genotype))) + 
  scale_color_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  # geom_text(aes(label=node), size = 3)+ ## to check for rotation
  theme_tree2()
plot_tree_Genotypes_nonSTN = rotate(plot_tree_Genotypes_nonSTN, 1656) ## to put yellow clade up
plot_tree_Genotypes_nonSTN = rotate(plot_tree_Genotypes_nonSTN, 1657)
plot_tree_Genotypes_nonSTN = rotate(plot_tree_Genotypes_nonSTN, 2435)
plot_tree_Genotypes_nonSTN = gheatmap(plot_tree_Genotypes_nonSTN, Genotypes_nonSTN, offset=0.1, width=0.10, 
                                   colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

groups_nonSTN = matrix(dataset_with_nodes_nonSTN$groups, ncol = 1)
colnames(groups_nonSTN) = 'groups'
rownames(groups_nonSTN) = dataset_with_nodes_nonSTN$name_seq
cols = as.character(colors_groups_nonSTN)
names(cols) = as.character(1:max(as.numeric(name_groups_nonSTN)))
plot_tree_groups_nonSTN <- ggtree(tree_nonSTN, mrsd=lubridate::date_decimal(max(dataset_with_nodes_nonSTN$time)), 
                               size = 0.15,
                               aes(color = as.character(split_nonSTN$groups))) + 
  scale_color_manual(values = cols, na.value = 'grey90')+
  # geom_text(aes(label=node), size = 3)+ ## to check for rotation
  theme_tree2()
plot_tree_groups_nonSTN = rotate(plot_tree_groups_nonSTN, 1656) ## to put yellow clade up
plot_tree_groups_nonSTN = rotate(plot_tree_groups_nonSTN, 1657)
plot_tree_groups_nonSTN = rotate(plot_tree_groups_nonSTN, 2435)
plot_tree_groups_nonSTN = gheatmap(plot_tree_groups_nonSTN, groups_nonSTN, offset=0.1, width=0.10, 
                                colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = cols, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

comparison = plot_grid(plot_tree_groups_nonSTN, plot_tree_Genotypes_nonSTN, 
                       rel_widths = c(1, 1), labels = '', ncol = 2)

ggsave(filename = '2_analysis_index/2_find_index_groups/Comparison_groups_genotypes_nonSTN_20250405.pdf', plot = comparison, device = 'pdf', scale = 1,
       width = 7, height = 7, units = 'cm')
########################################################################################################################################

########################################################################################################################################
## Plot comparison - MSM pandemic
########################################################################################################################################
Genotypes_MSM_pandemic = matrix(as.character(dataset_with_nodes_MSM_pandemic$Genotype), ncol = 1)
colnames(Genotypes_MSM_pandemic) = 'groups'
rownames(Genotypes_MSM_pandemic) = dataset_with_nodes_MSM_pandemic$name_seq
plot_tree_Genotypes_MSM_pandemic <- ggtree(tree_MSM_pandemic, mrsd=lubridate::date_decimal(max(dataset_with_nodes_MSM_pandemic$time)), 
                                     size = 0.15,
                                     aes(color = as.character(dataset_with_nodes_MSM_pandemic$Genotype))) + 
  scale_color_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  # geom_text(aes(label=node), size = 3)+ ## to check for rotation
  theme_tree2()
plot_tree_Genotypes_MSM_pandemic = gheatmap(plot_tree_Genotypes_MSM_pandemic, Genotypes_MSM_pandemic, offset=0.1, width=0.10, 
                                      colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = colors_groups_full_Genotype, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

groups_MSM_pandemic = matrix(dataset_with_nodes_MSM_pandemic$groups, ncol = 1)
colnames(groups_MSM_pandemic) = 'groups'
rownames(groups_MSM_pandemic) = dataset_with_nodes_MSM_pandemic$name_seq
cols = as.character(colors_groups_MSM_pandemic)
names(cols) = as.character(1:max(as.numeric(name_groups_MSM_pandemic)))
plot_tree_groups_MSM_pandemic <- ggtree(tree_MSM_pandemic, mrsd=lubridate::date_decimal(max(dataset_with_nodes_MSM_pandemic$time)), 
                                  size = 0.15,
                                  aes(color = as.character(split_MSM_pandemic$groups))) + 
  scale_color_manual(values = cols, na.value = 'grey90')+
  # geom_text(aes(label=node), size = 3)+ ## to check for rotation
  theme_tree2()
plot_tree_groups_MSM_pandemic = gheatmap(plot_tree_groups_MSM_pandemic, groups_MSM_pandemic, offset=0.1, width=0.10, 
                                   colnames=FALSE, legend_title="Group", color=NA) +
  scale_fill_manual(values = cols, na.value = 'grey90')+
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3))+
  theme(legend.position = 'none')

comparison = plot_grid(plot_tree_groups_MSM_pandemic, plot_tree_Genotypes_MSM_pandemic, 
                       rel_widths = c(1, 1), labels = '', ncol = 2)

ggsave(filename = '2_analysis_index/2_find_index_groups/Comparison_groups_genotypes_MSM_pandemic_20250405.pdf', plot = comparison, device = 'pdf', scale = 1,
       width = 7, height = 7, units = 'cm')
########################################################################################################################################




## HEATMAPS FULL DATASET
#######################################################################################################################################
## Plot heatmap groups vs ESBLR
#######################################################################################################################################
dataset_with_nodes$ESBLR = as.character(dataset_with_nodes$ESBLR)
dataset_with_nodes$ESBLR = factor(dataset_with_nodes$ESBLR)

grouping1 = levels(as.factor(dataset_with_nodes$ESBLR))
grouping2 = levels(as.factor(dataset_with_nodes$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes$ESBLR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,  #TRUE,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_ESBLR_20250405.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('ESBLR-1', 'ESBLR-0'))
dev.off()

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_ESBLR_20250405.csv')

auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
correspondance_groups_N_EBSL = correspondance_groups_N
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_ESBLR_20250405_row.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('ESBLR-1', 'ESBLR-0'))
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs CiprofloxacinR
#######################################################################################################################################
dataset_with_nodes$CiprofloxacinR = as.character(dataset_with_nodes$CiprofloxacinR)
dataset_with_nodes$CiprofloxacinR = factor(dataset_with_nodes$CiprofloxacinR)

grouping1 = levels(as.factor(dataset_with_nodes$CiprofloxacinR))
grouping2 = levels(as.factor(dataset_with_nodes$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes$CiprofloxacinR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

correspondance_groups_N_CIPRO = correspondance_groups_N
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Cipro_20250405.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('CIPR-1', 'CIPR-0'))
dev.off()

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_Cipro_20250405.csv')

auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
correspondance_groups_N_CIPRO = correspondance_groups_N
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Cipro_20250405_row.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('CIPR-1', 'CIPR-0'))
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs Az
#######################################################################################################################################
dataset_with_nodes$AzithromycinR = as.character(dataset_with_nodes$AzithromycinR)
dataset_with_nodes$AzithromycinR = factor(dataset_with_nodes$AzithromycinR)

grouping1 = levels(as.factor(dataset_with_nodes$AzithromycinR))
grouping2 = levels(as.factor(dataset_with_nodes$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes$AzithromycinR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  # t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}


pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Az_20250405.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('AzR-1', 'AzR-0'))
dev.off()

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_Az_20250405.csv')

auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
correspondance_groups_N_Az = correspondance_groups_N
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Az_20250405_row.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('AzR-1', 'AzR-0'))
dev.off()
#######################################################################################################################################


#######################################################################################################################################
## Plot heatmap STN vs pMSM
#######################################################################################################################################
dataset_with_nodes$pmsmdemdef = as.character(dataset_with_nodes$pmsmdemdef)
dataset_with_nodes$pmsmdemdef = factor(dataset_with_nodes$pmsmdemdef)

grouping1 = levels(as.factor(dataset_with_nodes$pmsmdemdef))
grouping2 = levels(as.factor(dataset_with_nodes$STN))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes$pmsmdemdef == grouping1[i])
  t = table((as.character(dataset_with_nodes$STN[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  # t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp

tmp = colSums(correspondance_groups_N)
correspondance_groups_N[1,] = correspondance_groups_N[1,]/tmp
correspondance_groups_N[2,] = correspondance_groups_N[2,]/tmp
correspondance_groups_N[3,] = correspondance_groups_N[3,]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)

correspondance_groups_N = correspondance_groups_N[nrow(correspondance_groups_N):1,]
colnames(correspondance_groups_N) = c('pMSM', 'non-pMSM')
rownames(correspondance_groups_N) = c('Travel', 'STN', 'non-STN')
par(oma = c(2,2,2,2), mar = c(2,2,2,2), mgp = c(2,0.5,0), mfrow = c(4,4),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
Heatmap((correspondance_groups_N), col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c('white', RColorBrewer::brewer.pal(name = 'Reds', n = 9)[c(2,3,6,9)])),
          heatmap_legend_param = list(at = c(0, 0.5, 1), 
                                      title_gp = gpar(fontsize = 3),
                                      labels_gp = gpar(fontsize = 2),
                                      legend_height = unit(0.25, "cm"),
                                      legend_width = unit(0.025, "cm")), 
          show_row_dend = F, show_column_dend = F, 
          cluster_rows = F, cluster_columns = F,
          show_heatmap_legend = T, name = 'Proportion', column_names_rot = 90, 
          column_names_gp = grid::gpar(fontsize = 8),
          row_names_gp = grid::gpar(fontsize = 8),
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            if(!is.na(correspondance_groups_N[i, j])){
              grid.text(round(correspondance_groups_N[i, j], digits = 2), x, y, gp = gpar(fontsize = 8))
            }
          })
#######################################################################################################################################









## HEATMAPS STN
dataset_with_nodes_heatmap = dataset_with_nodes_reconstructed[which(dataset_with_nodes_reconstructed$STN == "pMSM"),]
#######################################################################################################################################
## Plot heatmap groups vs ESBLR
#######################################################################################################################################
dataset_with_nodes_heatmap$ESBLR = as.character(dataset_with_nodes_heatmap$ESBLR)
dataset_with_nodes_heatmap$ESBLR = factor(dataset_with_nodes_heatmap$ESBLR)

grouping1 = levels(as.factor(dataset_with_nodes_heatmap$ESBLR))
grouping2 = levels(as.factor(dataset_with_nodes_heatmap$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_heatmap$ESBLR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_heatmap$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,  #TRUE,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_ESBLR_20250113_STN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('ESBLR-1', 'ESBLR-0'))
dev.off()

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_ESBLR_20250113_STN.csv')


auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_ESBLR_20250113_row_STN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('ESBLR-1', 'ESBLR-0'))
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs CiprofloxacinR
#######################################################################################################################################
dataset_with_nodes_heatmap$CiprofloxacinR = as.character(dataset_with_nodes_heatmap$CiprofloxacinR)
dataset_with_nodes_heatmap$CiprofloxacinR = factor(dataset_with_nodes_heatmap$CiprofloxacinR)

grouping1 = levels(as.factor(dataset_with_nodes_heatmap$CiprofloxacinR))
grouping2 = levels(as.factor(dataset_with_nodes_heatmap$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_heatmap$CiprofloxacinR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_heatmap$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Cipro_20250113_STN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('CIPR-1', 'CIPR-0'))
dev.off()

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_Cipro_20250113_STN.csv')


auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Cipro_20250113_row_STN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('CIPR-1', 'CIPR-0'))
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs Az
#######################################################################################################################################
dataset_with_nodes_heatmap$AzithromycinR = as.character(dataset_with_nodes_heatmap$AzithromycinR)
dataset_with_nodes_heatmap$AzithromycinR = factor(dataset_with_nodes_heatmap$AzithromycinR)

grouping1 = levels(as.factor(dataset_with_nodes_heatmap$AzithromycinR))
grouping2 = levels(as.factor(dataset_with_nodes_heatmap$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_heatmap$AzithromycinR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_heatmap$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Az_20250113_STN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('AzR-1', 'AzR-0'))
dev.off()

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_Az_20250113_STN.csv')


auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Az_20250113_row_STN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('AzR-1', 'AzR-0'))
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs Genotypes
#######################################################################################################################################
dataset_with_nodes_STN$Genotype = as.character(dataset_with_nodes_STN$Genotype)
dataset_with_nodes_STN$Genotype = factor(dataset_with_nodes_STN$Genotype)

grouping1 = levels(as.factor(dataset_with_nodes_STN$Genotype))
grouping2 = levels(as.factor(dataset_with_nodes_STN$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_STN$Genotype == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_STN$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = TRUE,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_genotypes_20250405_STN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5)
dev.off()

colnames(correspondance_groups_N) = paste0('_', colnames(correspondance_groups_N))
write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_genotypes_20250405_STN.csv')
#######################################################################################################################################



## HEATMAPS non-STN
dataset_with_nodes_heatmap = dataset_with_nodes_reconstructed[which(dataset_with_nodes_reconstructed$STN == "non-pMSM"),]
#######################################################################################################################################
## Plot heatmap groups vs ESBLR
#######################################################################################################################################
dataset_with_nodes_heatmap$ESBLR = as.character(dataset_with_nodes_heatmap$ESBLR)
dataset_with_nodes_heatmap$ESBLR = factor(dataset_with_nodes_heatmap$ESBLR)

grouping1 = levels(as.factor(dataset_with_nodes_heatmap$ESBLR))
grouping2 = levels(as.factor(dataset_with_nodes_heatmap$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_heatmap$ESBLR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_heatmap$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,  #TRUE,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_ESBLR_20250113_nonSTN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('ESBLR-1', 'ESBLR-0'))
dev.off()

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_ESBLR_20250113_nonSTN.csv')


auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_ESBLR_20250113_row_nonSTN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('ESBLR-1', 'ESBLR-0'))
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs CiprofloxacinR
#######################################################################################################################################
dataset_with_nodes_heatmap$CiprofloxacinR = as.character(dataset_with_nodes_heatmap$CiprofloxacinR)
dataset_with_nodes_heatmap$CiprofloxacinR = factor(dataset_with_nodes_heatmap$CiprofloxacinR)

grouping1 = levels(as.factor(dataset_with_nodes_heatmap$CiprofloxacinR))
grouping2 = levels(as.factor(dataset_with_nodes_heatmap$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_heatmap$CiprofloxacinR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_heatmap$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Cipro_20250113_nonSTN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('CIPR-1', 'CIPR-0'))
dev.off()

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_Cipro_20250113_nonSTN.csv')


auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Cipro_20250113_row_nonSTN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('CIPR-1', 'CIPR-0'))
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs Az
#######################################################################################################################################
dataset_with_nodes_heatmap$AzithromycinR = as.character(dataset_with_nodes_heatmap$AzithromycinR)
dataset_with_nodes_heatmap$AzithromycinR = factor(dataset_with_nodes_heatmap$AzithromycinR)

grouping1 = levels(as.factor(dataset_with_nodes_heatmap$AzithromycinR))
grouping2 = levels(as.factor(dataset_with_nodes_heatmap$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_heatmap$AzithromycinR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_heatmap$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Az_20250113_nonSTN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('AzR-1', 'AzR-0'))
dev.off()

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_Az_20250113_nonSTN.csv')


auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Az_20250113_row_nonSTN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('AzR-1', 'AzR-0'))
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs Genotypes
#######################################################################################################################################
dataset_with_nodes_nonSTN$Genotype = as.character(dataset_with_nodes_nonSTN$Genotype)
dataset_with_nodes_nonSTN$Genotype = factor(dataset_with_nodes_nonSTN$Genotype)

grouping1 = levels(as.factor(dataset_with_nodes_nonSTN$Genotype))
grouping2 = levels(as.factor(dataset_with_nodes_nonSTN$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_nonSTN$Genotype == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_nonSTN$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = TRUE,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_genotypes_20250405_nonSTN.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5)
dev.off()

colnames(correspondance_groups_N) = paste0('_', colnames(correspondance_groups_N))
write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_genotypes_20250405_nonSTN.csv')
#######################################################################################################################################




## HEATMAPS MSM_pandemic
#######################################################################################################################################
## Plot heatmap groups vs ESBLR
#######################################################################################################################################
dataset_with_nodes_MSM_pandemic$ESBLR = as.character(dataset_with_nodes_MSM_pandemic$ESBLR)
dataset_with_nodes_MSM_pandemic$ESBLR = factor(dataset_with_nodes_MSM_pandemic$ESBLR)

grouping1 = levels(as.factor(dataset_with_nodes_MSM_pandemic$ESBLR))
grouping2 = levels(as.factor(dataset_with_nodes_MSM_pandemic$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_MSM_pandemic$ESBLR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_MSM_pandemic$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,  #TRUE,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_ESBLR_20250113_MSM_pandemic.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('ESBLR-1', 'ESBLR-0'))
dev.off()

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_ESBLR_20250113_MSM_pandemic.csv')


auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_ESBLR_20250113_row_MSM_pandemic.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('ESBLR-1', 'ESBLR-0'))
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs CiprofloxacinR
#######################################################################################################################################
dataset_with_nodes_MSM_pandemic$CiprofloxacinR = as.character(dataset_with_nodes_MSM_pandemic$CiprofloxacinR)
dataset_with_nodes_MSM_pandemic$CiprofloxacinR = factor(dataset_with_nodes_MSM_pandemic$CiprofloxacinR)

grouping1 = levels(as.factor(dataset_with_nodes_MSM_pandemic$CiprofloxacinR))
grouping1 = c("0", "1")
grouping2 = levels(as.factor(dataset_with_nodes_MSM_pandemic$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_MSM_pandemic$CiprofloxacinR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_MSM_pandemic$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Cipro_20250113_MSM_pandemic.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('CIPR-1', 'CIPR-0'))
dev.off()

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_Cipro_20250113_MSM_pandemic.csv')


auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Cipro_20250113_row_MSM_pandemic.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('CIPR-1', 'CIPR-0'))
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs Az
#######################################################################################################################################
dataset_with_nodes_MSM_pandemic$AzithromycinR = as.character(dataset_with_nodes_MSM_pandemic$AzithromycinR)
dataset_with_nodes_MSM_pandemic$AzithromycinR = factor(dataset_with_nodes_MSM_pandemic$AzithromycinR)

grouping1 = levels(as.factor(dataset_with_nodes_MSM_pandemic$AzithromycinR))
grouping2 = levels(as.factor(dataset_with_nodes_MSM_pandemic$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_MSM_pandemic$AzithromycinR == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_MSM_pandemic$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = F,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Az_20250113_MSM_pandemic.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('AzR-1', 'AzR-0'))
dev.off()

write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_Az_20250113_MSM_pandemic.csv')


auto_groups = function(){
  heatmap(correspondance_groups_N, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}
pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_Az_20250113_row_MSM_pandemic.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_N, 
        Colv = NA, Rowv = NA, 
        scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5,
        labCol = c('AzR-1', 'AzR-0'))
dev.off()
#######################################################################################################################################

#######################################################################################################################################
## Plot heatmap groups vs Genotypes
#######################################################################################################################################
dataset_with_nodes_MSM_pandemic$Genotype = as.character(dataset_with_nodes_MSM_pandemic$Genotype)
dataset_with_nodes_MSM_pandemic$Genotype = factor(dataset_with_nodes_MSM_pandemic$Genotype)

grouping1 = levels(as.factor(dataset_with_nodes_MSM_pandemic$Genotype))
grouping2 = levels(as.factor(dataset_with_nodes_MSM_pandemic$groups))
correspondance_groups = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups) = grouping1
rownames(correspondance_groups) = grouping2
correspondance_groups_N = matrix(0, ncol = length(grouping1), nrow = length(grouping2))
colnames(correspondance_groups_N) = grouping1
rownames(correspondance_groups_N) = grouping2
for(i in 1:length(grouping1)){
  idx = which(dataset_with_nodes_MSM_pandemic$Genotype == grouping1[i])
  t = table(as.numeric(as.character(dataset_with_nodes_MSM_pandemic$groups[idx])))
  correspondance_groups_N[match(names(t), grouping2), i] = t
  t = t/sum(t)
  correspondance_groups[match(names(t), grouping2), i] = t
}
tmp = rowSums(correspondance_groups_N)
correspondance_groups_N[,1] = correspondance_groups_N[,1]/tmp
correspondance_groups_N[,2] = correspondance_groups_N[,2]/tmp
orders = slanter::slanted_orders(
  data = correspondance_groups,
  order_rows = F,
  order_cols = TRUE,
  squared_order = TRUE,
  same_order = FALSE,
  discount_outliers = F,
  max_spin_count = 100
)
correspondance_groups_auto = correspondance_groups[orders$rows,rev(orders$cols)]
correspondance_groups_N = correspondance_groups_N[orders$rows,rev(orders$cols)]
coul <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)
auto_groups = function(){
  heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)))
}

pdf('2_analysis_index/2_find_index_groups/Heatmap_SSonnei_genotypes_20250113_MSM_pandemic.pdf', width = 5/2.54, height = 5/2.54)
par(oma = c(1,0,1,0), mar = c(1,0,1,0), mgp = c(2,0.5,0), mfrow = c(1,1),
    cex.axis = 0.01, cex.lab = 0.01, cex.main = 0.05)
heatmap(correspondance_groups_auto, Colv = NA, Rowv = NA, scale="none", col = c('white' , colorRampPalette(RColorBrewer::brewer.pal(name = 'Reds', n = 9))(25)), 
        cex.main = 0.05, margins = c(1,1), cexRow = 0.5, cexCol = 0.5)
dev.off()

colnames(correspondance_groups_N) = paste0('_', colnames(correspondance_groups_N))
write.csv(correspondance_groups_N, file = '2_analysis_index/2_find_index_groups/Data_heatmap_SSonnei_genotypes_20250113_MSM_pandemic.csv')
#######################################################################################################################################






#######################################################################################################################################
library(fossil)
idx = which(!is.na(dataset_with_nodes$Genotype_simplified))
group1 = as.numeric(dataset_with_nodes$Genotype_simplified[idx])
group2 = as.numeric(dataset_with_nodes$groups[idx])
rand.treeannotator_partitions = fossil::adj.rand.index(group1, group2)
rand.treeannotator_partitions
# 0.6496448
#######################################################################################################################################

idx = which(!is.na(dataset_with_nodes_nonSTN$Genotype))
group1 = as.numeric(dataset_with_nodes_nonSTN$Genotype[idx])
group2 = as.numeric(dataset_with_nodes_nonSTN$groups[idx])
rand.treeannotator_partitions = fossil::adj.rand.index(group1, group2)
rand.treeannotator_partitions
# 0.5123906

idx = which(!is.na(dataset_with_nodes_STN$Genotype))
group1 = as.numeric(dataset_with_nodes_STN$Genotype[idx])
group2 = as.numeric(dataset_with_nodes_STN$groups[idx])
rand.treeannotator_partitions = fossil::adj.rand.index(group1, group2)
rand.treeannotator_partitions
# 0.9699635
#######################################################################################################################################


