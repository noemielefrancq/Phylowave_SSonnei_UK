########################################################################################################################################
## 1.1: Sort datasets and compute index, output: Initial_index_computation_and_parameters_20251218.Rdata
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

########################################################################################################################################
## Useful functions
########################################################################################################################################
## Load index functions (from phylowave githhub https://github.com/noemielefrancq/Phylowave_Learning-fitness-dynamics-pathogens-in-phylogenies)
source('~/Dropbox/Projects/202312_Index_paper_methods_Sc2_H3N2_BP_TB/Phylowave_Learning-fitness-dynamics-pathogens-in-phylogenies/2_Functions/2_1_Index_computation_20240909.R')
source('~/Dropbox/Projects/202312_Index_paper_methods_Sc2_H3N2_BP_TB/Phylowave_Learning-fitness-dynamics-pathogens-in-phylogenies/2_Functions/2_2_Lineage_detection_20240909.R')

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

########################################################################################################################################
## Load data & set parameters
########################################################################################################################################
## Set directory
setwd('~/Dropbox/Projects/2025_Phylowave_SSonnei/Phylowave_SSonnei/')

## Tree
bactres = readRDS('1_Data/ssonneigeospatialbacdate030524.rds')
tree = ladderize(bactres$tree, right = F)
tree$node.label = paste0('NODE', 1:tree$Nnode)
dates = bactres$record[1:length(tree$tip.label)]
metadata = read.csv('1_Data/datforphylowave121625.csv') ## Contains only sequences up to Feb 2020

## Names all sequences
names_seqs = tree$tip.label
names_seqs_pMSM = metadata$Accession[which(metadata$Group == 'pMSM')]
names_seqs_nonpMSM = metadata$Accession[which(metadata$Group == 'non-pMSM')]
names_seqs_travelpMSM = metadata$Accession[which(metadata$Group == 'Travel')]

## Trees all sequences
tree = drop.tip(tree, tip = tree$tip.label[which(is.na(match(tree$tip.label, metadata$Accession)))])
tree_pMSM = drop.tip(tree, tip = names_seqs[which(is.na(match(names_seqs, names_seqs_pMSM)))])
tree_nonpMSM = drop.tip(tree, tip = names_seqs[which(is.na(match(names_seqs, names_seqs_nonpMSM)))])
tree_travelpMSM = drop.tip(tree, tip = names_seqs[which(is.na(match(names_seqs, names_seqs_travelpMSM)))])
names_seqs = tree$tip.label
names_seqs_pMSM = tree_pMSM$tip.label
names_seqs_nonpMSM = tree_nonpMSM$tip.label
names_seqs_travelpMSM = tree_travelpMSM$tip.label

## Mutation rate
mu = 6.39E-7 ## in mutations/site/year

## Length genome 
genome_length = 4971117 ## Length genome

## Clean genotype column 
metadata$Genotype[which(metadata$Genotype == 'Uncalled')] = NA
metadata$Genotype[which(metadata$Genotype == 'uncalled')] = NA
metadata$Genotype[which(metadata$Genotype == '2.2999999999999998')] = "2.3"
########################################################################################################################################

########################################################################################################################################
## Preparation data tips
########################################################################################################################################
## Create dataset with names of each sequence, time sampling, prn type etc
dataset_tips = data.frame('ID' = 1:length(names_seqs),
                                 'name_seq' = names_seqs,
                                 'time' = NA)
dataset_tips$time = as.numeric(dataset_tips$time)

## pMSM
dataset_tips_pMSM = data.frame('ID' = 1:length(names_seqs_pMSM),
                              'name_seq' = names_seqs_pMSM,
                              'time' = NA)
dataset_tips_pMSM$time = as.numeric(dataset_tips_pMSM$time)

## Non-pMSM
dataset_tips_nonpMSM = data.frame('ID' = 1:length(names_seqs_nonpMSM),
                                 'name_seq' = names_seqs_nonpMSM,
                                 'time' = NA)
dataset_tips_nonpMSM$time = as.numeric(dataset_tips_nonpMSM$time)

## Travel pMSM
dataset_tips_travelpMSM = data.frame('ID' = 1:length(names_seqs_travelpMSM),
                               'name_seq' = names_seqs_travelpMSM,
                               'time' = NA)
dataset_tips_travelpMSM$time = as.numeric(dataset_tips_travelpMSM$time)
########################################################################################################################################

########################################################################################################################################
## Preparation data nodes
########################################################################################################################################
## Compute distance between each pair of sequences AND NODES in the tree
genetic_distance_mat = dist.nodes.with.names(tree)

## Get the time each node
nroot = length(tree$tip.label)+1 ## Checked and it's the root 
distance_to_root = genetic_distance_mat[nroot,]
root_height = tree$root.time
nodes_height = root_height + distance_to_root[length(names_seqs)+(1:(length(names_seqs)-1))]
all_nodes_heights = root_height + distance_to_root
  
# Meta-data with nodes 
dataset_with_nodes = data.frame('ID' = c(1:length(names_seqs), length(names_seqs)+(1:(length(names_seqs)-1))),
                                'name_seq' = c(names_seqs, length(names_seqs)+(1:(length(names_seqs)-1))),
                                'name_seq2' = c(names_seqs, tree$node.label),
                                'time' = all_nodes_heights,
                                'is.node' = c(rep('no', length(names_seqs)), rep('yes', (length(names_seqs)-1))),
                                'AMRany' = c(metadata$AMRany[match(names_seqs, metadata$Accession)], rep(NA, (length(names_seqs)-1))),
                                'ESBLR' = c(metadata$ESBLR[match(names_seqs, metadata$Accession)], rep(NA, (length(names_seqs)-1))),
                                'CiprofloxacinR' = c(metadata$CiprofloxacinR[match(names_seqs, metadata$Accession)], rep(NA, (length(names_seqs)-1))),
                                'AzithromycinR' = c(metadata$AzithromycinR[match(names_seqs, metadata$Accession)], rep(NA, (length(names_seqs)-1))),
                                'Genotype' = c(metadata$Genotype[match(names_seqs, metadata$Accession)], rep(NA, (length(names_seqs)-1))),
                                'pmsmdemdef' = c(metadata$Group[match(names_seqs, metadata$Accession)], rep(NA, (length(names_seqs)-1))))

## pMSM
## Compute distance between each pair of sequences AND NODES in the tree
genetic_distance_mat_pMSM = dist.nodes.with.names(tree_pMSM)

# Meta-data with nodes 
dataset_with_nodes_pMSM = data.frame('ID' = c(1:length(names_seqs_pMSM), length(names_seqs_pMSM)+(1:(length(names_seqs_pMSM)-1))),
                                    'name_seq' = c(names_seqs_pMSM, length(names_seqs_pMSM)+(1:(length(names_seqs_pMSM)-1))),
                                    'name_seq2' = c(names_seqs_pMSM, tree_pMSM$node.label),
                                    'time' = dataset_with_nodes$time[match(c(names_seqs_pMSM, tree_pMSM$node.label), dataset_with_nodes$name_seq2)],
                                    'is.node' = c(rep('no', length(names_seqs_pMSM)), rep('yes', (length(names_seqs_pMSM)-1))),
                                    'pMSM' = c(metadata$Group[match(names_seqs_pMSM, metadata$Accession)], rep(NA, (length(names_seqs_pMSM)-1))),
                                    'AMRany' = c(metadata$AMRany[match(names_seqs_pMSM, metadata$Accession)], rep(NA, (length(names_seqs_pMSM)-1))),
                                    'ESBLR' = c(metadata$ESBLR[match(names_seqs_pMSM, metadata$Accession)], rep(NA, (length(names_seqs_pMSM)-1))),
                                    'CiprofloxacinR' = c(metadata$CiprofloxacinR[match(names_seqs_pMSM, metadata$Accession)], rep(NA, (length(names_seqs_pMSM)-1))),
                                    'AzithromycinR' = c(metadata$AzithromycinR[match(names_seqs_pMSM, metadata$Accession)], rep(NA, (length(names_seqs_pMSM)-1))),
                                    'Genotype' = c(metadata$Genotype[match(names_seqs_pMSM, metadata$Accession)], rep(NA, (length(names_seqs_pMSM)-1))),
                                    'pmsmdemdef' = c(metadata$Group[match(names_seqs_pMSM, metadata$Accession)], rep(NA, (length(names_seqs_pMSM)-1))))

## Non-pMSM
## Compute distance between each pair of sequences AND NODES in the tree
genetic_distance_mat_nonpMSM = dist.nodes.with.names(tree_nonpMSM)

# Meta-data with nodes 
dataset_with_nodes_nonpMSM = data.frame('ID' = c(1:length(names_seqs_nonpMSM), length(names_seqs_nonpMSM)+(1:(length(names_seqs_nonpMSM)-1))),
                                       'name_seq' = c(names_seqs_nonpMSM, length(names_seqs_nonpMSM)+(1:(length(names_seqs_nonpMSM)-1))),
                                       'name_seq2' = c(names_seqs_nonpMSM, tree_nonpMSM$node.label),
                                       'time' = dataset_with_nodes$time[match(c(names_seqs_nonpMSM, tree_nonpMSM$node.label), dataset_with_nodes$name_seq2)],
                                       'is.node' = c(rep('no', length(names_seqs_nonpMSM)), rep('yes', (length(names_seqs_nonpMSM)-1))),
                                       'pMSM' = c(metadata$Group[match(names_seqs_nonpMSM, metadata$Accession)], rep(NA, (length(names_seqs_nonpMSM)-1))),
                                       'AMRany' = c(metadata$AMRany[match(names_seqs_nonpMSM, metadata$Accession)], rep(NA, (length(names_seqs_nonpMSM)-1))),
                                       'ESBLR' = c(metadata$ESBLR[match(names_seqs_nonpMSM, metadata$Accession)], rep(NA, (length(names_seqs_nonpMSM)-1))),
                                       'CiprofloxacinR' = c(metadata$CiprofloxacinR[match(names_seqs_nonpMSM, metadata$Accession)], rep(NA, (length(names_seqs_nonpMSM)-1))),
                                       'AzithromycinR' = c(metadata$AzithromycinR[match(names_seqs_nonpMSM, metadata$Accession)], rep(NA, (length(names_seqs_nonpMSM)-1))),
                                       'Genotype' = c(metadata$Genotype[match(names_seqs_nonpMSM, metadata$Accession)], rep(NA, (length(names_seqs_nonpMSM)-1))),
                                       'pmsmdemdef' = c(metadata$Group[match(names_seqs_nonpMSM, metadata$Accession)], rep(NA, (length(names_seqs_nonpMSM)-1))))

## Travel travelpMSM
## Compute distance between each pair of sequences AND NODES in the tree
genetic_distance_mat_travelpMSM = dist.nodes.with.names(tree_travelpMSM)

# Meta-data with nodes 
dataset_with_nodes_travelpMSM = data.frame('ID' = c(1:length(names_seqs_travelpMSM), length(names_seqs_travelpMSM)+(1:(length(names_seqs_travelpMSM)-1))),
                                     'name_seq' = c(names_seqs_travelpMSM, length(names_seqs_travelpMSM)+(1:(length(names_seqs_travelpMSM)-1))),
                                     'name_seq2' = c(names_seqs_travelpMSM, tree_travelpMSM$node.label),
                                     'time' = dataset_with_nodes$time[match(c(names_seqs_travelpMSM, tree_travelpMSM$node.label), dataset_with_nodes$name_seq2)],
                                     'is.node' = c(rep('no', length(names_seqs_travelpMSM)), rep('yes', (length(names_seqs_travelpMSM)-1))),
                                     'travelpMSM' = c(metadata$Group[match(names_seqs_travelpMSM, metadata$Accession)], rep(NA, (length(names_seqs_travelpMSM)-1))),
                                     'AMRany' = c(metadata$AMRany[match(names_seqs_travelpMSM, metadata$Accession)], rep(NA, (length(names_seqs_travelpMSM)-1))),
                                     'ESBLR' = c(metadata$ESBLR[match(names_seqs_travelpMSM, metadata$Accession)], rep(NA, (length(names_seqs_travelpMSM)-1))),
                                     'CiprofloxacinR' = c(metadata$CiprofloxacinR[match(names_seqs_travelpMSM, metadata$Accession)], rep(NA, (length(names_seqs_travelpMSM)-1))),
                                     'AzithromycinR' = c(metadata$AzithromycinR[match(names_seqs_travelpMSM, metadata$Accession)], rep(NA, (length(names_seqs_travelpMSM)-1))),
                                     'Genotype' = c(metadata$Genotype[match(names_seqs_travelpMSM, metadata$Accession)], rep(NA, (length(names_seqs_travelpMSM)-1))),
                                     'travelpMSMdemdef' = c(metadata$Group[match(names_seqs_travelpMSM, metadata$Accession)], rep(NA, (length(names_seqs_travelpMSM)-1))))
########################################################################################################################################

########################################################################################################################################
## Reconstruction genotypes on all nodes
########################################################################################################################################
## full data
dataset_with_nodes$Gen = dataset_with_nodes$Genotype
snp_data = dataset_with_nodes$Genotype[which(dataset_with_nodes$is.node == 'no')]
snp_data = as.numeric(as.factor(snp_data))
tree_tmp = tree
tree_tmp$node.label = NULL
rec = ace(snp_data, tree_tmp, type = 'discrete')
rec_nodes = apply(rec$lik.anc, MAR = 1, function(x)which.max(x))
rec_nodes_score = apply(rec$lik.anc, MAR = 1, function(x)max(x))
rec_nodes[which(rec_nodes_score < 0.75)] = NA
rec_all = c(as.numeric(as.character(snp_data)), rec_nodes) ## List of all states: first all tips, then all nodes
dataset_with_nodes$Genotype = rec_all

## non-pMSM
dataset_with_nodes_nonpMSM$Genotype = dataset_with_nodes$Genotype[match(dataset_with_nodes_nonpMSM$name_seq2, 
                                                                       dataset_with_nodes$name_seq2)]

## pMSM
dataset_with_nodes_pMSM$Genotype = dataset_with_nodes$Genotype[match(dataset_with_nodes_pMSM$name_seq2, 
                                                                    dataset_with_nodes$name_seq2)]

## Travel pMSM
dataset_with_nodes_travelpMSM$Genotype = dataset_with_nodes$Genotype[match(dataset_with_nodes_travelpMSM$name_seq2, 
                                                                     dataset_with_nodes$name_seq2)]
########################################################################################################################################

########################################################################################################################################
## Genotype simplified
########################################################################################################################################
eq_Gen_Genotypes = matrix(NA, ncol = 2, nrow = length(unique(dataset_with_nodes$Genotype)))
eq_Gen_Genotypes[,1] = unique(dataset_with_nodes$Genotype)
eq_Gen_Genotypes[,2] = dataset_with_nodes$Gen[match(unique(dataset_with_nodes$Genotype), dataset_with_nodes$Genotype)]

dataset_with_nodes$Genotype_simplified = eq_Gen_Genotypes[match(dataset_with_nodes$Genotype, eq_Gen_Genotypes[,1]),2]
dataset_with_nodes$Genotype_simplified = as.factor(dataset_with_nodes$Genotype_simplified )

# levels(dataset_with_nodes$Genotype_simplified) = c("1.5", "1.5.1",
#                                                     "2",
#                                                     "2.1", "2.1.1", "2.1.3", "2.1.7",
#                                                     "2.10.1", "2.10.2",
#                                                     "2.11.2", "2.11.4", "2.11.5",
#                                                     "2.12.4",
#                                                     "2.3",
#                                                     "2.4.1", "2.4.2", "2.4.3",
#                                                     "2.5.1",
#                                                     "2.6",
#                                                     "2.7.1", "2.7.3", "2.7.4",
#                                                     "3",
#                                                     "3.1",
#                                                     "3.4.1",
#                                                     "3.6",
#                                                     "3.6.1", "3.6.1.1", "3.6.1.1.1", "3.6.1.1.2", "3.6.1.1.3", "3.6.1.1.3.1",
#                                                     "3.6.2",
#                                                     "3.6.3",
#                                                     "3.6.4",
#                                                     "3.7",
#                                                     "3.7.10", "3.7.11", "3.7.12", "3.7.15", "3.7.16", "3.7.17", "3.7.18", "3.7.19",
#                                                     "3.7.20", "3.7.21", "3.7.25", "3.7.26", "3.7.27", "3.7.28",
#                                                     "3.7.29.1", "3.7.29.1.2", "3.7.29.1.2.1", "3.7.29.1.4", "3.7.29.1.4.1",
#                                                     "3.7.3",
#                                                     "3.7.30.1", "3.7.30.4", "3.7.30.4.1",
#                                                     "3.7.4",
#                                                     "3.7.6",
#                                                     "3.7.7",
#                                                     "3.7.8",
#                                                     "3.7.9",
#                                                     "5.1.6")

levels(dataset_with_nodes$Genotype_simplified) = c("ZOthers", "ZOthers",
                                                    "2",
                                                    "2", "2", "2", "2",
                                                    "2", "2",
                                                    "2", "2", "2",
                                                    "2",
                                                    "2",
                                                    "2", "2", "2",
                                                    "2",
                                                    "2",
                                                    "2", "2", "2",
                                                    "ZOthers",
                                                    "ZOthers",
                                                    "ZOthers",
                                                    "ZOthers",
                                                    "3.6.1", "3.6.1.1", "ZOthers", "3.6.1.1.2", "ZOthers", "ZOthers",
                                                    "3.6.2",
                                                    "3.6.3",
                                                    "3.6.4",
                                                    "3.7",
                                                    "ZOthers", "ZOthers", "ZOthers", "ZOthers", "3.7.16", "ZOthers", "3.7.18", "ZOthers",
                                                    "ZOthers", "ZOthers", "3.7.25", "ZOthers", "ZOthers", "ZOthers",
                                                    "3.7.29.1", "3.7.29.1", "3.7.29.1", "3.7.29.1", "3.7.29.1",
                                                    "ZOthers",
                                                    "ZOthers", "ZOthers", "3.7.30.4.1",
                                                    "ZOthers",
                                                    "ZOthers",
                                                    "ZOthers",
                                                    "ZOthers",
                                                    "ZOthers",
                                                    "ZOthers")

## non-pMSM
dataset_with_nodes_nonpMSM$Genotype_simplified = dataset_with_nodes$Genotype_simplified[match(dataset_with_nodes_nonpMSM$name_seq2, 
                                                                                             dataset_with_nodes$name_seq2)]

## pMSM
dataset_with_nodes_pMSM$Genotype_simplified = dataset_with_nodes$Genotype_simplified[match(dataset_with_nodes_pMSM$name_seq2, 
                                                                                          dataset_with_nodes$name_seq2)]

## Travel pMSM
dataset_with_nodes_travelpMSM$Genotype_simplified = dataset_with_nodes$Genotype_simplified[match(dataset_with_nodes_travelpMSM$name_seq2, 
                                                                                           dataset_with_nodes$name_seq2)]
########################################################################################################################################

########################################################################################################################################
## Compute index of every tip and node
########################################################################################################################################
wind = 1
timescale = 2 ## Timescale
dataset_with_nodes$index = compute.index(time_distance_mat = genetic_distance_mat, 
                                                timed_tree = tree, 
                                                time_window = wind,
                                                metadata = dataset_with_nodes, 
                                                mutation_rate = mu,
                                                timescale = timescale,
                                                genome_length = genome_length)

dataset_with_nodes_pMSM$index = compute.index(time_distance_mat = genetic_distance_mat_pMSM, 
                                             timed_tree = tree_pMSM, 
                                             time_window = wind,
                                             metadata = dataset_with_nodes_pMSM, 
                                             mutation_rate = mu,
                                             timescale = timescale,
                                             genome_length = genome_length)
dataset_with_nodes_nonpMSM$index = compute.index(time_distance_mat = genetic_distance_mat_nonpMSM, 
                                                timed_tree = tree_nonpMSM, 
                                                time_window = wind,
                                                metadata = dataset_with_nodes_nonpMSM, 
                                                mutation_rate = mu,
                                                timescale = timescale,
                                                genome_length = genome_length)
dataset_with_nodes_travelpMSM$index = compute.index(time_distance_mat = genetic_distance_mat_travelpMSM, 
                                                   timed_tree = tree_travelpMSM, 
                                                   time_window = wind,
                                                   metadata = dataset_with_nodes_travelpMSM, 
                                                   mutation_rate = mu,
                                                   timescale = timescale,
                                                   genome_length = genome_length)
########################################################################################################################################

########################################################################################################################################
## Colors datastest
########################################################################################################################################
library(MetBrewer)
## Colours based on full dataset
tmp = unlist(lapply(dataset_with_nodes$Genotype, function(x)paste0(str_split(x, '\\.')[[1]][1:min(c(3,length(str_split(x, '\\.')[[1]])))], collapse = '.')))
tmp[which(tmp == 'NA')] = NA
clades = c(levels(as.factor(tmp)))
colors_clades = c(met.brewer(name="Cross", n=length(clades), type="continuous", direction = 1))

## Full
tmp = unlist(lapply(dataset_with_nodes$Genotype, function(x)paste0(str_split(x, '\\.')[[1]][1:min(c(3,length(str_split(x, '\\.')[[1]])))], collapse = '.')))
dataset_with_nodes$Genotype_color = factor(tmp, levels = clades)
levels(dataset_with_nodes$Genotype_color) = colors_clades
dataset_with_nodes$Genotype_color = as.character(dataset_with_nodes$Genotype_color)
dataset_with_nodes$Genotype_color[which(is.na(dataset_with_nodes$Genotype_color))] = 'grey'

## pMSM
tmp = unlist(lapply(dataset_with_nodes_pMSM$Genotype, function(x)paste0(str_split(x, '\\.')[[1]][1:min(c(3,length(str_split(x, '\\.')[[1]])))], collapse = '.')))
dataset_with_nodes_pMSM$Genotype_color = factor(tmp, levels = clades)
levels(dataset_with_nodes_pMSM$Genotype_color) = colors_clades
dataset_with_nodes_pMSM$Genotype_color = as.character(dataset_with_nodes_pMSM$Genotype_color)
dataset_with_nodes_pMSM$Genotype_color[which(is.na(dataset_with_nodes_pMSM$Genotype_color))] = 'grey'

## Non-pMSM
tmp = unlist(lapply(dataset_with_nodes_nonpMSM$Genotype, function(x)paste0(str_split(x, '\\.')[[1]][1:min(c(3,length(str_split(x, '\\.')[[1]])))], collapse = '.')))
dataset_with_nodes_nonpMSM$Genotype_color = factor(tmp, levels = clades)
levels(dataset_with_nodes_nonpMSM$Genotype_color) = colors_clades
dataset_with_nodes_nonpMSM$Genotype_color = as.character(dataset_with_nodes_nonpMSM$Genotype_color)
dataset_with_nodes_nonpMSM$Genotype_color[which(is.na(dataset_with_nodes_nonpMSM$Genotype_color))] = 'grey'

## Travel pMSM
tmp = unlist(lapply(dataset_with_nodes_travelpMSM$Genotype, function(x)paste0(str_split(x, '\\.')[[1]][1:min(c(3,length(str_split(x, '\\.')[[1]])))], collapse = '.')))
dataset_with_nodes_travelpMSM$Genotype_color = factor(tmp, levels = clades)
levels(dataset_with_nodes_travelpMSM$Genotype_color) = colors_clades
dataset_with_nodes_travelpMSM$Genotype_color = as.character(dataset_with_nodes_travelpMSM$Genotype_color)
dataset_with_nodes_travelpMSM$Genotype_color[which(is.na(dataset_with_nodes_travelpMSM$Genotype_color))] = 'grey'
########################################################################################################################################

########################################################################################################################################
## Colors Genotype_simplified
########################################################################################################################################
library(MetBrewer)
## Colours based on full dataset
tmp = unlist(lapply(dataset_with_nodes$Genotype_simplified, function(x)paste0(str_split(x, '\\.')[[1]][1:min(c(10,length(str_split(x, '\\.')[[1]])))], collapse = '.')))
tmp[which(tmp == 'NA' | tmp == 'ZOthers')] = NA
clades_simplified = c(levels(as.factor(tmp)))
colors_clades_simplified = c(met.brewer(name="Cross", n=length(clades_simplified), type="continuous", direction = 1))

## Full
tmp = unlist(lapply(dataset_with_nodes$Genotype_simplified, function(x)paste0(str_split(x, '\\.')[[1]][1:min(c(10,length(str_split(x, '\\.')[[1]])))], collapse = '.')))
dataset_with_nodes$Genotype_simplified_color = factor(tmp, levels = clades_simplified)
levels(dataset_with_nodes$Genotype_simplified_color) = colors_clades_simplified
dataset_with_nodes$Genotype_simplified_color = as.character(dataset_with_nodes$Genotype_simplified_color)
dataset_with_nodes$Genotype_simplified_color[which(is.na(dataset_with_nodes$Genotype_simplified_color))] = 'grey'

## pMSM
tmp = unlist(lapply(dataset_with_nodes_pMSM$Genotype_simplified, function(x)paste0(str_split(x, '\\.')[[1]][1:min(c(10,length(str_split(x, '\\.')[[1]])))], collapse = '.')))
dataset_with_nodes_pMSM$Genotype_simplified_color = factor(tmp, levels = clades_simplified)
levels(dataset_with_nodes_pMSM$Genotype_simplified_color) = colors_clades_simplified
dataset_with_nodes_pMSM$Genotype_simplified_color = as.character(dataset_with_nodes_pMSM$Genotype_simplified_color)
dataset_with_nodes_pMSM$Genotype_simplified_color[which(is.na(dataset_with_nodes_pMSM$Genotype_simplified_color))] = 'grey'

## Non-pMSM
tmp = unlist(lapply(dataset_with_nodes_nonpMSM$Genotype_simplified, function(x)paste0(str_split(x, '\\.')[[1]][1:min(c(10,length(str_split(x, '\\.')[[1]])))], collapse = '.')))
dataset_with_nodes_nonpMSM$Genotype_simplified_color = factor(tmp, levels = clades_simplified)
levels(dataset_with_nodes_nonpMSM$Genotype_simplified_color) = colors_clades_simplified
dataset_with_nodes_nonpMSM$Genotype_simplified_color = as.character(dataset_with_nodes_nonpMSM$Genotype_simplified_color)
dataset_with_nodes_nonpMSM$Genotype_simplified_color[which(is.na(dataset_with_nodes_nonpMSM$Genotype_simplified_color))] = 'grey'

## Travel pMSM
tmp = unlist(lapply(dataset_with_nodes_travelpMSM$Genotype_simplified, function(x)paste0(str_split(x, '\\.')[[1]][1:min(c(10,length(str_split(x, '\\.')[[1]])))], collapse = '.')))
dataset_with_nodes_travelpMSM$Genotype_simplified_color = factor(tmp, levels = clades_simplified)
levels(dataset_with_nodes_travelpMSM$Genotype_simplified_color) = colors_clades_simplified
dataset_with_nodes_travelpMSM$Genotype_simplified_color = as.character(dataset_with_nodes_travelpMSM$Genotype_simplified_color)
dataset_with_nodes_travelpMSM$Genotype_simplified_color[which(is.na(dataset_with_nodes_travelpMSM$Genotype_simplified_color))] = 'grey'
########################################################################################################################################

########################################################################################################################################
## Save results
########################################################################################################################################
save.image(file='2_analysis_index/1_index_computations/Initial_index_computation_and_parameters_20251218.Rdata')
########################################################################################################################################

