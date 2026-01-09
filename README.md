# Code supporting to replicate the Phylowave analysi of the paper "The spread of sexually transmissible drug-resistant shigellosis"

## Pre-requisites
Before starting, make sure you cloned the Phylowave repo: https://github.com/noemielefrancq/Phylowave_Learning-fitness-dynamics-pathogens-in-phylogenies

## Data
You can find the data used to perform the analysis in the "1_Data" folder.

## Run analaysis
To replicate the analysis, all the R codes are procided in the folder "2_analysis_index". 
The analysis is organised in 2 steps:
- 1: Compute the index dynamics "1_index_computations"
  - 1.1: 1_1_initial_index_computation_and_parameters.R
  - 1.2: 1_2_plot_initial_index.R
- 2: Detect lineages and explore them "2_find_index_groups"
  - 2_1_exploration_lineage_detection_122024.R
  - 2_2_plot_exploration.R
  - 2_3_detect_lienages_chosen.R
  - 2_4_plot_detected_lineages.R
  - 2_5_comparison_groups_clades
Corresponding folders are included in the folder "2_analysis_index".
