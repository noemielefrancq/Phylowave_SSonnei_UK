# Code supporting to replicate the Phylowave analysi of the paper "The spread of sexually transmissible drug-resistant shigellosis"

## Prerequisites
Before starting, make sure you cloned the *Phylowave repo*: https://github.com/noemielefrancq/Phylowave_Learning-fitness-dynamics-pathogens-in-phylogenies, and have installed all the required libraries listed there.

## Data
You can find the data used to perform the analysis in the "1_Data" folder.

## Run analaysis
To replicate the analysis, all the R codes are procided in the folder "2_analysis_index". 
The analysis is organised in 2 steps:
- **Compute the index dynamics, folder "1_index_computations"**

  The following codes have to be run *sequentially*:
  - 1.1: 1_1_initial_index_computation_and_parameters.R
  - 1.2: 1_2_plot_initial_index.R
- **Detect lineages and explore them, folder "2_find_index_groups"**

  After the script 1.1, the following codes have to be run *sequentially*:
  - 2.1: 2_1_exploration_lineage_detection_122024.R
  - 2.2: 2_2_plot_exploration.R
  - 2.3: 2_3_detect_lienages_chosen.R
  - 2.4: 2_4_plot_detected_lineages.R
  - 2.5: 2_5_comparison_groups_clades.R
 
Raw figures in PDFs are included in the repository.
