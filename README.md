# Linking pattern to process using phylogenies—insights from a comparison of eight simulation models

Allen H. Hurlbert, Florian Hartig, Juliano Sarmento Cabral, Rampal S. Etienne, Oskar Hagen, Shan Huang, Yun-Ting Jang, Thijs Janzen, Antonin Machac, Pedro Santos Neves, Mikael Pontarp, David Storch, Liang Xu  

This repository contains code underlying the comparison of 8 eco-evolutionary simulation models and the evaluation of links between eco-evolutionary processes and phylogenetic patterns.  

## Organization  

### Folders  
`trees`  
--`uniform_sampling_experiment`: all simulated phylogenetic trees used in the study, with file names in the form *mm_xxxxx.tre* where *mm* is the simulation model abbreviation, and *xxxx* is a unique simulation ID.  
--`empirical`: 16 empirical trees for 8 pairs of sister clades (4 pairs of mammals and 4 pairs of birds)

`derived_tree_data`
--`simulation_parameter_keys.csv`: key to the parameters used across the 8 different models, specifying which process the parameter was related to, its name, its description, and the direction (sign) of its relation to the stated process  
--`USE_treeStats.txt`: 70+ tree metrics calculated from the `treestats` R package for all simulated trees  
--`process_parameter_values_and_tree_metrics_sign_corrected.csv`: all simulated trees, their parameter settings grouped by relevant process and sign-corrected as necessary to standardize the direction of effect, and their 70+ tree metrics  
--`empiricalSisterClades_treeStats.txt`: 70+ tree metrics calculated from the `treestats` R package for each empirical tree  
--`scaled_predictions_empirical_trees.csv`: for each combination of empirical tree, simulation model, and eco-evolutionary process, the inferred parameter value (`value`) related to that process based on a random forest trained on trees simulated by that model, as well as the value rescaled to vary between 0 and 1 (`scaledValue`) based on the minimum and maximum parameter values actually used across all simulations.  

`code`  
--`analysis_functions.r`: functions used for calculating tree stats and for aligning parameters with processes.  

`figures`  
--manuscript figures

### Analysis workflow
`Analysis.qmd`: Complete analysis workflow given the `USE_treeStats.txt` and `empiricalSisterClades_treeStats.txt` output  
