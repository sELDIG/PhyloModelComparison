###########################################################################

library(devtools)
#install_github("thijsjanzen/treestats")
library(treestats)
library(ape)
library(geiger)
library(stringr)
library(dplyr)

# Calculating various attributes on a phylogeny in Newick format using
# Thijs Janzen's treestats package
# https://github.com/thijsjanzen/treestats
#
#   treeInput: a phylogenetic tree of class phylo

treeStats = function(treeInput) {
  
  if(class(treeInput) != "phylo") {
    stop("The object passed to treeMetrics() is not of class 'phylo'")
  }
  
  require(treestats)
  
  # Drop root edge
  treeInput$root.edge = 0
  
  # Prune out extinct species
  treeExtant = tryCatch({
    drop.extinct(treeInput, tol = 1e-7)
  }, error = function(e) {
    treeExtant = drop.fossil(treeInput)
  })
  
  # Tree must have at least 2 extant tips to run calc_all_stats below. If it does not, return NA
  
  if (length(treeExtant$tip.label) < 2) {
    
    treestats = length(treeExtant$tip.label)
    names(treestats) = "number_of_lineages"
    return(treestats)
  }
  
  # Randomly resolve any existing polytomies
  tree = multi2di(treeExtant)
  
  # Absolute tree length, which is the diagonal of the vcv matrix
  # (distance from tip to root; https://bio.libretexts.org/Bookshelves/Evolutionary_Developmental_Biology/Phylogenetic_Comparative_Methods_(Harmon)/03%3A_Introduction_to_Brownian_Motion/3.04%3A_Brownian_Motion_on_a_Phylogenetic_Tree)
  v.matrix = vcv(tree, corr=F)
  tree.length = diag(v.matrix)[1]
  
  # Tree with branch lengths scaled by total tree length
  tree.scaled = tree
  tree.scaled$edge.length = tree$edge.length/tree.length
  
  treestats = calc_all_stats(tree.scaled)
  
  return(treestats)
}



# Function for providing a list of tree filenames over which to calculate all tree stats
metricsForManyTrees = function(treefiles = NULL, 
                               minimumTreeSize = 20, 
                               fileOut, 
                               treedir = 'trees', 
                               append = TRUE) {
  
  if(is.null(treefiles)) {
    treefiles = list.files(treedir)[grepl(".tre", list.files(treedir))]
  }  
  
  # If append == FALSE, then create a new file with tree metric headers
  if (!append) {
    sink(fileOut)
    cat(paste("model", "simID", "number_of_lineages", "tree_height", "phylogenetic_div", 
              "gamma", "beta", "colless", "sackin", "var_depth", "psv", "i_stat", "mpd", 
              "vpd", "laplac_spectrum_e", "laplac_spectrum_a", "laplac_spectrum_p", 
              "laplac_spectrum_g", "nltt_base", "blum", "crown_age", "pigot_rho", 
              "avg_ladder", "max_ladder", "cherries", "il_number", "pitchforks", 
              "stairs", "imbalance_steps", "j_one", "b1", "b2", "area_per_pair", 
              "average_leaf_depth", "ew_colless", "max_del_width", "max_depth", 
              "max_width", "rogers", "stairs2", "tot_coph", "symmetry_nodes", "mntd", 
              "j_stat", "rquartet", "wiener", "max_betweenness", "max_closeness", 
              "diameter", "eigenvector", "mean_branch_length", "var_branch_length", 
              "mean_branch_length_int", "mean_branch_length_ext", "var_branch_length_int", 
              "var_branch_length_ext", "\n", sep = "\t"))
    sink()
    
  }

    for (treefile in treefiles) {
    
    tree = read.tree(paste(treedir, "/", treefile, sep = ""))
    
    if(tree$Nnode + 1 >= minimumTreeSize) {
      print(paste(treefile, Sys.time()))
      
      metrics = treeStats(tree)

      model = str_extract(treefile, "^[A-Za-z]*")
      simID = str_extract(treefile, "[0-9]+")
      
      

      sink(fileOut, append = TRUE)
      cat(paste(model, 
                simID,
                metrics[names(metrics) == "number_of_lineages"], 
                metrics[names(metrics) == "tree_height"],
                metrics[names(metrics) == "phylogenetic_div"], 
                metrics[names(metrics) == "gamma"],
                metrics[names(metrics) == "beta"], 
                metrics[names(metrics) == "colless"], 
                metrics[names(metrics) == "sackin"],
                metrics[names(metrics) == "var_depth"], #VRD
                metrics[names(metrics) == "psv"], 
                metrics[names(metrics) == "i_stat"],    #mean.Iprime
                metrics[names(metrics) == "mpd"], 
                metrics[names(metrics) == "vpd"],
                metrics[names(metrics) == "laplac_spectrum_e"], 
                metrics[names(metrics) == "laplac_spectrum_a"], 
                metrics[names(metrics) == "laplac_spectrum_p"], 
                metrics[names(metrics) == "laplac_spectrum_g"],
                metrics[names(metrics) == "nltt_base"],
                metrics[names(metrics) == "blum"],
                metrics[names(metrics) == "crown_age"],
                metrics[names(metrics) == "pigot_rho"],
                metrics[names(metrics) == "avg_ladder"],
                metrics[names(metrics) == "max_ladder"],
                metrics[names(metrics) == "cherries"],
                metrics[names(metrics) == "il_number"],
                metrics[names(metrics) == "pitchforks"],
                metrics[names(metrics) == "stairs"],
                metrics[names(metrics) == "imbalance_steps"],
                metrics[names(metrics) == "j_one"],
                metrics[names(metrics) == "b1"],
                metrics[names(metrics) == "b2"],
                metrics[names(metrics) == "area_per_pair"],
                metrics[names(metrics) == "average_leaf_depth"], #sounds like MRD"], but definition in Table doesn't match
                metrics[names(metrics) == "ew_colless"],
                metrics[names(metrics) == "max_del_width"],
                metrics[names(metrics) == "max_depth"],
                metrics[names(metrics) == "max_width"],
                metrics[names(metrics) == "rogers"],
                metrics[names(metrics) == "stairs2"],
                metrics[names(metrics) == "tot_coph"],
                metrics[names(metrics) == "symmetry_nodes"],
                metrics[names(metrics) == "mntd"],
                metrics[names(metrics) == "j_stat"],
                metrics[names(metrics) == "rquartet"],
                metrics[names(metrics) == "wiener"],
                metrics[names(metrics) == "max_betweenness"],
                metrics[names(metrics) == "max_closeness"],
                metrics[names(metrics) == "diameter"],
                metrics[names(metrics) == "eigenvector"],
                metrics[names(metrics) == "mean_branch_length"],
                metrics[names(metrics) == "var_branch_length"],
                metrics[names(metrics) == "mean_branch_length_int"],
                metrics[names(metrics) == "mean_branch_length_ext"],
                metrics[names(metrics) == "var_branch_length_int"],
                metrics[names(metrics) == "var_branch_length_ext"],
                '\n',
                sep = '\t'))
      sink()
      
      
    } else {
      print(paste(treefile, "skipped -- not enough species"))
    }
    
  } # end for loop  
  
} # end function



# Function for aligning simulation parameters across models according to the process each parameter is associated with.
# I.e., create columns for each process (in some cases, two columns for a process because some models have two parameters
# associated with that process) in which the relevant parameter value is stored.

# Associations between parameters and processes is originally given here: 
# https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=1171496897
# but has been written to 'experiments/uniform_sampling_experiment/simulation_parameters_key.csv'

# Parameters are multiplied by the sign specified in the above spreadsheet to ensure that the strength of the process
# increases with an increase in the parameter value.

alignParametersWithProcesses = function(modelAbbrev, paramKeyPath) {
  
  params = read.csv(paste(paramKeyPath, "/", modelAbbrev, "_USE_parameters.csv", sep = ""), header = T)
  
  if ("scenario" %in% names(params)) {
    params$model2 = paste(modelAbbrev, ".", params$scenario, sep = "")
  } else {
    params$model2 = params$model
  }
  
  # Here we drop the scenario description, assuming that the process-parameter association is not scenario-dependent
  paramKey <- read.csv(paste(paramKeyPath, "/", "simulation_parameters_key.csv", sep = ""), header = T) %>% 
    mutate(model = word(model, sep = "\\.")) %>%
    dplyr::filter(model == modelAbbrev) %>%
    distinct()
  
  # parameter names associated with each process
  env1name = ifelse("env" %in% paramKey$experiment, paramKey$parameterName[paramKey$experiment == "env"][1], NA)
  env2name = ifelse(length(unique(paramKey$parameterName[paramKey$experiment == "env"])) == 2, 
                    paramKey$parameterName[paramKey$experiment == "env"][2], NA)
  dis1name = ifelse("dis" %in% paramKey$experiment, paramKey$parameterName[paramKey$experiment == "dis"][1], NA)
  dis2name = ifelse(length(unique(paramKey$parameterName[paramKey$experiment == "dis"])) == 2, 
                    paramKey$parameterName[paramKey$experiment == "dis"][2], NA)
  nic1name = ifelse("nic" %in% paramKey$experiment, paramKey$parameterName[paramKey$experiment == "nic"][1], NA)
  nic2name = ifelse(length(unique(paramKey$parameterName[paramKey$experiment == "nic"])) == 2, 
                    paramKey$parameterName[paramKey$experiment == "nic"][2], NA)
  mut1name = ifelse("mut" %in% paramKey$experiment, paramKey$parameterName[paramKey$experiment == "mut"][1], NA)
  mut2name = ifelse(length(unique(paramKey$parameterName[paramKey$experiment == "mut"])) == 2, 
                    paramKey$parameterName[paramKey$experiment == "mut"][2], NA)
  com1name = ifelse("com" %in% paramKey$experiment, paramKey$parameterName[paramKey$experiment == "com"][1], NA)
  com2name = ifelse(length(unique(paramKey$parameterName[paramKey$experiment == "com"])) == 2, 
                    paramKey$parameterName[paramKey$experiment == "com"][2], NA)
  
  
  outputDF = params %>%
    mutate(
      env1Name = env1name,
      env2Name = env2name,
      dis1Name = dis1name,
      dis2Name = dis2name,
      nic1Name = nic1name,
      nic2Name = nic2name,
      mut1Name = mut1name,
      mut2Name = mut2name,
      com1Name = com1name,
      com2Name = com2name
    )
  
  outputDF$env1 = ifelse(!is.na(outputDF$env1Name), paramKey$sign[paramKey$parameterName == env1name] * outputDF[, env1name], NA)
  outputDF$env2 = ifelse(!is.na(outputDF$env2Name), paramKey$sign[paramKey$parameterName == env2name] * outputDF[, env2name], NA)
  outputDF$dis1 = ifelse(!is.na(outputDF$dis1Name), paramKey$sign[paramKey$parameterName == dis1name] * outputDF[, dis1name], NA)
  outputDF$dis2 = ifelse(!is.na(outputDF$dis2Name), paramKey$sign[paramKey$parameterName == dis2name] * outputDF[, dis2name], NA)
  outputDF$nic1 = ifelse(!is.na(outputDF$nic1Name), paramKey$sign[paramKey$parameterName == nic1name] * outputDF[, nic1name], NA)
  outputDF$nic2 = ifelse(!is.na(outputDF$nic2Name), paramKey$sign[paramKey$parameterName == nic2name] * outputDF[, nic2name], NA)
  outputDF$mut1 = ifelse(!is.na(outputDF$mut1Name), paramKey$sign[paramKey$parameterName == mut1name] * outputDF[, mut1name], NA)
  outputDF$mut2 = ifelse(!is.na(outputDF$mut2Name), paramKey$sign[paramKey$parameterName == mut2name] * outputDF[, mut2name], NA)
  outputDF$com1 = ifelse(!is.na(outputDF$com1Name), paramKey$sign[paramKey$parameterName == com1name] * outputDF[, com1name], NA)
  outputDF$com2 = ifelse(!is.na(outputDF$com2Name), paramKey$sign[paramKey$parameterName == com2name] * outputDF[, com2name], NA)
  
  output = outputDF %>%
    dplyr::select(model, model2, simID, env1Name:com2)
  
  return(output)  
}







