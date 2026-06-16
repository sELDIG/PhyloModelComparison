# Linking pattern to process using phylogenies—insights from a comparison of eight simulation models

Allen H. Hurlbert\*, Juliano Sarmento Cabral, Rampal S. Etienne, Oskar Hagen, Shan Huang,
Yun-Ting Jang, Thijs Janzen, Antonin Machac, Pedro Santos Neves, Mikael Pontarp, David Storch,
Liang Xu, Florian Hartig\*

\*Equally contributing lead authors

This repository contains code underlying the comparison of 8 eco-evolutionary simulation models
and the evaluation of links between eco-evolutionary processes and phylogenetic patterns.

---

## Description of the data and file structure

All folders and files below are in the zipped `PhyloModelComparison-1.0.0.zip`.

### Folders

#### `trees/`
- **`uniform_sampling_experiment/`**: All simulated phylogenetic trees used in the study, with
  file names in the form `mm_xxxxx.tre` where `mm` is the simulation model abbreviation (see
  [Model Abbreviations](#model-abbreviations) below) and `xxxxx` is a unique simulation ID.
- **`empirical/`**: 16 empirical trees for 8 pairs of sister clades (4 pairs of mammals and 4
  pairs of birds).

#### `derived_tree_data/`
- **`simulation_parameter_keys.csv`**: Key to the parameters used across the 8 different models,
  specifying which eco-evolutionary process each parameter is related to, its name, its
  description, and the direction (sign) of its relation to the stated process. See
  [Data Dictionary: simulation_parameter_keys.csv](#simulation_parameter_keyscsv).
- **`USE_treeStats.txt`**: 70+ tree shape and diversification metrics calculated using the
  `treestats` R package for all simulated trees. See
  [Data Dictionary: USE_treeStats.txt](#use_treestatstxt).
- **`process_parameter_values_and_tree_metrics_sign_corrected.csv`**: All simulated trees, their
  parameter settings grouped by relevant eco-evolutionary process and sign-corrected as necessary
  to standardize the direction of effect, alongside their 70+ tree metrics. See
  [Data Dictionary: process_parameter_values_and_tree_metrics_sign_corrected.csv](#process_parameter_values_and_tree_metrics_sign_correctedcsv).
- **`empiricalSisterClades_treeStats.txt`**: 70+ tree metrics calculated using the `treestats` R
  package for each empirical tree. Shares the same metric columns as `USE_treeStats.txt`, but
  with different identifier columns. See
  [Data Dictionary: empiricalSisterClades_treeStats.txt](#empiricalsistercladestreestatstxt).
- **`scaled_predictions_empirical_trees.csv`**: For each combination of empirical tree,
  simulation model, and eco-evolutionary process, the inferred parameter value (`value`) related
  to that process based on a random forest trained on trees simulated by that model, as well as
  the value rescaled to vary between 0 and 1 (`scaledValue`) based on the minimum and maximum
  parameter values actually used across all simulations. See
  [Data Dictionary: scaled_predictions_empirical_trees.csv](#scaled_predictions_empirical_treescsv).

#### `code/`
- **`analysis_functions.r`**: Functions used for calculating tree statistics and for aligning
  parameters with processes.
- **`plot_functions.r`**: Functions used for making figures.

#### `bootstrap/`
- Results from the bootstrap analysis in which simulation results were resampled with
  replacement. Each of 250 bootstrap runs is presented as its own Rdata file.

---

## Analysis workflow

**First, run:**

**`Calculations.qmd`**: Workflow that runs the analysis and saves results (e.g., average
correlations, cross-validation) separately for the dataset as a whole or for a single bootstrap
run.

**Then run:**

**`bootstrapLoopParallel.qmd`**: Loop for running the analysis across all bootstrap iterations.

**Finally, run:**

**`Analysis.qmd`**: Complete analysis workflow given the `USE_treeStats.txt` and
`empiricalSisterClades_treeStats.txt` output, or reading in the individual bootstrap run Rdata
files, that produces manuscript figures.

---

## Sharing/Access information

These files are available in the GitHub repository:
https://github.com/sELDIG/PhyloModelComparison

---

## Data Dictionary

This section provides detailed descriptions of all variables (column names) in each tabular data
file in the `derived_tree_data/` folder. It includes the full meaning of all abbreviations used,
units of measurement where applicable, and explanations of missing (`NA`) values.

---

### General Notes on Units

All phylogenetic trees in this dataset — both simulated and empirical — are scaled so that
`tree_height = 1`. Consequently:

- All branch length–based metrics are expressed relative to tree height (i.e., as a fraction of
  total tree height) and are therefore dimensionless on this normalized scale.
- All time-related statistics derived from branch lengths are similarly dimensionless.
- Statistics based on counts of nodes, tips, or lineages are unitless integers.
- Statistics that are ratios or indices are dimensionless.

For full mathematical definitions and original references for each tree metric, consult the
[`treestats` R package documentation](https://cran.r-project.org/package=treestats) and Table S3.

---

### Model Abbreviations

The following abbreviations identify the eight simulation models used throughout the data files.
Some models have multiple variants (indicated by suffixes) reflecting different experimental
conditions or parameterizations within that model. Refer to the associated manuscript for full
model descriptions, biological rationale, and references.

| Abbreviation | Variants in data files             |
|--------------|------------------------------------|
| `ca`         | `ca`                               |
| `hs`         | `hs.K`, `hs.S`, `hs.D`, `hs.N`    |
| `ve`         | `ve`                               |
| `xe`         | `xe`                               |
| `fh`         | `fh.1`, `fh.2`, `fh.3`            |
| `gen`        | `gen.NULL`, `gen.TEMPSP`, `gen.K`  |
| `la`         | `la`                               |
| `pontarp`    | `pontarp`                          |

---

### Process/Experiment Category Abbreviations

The following abbreviations are used in the `experiment` column of
`simulation_parameter_keys.csv` and as prefixes of the parameter value columns (e.g., `env1`,
`dis2`) in `process_parameter_values_and_tree_metrics_sign_corrected.csv`:

| Abbreviation | Full name                        | Description |
|--------------|----------------------------------|-------------|
| `env`        | Environmental filtering          | Parameters controlling how strongly the abiotic environment (e.g., temperature) filters species persistence or fitness. |
| `dis`        | Dispersal                        | Parameters controlling the rate and/or distance of species dispersal. |
| `nic`        | Niche evolution                  | Parameters controlling the rate or magnitude of niche trait evolution (niche conservatism vs. lability). |
| `mut`        | Mutation / Speciation            | Parameters controlling the rate of point mutation, speciation, cladogenesis, or anagenesis. |
| `com`        | Competition / Community interactions | Parameters controlling the strength of interspecific competition or other community-level biotic interactions. |
| `tim`        | Time                             | Parameters controlling time-related aspects of the simulation (e.g., rate of environmental change, timing of tree recording). This category appears only in `simulation_parameter_keys.csv` for the `ca` model and has no corresponding columns in `process_parameter_values_and_tree_metrics_sign_corrected.csv`. |

---

### `simulation_parameter_keys.csv`

This file provides a cross-model key linking each simulation parameter to the eco-evolutionary
process it represents, along with the sign of its effect. This information is used to
sign-correct parameter values in `process_parameter_values_and_tree_metrics_sign_corrected.csv`.

| Column | Full name | Description | Values |
|--------|-----------|-------------|--------|
| `model` | Model abbreviation | Abbreviation of the simulation model (see [Model Abbreviations](#model-abbreviations)). | Character string (e.g., `ca`, `hs.K`, `pontarp`) |
| `experiment` | Experiment / process category | The eco-evolutionary process category associated with this parameter (see [Process/Experiment Category Abbreviations](#processexperiment-category-abbreviations)). | One of: `env`, `dis`, `nic`, `mut`, `com`, `tim` |
| `parameterName` | Parameter name | The name of the parameter as used in the simulation model's code or documentation. | Character string (e.g., `Mutation_rate`, `dispersal`, `w`) |
| `sign` | Direction of effect | The direction of the parameter's effect on the stated process. `1` means higher parameter values indicate a stronger process effect. `-1` means higher parameter values indicate a weaker process effect (inverse relationship). This value is used to multiply raw parameter values in `process_parameter_values_and_tree_metrics_sign_corrected.csv` so that, after correction, higher values always indicate stronger process influence. | Integer: `1` or `-1` |
| `explanationOfEffect` | Explanation of effect | A plain-language description of what the parameter controls and how its value relates to the process category. | Character string |


---

### `USE_treeStats.txt`

This tab-delimited file contains 70+ phylogenetic tree shape and diversification metrics
computed using the `treestats` R package for all simulated trees.

#### Identifier columns

| Column | Full name | Description | Values |
|--------|-----------|-------------|--------|
| `model` | Model abbreviation | Abbreviation of the simulation model that generated this tree (see [Model Abbreviations](#model-abbreviations)). | Character string |
| `simID` | Simulation ID | Unique integer identifier for each simulation run within a given model. Corresponds to the numeric portion of the tree filename (`mm_xxxxx.tre`) in `trees/uniform_sampling_experiment/`. | Integer |

#### Tree metric columns

All metrics below are computed by the `treestats` R package. Unless otherwise noted, metrics are
dimensionless because all trees are normalized to `tree_height = 1`.

| Column | Full name | Description | Units |
|--------|-----------|-------------|-------|
| `number_of_lineages` | Number of lineages | Total number of tips (extant species) in the phylogenetic tree. | Count (integer) |
| `tree_height` | Tree height | Total height of the tree from root to tips. All trees in this dataset are scaled to a height of 1 prior to metric computation. | Relative (= 1 for all trees in this dataset) |
| `phylogenetic_div` | Phylogenetic diversity (PD) | Faith's phylogenetic diversity; the sum of all branch lengths in the tree. | Branch length units (relative to `tree_height` = 1) |
| `gamma` | Gamma statistic (γ) | Measures whether internal nodes are concentrated toward the root (negative: early burst of diversification) or toward the tips (positive: recent radiation). Under a constant-rate pure-birth process, γ approximately follows N(0,1). (Pybus & Harvey 2000) | Dimensionless |
| `beta` | Aldous' beta (β) | Measures tree shape/balance. β = −2 for a fully pectinate (caterpillar) tree; β = 0 under a Yule model; larger positive values indicate more balanced trees. (Aldous 2001) | Dimensionless |
| `colless` | Colless imbalance index | Sum of the absolute differences in tip-descendant counts between sister clades across all internal nodes. Higher values indicate more imbalanced trees. (Colless 1982) | Count-based (integer) |
| `sackin` | Sackin imbalance index | Sum of the number of ancestors for each tip across the tree. Higher values indicate more imbalanced (pectinate) trees. (Sackin 1972) | Count-based (integer) |
| `var_depth` | Variance in node depth | Variance in the number of nodes (edges) from the root to each tip. Low values indicate a balanced tree; high values indicate uneven tip depths. | Dimensionless |
| `psv` | Phylogenetic species variability | Ranges from 0 to 1. Values near 1 indicate species are maximally distantly related (star phylogeny); values near 0 indicate high phylogenetic clustering or severe imbalance. (Helmus et al. 2007) | Dimensionless (0–1) |
| `i_stat` | Mean I' statistic (Ī') | Mean of the I' balance index across all internal nodes. Ranges from 0 (completely imbalanced) to 1 (perfectly balanced). | Dimensionless (0–1) |
| `mpd` | Mean pairwise distance (MPD) | Mean phylogenetic distance (sum of connecting branch lengths) between all pairs of tips. | Branch length units (relative to `tree_height` = 1) |
| `vpd` | Variance in pairwise distances (VPD) | Variance in phylogenetic distances across all pairs of tips. | Squared branch length units |
| `laplac_spectrum_e` | Laplacian spectrum — principal eigenvalue | The largest eigenvalue of the graph Laplacian matrix of the tree. Summarizes overall tree connectivity and structure. | Dimensionless |
| `laplac_spectrum_a` | Laplacian spectrum — asymmetry | Asymmetry (skewness) of the distribution of Laplacian eigenvalues. | Dimensionless |
| `laplac_spectrum_p` | Laplacian spectrum — peakedness | Peakedness (kurtosis) of the distribution of Laplacian eigenvalues. | Dimensionless |
| `laplac_spectrum_g` | Laplacian spectrum — eigengap | Eigengap of the Laplacian spectrum: the difference between the largest and second-largest eigenvalues. | Dimensionless |
| `nltt_base` | Normalized lineage-through-time statistic (nLTT) | Area between the normalized lineage-through-time (LTT) curve of the tree and that of a reference pure-birth process. A value of 0 indicates a perfect match to the pure-birth reference. (Janzen et al. 2015) | Dimensionless (≥ 0) |
| `blum` | Blum & François' N̄ statistic | Measures tree balance by comparing the tree to expectations under a random branching model. Higher values indicate more balanced trees. (Blum & François 2006) | Dimensionless |
| `crown_age` | Crown age | Age of the most recent common ancestor (MRCA) of all extant tips. Because trees are ultrametric and scaled to `tree_height = 1`, this value equals 1 for all trees and is therefore uninformative on its own in this dataset. | Relative (= 1 for all trees in this dataset) |
| `pigot_rho` | Pigot's rho (ρ) | Measures temporal acceleration or deceleration of lineage accumulation. Positive values indicate acceleration toward the present; negative values indicate deceleration. (Pigot et al. 2010) | Dimensionless |
| `avg_ladder` | Average ladder length | Average length of "ladder" sequences — consecutive runs of internal nodes each with exactly one immediate tip descendant. Longer average ladders indicate more pectinate (asymmetric) tree sections. | Average count of nodes |
| `max_ladder` | Maximum ladder length | Maximum length of any single ladder sequence in the tree. | Count (integer) |
| `cherries` | Number of cherries | Count of "cherries" — pairs of tips sharing an immediate common ancestor with no other tips. More cherries indicate a more balanced, symmetric tree. | Count (integer) |
| `il_number` | IL node count | Number of internal nodes that have exactly one immediate tip descendant (also called "ladder nodes" or IL nodes — Internal-Leaf nodes). | Count (integer) |
| `pitchforks` | Number of pitchforks | Count of internal nodes with exactly three tip descendants. | Count (integer) |
| `stairs` | Staircase-ness (metric 1) | Proportion of internal nodes in the tree that are ladder (IL) nodes. Ranges from 0 to 1; higher values indicate more pectinate trees. | Proportion (0–1) |
| `imbalance_steps` | Imbalance steps | Number of node-resolution steps required to transform the observed tree into a fully balanced binary tree. | Count (integer) |
| `j_one` | J1 balance statistic | Measures phylogenetic imbalance; lower values indicate more balanced trees. | Dimensionless |
| `b1` | B1 balance statistic | For each internal node, the inverse of the maximum path length (in nodes) to any tip below it, summed across all internal nodes. Higher values indicate more balanced trees. (Shao & Sokal 1990) | Dimensionless |
| `b2` | B2 balance statistic | Entropy-based measure of tree balance. Higher values indicate more balanced (even) branching. | Dimensionless |
| `area_per_pair` | Area per lineage pair | Area under the lineage-through-time (LTT) curve normalized per pair of co-existing lineages. Larger values suggest more rapid early diversification relative to late diversification. | Relative units |
| `average_leaf_depth` | Average leaf depth | Mean number of nodes (edges) separating all tips from the root. | Average node count |
| `ew_colless` | Equal-weights Colless index | Variant of the Colless index in which each internal node is weighted equally regardless of clade size, reducing the influence of large clades. | Dimensionless |
| `max_del_width` | Maximum delta width | Maximum difference in the number of co-existing lineages between any two consecutive time slices in the LTT plot. | Count (integer) |
| `max_depth` | Maximum depth | Maximum number of nodes (edges) from the root to any node in the tree. | Count (integer) |
| `max_width` | Maximum width | Maximum number of lineages co-existing at any single point in time (peak of the LTT curve). | Count (integer) |
| `rogers` | Rogers' J statistic | Measures tree imbalance; higher values indicate more asymmetric trees. (Rogers 1996) | Dimensionless |
| `stairs2` | Staircase-ness (metric 2) | A second measure of staircase-ness, computed differently from `stairs`. Higher values indicate more pectinate trees. | Dimensionless |
| `tot_coph` | Total cophenetic index | Sum of the depths (distances from the root) of the most recent common ancestors for all pairs of tips. Higher values indicate deeper internal branching events. (Mir et al. 2013) | Count-based (integer) |
| `symmetry_nodes` | Number of symmetric nodes | Count of internal nodes where both daughter clades subtend exactly the same number of tips. Higher values indicate a more balanced tree. | Count (integer) |
| `mntd` | Mean nearest taxon distance (MNTD) | For each tip, the minimum phylogenetic distance to any other tip; averaged across all tips. | Branch length units (relative to `tree_height` = 1) |
| `j_stat` | J statistic | Measures phylogenetic imbalance related to the distribution of branching times. | Dimensionless |
| `rquartet` | r\* quartet index | Measures tree shape/balance by comparing the frequency of resolved quartet topologies to expectations under a null model. | Dimensionless |
| `wiener` | Wiener index | Sum of all pairwise phylogenetic distances between all nodes (tips and internal nodes) in the tree. | Branch length units (relative to `tree_height` = 1) |
| `max_betweenness` | Maximum betweenness centrality | Maximum betweenness centrality of any single node — measures how often a node lies on the shortest path between other pairs of nodes. | Dimensionless |
| `max_closeness` | Maximum closeness centrality | Maximum closeness centrality of any single node — the inverse of the mean distance from that node to all other nodes. | Dimensionless |
| `diameter` | Tree diameter | Maximum path length (sum of branch lengths) between any two nodes in the tree. | Branch length units (relative to `tree_height` = 1) |
| `eigenvector` | Eigenvector centrality | Eigenvector centrality of the tree, measuring nodal influence based on the centrality of neighboring nodes. | Dimensionless |
| `mean_branch_length` | Mean branch length | Mean length of all branches (both internal and external/tip) in the tree. | Branch length units (relative to `tree_height` = 1) |
| `var_branch_length` | Variance in branch length | Variance in lengths across all branches (internal and external) in the tree. | Squared branch length units |
| `mean_branch_length_int` | Mean internal branch length | Mean length of internal branches (branches connecting two internal nodes) only. | Branch length units (relative to `tree_height` = 1) |
| `mean_branch_length_ext` | Mean external branch length | Mean length of external (terminal/tip) branches only. | Branch length units (relative to `tree_height` = 1) |
| `var_branch_length_int` | Variance in internal branch length | Variance in lengths of internal branches only. | Squared branch length units |
| `var_branch_length_ext` | Variance in external branch length | Variance in lengths of external (tip) branches only. | Squared branch length units |

**Missing values (`NA`):** No `NA` values are expected in tree metric columns under normal
conditions, as all metrics should be computable for any valid phylogenetic tree. If `NA` values
are present in metric columns, they indicate that a particular statistic is mathematically
undefined for the topology of that specific tree (e.g., some balance statistics require a minimum
number of tips).

---

### `empiricalSisterClades_treeStats.txt`

This tab-delimited file contains the same 70+ tree metrics as `USE_treeStats.txt` (see
descriptions above), computed for the 16 empirical phylogenetic trees. All tree metric columns
are identical in meaning; only the two identifier columns differ:

| Column | Full name | Description | Values |
|--------|-----------|-------------|--------|
| `tree` | Tree filename | Filename of the empirical phylogenetic tree (e.g., `pair1_afrotheria.tre` or `pair1_afrotheria sister clade xenarthra.tre`), identifying the taxonomic group and its role within the sister clade pair. Corresponds to tree files in `trees/empirical/`. | Character string (tree filename) |
| `simID` | Sister clade pair number | An integer identifying the sister clade pair. **Note:** In this file, `simID` does not refer to a simulation run ID. It identifies the empirical pair number (1–8), where each pair comprises two sister clades from the same taxonomic group (4 mammal pairs and 4 bird pairs). | Integer (1–8) |

All remaining columns follow the same descriptions as those given for `USE_treeStats.txt` above.

**Missing values (`NA`):** As in `USE_treeStats.txt`, any `NA` in a tree metric column indicates
that the statistic is undefined for the topology of that particular empirical tree.

---

### `process_parameter_values_and_tree_metrics_sign_corrected.csv`

This comma-delimited file combines simulation parameter values — grouped by eco-evolutionary
process and sign-corrected to standardize the direction of effect — with the full set of tree
metrics for each simulated tree. Sign correction is performed by multiplying each raw parameter
value by the corresponding `sign` value from `simulation_parameter_keys.csv`, so that after
correction, **higher values consistently indicate stronger influence of the given process**,
regardless of the original direction of the parameter in the model.

#### Identifier columns

| Column | Full name | Description | Values |
|--------|-----------|-------------|--------|
| `model` | Model abbreviation | Abbreviation of the simulation model (see [Model Abbreviations](#model-abbreviations)). | Character string |
| `model2` | Secondary model identifier | In most cases, duplicates `model`. May be used for internal grouping or figure plotting. | Character string |
| `simID` | Simulation ID | Unique integer identifier for each simulation run within a given model. | Integer |

#### Parameter name columns

These columns record the original name of each parameter (as used in the simulation model,
matching `parameterName` in `simulation_parameter_keys.csv`) for each process category. The
suffixes `1` and `2` denote the first and second parameter for that category within a given model.

| Column | Full name | Description |
|--------|-----------|-------------|
| `env1Name` | Environmental parameter 1 name | Name of the first environmental filtering parameter for this model. |
| `env2Name` | Environmental parameter 2 name | Name of the second environmental filtering parameter, if the model has one. `NA` otherwise. |
| `dis1Name` | Dispersal parameter 1 name | Name of the first dispersal parameter for this model. |
| `dis2Name` | Dispersal parameter 2 name | Name of the second dispersal parameter, if the model has one. `NA` otherwise. |
| `nic1Name` | Niche evolution parameter 1 name | Name of the first niche evolution parameter for this model. |
| `nic2Name` | Niche evolution parameter 2 name | Name of the second niche evolution parameter, if the model has one. `NA` otherwise. |
| `mut1Name` | Mutation/speciation parameter 1 name | Name of the first mutation/speciation parameter for this model. |
| `mut2Name` | Mutation/speciation parameter 2 name | Name of the second mutation/speciation parameter, if the model has one. `NA` otherwise. |
| `com1Name` | Competition parameter 1 name | Name of the first competition/community parameter for this model. |
| `com2Name` | Competition parameter 2 name | Name of the second competition/community parameter, if the model has one. `NA` otherwise. |

#### Sign-corrected parameter value columns

These columns contain the numerical values of each parameter after sign correction (raw value ×
`sign` from `simulation_parameter_keys.csv`). After correction, higher values indicate stronger
process influence.

| Column | Full name | Description | Units |
|--------|-----------|-------------|-------|
| `env1` | Environmental parameter 1 value | Sign-corrected value of the first environmental filtering parameter. | Original model parameter units (sign-corrected) |
| `env2` | Environmental parameter 2 value | Sign-corrected value of the second environmental filtering parameter. `NA` if the model has only one or no environmental parameter. | Original model parameter units (sign-corrected) |
| `dis1` | Dispersal parameter 1 value | Sign-corrected value of the first dispersal parameter. | Original model parameter units (sign-corrected) |
| `dis2` | Dispersal parameter 2 value | Sign-corrected value of the second dispersal parameter. `NA` if not applicable. | Original model parameter units (sign-corrected) |
| `nic1` | Niche evolution parameter 1 value | Sign-corrected value of the first niche evolution parameter. | Original model parameter units (sign-corrected) |
| `nic2` | Niche evolution parameter 2 value | Sign-corrected value of the second niche evolution parameter. `NA` if not applicable. | Original model parameter units (sign-corrected) |
| `mut1` | Mutation/speciation parameter 1 value | Sign-corrected value of the first mutation/speciation parameter. | Original model parameter units (sign-corrected) |
| `mut2` | Mutation/speciation parameter 2 value | Sign-corrected value of the second mutation/speciation parameter. `NA` if not applicable. | Original model parameter units (sign-corrected) |
| `com1` | Competition parameter 1 value | Sign-corrected value of the first competition/community parameter. | Original model parameter units (sign-corrected) |
| `com2` | Competition parameter 2 value | Sign-corrected value of the second competition/community parameter. `NA` if not applicable. | Original model parameter units (sign-corrected) |

#### Tree metric columns

All columns from `number_of_lineages` onward are identical in meaning to those described in
[`USE_treeStats.txt`](#use_treestatstxt) above.

**Missing values (`NA`):** `NA` in any parameter name or value column indicates one of the
following: (1) the model has no parameter assigned to that process category; or (2) the model
has only one parameter for that category, making the second (`...2`) column `NA`. Note that `NA`
does not necessarily mean a process was absent from the model — it may reflect that no single
model parameter maps cleanly to that process category in the standardized framework. For tree
metric columns, `NA` indicates the metric could not be computed for that tree topology.

---

### `scaled_predictions_empirical_trees.csv`

This tab-delimited file contains random forest–based inferences of eco-evolutionary process
strength for each empirical tree, derived from models trained on simulated trees. For each
combination of empirical tree, simulation model, and process category, a predicted parameter
value and a rescaled version of that value are provided.

| Column | Full name | Description | Units / Values |
|--------|-----------|-------------|----------------|
| `tree` | Tree name | Name of the empirical tree and taxon group (e.g., `pair1_afrotheria`). Corresponds to tree filenames in `trees/empirical/` (without the `.tre` extension). | Character string |
| `model` | Model abbreviation | Abbreviation of the simulation model used to train the random forest (see [Model Abbreviations](#model-abbreviations)). | Character string |
| `predictor` | Process-parameter predictor | The process-parameter identifier for which a value was inferred (e.g., `com1`, `env1`, `dis1`). Corresponds directly to the sign-corrected parameter value column names in `process_parameter_values_and_tree_metrics_sign_corrected.csv`. | Character string |
| `value` | Predicted parameter value | The parameter value inferred for the empirical tree by the random forest model trained on that model's simulations. | Scale varies by model and process. |
| `minVal` | Minimum simulation parameter value | The minimum value of the parameter (on the same scale as `value`) used across all simulations for this model and process. Serves as the lower bound for computing `scaledValue`. | Same units as `value` |
| `maxVal` | Maximum simulation parameter value | The maximum value of the parameter (on the same scale as `value`) used across all simulations for this model and process. Serves as the upper bound for computing `scaledValue`. | Same units as `value` |
| `pair` | Sister clade pair number | Integer (1–8) identifying the empirical sister clade pair. | Integer (1–8) |
| `clade` | Clade within pair | Integer (1 or 2) identifying which of the two clades within a sister pair this row describes. | Integer (1 or 2) |
| `scaledValue` | Rescaled predicted parameter value | The predicted parameter value rescaled to [0, 1] using the formula: (`value` − `minVal`) / (`maxVal` − `minVal`). A value of 0 corresponds to the weakest inferred process strength relative to the simulation parameter space; 1 corresponds to the strongest. Values slightly outside [0, 1] are possible if the empirical tree's inferred value falls outside the range used in simulations. | Dimensionless (nominally 0–1) |
| `newLabel` | Plot label | Model label used for figure plotting. In most cases identical to `model`. | Character string |
| `color` | Plot color | Hexadecimal color code (e.g., `#000000` for black) assigned to the model for use in manuscript figures. | Hex color string (e.g., `#RRGGBB`) |

