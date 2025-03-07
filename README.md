# partial-clustering-deffs
R code for a simulation study to validate GEE-based design effects for partially clustered trials. The accompanying manuscript is currently in preparation.

# Files
**R/**: R functions and simulation scenarios

  - partial_clus_deffs.R: Functions for calculating the GEE-based design effects that are derived in the manuscript
  - partial_clus_functions.R: Functions for generating data, analysing the datasets, and generating the tables of results
  - scenarios_cont.xlsx: Parameters for the simulation scenarios for continuous outcomes
  - scenarios_bin.xlsx: Parameters for the simulation scenarios for binary outcomes

**scripts/**: The R scripts for running the simulation study

  - runsim_cont.R: Simulation code for continuous outcomes
  - runsim_bin.R: Simulation code for binary outcomes
  - figures.R: Figures of the theoretical design effect functions
  - runsim_cont_multinom.R: Simulation code for sensitivity analysis 2 for continuous outcomes
  - runsim_bin_multinom.R: Simulation code for sensitivity analysis 2 for binary outcomes
