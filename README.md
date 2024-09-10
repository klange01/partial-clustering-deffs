# partial-clustering-deffs
R code for a simulation study to validate GEE-based design effects for partially clustered trials. The accompanying manuscript is currently in preparation.

# Files
**R/**: R functions and simulation scenarios

  - partial_clus_deffs.R: Functions for calculating the GEE-based design effects that are derived in the manuscript
  - partial_clus_functions.R: Functions for generating simulated data, analysiing the datasets, and generating the tables of results
  - scenarios_cont.xlsx: Parameters for the simulation scenarios for continuous outcomes
  - scenarios_bin.xlsx: Parameters for the simulation scenarios for binary outcomes

**scripts/**: The R scripts for running the simulation study

  - runsim_cont.R: Simulation study for continuous outcomes
  - runsim_bin.R: Simulation study for binary outcomes
  - figures.R: Figures of the theoretical design effect functions



