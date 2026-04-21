# BA

This repository contains accompanying and data for the manuscript

> Coloring for dispersion: A polynomial-time algorithm for cardinality-constrained 2-anticlustering

Directory `R/` has the R code to reproduce the Figures 1-4 illustrated in the paper. The C++ code implementing the 2-COL-CC algorithm (and which is called from an R interface function in `R/sc2col.R` and `R sc2colheap.R`) resides in `src`. It was written by Lin Mu. Directory `data/` has the data that was generated for the paper by running the scripts in the Directory `R/`. To run the scripts, the R working directory must reside in the uppermost directory of this repository (i.e., where the README.md is).

## Dependencies

Data was generated using R version 4.5.3 and the anticlust package version 0.8.13. To generate Figure 1 and Figure 3, we also used [Gurobi](https://www.gurobi.com/) (version 13.01), which distributes an `R` package gurobi, which we used as backend solver of the ILP model (version 13.0-1).
