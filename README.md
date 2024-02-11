# Latitudinal CNDD analysis by Hülsmann et al.

This repository contains the code to reproduce the analyses by Hülsmann et al. "Latitudinal patterns in stabilizing density dependence of forest communities", which uses repeated census data from 23 large forest sites from the [ForestGEO network](https://forestgeo.si.edu/) to analyze latitudinal patterns in stabilizing conspecific negative density dependence (CNDD).

The repository holds the folders [**code**](/code) (all Rscripts), [**data_prep**](/data_prep) (input and derived output data to run the analyses) and [**reproducibility_exports**](/reproducibility_exports) (CNDD estimates to reproduce the main meta-regressions). Moreover, when running the analysis pipeline, a folder **out** will be created with the results of mortality models, meta-regressions, and tables for the Extended Data.

Because tree data from most ForestGEO sites cannot be made freely available, we generated 'fake sites' based on small samples of the dataset from Barro Colorado Island (BCI, Panama) available [here](https://datadryad.org/stash/dataset/doi:10.15146/5xcp-0d46). These fake sites were artificially aligned along a latitudinal gradient to be able to run the full analysis pipeline.

## Methodology

We used repeated census data from 23 large forest sites around the globe to analyze latitudinal patterns in stabilizing conspecific negative density dependence (CNDD) following a three-step approach: First, we fitted species-site-specific mortality models from repeated observations of individual trees. Second, we used these models to quantify CNDD for each species and site using an estimator designed to maximize robustness, comparability, and relevance for fitness and stabilization. Third, we used meta-regressions to explore three distinct latitudinal patterns in CNDD derived from the hypothesis that CNDD is more influential for maintaining local tree species diversity in the tropics. Robustness of the analysis pipeline was validated by model diagnostics and randomization.

This approach is based on recently developed best-practice statistical methods for estimating CNDD. Crucially, the use of dynamic mortality data allowed us to avoid the statistical pitfalls of previous CNDD studies, in particular analyses of the static relationship of number of saplings to number of adults, where the null hypothesis is a positive linear relationship but regression dilution flattens this relationship and thus biases analyses towards finding CNDD, especially for rare species. By fitting mortality models where the null hypothesis is no relationship between survival and number of conspecific neighbors, we ensure that any regression dilution has a conservative effect by reducing CNDD estimates. We also addressed other recently identified limitations of CNDD analyses, namely non-linear and saturating CNDD, the comparability of CNDD among species and sites, and the extent to which CNDD estimates are meaningful for stabilization and species coexistence.

## System and R Packages

All analyses were conducted in R version 4.2.1 on Ubuntu 20.04.4 LTS, but the code has also been tested on macOS Monterey 12.6.

Major steps of the analyses were carried out with the following R packages:

-   mgcv (Version 1.8-40)

-   metafor (Version 3.4-0)

-   DHARMa (Version 0.4.6)

The expected runtime for the demo dataset on a normal computer is approx. a few minutes for the data preparation, approx. 24h for the mortality models (the grid search is computationally quite expensive), 20min for the meta-regressions, and a few seconds for the tables. Runtime can be reduced when using more cores for the mortality models and meta-regressions.

## Contact

Lisa Hülsmann, EASI lab, University of Bayreuth, [lisa.huelsmann\@uni-bayreuth.de](mailto:lisa.huelsmann@uni-bayreuth.de){.email}
