# Latitudinal CNDD analysis by Hülsmann et al.

This repository contains the code to reproduce the analyses by Hülsmann et al., which uses repeated census data from twenty-three large forest sites around the globe to analyze latitudinal patterns in conspecific negative density dependence (CNDD).

The repository holds the folders [code](/code) (all Rscripts), [data_prep](/data_prep) (input and derived output data to run the analyses), and [out](/out) (results of mortality models, meta-regressions, and more general tables).

Because tree data from most ForestGEO sites cannot be made freely available, we generated 'fake sites' based on small samples of the dataset from Barro Colorado Island (BCI, Panama) available [here](https://datadryad.org/stash/dataset/doi:10.15146/5xcp-0d46). These fake sites were artificially aligned along a latitudinal gradient to be able to run the full analysis pipeline.

## Methodology

We analyzed latitudinal patterns in conspecific negative density dependence (CNDD) following a three-step approach: First, we fitted species-site-specific mortality models from repeated observations of individual trees. Second, we used these models to quantify CNDD for each species and site using an estimator designed to maximize robustness, comparability, and relevance for fitness and stabilization. Third, we used meta-regressions to explore latitudinal patterns in CNDD. Robustness of the analysis pipeline was validated by model diagnostics and randomization.

## Contact

Lisa Hülsmann, EASI lab, University of Bayreuth, [lisa.huelsmann\@uni-bayreuth.de](mailto:lisa.huelsmann@uni-bayreuth.de){.email}
