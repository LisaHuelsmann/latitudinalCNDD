# Shared data for reproducibility

This folder contains an Rdata file that includes the R objects necessary to calibrate the main meta-regressions with the species-specific CNDD estimates (global_mortality.Rdata). The file includes the following objects

-   nsp_global: species-site-level information on number of observations, unique values and range of conspecific density (includes all species available per site)

-   sums_global: summaries of the generalized additive models for mortality (includes only species for which species-specific models could be fitted, and rare species groups)

-   AMEsums_global: estimates and standard errors for *absolute* average marginal effects for three different changes in conspecific density (equilbrium, invasion, and interquantile range (see methods of Hülsmann et al.) (includes only species for which species-specific models could be fitted, and rare species groups)

-   rAMEsums_global: estimates and standard errors for *relative* average marginal effects for three different changes in conspecific density (equilbrium, invasion, and interquantile range (see methods of Hülsmann et al.) (includes only species for which species-specific models could be fitted, and rare species groups)

The objects can be used to run the code in the Rscript source_meta_models.R starting L141. Species names were anonymized.
