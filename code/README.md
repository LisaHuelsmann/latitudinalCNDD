# Code to run all steps of the latitudinal CNDD analysis

It includes the following steps, which can be completed by executing the respective R script.

1)  Data preparation [1_data_prep.R](/code/1_data_prep.R)
2)  Mortality models [2_mortality_models.R](/code/2_mortality_models.R)
3)  Meta models [3_meta_models.R](/code/3_meta_models.R)
4)  Tables [4_tables.R](/code/4_tables.R)

Most of these scripts source other lower level scripts for additional functions and analysis steps, which can be found in the respective [code folder](/code/).

The results that are produced during each step of the analysis are stored either in the folder [data_prep](/data_prep) or in the folder [out](/out) (mortality models, meta-regressions,d and tables).
