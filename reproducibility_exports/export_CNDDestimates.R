

# Export species-annonymized CNDD estimates

library(dplyr)




# Define location of input and output -------------------------------------

# if not stated otherwise, load and save from/in repo folders


# Data prep input
if (!exists("path_input")) path_input = "data_prep/input/"

# Data prep output
if (!exists("path_output")) path_output = "data_prep/output/"

# Analysis outputs
if (!exists("path_mortality")) path_mortality = "out/mortality_models/"
if (!exists("path_meta")) path_meta = "out/meta_models/"

# Tables
if (!exists("path_tables")) path_meta = "out/tables/"







# DERIVED CNDD ESTIMATES --------------------------------------------------



run = "main"




# Load global mortality runs ----------------------------------------------


# load combined results from run
load(paste0(path_mortality, run, "/global_mortality.Rdata"))



# Resemble data prep from source_meta_models.R ----------------------------

# chose available results
terms = unique(AMEsamples_global$term)[startsWith(unique(AMEsamples_global$term), "con_")]
types = c("AME", "rAME")
changes = unique(AMEsamples_global$change)
sites = sort(unique(AMEsamples_global$site))





# Transformations ---------------------------------------------------------


# functions to transform response variables

# for AME, based on residuals
# trans_AME = function(x) sign(x)*(abs(x)^(1/3))
# backtrans_AME = function(x) sign(x)*((abs(x))^3)
trans_AME = function(x) x
backtrans_AME = function(x) x


# for rAME, based on common practice for RR
trans_rAME = function(x) log(x+1)
backtrans_rAME = function(x) exp(x)-1






# Uncertainty and significance --------------------------------------------


# point estimates from MLE
# significance and std.errors from samples (transformed)

for (type in types) {
  
  # get sample and aggregated object with MLE
  temp = get(paste0(type, "samples_global"))
  
  # get transformation
  trans = get(paste0("trans_", type))
  
  temp %>% 
    group_by(across(-c(estimate, MLE))) %>%     # group by all except estimate and MLE, results in grouping per site_sp x change with n = iter
    mutate(q2.5 = stats::quantile(estimate, 0.025),    #  quantiles for significance
           q97.5 = stats::quantile(estimate, 0.975),
           transSamples = trans(estimate),
           transMLE = trans(MLE)) %>% 
    summarise(std.error = diff(stats::quantile(transSamples, c(pnorm(-1), pnorm(1))))/2,
              estimate = unique(transMLE),
              significant = unique(q2.5 > 0 | q97.5 < 0)) %>% 
    as.data.frame() -> temp
  
  
  # store site as factor
  temp$site = as.factor(temp$site)
  
  
  assign(paste0(type, "sums_global"), temp)
  
}




# Add meta information -----------------------------------------------------


output_objects = gsub("_global", "", ls()[grepl("_global", ls())])

source("code/meta_models/source_add_meta_information.R", local = T)




# Anonymization -----------------------------------------------------------

# remove row.names
# change column sp
# remove site_sp

# generate fake letters
sp_names = expand.grid(letters, letters, letters)
sp_names = sp_names[, c(3, 2, 1)]
sp_names = apply(sp_names, 1, paste, collapse = "")

# use nsp_global to make species names anonymous
nsp_global %>% 
  group_by(site) %>% 
  mutate(sp_new = sp_names[1:n()],
         site_sp = paste(site, sp, sep = "_")) -> nsp_new


# change in all output objects
output_objects = ls()[grepl("_global", ls())]

for (i in output_objects) {
  
  temp = get(i)
  
  # remove rownames
  rownames(temp) = NULL
  
  # change sp
  temp$site_sp = paste(temp$site, temp$sp, sep = "_")
  sel = !grepl("Rare", temp$sp) # non rare
  temp$sp[sel] = nsp_new$sp_new[match(temp$site_sp[sel], nsp_new$site_sp)]
  temp$site_sp = NULL
  
  # overwrite old object
  assign(i, temp)
  
}




# Save result -------------------------------------------------------------


save(list = c("AMEsums_global", "rAMEsums_global", "sums_global", "nsp_global")
     , file = "reproducibility_exports/global_mortality.Rdata")

