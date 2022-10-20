

# Script restructures tree data to be used for mortality analyses


### includes the steps
# import tree data
# plausibility of sequences
# restructure of tree data to data.frame 
# select observations for mortality analyses
# exports tree data


library(lubridate)





# add site to each list element
tree = Map(cbind, tree, site = site)

# add census to each list element
tree = Map(cbind, tree, census = names(tree))









# Restructure data for mortality analyses ---------------------------------


for (census in 1:(length(tree)-1)) {
  tree[[census]]$surv = ifelse(grepl("A", tree[[census]]$status), 1, 0)            # status at t=0
  tree[[census]]$surv_next = ifelse(grepl("A", tree[[census+1]]$status), 1, 0)     # status at t=1
  
  tree[[census]]$mort = ifelse(grepl("A", tree[[census]]$status), 0, 1)            # status at t=0
  tree[[census]]$mort_next = ifelse(grepl("A", tree[[census+1]]$status), 0, 1)     # status at t=1
  
  tree[[census]]$interval = time_length(difftime(tree[[census+1]]$date
                                                 , tree[[census]]$date), "years")  # length of interval
}



# generate dataframe with all censuses to be used for mortality analyses
tree_mort = as.data.frame(data.table::rbindlist(tree, fill = T))




# select only observations of interest ------------------------------------


# 1) living and then living or dead
tree_mort = tree_mort[tree_mort$surv == 1 
                      & !is.na(tree_mort$surv)
                      & !is.na(tree_mort$surv_next), ]

# 2) observations inside buffer
outside = tree_mort$gx < plotrange[1, 1]+max(dist_def) |
  tree_mort$gx > plotrange[1, 2]-max(dist_def) |
  tree_mort$gy < plotrange[2, 1]+max(dist_def) |
  tree_mort$gy > plotrange[2, 2]-max(dist_def)

# 2a) additional outside condition for Ituri 30m around x=200m
if (site %in% c("edo", "len")) {
  outside = outside |
    (tree_mort$gx < (200 + max(dist_def)) & 
       tree_mort$gy > (200 - max(dist_def)))
}
tree_mort = tree_mort[!outside, ]


# 3) observations with dbh>10mm at t=0
tree_mort = tree_mort[!is.na(tree_mort$dbh), ]
tree_mort = tree_mort[tree_mort$dbh >= 10, ]




# Chose columns to make a separate file for exp and expn dens ------------


exp_cols = colnames(tree_mort)[grepl("exp_", colnames(tree_mort))]
expn_cols = colnames(tree_mort)[grepl("expn_", colnames(tree_mort))]
others = colnames(tree_mort)[!colnames(tree_mort) %in% c(exp_cols, expn_cols)]

# function that saves selected columns unter the name tree_mort
save_tree_mort = function(dat, style) {
  tree_mort = dat
  save(tree_mort, file = paste0(path_output, "data_3a_mortality/", site, "_tree_3a_mortality_", style, ".Rdata"))
}





# Save files for site -----------------------------------------------------


save_tree_mort(tree_mort[, c(others, exp_cols)]
               , style = "exp")
save_tree_mort(tree_mort[, c(others, expn_cols)]
               , style = "expn")


rm(list = ls()[!grepl("site", ls()) & !grepl("path", ls())])




