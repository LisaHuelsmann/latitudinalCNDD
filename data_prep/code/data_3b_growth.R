

# Script restructures tree data to be used for growth analyses


### includes the steps
# import tree data
# plausibility of sequences
# restructure of tree data to data.frame 
# select observations for growth analyses
# exports tree data


library(lubridate)





# add site to each list element
tree = Map(cbind, tree, site = site)

# add census to each list element
tree = Map(cbind, tree, census = names(tree))









# Prepare data for growth analyses ----------------------------------------


for (census in 1:(length(tree)-1)) {
  
  tree[[census]]$surv = ifelse(grepl("A", tree[[census]]$status), 1, 0)            # status at t=0
  tree[[census]]$surv_next = ifelse(grepl("A", tree[[census+1]]$status), 1, 0)     # status at t=1

  tree[[census]]$incr = tree[[census+1]]$dbh-tree[[census]]$dbh
  
  tree[[census]]$interval = time_length(difftime(tree[[census+1]]$date
                                                 , tree[[census]]$date), "years")  # length of interval
}



# generate dataframe with all censuses
tree_growth = as.data.frame(data.table::rbindlist(tree, fill = T))




# select only observations of interest ------------------------------------


# 1) living and then again living
tree_growth = tree_growth[tree_growth$surv == 1 
                      & tree_growth$surv_next == 1
                      & !is.na(tree_growth$surv)
                      & !is.na(tree_growth$surv_next), ]


# # 2) observations inside buffer
# outside = tree_growth$gx < plotrange[1, 1]+max(dist_def) |
#   tree_growth$gx > plotrange[1, 2]-max(dist_def) |
#   tree_growth$gy < plotrange[2, 1]+max(dist_def) |
#   tree_growth$gy > plotrange[2, 2]-max(dist_def)
# 
# # 2a) additional outside condition for Ituri 30m around x=200m
# if (site %in% c("edo", "len")) {
#   outside = outside |
#     (tree_growth$gx < (200 + max(dist_def)) & 
#        tree_growth$gy > (200 - max(dist_def)))
# }
# tree_growth = tree_growth[!outside, ]


# 3) observations with dbh>10mm at t=0
tree_growth = tree_growth[!is.na(tree_growth$dbh), ]
tree_growth = tree_growth[tree_growth$dbh >= 10, ]
# summary(tree_growth)
# table(tree_growth$mort_next)


# 4) observations with growth 
tree_growth = tree_growth[!is.na(tree_growth$incr), ]




# Save file for site ------------------------------------------------------


save(tree_growth, file = paste0(path_output, "data_3b_growth/", site, "_tree_3b_growth.Rdata"))


rm(list = ls()[!grepl("site", ls()) & !grepl("path", ls())])




