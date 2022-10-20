


# Code to add meta information at the site and species level
# information can be added to a single dataframe or several dataframes 
# this is done separately per site!

# requires the following objects
#    site
#    output_objects
#    and the _global objects to which the output_objects are assigned

# add the following variables
#    site
#    latitude
#    abundance
#    abundance for rare species
#    growth form / stature



# Load required meta data objects -----------------------------------------


# Meta info sites
Site_table = readxl::read_xlsx("data_prep/plot_sites_information.xlsx", sheet = 1)

# Abundances
load(paste0("data_prep/meta_abundances/", site, "_abundances.Rdata"))

# Clusters
load(paste0("data_prep/meta_clusters/", site, "_clusters.Rdata"))

# Stature
load(paste0("data_prep/meta_stature/", site, "_stature.Rdata"))
stature = as.data.frame(stature)

# Growth
load(paste0("data_prep/meta_growth/", site, "_growth.Rdata"))

# Mortality
load(paste0("data_prep/meta_mortality/", site, "_mortality.Rdata"))




# add additional info to all output style objects
for (i in output_objects) {
  
  temp = get(i)
  temp$site = site
  
  
  ### Site specific variables
  
  # add latitude
  temp$latitude = Site_table$lat[match(temp$site, Site_table$ID)]
  
  # add climate variables
  temp$MAT = Site_table$mat[match(temp$site, Site_table$ID)]
  temp$MAP = Site_table$map[match(temp$site, Site_table$ID)]
  temp$PET = Site_table$PET_annual[match(temp$site, Site_table$ID)]
  
  
  ### Species specific variables
  
  # add abundances per species (N)
  temp$abundance = abundances$Nha[match(temp$sp, abundances$sp)]  
  
  # add mean abundance for rare species
  temp$abundance[temp$sp == "Rare_tree"] = mean(abundances$Nha[abundances$sp %in% 
                                                                 nsp$sp[which(nsp$rare_stature == "Rare_tree")]])
  temp$abundance[temp$sp == "Rare_shrub"] = mean(abundances$Nha[abundances$sp %in% 
                                                                  nsp$sp[which(nsp$rare_stature == "Rare_shrub")]])
  
  # add abundances also to stature, growth and mort for abundance weighted means
  stature$abundance = abundances$Nha[match(stature$sp, abundances$sp)]  
  growth$abundance = abundances$Nha[match(growth$sp, abundances$sp)]  
  mort$abundance = abundances$Nha[match(mort$sp, abundances$sp)]  
  
  
  # add abundances per species (BA)
  temp$abundance_BA = abundances$BAha[match(temp$sp, abundances$sp)]  
  
  # add mean abundance for rare species
  temp$abundance_BA[temp$sp == "Rare_tree"] = mean(abundances$BAha[abundances$sp %in% 
                                                                 nsp$sp[which(nsp$rare_stature == "Rare_tree")]])
  temp$abundance_BA[temp$sp == "Rare_shrub"] = mean(abundances$BAha[abundances$sp %in% 
                                                                  nsp$sp[which(nsp$rare_stature == "Rare_shrub")]])
  
  # add abundances also to stature, growth and mort for abundance weighted means
  stature$abundance_BA = abundances$BAha[match(stature$sp, abundances$sp)]  
  growth$abundance_BA = abundances$BAha[match(growth$sp, abundances$sp)]  
  mort$abundance_BA = abundances$BAha[match(mort$sp, abundances$sp)]  
  
  
  
  # add cluster per species
  temp$cluster = clusters$clust_hier[match(temp$sp, clusters$sp)]  
  
  # add mean abundance for rare species
  temp$cluster[temp$sp == "Rare_tree"] = names(sort(table(clusters$clust_hier[clusters$sp %in% 
                                                                                nsp$sp[which(nsp$rare_stature == "Rare_tree")]]), decreasing = T)[1])
  temp$cluster[temp$sp == "Rare_shrub"] = names(sort(table(clusters$clust_hier[clusters$sp %in% 
                                                                                 nsp$sp[which(nsp$rare_stature == "Rare_shrub")]]), decreasing = T)[1])
  
  
  # add max dbh variables
  for (vars in names(stature)[grepl("dbh_", names(stature))]) {
    
    # add stature variables per species
    temp[, vars] = stature[match(temp$sp, stature$sp), vars]  
    
    # add mean max dbh variables for rare species
    temp[temp$sp == "Rare_tree", vars] = weighted.mean(stature[stature$sp %in% nsp$sp[which(nsp$rare_stature == "Rare_tree")], vars]
                                                       , w = stature[stature$sp %in% nsp$sp[which(nsp$rare_stature == "Rare_tree")], "abundance"]
                                                       , na.rm = T)
    temp[temp$sp == "Rare_shrub", vars] = weighted.mean(stature[stature$sp %in% nsp$sp[which(nsp$rare_stature == "Rare_shrub")], vars]
                                                        , w = stature[stature$sp %in% nsp$sp[which(nsp$rare_stature == "Rare_shrub")], "abundance"]
                                                        , na.rm = T)
    
  }
  
  # add stature per species
  temp$stature = stature$stature[match(temp$sp, stature$sp)]  
  
  # add stature per species group
  temp$stature[grepl("Rare_", temp$sp)] = gsub("Rare_", "", temp$sp[grepl("Rare_", temp$sp)])
  

  
  
  # add growth variables
  for (incr in names(growth)[grepl("incr_", names(growth))]) {
    
    # add growth per species
    temp[, incr] = growth[match(temp$sp, growth$sp), incr]  
    
    # add mean growth for rare species
    temp[temp$sp == "Rare_tree", incr] = weighted.mean(growth[growth$sp %in% nsp$sp[which(nsp$rare_stature == "Rare_tree")], incr]
                                                       , w = growth[growth$sp %in% nsp$sp[which(nsp$rare_stature == "Rare_tree")], "abundance"]
                                                       , na.rm = T)
    temp[temp$sp == "Rare_shrub", incr] = weighted.mean(growth[growth$sp %in% nsp$sp[which(nsp$rare_stature == "Rare_shrub")], incr]
                                                        , w = growth[growth$sp %in% nsp$sp[which(nsp$rare_stature == "Rare_shrub")], "abundance"]
                                                        , na.rm = T)
    
  }
  
  # add mortality variables
  for (vars in names(mort)[grepl("mort_", names(mort))]) {
    
    # add mort per species
    temp[, vars] = mort[match(temp$sp, mort$sp), vars]  
    
    # add mean mort for rare species
    temp[temp$sp == "Rare_tree", vars] = weighted.mean(mort[mort$sp %in% nsp$sp[which(nsp$rare_stature == "Rare_tree")], vars]
                                                       , w = mort[mort$sp %in% nsp$sp[which(nsp$rare_stature == "Rare_tree")], "abundance"]
                                                       , na.rm = T)
    temp[temp$sp == "Rare_shrub", vars] = weighted.mean(mort[mort$sp %in% nsp$sp[which(nsp$rare_stature == "Rare_shrub")], vars]
                                                        , w = mort[mort$sp %in% nsp$sp[which(nsp$rare_stature == "Rare_shrub")], "abundance"]
                                                        , na.rm = T)
    
  }
  

  assign(paste0(i, "_global"), rbind(get(paste0(i, "_global")), temp))

}



