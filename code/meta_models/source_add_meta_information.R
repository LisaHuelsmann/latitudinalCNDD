


# Code to add meta information at the site and species level
# information can be added to a single data.frame or several data.frames 


# requires that the global output objects are indicated by _global


# adds the following variables
#    latitude
#    climate
#    abundance
#    stature
#    growth rate
#    mortality rate





# Load required meta data objects -----------------------------------------


# Meta info at the site level
Site_table = readxl::read_xlsx(paste0(path_input, "global_metainfo/plot_sites_information.xlsx"), sheet = 1)


# Meta info at the species level
for (type in c("abundances", "growth", "mortality", "stature")) {
  
  files <- data.frame(dir = dir(paste0(path_output, "meta_", type, "/"), 
                                pattern = "*.Rdata", full.names = T))
  files$site = dir(paste0(path_output, "meta_", type, "/"), pattern = "*.Rdata")
  files$site = gsub(paste0("_", type, ".Rdata"), "", files$site)
  
  assign(paste0(type, "_meta"), data.frame())
  
  for (i in 1:nrow(files)) {
    load(files$dir[i])
    if (type == "mortality") temp = get("mort") else temp = get(type) # mortality objects have non-standard name
    temp$site = files$site[i] # add site
    temp$site_sp = paste(temp$site, temp$sp, sep = "_") # add site species identifier
    assign(paste0(type, "_meta"), rbind(get(paste0(type, "_meta")), temp))
  } 
}



# Add meta information to output objects -----------------------------------


for (i in output_objects) {
  
  temp = get(paste0(i, "_global"))
  
  
  ### Site specific variables
  
  # add latitude
  temp$latitude = Site_table$lat[match(temp$site, Site_table$ID)]
  
  # add climate variables
  temp$MAT = Site_table$mat[match(temp$site, Site_table$ID)]
  temp$MAP = Site_table$map[match(temp$site, Site_table$ID)]
  temp$PET = Site_table$PET_annual[match(temp$site, Site_table$ID)]
  
  
  ### Species specific variables
  
  # add site species identifier
  temp$site_sp = paste(temp$site, temp$sp, sep = "_")
  nsp_global$site_sp = paste(nsp_global$site, nsp_global$sp, sep = "_")
  
  
  ## Abundances
  
  # add abundances per species (N and BA)
  temp$abundance = abundances_meta$Nha[match(temp$site_sp, abundances_meta$site_sp)]
  temp$abundance_BA = abundances_meta$BAha[match(temp$site_sp, abundances_meta$site_sp)]  
  
  # add average abundances for rare species
  abundances_meta$rare = nsp_global$rare_stature[match(abundances_meta$site_sp, nsp_global$site_sp)]
  abundances_meta %>% 
    filter(!is.na(rare)) %>%
    mutate(site_sp = paste(site, rare, sep = "_")) %>% 
    group_by(site_sp) %>% 
    summarise(Nha = mean(Nha),
              BAha = mean(BAha)) -> abundances_rare
  temp$abundance[grepl("Rare", temp$sp)] = abundances_rare$Nha[match(temp$site_sp[grepl("Rare", temp$sp)], abundances_rare$site_sp)]
  temp$abundance_BA[grepl("Rare", temp$sp)] = abundances_rare$BAha[match(temp$site_sp[grepl("Rare", temp$sp)], abundances_rare$site_sp)]
  
  # add abundances also to stature, growth and mortality for abundance weighted means
  stature_meta$abundance = abundances_meta$Nha[match(stature_meta$sp, abundances_meta$sp)]  
  growth_meta$abundance = abundances_meta$Nha[match(growth_meta$sp, abundances_meta$sp)]  
  mortality_meta$abundance = abundances_meta$Nha[match(mortality_meta$sp, abundances_meta$sp)]  
  stature_meta$abundance_BA = abundances_meta$BAha[match(stature_meta$sp, abundances_meta$sp)]  
  growth_meta$abundance_BA = abundances_meta$BAha[match(growth_meta$sp, abundances_meta$sp)]  
  mortality_meta$abundance_BA = abundances_meta$BAha[match(mortality_meta$sp, abundances_meta$sp)]  
  
  
  ## Stature
  
  # add max dbh per species
  temp$dbh_q90 = stature_meta$dbh_q90[match(temp$site_sp, stature_meta$site_sp)]
  
  # add max dbh for rare species (weighted mean)
  stature_meta$rare = nsp_global$rare_stature[match(stature_meta$site_sp, nsp_global$site_sp)]
  stature_meta %>% 
    filter(!is.na(rare)) %>%
    mutate(site_sp = paste(site, rare, sep = "_")) %>% 
    group_by(site_sp) %>% 
    summarise(dbh_q90 = weighted.mean(x = dbh_q90, w = abundance, na.rm = T)) -> stature_rare
  temp$dbh_q90[grepl("Rare", temp$sp)] = stature_rare$dbh_q90[match(temp$site_sp[grepl("Rare", temp$sp)], stature_rare$site_sp)]
  
  
  # add stature per species
  temp$stature = stature_meta$stature[match(temp$site_sp, stature_meta$site_sp)]
  
  # add stature per species group
  temp$stature[grepl("Rare_", temp$sp)] = gsub("Rare_", "", temp$sp[grepl("Rare_", temp$sp)])
  
  
  ## Growth
  
  # add growth per species
  temp$incr_median = growth_meta$incr_median[match(temp$site_sp, growth_meta$site_sp)]
  
  # add growth for rare species (weighted mean)
  growth_meta$rare = nsp_global$rare_stature[match(growth_meta$site_sp, nsp_global$site_sp)]
  growth_meta %>% 
    filter(!is.na(rare)) %>%
    mutate(site_sp = paste(site, rare, sep = "_")) %>% 
    group_by(site_sp) %>% 
    summarise(incr_median = weighted.mean(x = incr_median, w = abundance, na.rm = T)) -> growth_rare
  temp$incr_median[grepl("Rare", temp$sp)] = growth_rare$incr_median[match(temp$site_sp[grepl("Rare", temp$sp)], stature_rare$site_sp)]
  
  
  ## Mortality
  
  # add mortality per species
  temp$mort_rate = mortality_meta$mort_rate[match(temp$site_sp, mortality_meta$site_sp)]
  
  # add mortality for rare species (weighted mean)
  mortality_meta$rare = nsp_global$rare_stature[match(mortality_meta$site_sp, nsp_global$site_sp)]
  mortality_meta %>% 
    filter(!is.na(rare)) %>%
    mutate(site_sp = paste(site, rare, sep = "_")) %>% 
    group_by(site_sp) %>% 
    summarise(mort_rate = weighted.mean(x = mort_rate, w = abundance, na.rm = T)) -> mortality_rare
  temp$mort_rate[grepl("Rare", temp$sp)] = mortality_rare$mort_rate[match(temp$site_sp[grepl("Rare", temp$sp)], stature_rare$site_sp)]
  
  # remove site_sp column
  temp$site_sp = NULL
  
  assign(paste0(i, "_global"), temp)
  
}



