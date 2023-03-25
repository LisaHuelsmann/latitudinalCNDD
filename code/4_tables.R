


# Generate tables for Extended Data


library(dplyr)
library(flextable)
library(officer)


# create folder for output
dir.create(path_tables)




# Load data ---------------------------------------------------------------


# select run
run = "main"

# create folder
dir.create(paste0(path_tables, run))

# mortality modeling results
load(paste0(path_mortality, run, "/global_mortality.Rdata"))

# sites
sites = sort(unique(sums_global$site))

# auxiliary information
Site_table = readxl::read_xlsx(paste0(path_input, "global_metainfo/plot_sites_information.xlsx"), sheet = 1)
Site_table = Site_table %>%
  dplyr::filter(ID %in% sites)







# Extended Data Table 1 ---------------------------------------------------



# Variables: 
# Site name, latitude, longitude, mean temperature, mean precipitation, 
# plot size, number of censuses, census years


# generated based on Site_table
# with additional info from 
#   _tree_1_metainfo.Rdata (censuses)

# select columns
ext_data_table_1 = Site_table[, c("site", "lat", "long", "mat", "map", "size_ha", "ID")]
ext_data_table_1$size_ha[ext_data_table_1$ID == "ucsc"] = 6 # subset was used

# sort by full site name
ext_data_table_1 = ext_data_table_1[order(ext_data_table_1$site), ]

# round lat and long
ext_data_table_1$lat = format(round(ext_data_table_1$lat, digits=2), nsmall = 2) 
ext_data_table_1$long = format(round(ext_data_table_1$long, digits=2), nsmall = 2) 


# add number of censuses and census years from _tree_1_metainfo.Rdata
ext_data_table_1$n_censuses = NA
ext_data_table_1$census_years = NA
for (site in sites) {
  load(paste0(path_output, "data_1_metainfo/", site, "_tree_1_metainfo.Rdata"))
  ext_data_table_1$n_censuses[ext_data_table_1$ID == site] = length(tree)
  get_year = function(x) names(sort(table(format(x$date, format="%Y")), decreasing = T))[1] # most common name per census
  ext_data_table_1$census_years[ext_data_table_1$ID == site] = paste(unlist(lapply(tree, get_year)), collapse = ", ")
}

ext_data_table_1 = ext_data_table_1 %>% 
  select(-ID)

# Formatting
ext_data_table_1 %>%
  flextable() %>%
  set_table_properties(layout = "autofit") %>%
  height_all(height = .21, part = "all") %>% 
  hrule(rule = "exact") %>% 
  font(fontname = "Arial", part = "all") %>%
  fontsize(part = "all", size = 8) %>% 
  set_header_labels(site = "Site",
                    lat = "Latitude (°)",
                    long = "Longitude (°)",
                    mat = "Mean temperature (°C)",
                    map = "Mean precipitation (mm/yr)",
                    size_ha = "Plot size (ha)",
                    n_censuses = "N censuses",
                    census_years = "Census years") %>% 
  theme_booktabs() -> var

template <- system.file(package = "officer",
                        "doc_examples", "landscape.docx")
doc_1 <- read_docx(path = template)
doc_1 %>%
  body_add_flextable(var) -> my_doc

print(my_doc, target = paste0(path_tables, run, "/extended_data_table_1.docx")) %>% invisible()






# Extended Data Table 2 ---------------------------------------------------


# Variables: 
# number of observations, number of species, number of species in species groups (shrubs/trees). 

# generated based on nsp_global
# with additional info from 
#   _tree_1_metainfo.Rdata (number of recorded trees)
#   rAME_global (rare species group convergence)

# add number of recorded trees from _tree_1_metainfo.Rdata (full list of tree censuses)
# nsp_global$n_trees = NA
# for (site in sites) {
#   load(paste0(path_output, "data_1_metainfo/", site, "_tree_1_metainfo.Rdata"))
#   nsp_global$n_trees[nsp_global$site == site] = nrow(tree[[1]])     # take first census but all have the same number of entries
# }



# calculate various n
nsp_global %>% 
  mutate(nobs = ndead + nsurv) %>% 
  group_by(site) %>% 
  summarise(# n_trees = unique(n_trees),
            n_mortality_observations = sum(nobs),
            n_species = n(),
            perc_species_fitted_individually = round(100*sum(!rare)/n_species, 1),
            n_species_fitted_as_rare_trees = sum(rare_stature == "Rare_tree", na.rm = T),
            n_species_fitted_as_rare_shrubs = sum(rare_stature == "Rare_shrub", na.rm = T),
            perc_mortality = round(100*sum(ndead)/sum(nobs), 1)) -> ext_data_table_2

# global numbers (for text)
nsp_global %>% 
  mutate(nobs = ndead + nsurv) %>% 
  summarise(n_mortality_observations = sum(nobs),
            n_species = n(),
            perc_species_fitted_individually = round(100-100*sum(!rare)/n_species, 1),
            perc_mortality = round(100*sum(ndead)/sum(nobs), 1))


# replace n_rare with 0 if model not converged
rare_convergence = rAME_global[grepl("Rare", rAME_global$sp) & rAME_global$change == "equilibrium", c("site", "sp")]
for (site in sites) {
  groups = c("Rare_shrub", "Rare_tree") 
  not_converged = groups[!groups %in% rare_convergence$sp[rare_convergence$site == site]]
  if (any(not_converged == "Rare_shrub")) ext_data_table_2$n_species_fitted_as_rare_shrubs[ext_data_table_2$site == site] = 0
  if (any(not_converged == "Rare_tree")) ext_data_table_2$n_species_fitted_as_rare_trees[ext_data_table_2$site == site] = 0
}

# replace with full site name
ext_data_table_2$site = Site_table$site[match(ext_data_table_2$site, Site_table$ID)]
ext_data_table_2 = ext_data_table_2[order(ext_data_table_2$site), ]

# Formatting
ext_data_table_2 %>%
  flextable() %>%
  set_table_properties(layout = "autofit") %>%
  height_all(height = .21, part = "all") %>% 
  hrule(rule = "exact") %>% 
  font(fontname = "Arial", part = "all") %>%
  fontsize(part = "all", size = 8) %>% 
  set_header_labels(site = "Site",
                    # n_trees = "N trees on entire plot", 
                    n_mortality_observations = "N status observations",
                    n_species = "N species for mortality analyses",
                    perc_species_fitted_individually = "% species fitted individually",
                    n_species_fitted_as_rare_trees = "N species fitted as rare trees",
                    n_species_fitted_as_rare_shrubs = "N species fitted as rare shrubs",
                    perc_mortality = "% dead status observations") %>% 
  theme_booktabs() -> var

template <- system.file(package = "officer",
                        "doc_examples", "landscape.docx")
doc_2 <- read_docx(path = template)
doc_2 %>%
  body_add_flextable(var) -> my_doc

print(my_doc, target = paste0(path_tables, run, "/extended_data_table_2.docx")) %>% invisible()




