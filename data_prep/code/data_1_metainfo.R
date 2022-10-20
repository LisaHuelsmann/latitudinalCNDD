

# Script imports all information for a site and merges meta data to the stem and tree table


### includes the steps
# import stem and tree data plus all other files for a site
# add elevation data
# add soil data
# exports stem and tree data



library(raster)




# Import data -------------------------------------------------------------


# import data

files = list.files(path_input, recursive = T, include.dirs = T)
files = files[grepl(site, files) & !grepl("data_raw", files) & grepl(".Rdata", files)]

for (i in files) {
  try(load(paste0(path_input, i)), silent = T)
}




# Create tree data from stem data if tree is missing ----------------------

# e.g. for Pasoh (pas)

if (!any(ls() == "tree")) {
  
  tree = stem
  names(tree) = gsub("stem", "tree", names(tree))
  
  tree = lapply(tree, transform, nostems = 1)
  
  
  # select final columns
  tree = lapply(tree, dplyr::select, c("treeID", "sp", "gx", "gy", "dbh",  
                                       "homchange", "ba", "nostems", "status", "date"))
  
}





# Census selection when incomplete ----------------------------------------


datasets_info = readxl::read_xlsx(paste0(path_input, "global_metainfo/datasets_info.xlsx"))
problems = datasets_info$censuses.intervals[datasets_info$ID == site & datasets_info$mortality == "no"]
tree[paste0("tree", problems)] = NULL




# Add elevation data ------------------------------------------------------

if (any(grepl("elevation", ls()))) {
  # stem = lapply(stem, transform, 
  #               elevation = raster::extract(elevation, cbind(gx, gy)))
  tree = lapply(tree, transform, 
                elevation = raster::extract(elevation, cbind(gx, gy)))
} else {
  # stem = lapply(stem, transform, 
  #               elevation = NA)
  tree = lapply(tree, transform, 
                elevation = NA)
}




# Add soil data -----------------------------------------------------------



if (any(grepl("soil", ls()))) {
  
  for (i in names(soil)) {
    
    # soil[[i]] in raster::extract() didn't work, unclear why
    raster = soil[[i]]
    
    # add soil variable
    tree = lapply(tree, transform,
                  temp = raster::extract(raster, cbind(gx, gy)))
    
    # correct column name from temp to value of i
    tree = lapply(tree, function(x) {
      colnames(x)[colnames(x)=="temp"] <- paste0("soil_", i) ; x
    })
  }
}
  
# BayesianTools::correlationPlot(tree[[1]][, names(tree[[1]]) %in% names(soil)])








# Save file for site ------------------------------------------------------


save(tree, file = paste0(path_output, "data_1_metainfo/", site, "_tree_1_metainfo.Rdata"))






