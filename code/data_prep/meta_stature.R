

# Assign species to shrub or tree by their stature



# Packages and data -------------------------------------------------------


library(dplyr)



files = list.files(path_input, recursive = T, include.dirs = T)
files = files[grepl(site, files) & !grepl("data_raw", files) & grepl(".Rdata", files)]

for (i in files) {
  try(load(paste0(path_input, i)), silent = T)
}


if (any(ls() == "tree")) {
  tree = Map(cbind, tree, census = names(tree))
  dat = do.call("rbind", tree)
} else {
  stem = Map(cbind, stem, census = names(stem))
  dat = do.call("rbind", stem)
}






# Calculate max size of large trees ---------------------------------------

 
dat %>% 
  filter(!is.na(dbh) & status == "A") %>%
  arrange(census, sp, desc(dbh)) %>% 
  group_by(sp, census) %>% 
  slice(1:6) %>% 
  summarise(n=n(), 
            mean.dbh=mean(dbh)) %>%
  group_by(sp) %>% 
  summarise(dbh_largetree = max(mean.dbh)) %>% 
  mutate(stature = ifelse(dbh_largetree < 100, "shrub", "tree"))     -> stature




# Calculate max size of all trees based on quantile -----------------------


dat %>% 
  filter(!is.na(dbh) & status == "A") %>%
  group_by(sp, census) %>% 
  summarise(n=n(), 
            dbh_q90 = quantile(dbh, 0.9)) %>% 
  group_by(sp) %>% 
  summarise(dbh_q90 = max(dbh_q90)) -> temp
  

# combine
stature$dbh_q90 = temp$dbh_q90[match(stature$sp, temp$sp)]

plot(stature$dbh_largetree, stature$dbh_q90)





# Save result -------------------------------------------------------------

save(stature, file = paste0(path_output, "meta_stature/", site, "_stature.Rdata"))



# Clean environment -------------------------------------------------------

rm(list = ls()[!grepl("site", ls()) & !grepl("path", ls())])



