


# Calculate species abundances


# Packages and data -------------------------------------------------------


library(dplyr)
library(tidyr)




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






# Calculate plot area -----------------------------------------------------

area = prod(plot$plotdimension[, 2]) / 1e4







# Abundances --------------------------------------------------------------


# average over all censuses 
# remember to add zeros!

dat %>% 
  filter(status == "A") %>% 
  group_by(sp, census) %>% 
  summarise(Nha = n()/area,
            BAha = sum(ba, na.rm = T)/area) %>% 
  pivot_wider(id_cols = sp, 
              names_from = census,
              values_from = c(Nha, BAha),
              values_fill = 0) %>% 
  pivot_longer(cols = -sp,
               names_to = c(".value", "census"),
               names_sep = "_") %>% 
  group_by(sp) %>% 
  summarise(Nha = mean(Nha),
            BAha = mean(BAha)) -> abundances


save(abundances, file = paste0(path_output, "meta_abundances/", site, "_abundances.Rdata"))






# Clean environment -------------------------------------------------------


rm(list = ls()[!grepl("site", ls()) & !grepl("path", ls())])
