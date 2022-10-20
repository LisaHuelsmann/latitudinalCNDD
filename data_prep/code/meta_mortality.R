



# Calculate species mortality rates





# Packages and data -------------------------------------------------------


library(dplyr)
library(tidyr)


load(paste0(path_output, "data_3a_mortality/", site, "_tree_3a_mortality_exp.Rdata"))






# Mean species growth rate ------------------------------------------------


# annualized via cloglog model with offset
mort_glm = function(mort_next, interval) {
  mod =  glm(mort_next ~ 1 + offset(log(interval)), family = binomial(link = "cloglog"))
  return(mod$family$linkinv(coef(mod)[1]))
}

tree_mort %>%  
  filter(!is.na(interval) & interval > 0) %>% 
  group_by(sp) %>%
  summarise(mort_rate = mort_glm(mort_next, interval)) %>% 
  as.data.frame() -> mort




# Size-dependent growth ---------------------------------------------------


mort_glm_site = function(mort_next, interval, dbh) {
  mod =  glm(mort_next ~ I(dbh-20) + offset(log(interval)), family = binomial(link = "cloglog"))
  return(mod$family$linkinv(coef(mod)[1]))
}

tree_mort %>% 
  filter(!is.na(interval) & interval > 0) %>% 
  group_by(sp) %>% 
  summarise(mort_rate = mort_glm_site(mort_next, interval, dbh)) %>% 
  as.data.frame() -> mort_size

mort$mort_rate_20 = mort_size$mort_rate[match(mort$sp, mort_size$sp)]

plot(mort$mort_rate, mort$mort_rate_20)




# Save --------------------------------------------------------------------

save(mort, file = paste0(path_output, "meta_mortality/", site, "_mortality.Rdata"))



# Clean environment -------------------------------------------------------

rm(list = ls()[!grepl("site", ls()) & !grepl("path", ls())])


