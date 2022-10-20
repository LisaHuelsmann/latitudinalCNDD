



# Calculate species growth rates





# Packages and data -------------------------------------------------------


library(dplyr)
library(tidyr)


load(paste0(path_output, "data_3b_growth/", site, "_tree_3b_growth.Rdata"))






# Mean species growth rate ------------------------------------------------


tree_growth %>% 
  # filter(interval > 0) %>% # remove 0 intervals, otherwise problem in linear model
  group_by(sp) %>% 
  summarise(n = n(), 
            incr_mean = mean(incr/interval, na.rm = T),
            incr_median = median(incr/interval, na.rm = T)) %>% 
  as.data.frame() -> growth


# hist(growth$incr_mean)
# hist(growth$incr_median)




# Size-dependent growth ---------------------------------------------------


tree_growth %>% 
  filter(interval > 0) %>% # remove 0 intervals, otherwise problem in linear model
  group_by(sp) %>% 
  mutate(n = n()) %>% 
  # filter(n>10) %>%
  summarise(incr_2cm = stats::coef(lm(I(incr/interval) ~ I(dbh-20)))[1],
            incr_5cm = stats::coef(lm(I(incr/interval) ~ I(dbh-50)))[1]) -> test

growth$incr_2cm = test$incr_2cm[match(growth$sp, test$sp)]
growth$incr_5cm = test$incr_5cm[match(growth$sp, test$sp)]

# plot(incr_5cm ~ incr_mean, growth); abline(0, 1)
# plot(incr_5cm ~ incr_median, growth); abline(0, 1)
# plot(incr_2cm ~ incr_median, growth); abline(0, 1)




# Save --------------------------------------------------------------------

save(growth, file = paste0(path_output, "meta_growth/", site, "_growth.Rdata"))



# Clean environment -------------------------------------------------------

rm(list = ls()[!grepl("site", ls()) & !grepl("path", ls())])

