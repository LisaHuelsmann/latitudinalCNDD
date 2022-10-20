


# Calculate species abundances


# Packages and data -------------------------------------------------------


library(dplyr)
library(FuzzyQ)
library(tidyr)
library(vegan)
library (cluster)





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





# Fuzzy abundances --------------------------------------------------------


# assign observations to quadrats
# make a matrix with species counts for fuzzyq
# do this per census

fuzzy_abundances_censuses = data.frame()
for (i in unique(dat$census)) {
  
  dat %>% 
    filter(census == i & status == "A") %>% 
    mutate(qx = cut(gx, breaks = seq(0, plot$plotdimension[1, 2], by = 10)), 
           qy = cut(gy, breaks = seq(0, plot$plotdimension[2, 2], by = 10)),
           quadrat = as.numeric(as.factor(paste(qx, qy)))) %>% 
    group_by(quadrat, sp) %>% 
    summarise(n = n()) %>% 
    pivot_wider(id_cols = quadrat,
                names_from = sp,
                values_from = n,
                values_fill = 0) %>% 
    ungroup() %>% 
    dplyr::select(-quadrat) %>% 
    as.matrix() -> M
  
  fuzzy = fuzzyq(M, sorting = TRUE)
  # plot(fuzzy$A_O, log = "xy", col = fuzzy$spp$cluster+1)
  
  fuzzy_abundances_censuses = rbind(fuzzy_abundances_censuses, cbind(fuzzy$spp, sp = rownames(fuzzy$spp), census = i))
  
}


# average over all censuses 
# remember to add zeros!
fuzzy_abundances_censuses %>% 
  pivot_wider(id_cols = sp, 
              names_from = census,
              values_from = c(Common.I, cluster),
              values_fill = 0) %>% 
  pivot_longer(cols = -sp,
               names_to = c(".value", "census"),
               names_sep = "_") %>% 
  group_by(sp) %>% 
  summarise(commonness = mean(Common.I),
            common = median(cluster)) -> fuzzy_abundances

# hist(fuzzy_abundances$commonness)
# table(fuzzy_abundances$common)


save(fuzzy_abundances, file = paste0(path_output, "meta_abundances/", site, "_fuzzy_abundances.Rdata"))




# Species clusters --------------------------------------------------------




# assign species to clusters based on there joint occurfence in quadrats
# censuses are ignored so that temporal co-occurrence is considered


dat %>% 
  filter(status == "A" & dbh < 100 & !is.na(sp)) %>% 
  mutate(qx = cut(gx, breaks = seq(0, plot$plotdimension[1, 2], by = 10)), 
         qy = cut(gy, breaks = seq(0, plot$plotdimension[2, 2], by = 10)),
         quadrat = as.numeric(as.factor(paste(qx, qy)))) %>% 
  group_by(sp) %>% 
  mutate(nobs = n()) %>% 
  filter(nobs > 20) %>% 
  group_by(quadrat, sp) %>% 
  summarise(n = n()) %>% 
  pivot_wider(id_cols = quadrat,
              names_from = sp,
              values_from = n,
              values_fill = 0) %>% 
  ungroup() %>% 
  dplyr::select(-quadrat) %>% 
  as.data.frame() -> M

M.log <- t(log1p (M))
bc.dist <- vegdist(M.log, method = 'bray')

# chose k, if less than 10 species per cluster (5), take 2
k = ifelse(ncol(M)/10 < 5, 2, 5)

clust <- agnes(sqrt(bc.dist), method = 'ward') 
clust.hclust <- as.hclust(clust)
# plot(clust.hclust)
groups <- cutree (clust, k = k)
# table(groups)
# group.order <- groups[clust.hclust$order]
# group.in.cluster <- unique (group.order)
# plot (clust.hclust)
# rect.hclust (clust.hclust, border = group.in.cluster, k = round(ncol(M)/10, 0)) 
# legend ('topleft', legend = paste ('Cluster', 1:round(ncol(M)/10, 0)), pch = 22, col = 1:round(ncol(M)/10, 0), bty = 'n')


ks = kmeans(x = as.matrix(M.log), centers = k)
# table(ks$cluster, groups)

clusters = data.frame(sp = clust.hclust$labels
                      , clust_hier = groups
                      , clust_kmeans = ks$cluster
                      , census = i)


save(clusters, file = paste0(path_output, "meta_clusters/", site, "_clusters.Rdata"))



# # check how clusters plot spatially
# census = "tree1"
# sel = dat$census == census & dat$status == "A" & dat$dbh < 100 & !is.na(dat$sp)
# res = dat[sel, ]
# cols = clusters$clust_hier[match(res$sp, clusters$sp)]
# table(cols)
# 
# plot(res$gx, res$gy
#      , xlim = c(0, 200), ylim = c(0, 200)
#      , pch = ".", col = 1+cols, cex = 2)
# abline(h = seq(0, 500, 10), v = seq(0, 1000, 10))






# Clean environment -------------------------------------------------------


rm(list = ls()[!grepl("site", ls()) & !grepl("path", ls())])
