# Libraries
library(tidyverse)
library(sf)
library(tmap)
library(parallel)
library(spdep)
library(spatialreg)
library(sp)
library(modelsummary)
library(stats)
library(osmdata)
library(nngeo)
library(tidygeocoder)


## DEFINE DIRECTORIES
## compute/kratos
if (Sys.getenv("USER")=="ckoenig") {
  data_raw <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/data_raw/"
  workdata <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/workdata/"
  tabs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/tables/"
  figs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/figures/"
}

#################################
##### READ & PREP CENTROIDS #####
#################################

## read pre-processed data on inhabited as well as all grid cells in cities
data <- readRDS(paste0(workdata,"data_core_inhabited_02.rds")) %>% 
  st_transform(3035) %>%
  st_centroid()

sample_city_list <- data$AGS %>% table() %>% names() %>% as.list()

## read municipality data
mun_data <- st_read(paste0(data_raw,"gdz_bkg/vg250_ew/VG250_GEM.shp")) %>% 
  st_transform(3035) %>%
  filter(EWZ>=50000)


sample_cities <- mun_data %>%
  filter(AGS %in% sample_city_list)

##################################
##### READ & PREP CITY HALLS #####
##################################
townhalls <- read.csv(paste0(workdata,"cityhalls.csv"),sep = ";")

townhalls_sf <- tidygeocoder::geocode(townhalls,
                                      address = townhall_address,
                                      return_addresses = TRUE) %>% 
  st_as_sf(coords = c("long","lat")) %>%
  st_set_crs(4326) %>%
  st_transform(3035)

townhalls_sf <- townhalls_sf %>% 
  inner_join(mun_data %>% filter(AGS %in% sample_city_list) %>% select(AGS,GEN) %>% st_drop_geometry(),by=c("city"="GEN"))

write_rds(townhalls_sf,paste0(workdata,"townhalls_sf_3035.rds"))
# townhalls_sf <- read_rds(paste0(workdata,"townhalls_sf_3035.rds"))

## visual check
tmap_mode("view")

p <- tm_shape(townhalls_sf)+
  tm_dots(fill = "red",fill_alpha = 0.75, size = 1)+
  tm_shape(mun_data %>% filter(GEN %in% townhalls_sf$city))+
  tm_polygons(fill = "orange",fill_alpha = .3)

tmap_save(p,paste0(figs,"townhalls_manual.html"))


x <- sample_city_list[[1]]

dist_to_cityhall <- mclapply(sample_city_list, function(x) {
  
  data_subset <- data %>% filter(AGS==x)
  cityhall <- townhalls_sf %>% filter(AGS==x)
  dist <- st_distance(data_subset,cityhall)
  
  data_subset$dist_cityhall <- dist %>% as.numeric()
  
  data_subset <- data_subset %>%
    mutate(dist_cityhall=dist_cityhall/1000)
  
  data_subset
  
},mc.cores = 20)


data_aux <- do.call(what = rbind, args = dist_to_cityhall)

city_level_centrality_cityhall <- data_aux %>%
  select(id,AGS,municipality,share_sgb_total,share_workingage_foreign,share_academic,dist_cityhall,grid_pop) %>%
  group_by(AGS,municipality) %>%
  summarise(
    corr_sgb2_distcityhall = cor(.data$share_sgb_total,.data$dist_cityhall),
    corr_foreign_distcityhall = cor(.data$share_workingage_foreign,.data$dist_cityhall),
    corr_acad_distcityhall = cor(.data$share_academic,.data$dist_cityhall),
    corr_gridpop_distcityhall = cor(.data$grid_pop,.data$dist_cityhall))

saveRDS(city_level_centrality_cityhall,paste0(workdata,"city_level_centrality_robustness.rds"))


