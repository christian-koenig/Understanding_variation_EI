library(sf)
library(tidyverse)

## compute/kratos
if (Sys.getenv("USER")=="ckoenig") {
  data_raw <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/data_raw/"
  workdata <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/workdata/"
  tabs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/tables/"
  figs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/figures/"
}

## read spatial grid (ETRS89-LAEA-1km)
grid_1km <- st_read(paste0(data_raw,"DE_Grid_ETRS89-LAEA_1km.gpkg/geogitter/DE_Grid_ETRS89-LAEA_1km.gpkg")) %>%
  select(OBJECTID,id,p_staat,f_land,p_land,p_wasser,Shape_Area,geom) %>%
  st_transform(32632)

## load FUA core shape files
fua_core_shape <- st_read(paste0(data_raw,"FUA/DEU_core.shp")) %>%
  st_zm(drop = TRUE) %>%
  st_transform(32632)

## Core of Dresden FUA is missspecified (including surrounding districts)
## Fix: Replace with Dresden municipality border
dresden_core <- st_read(paste0(data_raw,"gdz_bkg/VG250_GEM.shp")) %>%
  filter(GEN=="Dresden") %>%
  select(GEN,geometry) %>%
  rename(fuaname=GEN) %>%
  st_transform(32632)

fua_core_shape <- fua_core_shape %>%
  filter(fuaname !="Dresden") %>%
  select(-fuacode) %>%
  rbind(dresden_core) 

## save adjusted FUA Core dataset
saveRDS(fua_core_shape,paste0(workdata,"fua_core_shape_adjusted.rds"))

## filter on grid cells inside FUAs and calculate share of area located within FUA core
grid_1km_fuacore <- st_filter(grid_1km,fua_core_shape)

aux <- st_intersection(grid_1km_fuacore,fua_core_shape) %>%
              filter(f_land!=0) %>%
  mutate(area_in_fuacore=as.numeric(st_area(.)),
         share_area_fuacore=area_in_fuacore/Shape_Area) %>%
  st_drop_geometry() %>%
  group_by(id,p_staat,f_land,p_land,p_wasser)%>%
  summarise(share_area_fuacore=sum(share_area_fuacore)) %>%
  ungroup() %>%
  select(id,share_area_fuacore)

grid_1km_fuacore <- grid_1km_fuacore %>%
  inner_join(aux,by="id") %>% 
  filter(share_area_fuacore>=.5)
  
saveRDS(grid_1km_fuacore,paste0(workdata,"grid1km_fuacore.rds"))
