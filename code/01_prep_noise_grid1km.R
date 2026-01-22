library(tidyverse)
library(sf)
library(parallel)
library(purrr)

## DEFINE DIRECTORIES
## compute/kratos
if (Sys.getenv("USER")=="ckoenig") {
  data_raw <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/data_raw/"
  workdata <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/workdata/"
  tabs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/tables/"
  figs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/figures/"
}

######################
##### ROAD NOISE #####
######################

## generate list of raw roadnoise datasets
list_road <- list.files(path = paste0(data_raw,"noisedata/"),
                        pattern = "Aggroad_Lnight.shp",
                        recursive = TRUE,
                        full.names = TRUE)

## prepare road noise data in parallel
system.time(list_road_out <- mclapply(list_road, function(x) {
  
  df_sf <- st_read(x) %>%
    st_transform(32632) %>%
    st_make_valid() %>%
    filter(st_is_valid(.) & !st_is_empty(.)) %>%
    group_by(UnAggID, Agglomerat) %>% ## aggregate road noise polygons at city-level
    summarise(geometry = st_union(geometry)) %>%
    st_buffer(0) %>%
    st_simplify(dTolerance = 10, preserveTopology = TRUE) %>% ## simplify geometry (10m tolerance)
    st_make_valid() %>%
    filter(!st_is_empty(.))
  
  agg <- table(df_sf$Agglomerat) %>% names() %>% as.character()
  saveRDS(df_sf,file=paste0(workdata,"roadnoise_",agg,".rds")) ## save road noise polygons per city
  
  df_sf
},mc.cores = 20))

## combine polygons from different cities and save
road_simplified10 <- do.call(what=rbind, args=list_road_out)
saveRDS(road_simplified10,file=paste0(workdata,"road_simplified10.rds"))
road_simplified10 <- readRDS(paste0(workdata,"road_simplified10.rds"))

######################
##### RAIL NOISE #####
######################
## generate list of raw rail noise files
list_rail <- list.files(path = paste0(data_raw,"noisedata/envxrljlw_eisenbahn/"),
                        pattern = "Aggrail_Lnight.shp",
                        recursive = FALSE,
                        full.names = TRUE)

## prepare rail noise data
system.time(list_rail_out <- mclapply(list_rail, function(x) {
  
  df_sf <- st_read(x) %>%
    st_transform(32632) %>%
    st_make_valid() %>%
    st_simplify(dTolerance = 10, preserveTopology = FALSE) %>%
    st_make_valid() %>%
    filter(!st_is_empty(.))
  
  df_sf
  
},mc.cores = 10))

## combine and save
rail_simplified10 <- do.call(what=rbind, args=list_rail_out)
saveRDS(rail_simplified10,file=paste0(workdata,"rail_simplified10.rds"))
rail_simplified10 <- readRDS(paste0(workdata,"rail_simplified10.rds"))


##########################
##### AVIATION NOISE #####
##########################

## two sources of aviation noise data
## 1) Major airports
list_air1 <- list.files(path = paste0(data_raw,"noisedata/"),
                        pattern = "Mair_Lnight.shp",
                        recursive = TRUE,
                        full.names = TRUE)

system.time(list_air_out1 <- mclapply(list_air1, function(x) {
  
  df_sf <- st_read(x) %>%
    st_transform(32632) %>%
    st_make_valid() %>%
    st_simplify(dTolerance = 10, preserveTopology = FALSE) %>%
    st_make_valid() %>%
    filter(!st_is_empty(.))
  
  df_sf
  
},mc.cores = 20))

air1_simplified10 <- do.call(what=rbind, args=list_air_out1)

## 2) Air noise data Agglomerations (Airports of Dresden, Bremen, Dortmund)
list_air2 <- list.files(path = paste0(data_raw,"noisedata/"),
                       pattern = "Aggair_Lnight.shp",
                       recursive = TRUE,
                       full.names = TRUE)

list_air2 <- list_air2[c(1,3,4)] ## drop Info on Airport Leipzig/Halle (contained in first airport dataset above)

system.time(list_air_out2 <- mclapply(list_air2, function(x) {
  
  df_sf <- st_read(x) %>%
    st_transform(32632) %>%
    st_make_valid() %>%
    st_simplify(dTolerance = 10, preserveTopology = FALSE) %>%
    st_make_valid() %>%
    filter(!st_is_empty(.)) %>%
    select(Agglomerat,DB_Low,DB_High,geometry)
  
  df_sf
  
},mc.cores = 40))

air2_simplified10 <- do.call(what=rbind, args=list_air_out2)

## combine two aviation noise datasets
air_simplified10 <- air1_simplified10 %>%
  select(DB_Low,DB_High,geometry) %>%
  rbind(air2_simplified10 %>% select(DB_Low,DB_High,geometry))

saveRDS(air_simplified10,file=paste0(workdata,"air_simplified10.rds"))
air_simplified10 <- readRDS(paste0(workdata,"air_simplified10.rds"))


#####################################################
###### INTERSECTION WITH GRID RESIDENTIAL AREAS #####
#####################################################

## read spatial grid (ETRS89-LAEA-1km)
grid_core <- readRDS(paste0(workdata,"grid1km_fuacore.rds"))
st_crs(grid_core)

## read info about residential area per FUA/grid
residential <- readRDS(paste0(workdata,"residential_areas_core_valid_00.rds")) %>%
  mutate(area_residential = st_area(geometry))
st_crs(residential)

## intersect grid with residential areas
aux <- residential %>%
  inner_join(grid_core %>% select(id,Shape_Area) %>% st_drop_geometry(),by="id")

## intersect grid residential areas with aviation noise data
## calculate share of residential areas per grid affected by aviation noise
system.time(grid_air <- aux %>% st_intersection(air_simplified10, join=st_intersection,left=FALSE) %>%
  select(id,geometry,DB_Low,Shape_Area,area_residential) %>%
  mutate(area_airnoise=as.numeric(st_area(.))) %>%
  st_drop_geometry() %>%
  group_by(id,Shape_Area) %>%
  summarise(area_airnoise=sum(area_airnoise),
            area_residential=sum(area_residential)) %>%
  mutate(area_airnoise=(area_airnoise/area_residential)*100) %>%
  select(-Shape_Area))

## intersect grid residential areas with rail noise data
## calculate share of residential areas per grid affected by rail noise
system.time(grid_rail <- aux %>% st_intersection(rail_simplified10, join=st_intersection,left=FALSE) %>%
              select(id,geometry,DB_Low,Shape_Area,area_residential) %>%
              mutate(area_railnoise=as.numeric(st_area(.))) %>%
              st_drop_geometry() %>%
              group_by(id,Shape_Area) %>%
              summarise(area_railnoise=sum(area_railnoise),
                        area_residential=sum(area_residential)) %>%
              mutate(area_railnoise=(area_railnoise/area_residential)*100) %>%
              select(-Shape_Area))


## intersect grid residential areas with road noise data
## calculate share of residential areas per grid affected by road noise
aux <- aux %>%
  mutate(row_id=row_number())

## split into equal chunks; preparation of parallel computing to speed up 
system.time(aux_split_list <- aux %>% select(row_id) %>% st_drop_geometry() %>%
              split(factor(sort(rank(row.names(.))%%40))))

x <- aux_split_list[[1]] 

system.time(aux_road_list <- mclapply(aux_split_list, function(x){
  df <- inner_join(aux,x,by="row_id")
  df <- df %>% st_intersection(road_simplified10, join=st_intersection, left=FALSE) %>%
    select(id,geometry,Shape_Area,area_residential) %>%
    mutate(area_roadnoise=as.numeric(st_area(.))) %>%
    st_drop_geometry()
  
  df
},mc.cores = 21)) 

grid_road <- do.call(what = rbind,args = aux_road_list)

grid_road <- grid_road %>%
  group_by(id,Shape_Area) %>%
  summarise(area_roadnoise=sum(area_roadnoise),
            area_residential=sum(area_residential))%>%
  mutate(area_roadnoise=(area_roadnoise/area_residential)*100) %>%
  select(-Shape_Area)


## combine info on exposure to noise from different sources
grid_core_noise <- grid_road %>%
  select(-area_residential) %>%
  full_join(grid_rail %>% select(-area_residential),by="id") %>%
  full_join(grid_air %>% select(-area_residential),by="id") %>%
  inner_join(grid_core,by="id") %>% 
  st_as_sf() 


## read AGS5 shapefile (municipality boundaries)
## determine municipalities for which noise data is available; generate indicator
ags8 <- st_read(paste0(data_raw,"gdz_bkg/VG250_GEM.shp")) %>%
  st_transform(32632)

ags8_noise <- ags8 %>%
  select(GEN,AGS,geometry) %>%
  st_filter(grid_core_noise) %>%
  mutate(ags8_area=as.numeric(st_area(.))) %>%
  st_intersection(grid_core_noise) %>%
  mutate(share_area_intersection=as.numeric(st_area(.))/ags8_area) %>%
  ungroup() %>%
  group_by(GEN,AGS) %>%
  summarise(share_area_intersection=sum(share_area_intersection))%>%
  ungroup() %>%
  filter(share_area_intersection >.2)

noise_cities_name <- ags8_noise$GEN %>% table() %>% names() %>% as.list()
noise_cities_key <- ags8_noise$AGS %>% table() %>% names() %>% as.list()

ags8_noise <- ags8 %>% filter(GEN %in% noise_cities_name & AGS %in% noise_cities_key) %>% select(GEN,AGS)

saveRDS(ags8_noise,paste0(workdata,"ags8_Städte_noisedata.rds"))

## generate indicator of noise data availability
## grid cell considered to have noise data if more than 50% located within municipality with noise data available
grid_core_noise <- grid_core %>%
  left_join(grid_core_noise %>% st_drop_geometry() %>% select(id,area_roadnoise,area_railnoise,area_airnoise),by="id") %>%
  st_join(ags8_noise,left=FALSE,largest=TRUE) %>%
  mutate(across(.cols = c(area_roadnoise,area_railnoise,area_airnoise),
                .fns = ~as.numeric(.x) %>% replace_na(0)), ## handle NAs that are actually zeros
         noise_data_aval=1)

## identify grid cell IDs that have less than 50% overlap with noise data city
aux <- grid_core_noise %>% 
  st_intersection(ags8_noise %>% select(geometry)) %>%
  mutate(area_ags8_noise=as.numeric(st_area(.)),
         share_area_ags8_noise=area_ags8_noise/Shape_Area) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(share_area_ags8_noise=sum(share_area_ags8_noise)) %>%
  ungroup() %>%
  filter(share_area_ags8_noise>.5)

## drop grid cells that have less than 50% noise data coverage
grid_core_noise <- grid_core_noise %>% inner_join(aux)

saveRDS(grid_core_noise,paste0(workdata,"grid1km_noisedata.rds"))
