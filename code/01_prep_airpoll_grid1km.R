library(sf)
library(tidyverse)
library(parallel)
library(stats)

## Define directories (anonymized)
## compute/kratos
if (Sys.getenv("USER")=="ckoenig") {
  data_raw <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/data_raw/"
  workdata <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/workdata/"
  tabs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/tables/"
  figs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/figures/"
}

## read spatial grid (ETRS89-LAEA-1km)
grid_core <- readRDS(paste0(workdata,"grid1km_fuacore.rds"))

## read UBA air pollution data on 2x2km grid
uba_no2 <- readRDS(paste0(data_raw,"uba_2017/NO2_2017.rds")) %>%
  select(Index,MEAN,Shape) %>%
  rename(no2_mean=MEAN,geom=Shape)

uba_pm25 <- readRDS(paste0(data_raw,"uba_2017/PM25_2017.rds"))%>%
  select(Index,MEAN,Shape) %>%
  rename(pm25_mean=MEAN,geom=Shape)

uba_so2 <- readRDS(paste0(data_raw,"uba_2017/SO2_2017.rds"))%>%
  select(Index,MEAN,Shape) %>%
  rename(so2_mean=MEAN,geom=Shape)

## merge pollutant-specific data sets
## generate toxicity-weighted index of air pollution at the 2km x 2km grid
uba <- inner_join(uba_no2,uba_pm25 %>% st_drop_geometry(),by="Index") %>%
  inner_join(uba_so2 %>% st_drop_geometry(),by="Index") %>%
  st_transform(32632) %>%
  mutate(airpoll_w=pm25_mean*61.5+no2_mean*15.2+so2_mean*14.3)

## intersect air poll data with spatial grid
## assign pollution from 2km x 2km grids to 1km x 1km grid by spatial overlap
system.time(uba_grid <- grid_core %>% st_intersection(uba) %>%
  mutate(intersection=as.numeric(st_area(.))/Shape_Area) %>%
  mutate(across(
    .cols = c(no2_mean,pm25_mean,so2_mean,airpoll_w),
    .fns = ~.x*intersection))%>%
  st_drop_geometry() %>%
  group_by(id,Shape_Area) %>%
  summarise(across(.cols = c(no2_mean,pm25_mean,so2_mean,airpoll_w),
                   .fns = sum)))


saveRDS(uba_grid,paste0(workdata,"grid1km_uba_airpollution.rds"))
