# Libraries
library(tidyverse)
library(sf)
library(haven)
library(tmap)
library(parallel)

## DEFINE DIRECTORIES
## compute/kratos
if (Sys.getenv("USER")=="ckoenig") {
  data_raw <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/data_raw/"
  workdata <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/workdata/"
  tabs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/tables/"
  figs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/figures/"
}

ags_list_cities <- readRDS(paste0(workdata,"data_core_full_02.rds")) %>%
  select(AGS) %>%
  mutate(AGS=as.factor(AGS)) %>% 
  st_drop_geometry()

ags_list_cities <- ags_list_cities$AGS %>% 
  table() %>% names() %>% as.list()


## read AGS5 shapefile (Gemeindegrenzen)
## determine municipalities for which noise data is available; generate indicator
ags8 <- st_read(paste0(data_raw,"gdz_bkg/VG250_GEM.shp")) %>%
  st_transform(32632) %>%
  select(AGS,GEN) %>% 
  filter(AGS %in% ags_list_cities) %>%
  group_by(AGS,GEN) %>%
  summarise(geometry=st_union(geometry))


## read spatial grid
grid_core <- readRDS(paste0(workdata,"grid1km_fuacore.rds"))

system.time(grid_core <- grid_core %>% st_join(ags8, largest = TRUE) %>%
              filter(!is.na(AGS)))


## read BA data
ba_data <- read_dta(paste0(data_raw,"gesamt_final_kkz.dta"))

## AGS to character and add padding (leading zeros)
ba_data <- ba_data %>% mutate(
  ags=as.character(ags) %>% str_pad(8,pad = "0"))

## filter BA data to 2017
data <- ba_data %>% filter(jahr %in% c(2013:2017))



## merge spatial information to BA data
data <- data %>% inner_join(grid_core, by=c("Gitter1km"="id")) %>%
  st_as_sf() %>%
  st_set_crs(32632) %>%
  rename(id=Gitter1km) %>%
  select(id,jahr,bev_ges,akad,alo_ges,svb_ges,AGS,GEN) %>%
  filter(!is.na(akad) & !is.na(bev_ges) & !is.na(alo_ges) & !is.na(svb_ges))


data <- data %>%
  rename(municipality=GEN,
         grid_pop=bev_ges,
         year=jahr,
         sscemployed_total=svb_ges,
         unemployed_total=alo_ges,
         academic=akad) %>%
  filter(grid_pop!=0) %>% ## drop uninhabited grids
  filter(sscemployed_total>0|unemployed_total>0) ## filter on grids with at least 1 workingage resident

table(data$year,useNA = "always")

corr_acad_density <- data %>%
  filter(year==2017) %>%
  mutate(share_academic=academic*100/(sscemployed_total+unemployed_total)) %>% ## calculate grid-level share of academics (workingage)
  ungroup() %>%
  group_by(AGS,municipality) %>%
  summarise(corr_acad_density=cor(share_academic,grid_pop)) %>%
  select(AGS,municipality,corr_acad_density) %>%
  ungroup()

# corr_change_acad_density <- data %>%
#   mutate(share_academic=academic*100/(sscemployed_total+unemployed_total)) %>% ## calculate grid-level share of academics (workingage)
#   ungroup() %>%
#   group_by(id,AGS,municipality) %>%
#   arrange(year) %>%
#   mutate(share_academic_lag=lag(share_academic),
#          change_academic=share_academic-share_academic_lag) %>%
#   summarise(mean_change_academic=mean(change_academic,na.rm=TRUE),
#             mean_grid_pop=mean(grid_pop,na.rm=TRUE)) %>%
#   filter(!is.na(mean_change_academic)) %>%
#   ungroup() %>%
#   group_by(AGS,municipality) %>%
#   summarise(corr_change_acad_density=cor(mean_change_academic,mean_grid_pop)) %>%
#   select(AGS,municipality,corr_change_acad_density) %>%
#   ungroup()
# 
# 
# acad <- data %>%
#   filter(year==2017) %>%
#   ungroup() %>%
#   group_by(AGS,municipality) %>%
#   summarise(sscemployed_total=sum(sscemployed_total),
#             unemployed_total=sum(unemployed_total),
#             academic=sum(academic),
#             city_pop=sum(grid_pop)) %>%
#   mutate(share_academic=academic*100/(sscemployed_total+unemployed_total)) %>%
#   select(AGS,municipality,share_academic) %>%
#   ungroup()
# 
# change_acad <- data %>%
#   ungroup() %>%
#   group_by(year,AGS,municipality) %>%
#   summarise(sscemployed_total=sum(sscemployed_total),
#             unemployed_total=sum(unemployed_total),
#             academic=sum(academic),
#             city_pop=sum(grid_pop)) %>%
#   mutate(share_academic=academic*100/(sscemployed_total+unemployed_total)) %>%
#   ungroup() %>%
#   group_by(AGS,municipality) %>%
#   arrange(year) %>%
#   mutate(share_academic_lag=lag(share_academic),
#          change_academic=share_academic-share_academic_lag) %>%
#   summarise(mean_change_academic=mean(change_academic,na.rm=TRUE),
#             mean_grid_pop=mean(city_pop,na.rm=TRUE)) %>%
#   filter(!is.na(mean_change_academic)) %>%
#   select(AGS,municipality,mean_change_academic) %>%
#   ungroup()
  


# gentrification_cities <- acad %>%
#   inner_join(change_acad %>% st_drop_geometry(),by=c("AGS","municipality")) %>%
#   inner_join(corr_acad_density %>% st_drop_geometry(),by=c("AGS","municipality")) %>%
#   inner_join(corr_change_acad_density %>% st_drop_geometry(),by=c("AGS","municipality"))
# 
# class(gentrification_cities)


saveRDS(corr_acad_density,file = paste0(workdata,"gentrification_indices_cities.rds"))





