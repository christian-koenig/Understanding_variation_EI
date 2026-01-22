# Libraries
library(tidyverse)
library(sf)
library(sp)
library(spdep)
library(haven)
library(parallel)
library(stats)

## DEFINE DIRECTORIES
## compute/kratos
if (Sys.getenv("USER")=="ckoenig") {
  data_raw <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/data_raw/"
  workdata <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/workdata/"
  tabs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/tables/"
  figs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/figures/"
}

## read data from Federal Employment Agency (BA)
ba_data <- read_dta(paste0(data_raw,"gesamt_final_kkz.dta"))

## Fix municipality ID: AGS to character and add padding (leading zeros)
ba_data <- ba_data %>% mutate(
  ags=as.character(ags) %>% str_pad(8,pad = "0"))

## read spatial grid
grid_core <- readRDS(paste0(workdata,"grid1km_fuacore.rds"))

## merge spatial information to BA data
data <- ba_data %>% inner_join(grid_core, by=c("Gitter1km"="id")) %>%
  select(-ags,-BL,-Gemeindename,-Bev_ges,-dissim_sgb,-dissim_4800,-dissim_6600,-dissim_akad,-sgb_quote_stadt)


## read prepared environmental quality data at the grid level
uba <- readRDS(paste0(workdata,"grid1km_uba_airpollution.rds")) %>% select(-Shape_Area)

noise <- readRDS(paste0(workdata,"grid1km_noisedata.rds"))%>% 
  select(-OBJECTID,-share_area_fuacore,-p_staat,-f_land,-p_land,-p_wasser,-Shape_Area,-geom) %>%
  mutate(noise_exposure=area_roadnoise+area_railnoise+area_airnoise) %>%
  st_drop_geometry()

greens <- readRDS(paste0(workdata,"grid1km_urban_greenspace.rds"))

## filter data to 2017
data_2017 <- data %>% filter(jahr==2017) %>%
  st_as_sf() %>%
  st_set_crs(32632) %>%
  rename(id=Gitter1km)

## merge environmental quality data at 1x1km grid level
data_2017 <- data_2017 %>% left_join(uba,by="id")%>%
  left_join(greens,by="id") %>%
  left_join(noise,by="id")

## use data avalability indicators to set variables NA if area not covered by respective dataset
data_2017 <- data_2017 %>%
  mutate(across(.cols = c(greenspace_data_aval,noise_data_aval),
                .fns = ~replace_na(.,0)),
         across(.cols = c(parks18,greens18,lack_greens18),
                .fns = ~case_when(
                  greenspace_data_aval==1 ~ as.character(.x),
                  TRUE ~ as.character(NA)) %>% as.numeric()),
         across(.cols = c(area_roadnoise,area_railnoise,area_airnoise,noise_exposure),
                .fns = ~case_when(
                  noise_data_aval==1 ~ as.character(.x),
                  TRUE ~ as.character(NA)) %>% as.numeric()))

## read city-level contextual variables and fix city ID issues
city_variables <- haven::read_dta(paste0(data_raw,"city-level_indices.dta")) %>%
  mutate(ags=as.character(ags),
         ags=case_when(
           ags_name=="Hamburg" ~ "02000000",
           ags_name=="Berlin" ~ "11000000",
           TRUE ~ ags),
         ags=str_pad(ags,width=8,side="left",pad="0")) %>%
  select(ags,dissim_sgb,dissim_ausl,bev_log,Bev_ges,dichte_arm,dichte_ausl) %>%
  rename("mun_pop"="Bev_ges")
  
## restrict to those cities/grids with noise data available
data_2017 <- data_2017 %>% filter(noise_data_aval==1)
data_2017 <- left_join(data_2017,city_variables,by=c("AGS"="ags")) ## merge contextual variables by city ID

## standardize EQ vars (population-weighted) and generate index of multiple environmental burdens
weighted_sd <- function(x, w) {
  wm <- weighted.mean(x, w)
  variance <- sum(w * (x - wm)^2) / sum(w)
  sqrt(variance)
}

z_standardize_weighted <- function(x, w) {
  wm <- weighted.mean(x, w)
  wsd <- weighted_sd(x, w)
  (x - wm) / wsd
}

data_2017 <- data_2017 %>%
  rename(grid_pop=bev_ges) %>%
  mutate(
    across(
      .cols = c(airpoll_w,noise_exposure,greens18,lack_greens18),
      .fns = ~z_standardize_weighted(.x,w=grid_pop),
      .names = "{col}_std"),
    meb=airpoll_w_std+noise_exposure_std+lack_greens18_std,
    meb_std=z_standardize_weighted(meb,w=grid_pop))%>%
  select(id,Bev_km2,jahr,alo_ges,alo_sgb2,alo_ausl,alo_aka,sgb_pers,svb_ges,svb_aka,svb_ausl,svb_4800,grid_pop,unter_15,zw_15_65,ueber_65,unter65,akad,
         Shape_Area,share_area_fuacore,geom,airpoll_w,airpoll_w_std,greens18,greens18_std,lack_greens18,lack_greens18_std,greenspace_data_aval,GEN,AGS,
         noise_data_aval,share_area_ags8_noise,noise_exposure,noise_exposure_std,mun_pop,dissim_sgb,dissim_ausl,bev_log,dichte_arm,dichte_ausl,meb,meb_std) %>%
  mutate(rownum=row_number()) ## generate row number variable


## Check for grid cells with 0 neighbors (contiguity criterion)
data_core_full <- data_2017

## identify grids with no neighboring grids
listw <- poly2nb(data_core_full, queen = TRUE, snap = 10)
list_no_neighbors <- which(card(listw)==0) %>% as.list()

## drop 1 grid cell with no neighboring grids
data_core_full <- data_core_full %>%
  filter(!(rownum %in% list_no_neighbors))

## generate municipality level variables on EQ
municipality_eq <- data_core_full %>%
  select(airpoll_w,airpoll_w_std,greens18,greens18_std,lack_greens18,lack_greens18_std,noise_exposure,noise_exposure_std,meb,meb_std,grid_pop,GEN) %>%
  st_drop_geometry() %>%
  ungroup() %>%
  group_by(GEN) %>%
  summarise(airpoll_mun=mean(airpoll_w),
            greens_mun=mean(greens18),
            lack_greens_mun=mean(lack_greens18),
            noise_mun=mean(noise_exposure),
            meb_mun=mean(meb),
            meb_mun_w=weighted.mean(meb,grid_pop))

## save; re-load for city-level analyses later
saveRDS(municipality_eq,paste0(workdata,"municipality_eq_02.rds"))

## final data preparation steps; neighborhood composition variables
data <- data_core_full %>%
  select(id,AGS,GEN,grid_pop,unter_15,zw_15_65,ueber_65,unter65,sgb_pers,svb_ges,svb_ausl,svb_4800,alo_ges,alo_ausl,akad,
         airpoll_w,airpoll_w_std,greenspace_data_aval,greens18,greens18_std,lack_greens18,lack_greens18_std,meb,meb_std,
         noise_data_aval,noise_exposure,noise_exposure_std,share_area_fuacore,Shape_Area,
         dissim_sgb,dissim_ausl,bev_log,mun_pop,dichte_arm,dichte_ausl)%>%
  rename(municipality=GEN,
         age_0_14=unter_15,
         age_15_65=zw_15_65,
         age_over_65=ueber_65,
         age_0_65=unter65,
         sgb_total=sgb_pers,
         sscemployed_total=svb_ges,
         sscemployed_foreign=svb_ausl,
         unemployed_total=alo_ges,
         unemployed_foreign=alo_ausl,
         academic=akad) %>%
  filter(grid_pop!=0) %>%
  filter(sgb_total>0|sscemployed_total>0|unemployed_total>0)%>%
  filter(age_0_65>0) %>%
  mutate(mun_pop=mun_pop/100000,
         share_sgb_total=sgb_total*100/age_0_65,
         share_workingage_foreign=(unemployed_foreign+sscemployed_foreign)*100/(unemployed_total+sscemployed_total),
         share_workingage_4800=(svb_4800*100)/(unemployed_total+sscemployed_total),
         share_academic=academic*100/(sscemployed_total+unemployed_total)) %>%
  filter(complete.cases(share_sgb_total,share_workingage_foreign)) %>%
  mutate(across(.cols=c(age_0_14,age_0_65,age_15_65,age_over_65),
                .fns=~.x/grid_pop,
                .names="share_{.col}"),
         across(.cols=c(share_sgb_total,share_workingage_foreign,share_workingage_4800),
                .fns= ~ ifelse(is.nan(.x),NA,.x)),
         across(.cols=c(share_sgb_total,share_workingage_foreign,share_workingage_4800),
                .fns= ~ case_when(
                  .x>100 ~ "1",
                  between(.x,0,100) ~ "0",
                  TRUE ~ as.character(NA))%>%as.numeric(),
                .names = "flag_{.col}"),
         across(.cols=c(share_sgb_total,share_workingage_foreign,share_workingage_4800),
                .fns= ~ case_when(
                  .x>100 ~ "50",
                  TRUE ~ as.character(.x))%>%as.numeric()))


## add EQ vars at municipality level
data <- inner_join(data,municipality_eq,by=c("municipality"="GEN"))

## z-standardization of municipality-level variables
data <- data %>%
  mutate(across(
    .cols=c(bev_log,airpoll_mun,greens_mun,lack_greens_mun,noise_mun,meb_mun,meb_mun_w,dissim_sgb,dissim_ausl,dichte_arm,dichte_ausl),
    .fns=  ~ .x/sqrt(Hmisc::wtd.var(.x,weights = data$grid_pop)),
    .names="{.col}_std"))


saveRDS(data,paste0(workdata,"data_core_inhabited_02.rds"))
saveRDS(data_core_full,paste0(workdata,"data_core_full_02.rds"))
