library(sf)
library(tidyverse)
library(parallel)

## DEFINE DIRECTORIES
## compute/kratos
if (Sys.getenv("USER")=="ckoenig") {
  data_raw <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/data_raw/"
  workdata <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/workdata/"
  tabs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/tables/"
  figs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/figures/"
}

## read spatial grid (ETRS89-LAEA-1km)
grid_core <- readRDS(paste0(workdata,"grid1km_fuacore.rds"))

## this code creates a list of all files in the folder "urban atlas 2018" that are of type gpkg
urban_atlas_2018_list <- list.files(path = paste0(data_raw,"urban_atlas_2018/"), pattern = ".gpkg", recursive = TRUE) %>% as.list()

x <- urban_atlas_2018_list[[20]] 

system.time(UA2018_list <- mclapply(urban_atlas_2018_list, function(x){
  ## create filepath to the urban atlas data
  filepath <- paste0(data_raw,"urban_atlas_2018/",x)
  ## read in the urban atlas data
  df_UA <- st_read(filepath) 
  ## select only the columns that are needed for the analysis
  df_UA <- df_UA %>% select(-country,-prod_date,-perimeter,-area,-comment) 
  # transform the data 
  df_UA <- df_UA %>% st_transform(32632) 
  # create a new column that combines the code and the class
  df_UA <- df_UA %>% mutate(code_class=paste(code_2018,class_2018, sep = "_")) 
  # filter the data to only include urban green space, urban leisure space, and forests
  df_UA <- df_UA %>% filter(code_2018 %in% c("14100","14200","31000"))
  # create a new column that combines the classes into two categories (parks only, leisure space/forests)
  df_UA$class_agr <- case_when(df_UA$code_2018 %in% c("14100")~"urban_greens_2018",
                               df_UA$code_2018 %in% c("14200","31000")~"urban_leisure_forest_2018", 
                               TRUE~as.character(NA)) %>% as.factor()
  # select only the columns that are needed for the analysis
  df_UA <- df_UA %>% select(fua_name,geom,class_agr)
  
  # return df_UA to list
  df_UA
}, mc.cores=40))

system.time(UA_2018 <- do.call(what = rbind, args = UA2018_list))

#intersect the urban atlas shapefile with the grid
system.time(UA_grid_2018 <- UA_2018 %>% 
              st_intersection(grid_core, join = st_intersection,left=FALSE) %>%
              select(class_agr,id,geom) %>%
              mutate(class_agr_area=st_area(.)) %>%
              st_drop_geometry()) 
            
system.time(UA_grid_2018 <- UA_grid_2018 %>%
              group_by(id,class_agr)%>%
              #summarise the data by summing the area of each land use class
              summarise(class_agr_area=sum(as.numeric(class_agr_area))) %>%
              #pivot the data to get the land use classes as columns
              pivot_wider(id_cols = c(id),
                          names_from = class_agr,
                          values_from = class_agr_area))

## handle NAs, generate green space variables as share of grid surface area covered by green space
system.time(grid_joint <- left_join(grid_core,UA_grid_2018,by="id") %>%
  mutate(across(.cols=c(urban_leisure_forest_2018,urban_greens_2018),
           .fns = ~replace_na(.,0)),
    urban_green_leisure_forest_2018=urban_greens_2018+urban_leisure_forest_2018) %>%
  filter(f_land!=0) %>%
  mutate(area_fuacore=Shape_Area*share_area_fuacore,
         parks18=(urban_greens_2018/Shape_Area)*100,
         greens18=(urban_green_leisure_forest_2018/Shape_Area)*100,
         greenspace_data_aval=1,
         lack_greens18=100-greens18) %>%
  mutate_at(vars(parks18,greens18,lack_greens18),~case_when(
    .<0 ~0,
    .>100~100,
    TRUE~.)) %>%
  select(id,parks18,greens18,lack_greens18,greenspace_data_aval,geom) %>%
  st_drop_geometry())


saveRDS(grid_joint,paste0(workdata,"grid1km_urban_greenspace.rds"))


