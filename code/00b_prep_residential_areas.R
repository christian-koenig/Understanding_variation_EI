library(tidyverse)
library(sf)
library(parallel)
library(osmextract)

## compute/kratos
if (Sys.getenv("USER")=="ckoenig") {
  data_raw <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/data_raw/"
  workdata <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/workdata/"
  tabs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/tables/"
  figs <- "/wzb/samba/user/ckoenig/network_m/user/ckoenig/projects/github/Understanding_variation_EI/figures/"
}

## read spatial grid (ETRS89-LAEA-1km)
grid_core <- readRDS(paste0(workdata,"grid1km_fuacore.rds"))
st_crs(grid_core)

## read .pbf file once and transform into .gpkg
st_layers(paste0(data_raw,"germany-180101.osm.pbf"))
osm_full <- oe_read(file_path = paste0(data_raw,"germany-180101.osm.pbf"), layer = "multipolygons")

## read .gpkg
st_layers(paste0(data_raw,"germany-180101.gpkg")) ## "multipolygons" has 31915436 features

## split into equal chunks for parallel processing
total_rows <- 31915436
n_intervals <- 20

breaks <- round(seq(1, total_rows + 1, length.out = n_intervals + 1))

intervals <- vector("list", n_intervals)
for (i in seq_len(n_intervals)) {
  start <- breaks[i]
  end <- breaks[i + 1] - 1
  if (i == n_intervals) {
    end <- total_rows
  }
  intervals[[i]] <- c(start, end)
}

x <- intervals[[1]] ## list item for testing purposes

iter_read_list <- mclapply(intervals, function(x) {
  start <- x[[1]]
  end <- x[[2]]
  
  iter_query <- paste0("select * from multipolygons where rowid between ",eval(parse(text = start))," and ",eval(parse(text = end)),";")
  iter_read <- st_read(paste0(data_raw,"germany-180101.gpkg"), 
                       query = iter_query) %>%
    filter(landuse=="residential")
  
  ## repair invalid multipolygons
  iter_read_valid <- st_make_valid(iter_read)
  st_is_valid(iter_read_valid)
  
  iter_read_valid
}, mc.cores=20)

residential_areas <- do.call(what = rbind, args = iter_read_list) %>%
  st_transform(32632)

residential_areas_core <- st_filter(residential_areas,grid_core) %>%
  st_intersection(grid_core) %>%
  select(osm_id,id,landuse) %>%
  group_by(id,landuse) %>%
  summarise(geometry=st_union(geometry))

saveRDS(residential_areas %>% ungroup(), file = paste0(workdata,"residential_areas_valid_00.rds"))
saveRDS(residential_areas_core %>% ungroup(), file = paste0(workdata,"residential_areas_core_valid_00.rds"))
