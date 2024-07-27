## this script calculates the total area of 100m CORINE land cover
## in each 10km planning unit across the EU. It also calculates this
## for land cover within N2K protected areas in each planning unity
## Load packages
library(sf)
library(tidyverse)
library(terra)
library(raster)
library(fasterize)
library(exactextractr)

# Template for aligning rasters
LULC_template <- raster("data/landcover/CLC_MAES_100m/Corine_NonAgg_cl1-WoodlandForest_laea.tif")

# Natura 2000 data
natura <- st_read("data/Natura2000_end2021_Shapefile/Natura2000_end2021_epsg3035.shp") |>
  st_cast("POLYGON") |>
  dplyr::select(SITETYPE, geometry) |>
  mutate(n2k = 1)

# Create n2k raster for calculating LULC in PA network
if (!file.exists("data/Natura2000_end2021_Shapefile/natura_100m.tif")) {
  natura <- fasterize(natura, LULC_template, field = "n2k")
  natura <- classify(rast(natura), cbind(NA, 0))
  writeRaster(natura, "data/Natura2000_end2021_Shapefile/natura_100m.tif", overwrite = TRUE)
}
natura <- raster("data/Natura2000_end2021_Shapefile/natura_100m.tif")

plot(natura)
# load in all LULC layers
filelist_temp <- list.files("data/landcover/CLC_MAES_100m/", pattern = "*_laea.tif")
LULC <- rast(paste0("data/landcover/CLC_MAES_100m/", filelist_temp))
LULC <- stack(LULC)

woodland <- rast("data/landcover/CLC_MAES_100m/Corine_NonAgg_cl1-WoodlandForest_laea.tif")
WoodlandForest_n2k <- woodland*natura
plot(WoodlandForest_n2k)
#natura[is.na(natura[])] <- 0

# Create n2k rasters (only run if files doesnt exist- not a quick process)
if (!file.exists("data/landcover/100m_n2k/WoodlandForest_n2k.tif")) {
  #PrimaryForest_n2k <- primary_forest_raster*natura
  #writeRaster(PrimaryForest_n2k, "data/landcover/100m_n2k/PrimaryForest_n2k.tif")

  WoodlandForest_n2k <- rast("data/landcover/CLC_MAES_100m/Corine_NonAgg_cl1-WoodlandForest_laea.tif")
  WoodlandForest_n2k <- WoodlandForest_n2k*natura
  writeRaster(WoodlandForest_n2k, "data/landcover/100m_n2k/WoodlandForest_n2k.tif")

  HeathlandShrub_n2k <- rast("data/landcover/CLC_MAES_100m/Corine_NonAgg_cl2-HeathlandShrub_laea.tif")
  HeathlandShrub_n2k <- HeathlandShrub_n2k*natura
  writeRaster(HeathlandShrub_n2k, "data/landcover/100m_n2k/HeathlandShrub_n2k.tif")

  Grassland_n2k <- rast("data/landcover/CLC_MAES_100m/Corine_nonAgg_cl3-NaturalGrassland_laea.tif")
  Grassland_n2k <- HeathlandShrub_n2k*natura
  writeRaster(Grassland_n2k, "data/landcover/100m_n2k/Grassland_n2k.tif")

  Pasture_n2k <- LULC[["Corine_nonAgg_cl11.Pasture_laea"]]*natura
  writeRaster(Grassland_n2k, "data/landcover/100m_n2k/Pasture_n2k.tif")

  SparseVeg_n2k <- LULC[["Corine_NonAgg_cl4.SparseVeg_laea"]]*natura
  writeRaster(SparseVeg_n2k, "data/landcover/100m_n2k/SparseVeg_n2k.tif")

  Cropland_n2k <- LULC[["Corine_NonAgg_cl5.Cropland_laea"]]*natura
  writeRaster(Cropland_n2k, "data/landcover/100m_n2k/Cropland_n2k.tif")

  Urban_n2k <- LULC[["Corine_NonAgg_cl6.Urban_laea"]]*natura
  writeRaster(Urban_n2k, "data/landcover/100m_n2k/Urban_n2k.tif")

  Wetlands_n2k <- LULC[["Corine_NonAgg_cl7.Wetlands_laea"]]*natura
  writeRaster(Wetlands_n2k, "data/landcover/100m_n2k/Wetlands_n2k.tif")

  RiversLakes_n2k <- LULC[["Corine_NonAgg_cl8.RiversLakes_laea"]]*natura
  writeRaster(RiversLakes_n2k, "data/landcover/100m_n2k/RiversLakes_n2k.tif")

  MarineTransitional_n2k <- LULC[["Corine_NonAgg_cl9.MarineTransitional_laea"]]*natura
  writeRaster(MarineTransitional_n2k, "data/landcover/100m_n2k/MarineTransitional_n2k.tif")
}

# Make a planning unit (PU) template from one of the LU layers and give and ID
PU_template <- raster("data/landcover/10km/Corine_2018_cropland.tif") |>
  rasterToPolygons() |>
  st_as_sf() |>
  st_transform(crs = st_crs(natura)) |>
  dplyr::mutate(PUID = seq(1:length(geometry))) %>%
  dplyr::select(-Corine_2018_cropland)

# Make stack of 100m LC data sets (PU is the aggregate of these)
filelist_temp <- list.files("data/landcover/100m_n2k/", pattern = "*_n2k.tif")
LULC_n2k <- rast(paste0("data/landcover/100m_n2k/", filelist_temp))
LULC_n2k <- stack(LULC_n2k)

# woodlandforest primary

# woodlandforest not primary
n2k_NP<- overlay(LULC_n2k[["WoodlandForest_n2k"]],
                 primary_forest_raster,
                 fun=function(r1, r2){return(r1-r2)})

writeRaster(n2k_NP, "data/landcover/100m_n2k/WoodlandForest_n2k_NP.tif")
# extract 100m LC for each PU/natura
# takes a while so only if file does not exist
# note: exactextractr is more efficient w/ raster than raster stack...
#if (!file.exists("data/outputs/1-PU/PU_natura_lc.csv")) {
PU_natura_lc <- PU_template |>
  mutate(
    PrimaryForest = exactextractr::exact_extract(primary_forest_raster, PU_template, 'sum'),
    WoodlandForest_NP = exactextractr::exact_extract(n2k_NP, PU_template, 'sum'),
    WoodlandForest = exactextractr::exact_extract(LULC_n2k[["WoodlandForest_n2k"]], PU_template, 'sum'),
    HeathlandShrub = exactextractr::exact_extract(LULC_n2k[["HeathlandShrub_n2k"]], PU_template, 'sum'),
    Grassland = exactextractr::exact_extract(LULC_n2k[["Grassland_n2k"]], PU_template, 'sum'),
    Pasture = exactextractr::exact_extract(LULC_n2k[["Pasture_n2k"]], PU_template, 'sum'),
    SparseVeg = exactextractr::exact_extract(LULC_n2k[["SparseVeg_n2k"]], PU_template, 'sum'),
    Cropland = exactextractr::exact_extract(LULC_n2k[["Cropland_n2k"]], PU_template, 'sum'),
    Urban = exactextractr::exact_extract(LULC_n2k[["Urban_n2k"]], PU_template, 'sum'),
    Wetlands = exactextractr::exact_extract(LULC_n2k[["Wetlands_n2k"]], PU_template, 'sum'),
    RiversLakes = exactextractr::exact_extract(LULC_n2k[["RiversLakes_n2k"]], PU_template, 'sum'),
    MarineTransitional = exactextractr::exact_extract(LULC_n2k[["MarineTransitional_n2k"]], PU_template, 'sum')) |>
  mutate(area = st_area(geometry),
         Status = "N2K")

as_tibble(PU_natura_lc) |>
  dplyr::select(-geometry) |>
  write_csv("data/outputs/1-PU/PU_natura_lc2.csv")

z <- read_csv("data/outputs/1-PU/PU_lc.csv") |> pivot_longer(-c(PUID, area,Status))

PU_natura_lc_long <- PU_natura_lc |>
  pivot_longer(-c(PUID, geometry, area, Status))
raster_template <- raster("data/SpeciesData/CurrentSDMs/Amphibians/EnsembleNormThreshold__Alytes_dickhilleni.tif")

lc_n2k_raster_forest <- fasterize(sf = PU_natura_lc_long, raster = raster_template, field = "value", by = "name")

lc_n2k_raster_forest <- rast(lc_n2k_raster_forest)
writeRaster(lc_n2k_raster_forest, "data/lc_PA_10km.tif",  overwrite = TRUE)
rast_test <- rast("data/lc_PA_10km.tif")
names(rast_test)
plot(lc_n2k_raster_forest[[2]])

#}
#PU_natura_lc <- read_csv("data/outputs/1-PU/PU_natura_lc.csv")

# extract all LC for each PU
#if (!file.exists("data/outputs/1-PU/PU_lc.csv")) {
PU_lc <- PU_template |>
  mutate(WoodlandForest = exactextractr::exact_extract(LULC[["Corine_NonAgg_cl1.WoodlandForest_laea"]], PU_template, 'sum'),
         HeathlandShrub = exactextractr::exact_extract(LULC[["Corine_NonAgg_cl2.HeathlandShrub_laea"]], PU_template, 'sum'),
         Grassland = exactextractr::exact_extract(LULC[["Corine_nonAgg_cl3.NaturalGrassland_laea"]], PU_template, 'sum'),
         Pasture = exactextractr::exact_extract(LULC[["Corine_nonAgg_cl11.Pasture_laea"]], PU_template, 'sum'),
         SparseVeg = exactextractr::exact_extract(LULC[["Corine_NonAgg_cl4.SparseVeg_laea"]], PU_template, 'sum'),
         Cropland = exactextractr::exact_extract(LULC[["Corine_NonAgg_cl5.Cropland_laea"]], PU_template, 'sum'),
         Urban = exactextractr::exact_extract(LULC[["Corine_NonAgg_cl6.Urban_laea"]], PU_template, 'sum'),
         Wetlands = exactextractr::exact_extract(LULC[["Corine_NonAgg_cl7.Wetlands_laea"]], PU_template, 'sum'),
         RiversLakes = exactextractr::exact_extract(LULC[["Corine_NonAgg_cl8.RiversLakes_laea"]], PU_template, 'sum'),
         MarineTransitional = exactextractr::exact_extract(LULC[["Corine_NonAgg_cl9.MarineTransitional_laea"]], PU_template, 'sum')) |>
  mutate(area = st_area(geometry),
         Status = "all_PU")
as_tibble(PU_lc) |>
  dplyr::select(-geometry) |>
  write_csv("data/outputs/1-PU/PU_lc.csv")
#}
#PU_lc <- read_csv("data/outputs/1-PU/PU_lc.csv")

# check that total cover is mostly 10km (or 0 outside of land area)
PU_lc <- PU_lc |>
  mutate(total_cover = WoodlandForest + HeathlandShrub + Grassland + Pasture + SparseVeg +
           Cropland + Urban + Wetlands + RiversLakes + MarineTransitional)

# just area in natura - check that < 10km
PU_natura_lc <- PU_natura_lc |>
  mutate(total_cover = WoodlandForest + HeathlandShrub + Grassland + Pasture+ SparseVeg +
           Cropland + Urban + Wetlands + RiversLakes + MarineTransitional)

# Binary potential by PU
filelist_temp <- list.files("data/PotentialMAES/", pattern = "PotentialMAES__median_*")[1:6]
PotentialMAES <- rast(paste0("data/PotentialMAES/", filelist_temp))
PotentialMAES <- stack(PotentialMAES)
PotentialMAES <- projectRaster(PotentialMAES,
                               crs = crs(PU_template))


PU_potential_lc <- PU_template |>
  mutate(WoodlandForest = exactextractr::exact_extract(PotentialMAES[["PotentialMAES__median__Woodland.and.forest"]], PU_template, 'sum'),
         HeathlandShrub = exactextractr::exact_extract(PotentialMAES[["PotentialMAES__median__Heathland.and.shrub"]], PU_template, 'sum'),
         Grassland = exactextractr::exact_extract(PotentialMAES[["PotentialMAES__median__Grassland"]], PU_template, 'sum'),
         SparseVeg = exactextractr::exact_extract(PotentialMAES[["PotentialMAES__median__Sparsely.vegetated.areas"]], PU_template, 'sum'),
         Wetlands = exactextractr::exact_extract(PotentialMAES[["PotentialMAES__median__Wetlands"]], PU_template, 'sum'),
         MarineTransitional = exactextractr::exact_extract(PotentialMAES[["PotentialMAES__median__Marine.inlets.and.transitional.waters"]], PU_template, 'sum')) |>
  mutate(area = st_area(geometry),
         Status = "Potential")

# write csv
as_tibble(PU_potential_lc) |>
  dplyr::select(-geometry) |>
  write_csv("data/outputs/1-PU/PU_potential_lc.csv")

############## Planning unit costs per zone #############
# read in zone data
zone_id <- read_csv("data/formatted-data/zone_id.csv")
all_zones <- paste0("z", zone_id$id)

# Create a value of cost = 1 (area) for each PU
PU_template <- raster("data/landcover/10km/Corine_2018_cropland.tif") |>
  rasterToPolygons() |>
  st_as_sf() |>
  #st_transform(crs = st_crs(natura)) |>
  dplyr::mutate(PUID = seq(1:length(geometry))) |>
  dplyr::select(-Corine_2018_cropland)
st_write(PU_template, "data/formatted-data/PU_template.shp", overwrite = TRUE)

raster <- raster("data/landcover/10km/Corine_2018_cropland.tif")
# make it into a raster
PU_raster <- fasterize(PU_template, raster, field = "PUID")
writeRaster(PU_raster, "data/formatted-data/PU_raster.tif")

# PU in EU
# Get NUTS2 jurisdiction data per planning unit and make a binary stack
nuts2 <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nutsIDnum = seq(1:260)) |>
  st_transform(crs = crs(PU_template)) |>
  mutate(country = str_sub(NUTS_ID,1,2)) |>
  filter(country != "UK")
nuts2_raster <- fasterize(nuts2, PU_raster, field = "nutsIDnum", fun = "first")
names(nuts2_raster) <- "nuts2id"
writeRaster(nuts2_raster, "data/formatted-data/nuts2_raster.tif", overwrite = TRUE)
plot(nuts2_raster)
# filter out PU outside of the EU
pu_in_EU <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(nuts2_raster))) |>
  drop_na(nuts2id) |>
  rename(pu = layer) |>
  dplyr::select(-id)

# add a new id for problem formatting
pu_in_EU |>
  mutate(EU_id = seq(1:nrow(pu_in_EU))) |>
  dplyr::select(-cost) |>
  write_csv("data/formatted-data/pu_in_EU.csv")
