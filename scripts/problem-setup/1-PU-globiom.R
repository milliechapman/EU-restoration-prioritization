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
rm(list = ls())

globiom_lc <- stack("data/initial_lc_globiom/Protection_augmented_ETL2_shares_August_2023.tif")
perm_crop <- stack("data/initial_lc_globiom/Permanent_Cropland_split_share_total_crop_area_August_2023.tif")

PU_template_r <- raster("data/landcover/10km/Corine_2018_cropland.tif")
PU_template <- raster("data/landcover/10km/Corine_2018_cropland.tif") |>
  rasterToPolygons() |>
  st_as_sf() |>
  #st_transform(crs = st_crs(natura)) |>
  dplyr::mutate(PUID = seq(1:length(geometry))) |>
  dplyr::select(-Corine_2018_cropland)

PU_raster <- fasterize(PU_template, PU_template_r, field = "PUID")

globiom_lc <- projectRaster(globiom_lc,
                            crs = crs(PU_raster),
                            res = res(PU_raster))

globiom_lc_crop <- crop(globiom_lc, extent(PU_raster))
PU_raster <- crop(PU_raster, extent(globiom_lc_crop))
extent(globiom_lc_crop) <- extent(PU_raster)
stack <- stack(globiom_lc_crop, PU_raster)

nuts2 <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp")

perm_crop <- projectRaster(perm_crop,
                           crs = crs(PU_raster),
                           res = res(PU_raster))

perm_crop <- crop(perm_crop, extent(PU_raster))
extent(perm_crop) <- extent(PU_raster)
perm_crop <- stack(perm_crop, PU_raster)
perm_crop_df <- as.data.frame(perm_crop) |>
  rename(PUID = layer)
write_csv(perm_crop_df, "data/outputs/1-PU/PU_globiom_permcrop.csv")


z <- nuts2 |>
  mutate(value = exactextractr::exact_extract(stack, nuts2, fun = "sum"),
         value_IC = exact_extract(globiom_lc, nuts2, fun = "sum"))


PU_lc_globiom_df <- as.data.frame(stack) |>
  rename(PUID = layer)

glimpse(PU_lc_globiom_df)
hist(PU_lc_globiom_df$Grassland)
test <- PU_lc_globiom_df |>
  pivot_longer(-PUID) |>
  group_by(PUID) |> summarise(area = sum(value))

PU_natura_globiom_lc <- PU_lc_globiom_df |>
  mutate(WoodlandForest = Woodland.and.forest_protected,
         HeathlandShrub = Heathland.and.shrub_protected, #+ Transitional.woodland.shrub_protected,
         Grassland = Grassland_protected,
         Pasture = Pasture_protected,
         SparseVeg = Sparsely.vegetated.areas_protected,
         Cropland = Cropland_protected,
         Urban = Urban_protected,
         Wetlands = Wetlands_protected,
         RiversLakes = Rivers.and.lakes_protected,
         MarineTransitional = Marine.inlets.and.transitional.waters_protected + Marine,
         Status = "n2k") |>
  dplyr::select(PUID, WoodlandForest, HeathlandShrub, Grassland, Pasture, SparseVeg,
                Cropland, Urban, Wetlands, RiversLakes, MarineTransitional, Status)

write_csv(PU_natura_globiom_lc, "data/outputs/1-PU/PU_natura_globiom_lc.csv")

PU_globiom_lc <- PU_lc_globiom_df |>
  mutate(WoodlandForest = Woodland.and.forest_protected + Woodland.and.forest,
         HeathlandShrub = Heathland.and.shrub_protected + Heathland.and.shrub +
           Transitional.woodland.shrub, #+Transitional.woodland.shrub_protected,
         Grassland = Grassland_protected + Grassland,
         Pasture = Pasture_protected + Pasture,
         SparseVeg = Sparsely.vegetated.areas_protected + Sparsely.vegetated.areas,
         Cropland = Cropland_protected + Cropland,
         Urban = Urban_protected + Urban,
         Wetlands = Wetlands_protected +  Wetlands,
         RiversLakes = Rivers.and.lakes_protected + Rivers.and.lakes,
         MarineTransitional = Marine.inlets.and.transitional.waters_protected + Marine.inlets.and.transitional.waters + Marine,
         Status = "all_PU") |>
  dplyr::select(PUID, WoodlandForest, HeathlandShrub, Grassland, Pasture, SparseVeg,
                Cropland, Urban, Wetlands, RiversLakes, MarineTransitional, Status)

test <- PU_globiom_lc |>
  dplyr::select(-Status) |>
  drop_na() |>
  pivot_longer(-PUID) |>
  group_by(PUID) |>
  summarise(value = sum(value))
testn <- PU_natura_globiom_lc |>
  dplyr::select(-Status) |>
  drop_na() |>
  pivot_longer(-PUID) |>
  group_by(PUID) |>
  summarise(value = sum(value))

write_csv(PU_globiom_lc, "data/outputs/1-PU/PU_globiom_lc.csv")



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

