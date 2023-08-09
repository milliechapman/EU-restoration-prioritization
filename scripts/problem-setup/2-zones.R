## Calculating the boundary limits for each zones in planning units
## Load packages
library(tidyverse)
library(fasterize)
library(terra)
library(raster)
library(sf)
library(fst)
rm(list = ls())

# Apply initial globiom corrections
# Default should be FALSE
apply_initialglobiom <- TRUE

############# (0 - SETUP) pull in PU data from script 1-PU.R############
## LULC for PU overall (proportion from 100m data)
rast_template <- raster("data/landcover/10km/Corine_2018_cropland.tif")
rast_template[!is.na(rast_template)] <- 1

# Make a planning unit (PU) template from one of the LU layers and give and ID
PU_template <- st_read("data/formatted-data/PU_template.shp") |>
  st_transform(crs = st_crs(rast_template))

nuts2 <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nutsIDnum = seq(1:260)) |>
  mutate(country = str_sub(NUTS_ID,1,2)) |>
  filter(country != "UK") |>
  st_transform(crs = st_crs(rast_template))

# Calculate forest per PU
forests_multi <- rast("data/IntensityLayersPastureCrop/MultiFunctional_05deg__BAU_ver3.tif")#rast("data/BIOCLIMA_G4M/Summarized/MultiFunctional_05deg.tif")
forests_setaside <- rast("data/IntensityLayersPastureCrop/SetAside_05deg__BAU_ver3.tif")#rast("data/BIOCLIMA_G4M/Summarized/SetAside_05deg.tif")
forests_production <- rast("data/IntensityLayersPastureCrop/ProductionForest_05deg__BAU_ver3.tif")#rast("data/BIOCLIMA_G4M/Summarized/ProductionForest_05deg.tif")
#forests_production <- resample(forests_production, forests_setaside)
forests_total <- sum(forests_multi, forests_setaside ,forests_production, na.rm = TRUE)
forests_multi_perc <- forests_multi/forests_total
forests_setaside_perc <- forests_setaside/forests_total
forests_production_perc <- forests_production/forests_total
# plot(forests_production_perc)

# Calculate crop intensity per PU
crop_high <- rast("data/Dou_CropIntensity/DouEtAl_HighIntensityCropland.tif")
crop_mid <- rast("data/Dou_CropIntensity/DouEtAl_MediumIntensityCropland.tif")
crop_low <- rast("data/Dou_CropIntensity/DouEtAl_LowIntensityCropland.tif")
crop_total <- sum(crop_high,crop_low,crop_mid,na.rm = TRUE)
# crop_high_perc <- crop_high/crop_total
# crop_mid_perc <- crop_mid/crop_total
# crop_low_perc <- crop_low/crop_total
# # Correct for NA values (set to 0)
# crop_low_perc[is.na(crop_low_perc)] <- 0
# crop_mid_perc[is.na(crop_mid_perc)] <- 0
# crop_high_perc[is.na(crop_high_perc)] <- 0

# Take total cropland area directly from leo
pu_in_EU <- read_csv("data/formatted-data/pu_in_EU.csv")
PU_lc <- read_csv("data/outputs/1-PU/PU_globiom_lc.csv") |> drop_na(PUID) |>
  left_join(PU_template)

# left_join(pu_in_EU, by = c("PUID" = "pu")) |>
non_perm_crop <- read_csv("data/outputs/1-PU/PU_globiom_permcrop.csv") |> drop_na(PUID) |>
  left_join(PU_template) |> st_as_sf() |>
  st_transform(st_crs(PU_template))
# plot(non_perm_crop$Cropland)

PU_raster <- raster("data/landcover/10km/Corine_2018_cropland.tif")
non_perm_crop <- terra::rasterize(non_perm_crop, rast(rast_template), field = "Cropland", fun = "mean")
# plot(non_perm_crop)
# plot(rast(PU_raster))

Leo_cropland_rast <- terra::rasterize(PU_lc |> sf::st_as_sf() , rast(rast_template), field = "Cropland", fun = "sum")
#Leo_cropland_rast <- Leo_cropland_rast*non_perm_crop
plot(Leo_cropland_rast)

## Dou crop split
# Leo_cropland_low <- Leo_cropland_rast * crop_low_perc
# Leo_cropland_mid <- Leo_cropland_rast * crop_mid_perc
# Leo_cropland_high <- Leo_cropland_rast * crop_high_perc

# Induce Piero correction here!
if(apply_initialglobiom){
  # (1) globiom initial area
  initial <-  read.csv("data/ManagementInitialConditions/EUCropland__2020primes_ref_2020REFERENCE_ver3.csv") |>
    dplyr::select(NUTS2, Intensity_reclass, area_1000ha) |>
    dplyr::mutate(Intensity_reclass = forcats::fct_collapse(Intensity_reclass,
                                                            "Cropland_low_glob" = "MinimalCropland",
                                                            "Cropland_med_glob" = "LightCropland",
                                                            "Cropland_high_glob" = "IntenseCropland",
                                                            "PermanentCropland" = "PermanentCropland"
    )) |>
    dplyr::group_by(NUTS2,Intensity_reclass) |>
    dplyr::summarise(area_1000ha = sum(area_1000ha,na.rm=T)) |> dplyr::ungroup() |>
    # Attribute permanent cropland to the others
    tidyr::pivot_wider(names_from = "Intensity_reclass", values_from = "area_1000ha",values_fill = 0) |>
    dplyr::group_by(NUTS2) |>
    dplyr::mutate(Cropland_high_glob = Cropland_high_glob + (PermanentCropland/3),
                  Cropland_med_glob = Cropland_med_glob + (PermanentCropland/3),
                  Cropland_low_glob = Cropland_low_glob + (PermanentCropland/3)) |>
    dplyr::ungroup() |> dplyr::select(-PermanentCropland) |>
    tidyr::pivot_longer(cols = Cropland_high_glob:Cropland_low_glob,names_to = "Intensity_reclass",values_to = "area_1000ha") |>
    # Convert to shares
    dplyr::group_by(NUTS2) |>
    dplyr::mutate(propshare = area_1000ha/sum(area_1000ha)) |>
    dplyr::ungroup() |>
    dplyr::select(-area_1000ha) |>
    # To wide
    tidyr::pivot_wider(id_cols = NUTS2,names_from = Intensity_reclass, values_from = propshare,values_fill = 0)

  # Join in with NUTS
  nuts2_crop <- nuts2 |> dplyr::left_join(initial, by = c("NUTS_ID" = "NUTS2"))

  # Rasterize Globiom NUTS shares
  nuts2_globiom_crop_low <- terra::rasterize(nuts2_crop, rast_template, field = "Cropland_low_glob", fun = "mean")
  nuts2_globiom_crop_mid <- terra::rasterize(nuts2_crop, rast_template, field = "Cropland_med_glob", fun = "mean")
  nuts2_globiom_crop_high <- terra::rasterize(nuts2_crop, rast_template, field = "Cropland_high_glob", fun = "mean")
  # get total globiom nuts
  nuts2_globiom_crop_total <- nuts2_globiom_crop_low + nuts2_globiom_crop_mid + nuts2_globiom_crop_high
  assertthat::assert_that(
    round(terra::global(rast(nuts2_globiom_crop_total),fun = "max",na.rm = TRUE)[,1],3) == 1
  )
  # Assign flat shares
  crop_high_perc <- nuts2_globiom_crop_high
  crop_mid_perc <- nuts2_globiom_crop_mid
  crop_low_perc <- nuts2_globiom_crop_low
  assertthat::assert_that(
    all(cellStats(stack(crop_high_perc, crop_low_perc, crop_mid_perc), "max") <= 1)
  )
  crop_total <- sum( crop_low_perc, crop_mid_perc, crop_high_perc)
} else {
  crop_total <- sum(crop_low, crop_mid ,crop_high, na.rm = TRUE)
  crop_high_perc <- crop_high/crop_total
  crop_mid_perc <- crop_mid/crop_total
  crop_low_perc <- crop_low/crop_total
}

# same w/ pasture intensity
pasture_high <- rast("data/IntensityLayersPastureCrop/PastureCLC_HighintensityEPIC_ver2.tif") / 10000 #rast("data/ManagementIntensity/10000/PastureCLC_HighintensityEPIC.tif")
pasture_low <- rast("data/IntensityLayersPastureCrop/PastureCLC_LowintensityEPIC_ver2.tif")  / 10000 #rast("data/ManagementIntensity/10000/PastureCLC_LowintensityEPIC.tif")
if(apply_initialglobiom){
  # The same for cropland
  initial <-  read.csv("data/ManagementInitialConditions/EUPasture__2020primes_ref_2020REFERENCE_ver3.csv") |>
    dplyr::select(NUTS2, Intensity, area_1000ha) |>
    dplyr::mutate(Intensity = forcats::fct_collapse(Intensity,
                                                    "pasture_low_globiom" = "LowIntensityPasture",
                                                    "pasture_high_globiom" = "HighIntensityPasture"
    )) |>
    dplyr::group_by(NUTS2, Intensity) |>
    dplyr::summarise(area_1000ha = sum(area_1000ha,na.rm=T)) |> dplyr::ungroup() |>
    # Convert to shares
    dplyr::group_by(NUTS2) |>
    dplyr::mutate(propshare = area_1000ha/sum(area_1000ha)) |>
    dplyr::ungroup() |>
    # To wide
    dplyr::select(-area_1000ha) |>
    tidyr::pivot_wider(id_cols = NUTS2,names_from = Intensity, values_from = propshare,values_fill = 0)
  initial$pasture_high_globiom[is.nan(initial$pasture_high_globiom)] <- 0
  initial$pasture_low_globiom[is.nan(initial$pasture_low_globiom)] <- 0

  # Join in with NUTS
  nuts2_past <- nuts2 |> dplyr::left_join(initial, by = c("NUTS_ID" = "NUTS2"))
  # Rasterize Globiom NUTS shares
  nuts2_globiom_past_low <- terra::rasterize(nuts2_past, rast_template, field = "pasture_low_globiom", fun = "mean")
  nuts2_globiom_past_high <- terra::rasterize(nuts2_past, rast_template, field = "pasture_high_globiom", fun = "mean")

  # Summarize Dou per NUTS region and rasterize again as above
  pasture_total <- sum(nuts2_globiom_past_low, nuts2_globiom_past_high, na.rm = TRUE)
  pasture_high_perc <- nuts2_globiom_past_high
  pasture_low_perc <- nuts2_globiom_past_low

  assertthat::assert_that(
    raster::cellStats(pasture_high_perc, "max") <= 1,
    raster::cellStats(pasture_low_perc, "max") <= 1,
    round( raster::cellStats(pasture_total, "max"),2 ) == 1
  )
} else {
  pasture_total <- sum(pasture_high, pasture_low, na.rm = TRUE)
  pasture_high_perc <- pasture_high/pasture_total
  pasture_low_perc <- pasture_low/pasture_total
}

nuts2_pasture <- nuts2 |>
  mutate(pasture_high = exactextractr::exact_extract(pasture_high, nuts2, 'sum'),
         pasture_low = exactextractr::exact_extract(pasture_low, nuts2, 'sum')) |>
  as_tibble() |>
  mutate(pasture_high = replace_na(pasture_high,0),
         pasture_low = replace_na(pasture_low,0),
         pasture_total = pasture_high + pasture_low) |>
  mutate(pasture_high_perc = pasture_high/pasture_total,
         pasture_low_perc = pasture_low/pasture_total) |>
  dplyr::select(NUTS_ID, pasture_high_perc,pasture_low_perc)

write_csv(nuts2_pasture, "data/ag_intensity/nuts2_pasture_intensity_perc.csv")

# extract by PU
PU_intensity_forests <- PU_template |>
  st_transform(crs = st_crs(forests_production_perc)) |>
  mutate(forests_multi_perc = exactextractr::exact_extract(forests_multi_perc, PU_template, 'mean'),
         forests_setaside_perc = exactextractr::exact_extract(forests_setaside_perc, PU_template, 'mean'),
         forests_production_perc = exactextractr::exact_extract(forests_production_perc, PU_template, 'mean'))

PU_intensity_forests <- as_tibble(PU_intensity_forests) |>
  dplyr::select(-geometry)

PU_intensity_pasture <- PU_template |>
  mutate(pasture_high_perc = exactextractr::exact_extract(pasture_high_perc, PU_template, 'mean'),
         pasture_low_perc = exactextractr::exact_extract(pasture_low_perc, PU_template, 'mean'))

PU_intensity_pasture <- as_tibble(PU_intensity_pasture) |>
  dplyr::select(-geometry)

 PU_intensity_crop <- PU_template |>
  mutate(crop_high_perc = exactextractr::exact_extract(crop_high_perc, PU_template, 'mean'),
         crop_low_perc = exactextractr::exact_extract(crop_low_perc, PU_template, 'mean'),
         crop_mid_perc = exactextractr::exact_extract(crop_mid_perc, PU_template, 'mean'))

PU_intensity_crop <- as_tibble(PU_intensity_crop) |>
  dplyr::select(-geometry)

#PU_natura_lc <- read_csv("data/outputs/1-PU/PU_natura_lc.csv")
# PU_natura_lc <- read_csv("data/outputs/1-PU/PU_natura_lc2.csv")
PU_natura_lc <- read_csv("data/outputs/1-PU/PU_natura_globiom_lc.csv") |>
  drop_na(PUID)
PU_lc <- read_csv("data/outputs/1-PU/PU_globiom_lc.csv") |> drop_na(PUID)

## LULC for natura in PU (proportion from 100m)
PU_natura_lc <- PU_natura_lc |>
  left_join(PU_intensity_forests) |>
  left_join(PU_intensity_pasture) |>
  left_join(PU_intensity_crop) |>
  mutate(forests_multi_perc = replace_na(forests_multi_perc,0),
         forests_production_perc = replace_na(forests_production_perc,0),
         forests_setaside_perc = replace_na(forests_setaside_perc,0)) |>
  mutate(all_forests = forests_multi_perc+forests_production_perc+ forests_setaside_perc,
         all_pasture = pasture_high_perc +  pasture_low_perc,
         all_crop = crop_low_perc + crop_mid_perc + crop_high_perc) |>
  #mutate(perc_primary = PrimaryForest/WoodlandForest_NP) |>
  # mutate(forests_setaside_perc = ifelse(forests_setaside_perc-perc_primary>0,
  #                                       forests_setaside_perc-perc_primary, 0),
  #        forests_multi_perc = forests_multi_perc + perc_primary/2,
  #        forests_production_perc = forests_production_perc + perc_primary/2) |>
  # mutate(WoodlandForest = WoodlandForest_NP) |>
  mutate(#WoodlandForest_primary = PrimaryForest,
    WoodlandForest_setaside = ifelse(all_forests > 0,
                                     WoodlandForest*forests_setaside_perc,
                                     WoodlandForest*0.3),
    WoodlandForest_multi = ifelse(all_forests > 0,
                                  WoodlandForest*forests_multi_perc,
                                  WoodlandForest*0.3),
    WoodlandForest_prod = ifelse(all_forests > 0,
                                 WoodlandForest*forests_production_perc,
                                 WoodlandForest*0.4),
    Cropland_low = ifelse(all_crop > 0,
                          Cropland*crop_low_perc,
                          0),
    Cropland_med = ifelse(all_crop > 0,
                          Cropland*crop_mid_perc,
                          0),
    Cropland_high = ifelse(all_crop > 0,
                           Cropland*crop_high_perc,
                           Cropland),
    Pasture_low = ifelse(all_pasture > 0,
                         Pasture*pasture_low_perc,
                         Pasture*0.5),
    Pasture_high = ifelse(all_pasture > 0,
                          Pasture*pasture_high_perc,
                          Pasture*0.5)) |>
  # mutate(WoodlandForest_primary = WoodlandForest_primary + WoodlandForest_setaside) |>
  dplyr::select(-c(Pasture, Cropland, WoodlandForest, # WoodlandForest_NP,
                   #WoodlandForest_setaside,
                   #perc_primary, #PrimaryForest,
                   all_forests, all_pasture, all_crop)) |>
  rename(WoodlandForest_primary = WoodlandForest_setaside,
         HeathlandShrub_natural = HeathlandShrub,
         Grassland_natural = Grassland,
         SparseVeg_natural = SparseVeg,
         Urban_urban = Urban,
         Wetlands_natural = Wetlands,
         RiversLakes_natural = RiversLakes,
         MarineTransitional_natural = MarineTransitional) |>
  dplyr::select(-c(forests_setaside_perc,forests_multi_perc, forests_production_perc,
                   pasture_low_perc, pasture_high_perc,
                   crop_high_perc, crop_mid_perc, crop_low_perc,
                   #area
  ))

# Write output
write_csv(PU_natura_lc, "data/outputs/2-zones/PU_natura_lc_intensity.csv")


## LULC for PU overall (proportion from 100m data)
PU_lc <- PU_lc  |>
  left_join(PU_intensity_forests) |>
  left_join(PU_intensity_pasture) |>
  left_join(PU_intensity_crop) |>
  mutate(forests_multi_perc = replace_na(forests_multi_perc,0),
         forests_production_perc = replace_na(forests_production_perc,0),
         forests_setaside_perc = replace_na(forests_setaside_perc,0)) |>
  mutate(all_forests = forests_multi_perc+forests_production_perc+ forests_setaside_perc,
         all_pasture = pasture_high_perc +  pasture_low_perc,
         all_crop = crop_low_perc + crop_mid_perc + crop_high_perc) |>
  mutate(WoodlandForest_primary = ifelse(all_forests > 0,
                                         WoodlandForest*forests_setaside_perc,
                                         WoodlandForest*0.3),
         WoodlandForest_multi = ifelse(all_forests > 0,
                                       WoodlandForest*forests_multi_perc,
                                       WoodlandForest*0.3),
         WoodlandForest_prod = ifelse(all_forests > 0,
                                      WoodlandForest*forests_production_perc,
                                      WoodlandForest*0.4),
         Cropland_low = ifelse(all_crop > 0,
                               Cropland*crop_low_perc,
                               0),
         Cropland_med = ifelse(all_crop > 0,
                               Cropland*crop_mid_perc,
                               0),
         Cropland_high = ifelse(all_crop > 0,
                                Cropland*crop_high_perc,
                                Cropland),
         Pasture_low = ifelse(all_pasture > 0,
                              Pasture*pasture_low_perc,
                              Pasture*0.5),
         Pasture_high = ifelse(all_pasture > 0,
                               Pasture*pasture_high_perc,
                               Pasture*0.5)) |>
  dplyr::select(-c(Pasture, Cropland, WoodlandForest, all_forests, all_pasture, all_crop)) |>
  rename(HeathlandShrub_natural = HeathlandShrub,
         Grassland_natural = Grassland,
         SparseVeg_natural = SparseVeg,
         Urban_urban = Urban,
         Wetlands_natural = Wetlands,
         RiversLakes_natural = RiversLakes,
         MarineTransitional_natural = MarineTransitional) |>
  dplyr::select(-c(forests_setaside_perc,forests_multi_perc, forests_production_perc,
                   pasture_low_perc, pasture_high_perc,
                   crop_high_perc, crop_mid_perc, crop_low_perc#,
                   #area
  ))


write_csv(PU_lc, "data/outputs/2-zones/PU_lc_intensity.csv")


## potential LULC for PU (binary 10km)
PU_potential_lc <- read_csv("data/outputs/1-PU/PU_potential_lc.csv") |>
  dplyr::select(-c(area, Status)) |>
  mutate(WoodlandForest_multi = ifelse(WoodlandForest>0,1,0),
         WoodlandForest_prod = ifelse(WoodlandForest>0,1,0),
         WoodlandForest_primary = ifelse(WoodlandForest>0,1,0)) |>
  mutate(Cropland_low = 1,
         Cropland_med = 1,
         Pasture_low =1) |>
  dplyr::select(-c(WoodlandForest)) |>
  rename(HeathlandShrub_natural = HeathlandShrub,
         Grassland_natural = Grassland,
         SparseVeg_natural = SparseVeg,
         Wetlands_natural = Wetlands,
         MarineTransitional_natural = MarineTransitional)
write_csv(PU_potential_lc, "data/outputs/2-zones/PU_potential_lc_intensity.csv")

## restoration logic to restrict transitions
restoration_logic <- read_csv("data/restoration-transitions.csv")

############## (1) Conservation #################################
## lower - LC area already in n2k

PU_natura_lc <- read_csv("data/outputs/2-zones/PU_natura_lc_intensity.csv") |>
   drop_na(PUID)
PU_lc <- read_csv("data/outputs/2-zones/PU_lc_intensity.csv")

PU_potential_lc <- read_csv("data/outputs/2-zones/PU_potential_lc_intensity.csv")

assertthat::assert_that(sum(PU_natura_lc$WoodlandForest_primary, na.rm = TRUE)>100)

conservation_lower <- PU_natura_lc |>
  dplyr::select(-c(Status)) |>
  drop_na(PUID) |>
  pivot_longer(-PUID) |>
  rename(zone = name) |>
  mutate(zone = paste0(zone, "_conserve")) |>
  rename(lower = value) |>
  mutate(lower = ifelse(lower <0,0, lower)) #|>
# mutate(lower = lower/10000)
# low is the area in n2k

## LC area in total planning unit is upper limit
conservation_upper <- PU_lc |>
  dplyr::select(-c(Status)) |>
  drop_na(PUID) |>
  pivot_longer(-PUID) |>
  rename(zone = name) |>
  mutate(zone = paste0(zone, "_conserve")) |>
  rename(upper = value)# |>
# mutate(upper = upper/10000)

# Upper larger than lower?
assertthat::assert_that(sum(conservation_lower$lower, na.rm = TRUE)<=sum(conservation_upper$upper, na.rm = TRUE))
ggplot(conservation_lower) +
  geom_histogram(aes(x = lower)) +
  facet_wrap(~zone)## join these for bounding

conservation_bounds <- conservation_lower |>
  left_join(conservation_upper) |>
  filter(zone != "WoodlandForest_multi_conserve",
         zone != "WoodlandForest_prod_conserve",
         zone != "Cropland_low_conserve",
         zone != "Cropland_med_conserve",
         zone != "Cropland_high_conserve",
         zone != "Pasture_low_conserve",
         zone != "Pasture_high_conserve",
         zone != "Urban_urban_conserve")


conservation_bounds |>
  group_by(PUID) |>
  summarise(upper = sum(upper),
            lower = sum(lower),
            diff = upper-lower) |>
  arrange(-upper)


conservation_bounds |>
  group_by(PUID) |>
  summarise(upper = sum(upper),
            lower = sum(lower),
            diff = upper-lower) |>
  ggplot() +
  geom_histogram(aes(x = diff))

hist(PU_potential_lc$SparseVeg_natural)
# get total area of each LC
restore_all <- PU_lc |>
  drop_na(PUID) |>
  dplyr::select(-c(Status)) |>
  pivot_longer(-PUID) |>
  left_join(restoration_logic) |>
  rename(LC = name, areaLC = value) |>
  pivot_longer(-c(PUID, LC, areaLC)) |>
  rename(LCnew = name, transition = value) |>
  left_join(PU_potential_lc, by = "PUID") |>
  pivot_longer(-c(PUID:transition)) |>
  rename(zone = name,
         potentialMAES = value) |>
  mutate(include = ifelse(LCnew == zone, 1,0)) |> # clunky way to do this, but only include is zone = the new land cover
  group_by(PUID, zone) |>
  # area in land cover * transition logic * potential binary
  summarise(upper = sum(areaLC*transition*potentialMAES*include, na.rm = TRUE)) |>
  ungroup() |>
  mutate(upper = ifelse(upper>1,1, upper),
         lower = 0,
         zone = paste0(zone, "_restore"))

#restore_all <- restore_all |> mutate(upper = upper*10000)
restore_all |>
  ggplot(aes(x = upper)) + geom_histogram() +
  facet_wrap(~zone) + theme_classic()


############### Production ####################################
## LC inside PA (to exclude from upper limit)
lc_n2k <- PU_natura_lc |>
  dplyr::select(-c(Status)) |>
  pivot_longer(-PUID) |>
  mutate(place = "n2k") |>
  filter(name != "WoodlandForest_multi",
         name != "WoodlandForest_prod",
         name != "Cropland_low",
         name != "Cropland_med",
         name != "Cropland_high",
         name != "Pasture_low",
         name != "Pasture_high",
         name != "Urban_urban") |>
  group_by(PUID) |>
  summarize(area_protected = sum(value, na.rm = TRUE))# cropland and pastureland
hist(lc_n2k$area_protected)
assertthat::assert_that(all(lc_n2k$area_protected <= 1))
##not part of upper limit
## LC full PUT
## get 90% of current crop area
##
lc_all <- PU_lc |>
  dplyr::select(-c(Status)) |>
  pivot_longer(-PUID) |>
  mutate(place = "all") |>
  group_by(PUID) |>
  summarise(landarea = sum(value, na.rm = TRUE)) |>
  mutate(landarea = ifelse(landarea>1,1, landarea))

PU_lc_crop <- PU_lc |> dplyr::select(PUID, Cropland_high)
# get total area of each LC
hist(lc_all$landarea)
cropland_bounds_high <- lc_all |>
  left_join(lc_n2k) |>
  mutate(lc_np = landarea-area_protected) |>
  left_join(PU_lc_crop) |>
  mutate(upper = Cropland_high*1.5) |>
  mutate(upper = ifelse(upper > landarea, landarea, upper)) |>
  mutate(upper = ifelse(upper <0,0,upper)) |>
  mutate(upper = lc_np) |>
  mutate(lower = 0) |> dplyr::select(-c(landarea,Cropland_high,  area_protected)) |>
  mutate(zone = "Cropland_high_production") |> dplyr::select(-lc_np)
hist(cropland_bounds_high$upper)

cropland_bounds_low <- lc_all |>
  left_join(lc_n2k) |> left_join(PU_lc) |>
  mutate(upper = Cropland_low*1.2) |>
  mutate(upper = ifelse(upper > landarea, landarea, upper)) |>
  mutate(upper = ifelse(upper <0,0,upper)) |>
  mutate(lower = 0) |> dplyr::select(PUID, upper, lower) |>
  mutate(zone = "Cropland_low_production")
hist(cropland_bounds_low$upper)

cropland_bounds_med <- lc_all |>
  left_join(lc_n2k) |> left_join(PU_lc) |>
  mutate(upper = Cropland_med*1.2) |>
  mutate(upper = ifelse(upper > landarea, landarea, upper)) |>
  mutate(upper = ifelse(upper <0,0,upper)) |>
  mutate(lower = 0) |> dplyr::select(PUID, upper, lower) |>
  mutate(zone = "Cropland_med_production")
hist(cropland_bounds_med$upper)

pasture_bounds_low <- lc_all |>
  left_join(lc_n2k) |> left_join(PU_lc) |>
  mutate(upper = Pasture_low*1.2) |> #(Pasture_low + Pasture_high)*1.1) |>
  mutate(upper = ifelse(upper > landarea, landarea, upper)) |>
  mutate(upper = ifelse(upper <0,0,upper)) |>
  mutate(lower = 0) |> dplyr::select(PUID, upper, lower) |>
  mutate(zone = "Pasture_low_production")
hist(pasture_bounds_low$upper)

pasture_bounds_high <-lc_all |>
  left_join(lc_n2k) |> left_join(PU_lc) |>
  mutate(upper = Pasture_high*1.2) |> #(Pasture_low + Pasture_high)*1.1) |>
  mutate(upper = ifelse(upper > landarea, landarea, upper)) |>
  mutate(upper = ifelse(upper <0,0,upper)) |>
  mutate(lower = 0) |> dplyr::select(PUID, upper, lower) |>
  mutate(zone = "Pasture_high_production")
hist(pasture_bounds_high$upper)

forestry_bounds <- lc_all |>
  left_join(lc_n2k) |> left_join(PU_lc) |>
  mutate(upper = WoodlandForest_prod*1.5) |> #(Pasture_low + Pasture_high)*1.1) |>
  mutate(upper = ifelse(upper > landarea, landarea, upper)) |>
  mutate(upper = ifelse(upper <0,0,upper)) |>
  mutate(lower = 0) |> dplyr::select(PUID, upper, lower) |>
  mutate(zone = "WoodlandForest_prod_production")

hist(forestry_bounds$upper)

forestry_multi_bounds <- conservation_lower |>
  left_join(conservation_upper) |>
  filter(zone == "WoodlandForest_multi_conserve") |>
  mutate(zone = ifelse(zone ==  "WoodlandForest_multi_conserve", "WoodlandForest_multi_production", NA)) |>
  mutate(upper = ifelse(upper>1,1,upper)) |>
  mutate(upper = upper)
hist(forestry_multi_bounds$upper)

forestry_prod_bounds <- forestry_bounds |>
  mutate(zone = "WoodlandForest_prod_production")

production_bounds <- rbind(forestry_multi_bounds,
                           forestry_prod_bounds,
                           pasture_bounds_high,
                           pasture_bounds_low,
                           cropland_bounds_med,
                           cropland_bounds_low,
                           cropland_bounds_high) |>
  mutate(upper = upper,
         lower = lower)

production_bounds |> mutate(diff = upper-lower) |>
  ggplot(aes(x = upper)) + geom_histogram() +
  facet_wrap(~zone) + theme_classic()
################### Urban locked in ######################

urban_lockin <- PU_lc |>
  dplyr::select(-c(Status)) |>
  pivot_longer(-PUID) |>
  rename(zone = name) |>
  filter(zone == "Urban_urban") |>
  mutate(zone = paste0(zone, "_lockin")) |>
  rename(upper = value)  |>
  mutate(upper = upper) |>
  mutate(lower = upper) |>
  mutate(lower = ifelse(lower<0,0,lower)) |>
  mutate(lower = ifelse(lower >1, 0.99,lower))

hist(urban_lockin$lower)

###################### Manual Bounded Constraints ##############
manual_bounded_constraints <- rbind(restore_all,
                                    production_bounds,
                                    conservation_bounds,
                                    urban_lockin
) |>
  mutate(upper = ifelse(upper <0,0, upper)) |>
  mutate(upper = ifelse(upper >1,1, upper))

manual_bounded_constraints |>
  ggplot(aes(x = upper)) + geom_histogram() +
  facet_wrap(~zone) + theme_classic()


write_csv(manual_bounded_constraints, "data/formatted-data/manual_bounded_constraints_production_globiom.csv")

################## production flexibility ###############

nuts2 <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  st_transform(crs = st_crs(rast_template))

pu_in_EU <- read_csv("data/formatted-data/pu_in_EU.csv")
#
pu <- read_fst("data/formatted-data/pu_data.fst") |>
  left_join(pu_in_EU) |>
  rename(id = EU_id) |>
  dplyr::select(-c(pu, nuts2id)) |>
  drop_na(id)
#
zones <- read_csv("data/formatted-data/zone_id.csv") |>
  mutate(name = paste0("z", id))
#
of <- "data/formatted-data/manual_bounded_constraints_production_globiom.csv"

manual_bounded_constraints <- read_csv(of) |>
  rename(pu = PUID) |>
  left_join(pu_in_EU) |>
  drop_na(EU_id) |>
  mutate(pu = EU_id) |>
  left_join(zones) |>
  #dplyr::select(-c(zone)) |>
  #rename(zone=name) |>
  dplyr::select(pu, zone, lower, upper, nuts2id)

write_csv(manual_bounded_constraints, "data/formatted-data/manual_bounded_constraints_production_globiom_flex.csv")

########## Get names of all zones ################3
# write out for future use!
all_zones <- unique(manual_bounded_constraints$zone)
zone_ID <- as.data.frame(all_zones) |>
  rename(zone = all_zones) |>
  mutate(id = seq(1:length(all_zones)))
write_csv(zone_ID, "data/formatted-data/zone_id.csv")
