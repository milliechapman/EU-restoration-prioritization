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
pasture_high <- rast("data/IntensityLayersPastureCrop/PastureCLC_HighintensityEPIC_ver2.tif")#rast("data/ManagementIntensity/10000/PastureCLC_HighintensityEPIC.tif")
rast_template <- raster("data/landcover/10km/Corine_2018_cropland.tif")
rast_template[!is.na(rast_template)] <- 1

# Make a planning unit (PU) template from one of the LU layers and give and ID
PU_template <- st_read("data/formatted-data/PU_template.shp") |>
  st_transform(crs = st_crs(pasture_high))

nuts2 <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nutsIDnum = seq(1:260)) |>
  mutate(country = str_sub(NUTS_ID,1,2)) |>
  filter(country != "UK") |>
  st_transform(crs = st_crs(pasture_high))

# Calculate forest per PU
forests_multi <- rast("data/IntensityLayersPastureCrop/MultiFunctional_05deg__BAU_ver3.tif")#rast("data/BIOCLIMA_G4M/Summarized/MultiFunctional_05deg.tif")
forests_setaside <- rast("data/IntensityLayersPastureCrop/SetAside_05deg__BAU_ver3.tif")#rast("data/BIOCLIMA_G4M/Summarized/SetAside_05deg.tif")
forests_production <- rast("data/IntensityLayersPastureCrop/ProductionForest_05deg__BAU_ver3.tif")#rast("data/BIOCLIMA_G4M/Summarized/ProductionForest_05deg.tif")
#forests_production <- resample(forests_production, forests_setaside)
forests_total <- sum(forests_multi, forests_setaside ,forests_production, na.rm = TRUE)
forests_multi_perc <- forests_multi/forests_total
forests_setaside_perc <- forests_setaside/forests_total
forests_production_perc <- forests_production/forests_total


# Calculate crop intensity per PU
crop_high <- rast("data/Dou_CropIntensity/DouEtAl_HighIntensityCropland.tif")
crop_mid <- rast("data/Dou_CropIntensity/DouEtAl_MediumIntensityCropland.tif")
crop_low <- rast("data/Dou_CropIntensity/DouEtAl_LowIntensityCropland.tif")
crop_total <- sum(crop_low, crop_mid ,crop_high, na.rm = TRUE)
crop_high_perc <- crop_high/crop_total
crop_mid_perc <- crop_mid/crop_total
crop_low_perc <- crop_low/crop_total

# same w/ pasture intensity
pasture_high <- rast("data/IntensityLayersPastureCrop/PastureCLC_HighintensityEPIC_ver2.tif")#rast("data/ManagementIntensity/10000/PastureCLC_HighintensityEPIC.tif")
pasture_low <- rast("data/IntensityLayersPastureCrop/PastureCLC_LowintensityEPIC_ver2.tif") #rast("data/ManagementIntensity/10000/PastureCLC_LowintensityEPIC.tif")

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

pasture_total <- sum(pasture_high + pasture_low, na.rm = TRUE)
pasture_high_perc <- pasture_high/pasture_total
pasture_low_perc <- pasture_low/pasture_total

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

# Apply initial globiom correction
if(apply_initialglobiom){
  # Extract nuts2 id per PU
  nuts2_ras <- fasterize(nuts2, rast_template, field = "nutsIDnum", fun = "last")
  newpu <- PU_template
  newpu$nutsIDnum <- exactextractr::exact_extract(nuts2_ras, PU_template, "max")
  newpu <- newpu |> sf::st_drop_geometry() |> dplyr::filter(!is.na(nutsIDnum))
  # Get the actual NUTS2 id
  newpu <- left_join(newpu, nuts2 |> dplyr::select(NUTS2,nutsIDnum) |> sf::st_drop_geometry(),  "nutsIDnum")
  newpu <- newpu |> dplyr::select(-nutsIDnum)
  # Join in the nuts id with the PU_natura_lc
  PU_natura_lc <- left_join(PU_natura_lc, newpu, by = "PUID")
  rm(newpu)

  # Load in the initial targets for reference (assuming they are both relatively comparable)
  initial <- read.csv("data/ManagementInitialConditions/EUPasture__2020primes_ref_2020REFERENCE_ver3.csv") |>
    dplyr::select(NUTS2, Intensity, area_1000ha) |>
    dplyr::mutate(Intensity = forcats::fct_collapse(Intensity,
                                                    "pasture_high_perc_glob" = "HighIntensityPasture",
                                                    "pasture_low_perc_glob" = "LowIntensityPasture"
    )) |>
    # Convert to shares
    dplyr::group_by(NUTS2) |>
    dplyr::mutate(propshare = area_1000ha/sum(area_1000ha)) |>
    dplyr::ungroup() |>
    dplyr::select(-area_1000ha) |>
    # To wide
    tidyr::pivot_wider(id_cols = NUTS2,names_from = Intensity, values_from = propshare,values_fill = 0)
  initial$pasture_high_perc_glob[is.nan(initial$pasture_high_perc_glob)] <- 0
  initial$pasture_low_perc_glob[is.nan(initial$pasture_low_perc_glob)] <- 0

  # The same for cropland
  initial2 <-  read.csv("data/ManagementInitialConditions/EUCropland__2020primes_ref_2020REFERENCE_ver3.csv") |>
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

  # Combine the 2 initial shares
  initial_comb <- dplyr::left_join(initial,initial2)

  # Now join in new shares
  new <- dplyr::left_join(PU_natura_lc, initial_comb, by = "NUTS2")

  # --- #
  # Apply correction
  # Summarize per nuts2
  new <- dplyr::left_join(new,
                          new |> dplyr::group_by(NUTS2) |> dplyr::summarise(Pasture_low_total = sum(Pasture_low, na.rm = TRUE),
                                                                            Pasture_high_total = sum(Pasture_high, na.rm = TRUE),
                                                                            Cropland_low_total = sum(Cropland_low, na.rm = TRUE),
                                                                            Cropland_med_total = sum(Cropland_med, na.rm = TRUE),
                                                                            Cropland_high_total = sum(Cropland_high, na.rm = TRUE)
                                                                            )
  )
  o <- new$pasture_low_perc_glob * (new$Pasture_low/new$Pasture_low_total)
  # If sum of corrected shares should be equal to the globiom initial share (or at least not larger)
  assertthat::assert_that( sum(o,na.rm = T) <= sum(new$pasture_low_perc_glob,na.rm=TRUE) )
  new$Pasture_low <- o

  o <- new$pasture_high_perc_glob * (new$Pasture_high/new$Pasture_high_total)
  assertthat::assert_that( sum(o,na.rm = T) <= sum(new$pasture_high_perc_glob,na.rm=TRUE) )
  new$Pasture_high <- o
  # Cropland
  o <- new$Cropland_low_glob * (new$Cropland_low/new$Cropland_low_total)
  assertthat::assert_that( sum(o,na.rm = T) <= sum(new$Cropland_low_glob,na.rm=TRUE) )
  new$Cropland_low <- o
  o <- new$Cropland_med_glob * (new$Cropland_med/new$Cropland_med_total)
  assertthat::assert_that( sum(o,na.rm = T) <= sum(new$Cropland_med_glob,na.rm=TRUE) )
  new$Cropland_med <- o
  o <- new$Cropland_high_glob * (new$Cropland_high/new$Cropland_high_total)
  assertthat::assert_that( sum(o,na.rm = T) <= sum(new$Cropland_high_glob,na.rm=TRUE) )
  new$Cropland_high <- o

  new <- new |> dplyr::select(-Pasture_low_total,-Pasture_high_total,-pasture_low_perc_glob,-pasture_high_perc_glob,
                              -Cropland_low_glob,-Cropland_med_glob,-Cropland_high_glob,-Cropland_low_total,-Cropland_med_total,
                              -Cropland_high_total)
  assertthat::assert_that(!any(stringr::str_detect(names(new),"glob|total")))
  if("NUTS2" %in% names(new)) new <- new |> dplyr::select(-NUTS2)
  # --- #
  # Write output
  write_csv(new, "data/outputs/2-zones/PU_natura_lc_intensity_initialGLOBIOM.csv")
} else {
  # Write output
  write_csv(PU_natura_lc, "data/outputs/2-zones/PU_natura_lc_intensity.csv")
}

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

# Apply initial globiom correction
if(apply_initialglobiom){
  # Extract nuts2 id per PU
  nuts2_ras <- fasterize(nuts2, rast_template, field = "nutsIDnum", fun = "last")
  newpu <- PU_template
  newpu$nutsIDnum <- exactextractr::exact_extract(nuts2_ras, PU_template, "max")
  newpu <- newpu |> sf::st_drop_geometry() |> dplyr::filter(!is.na(nutsIDnum))
  # Get the actual NUTS2 id
  newpu <- left_join(newpu, nuts2 |> dplyr::select(NUTS2,nutsIDnum) |> sf::st_drop_geometry(),  "nutsIDnum")
  newpu <- newpu |> dplyr::select(-nutsIDnum)
  # Join in the nuts id with the PU_natura_lc
  PU_lc <- left_join(PU_lc, newpu, by = "PUID")
  rm(newpu)

  # Load in the initial targets for reference (assuming they are both relatively comparable)
  initial <- read.csv("data/ManagementInitialConditions/EUPasture__2020primes_ref_2020REFERENCE_ver3.csv") |>
    dplyr::select(NUTS2, Intensity, area_1000ha) |>
    dplyr::mutate(Intensity = forcats::fct_collapse(Intensity,
                                                    "pasture_high_perc_glob" = "HighIntensityPasture",
                                                    "pasture_low_perc_glob" = "LowIntensityPasture"
    )) |>
    # Convert to shares
    dplyr::group_by(NUTS2) |>
    dplyr::mutate(propshare = area_1000ha/sum(area_1000ha)) |>
    dplyr::ungroup() |>
    dplyr::select(-area_1000ha) |>
    # To wide
    tidyr::pivot_wider(id_cols = NUTS2,names_from = Intensity, values_from = propshare,values_fill = 0)
  initial$pasture_high_perc_glob[is.nan(initial$pasture_high_perc_glob)] <- 0
  initial$pasture_low_perc_glob[is.nan(initial$pasture_low_perc_glob)] <- 0

  # The same for cropland
  initial2 <-  read.csv("data/ManagementInitialConditions/EUCropland__2020primes_ref_2020REFERENCE_ver3.csv") |>
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

  # Combine the 2 initial shares
  initial_comb <- dplyr::left_join(initial,initial2)

  # Now join in new shares
  new <- dplyr::left_join(PU_lc, initial_comb, by = "NUTS2")

  # --- #
  # Apply correction
  # Summarize per nuts2
  new <- dplyr::left_join(new,
                          new |> dplyr::group_by(NUTS2) |> dplyr::summarise(Pasture_low_total = sum(Pasture_low, na.rm = TRUE),
                                                                            Pasture_high_total = sum(Pasture_high, na.rm = TRUE),
                                                                            Cropland_low_total = sum(Cropland_low, na.rm = TRUE),
                                                                            Cropland_med_total = sum(Cropland_med, na.rm = TRUE),
                                                                            Cropland_high_total = sum(Cropland_high, na.rm = TRUE)
                          )
  )
  o <- new$pasture_low_perc_glob * (new$Pasture_low/new$Pasture_low_total)
  # If sum of corrected shares should be equal to the globiom initial share (or at least not larger)
  assertthat::assert_that( sum(o,na.rm = T) <= sum(new$pasture_low_perc_glob,na.rm=TRUE) )
  new$Pasture_low <- o

  o <- new$pasture_high_perc_glob * (new$Pasture_high/new$Pasture_high_total)
  assertthat::assert_that( sum(o,na.rm = T) <= sum(new$pasture_high_perc_glob,na.rm=TRUE) )
  new$Pasture_high <- o
  # Cropland
  o <- new$Cropland_low_glob * (new$Cropland_low/new$Cropland_low_total)
  assertthat::assert_that( sum(o,na.rm = T) <= sum(new$Cropland_low_glob,na.rm=TRUE) )
  new$Cropland_low <- o
  o <- new$Cropland_med_glob * (new$Cropland_med/new$Cropland_med_total)
  assertthat::assert_that( sum(o,na.rm = T) <= sum(new$Cropland_med_glob,na.rm=TRUE) )
  new$Cropland_med <- o
  o <- new$Cropland_high_glob * (new$Cropland_high/new$Cropland_high_total)
  assertthat::assert_that( sum(o,na.rm = T) <= sum(new$Cropland_high_glob,na.rm=TRUE) )
  new$Cropland_high <- o

  new <- new |> dplyr::select(-Pasture_low_total,-Pasture_high_total,-pasture_low_perc_glob,-pasture_high_perc_glob,
                              -Cropland_low_glob,-Cropland_med_glob,-Cropland_high_glob,-Cropland_low_total,-Cropland_med_total,
                              -Cropland_high_total)
  assertthat::assert_that(!any(stringr::str_detect(names(new),"glob|total")))
  # --- #
  if("NUTS2" %in% names(new)) new <- new |> dplyr::select(-NUTS2)
  # Write output
  write_csv(new, "data/outputs/2-zones/PU_lc_intensity_initialGLOBIOM.csv")
} else {
  # Write output
  write_csv(PU_lc, "data/outputs/2-zones/PU_lc_intensity.csv")
}

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
if(apply_initialglobiom){
  PU_natura_lc <- read_csv("data/outputs/2-zones/PU_natura_lc_intensity_initialGLOBIOM.csv") |>
    drop_na(PUID)
  PU_lc <- read_csv("data/outputs/2-zones/PU_lc_intensity_initialGLOBIOM.csv")
} else {
  PU_natura_lc <- read_csv("data/outputs/2-zones/PU_natura_lc_intensity.csv") |>
    drop_na(PUID)
  PU_lc <- read_csv("data/outputs/2-zones/PU_lc_intensity.csv")
}
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

## LC inside PA (to exclude from upper limit)
# lc_n2k_forests <- PU_natura_lc |>
#   dplyr::select(-c(Status)) |>
#   pivot_longer(-PUID) |>
#   mutate(place = "n2k") |>
#   filter(name == "WoodlandForest_multi"|
#          name == "WoodlandForest_prod"|
#          name == "WoodlandForest_primary") |>
#   group_by(PUID) |>
#   summarize(area_protected = sum(value, na.rm = TRUE))# cropland and pastureland not part of upper limit
#
# hist(lc_n2k_forests$area_protected)

## LC full PUT
# lc_all_forests <- PU_lc |>
#   dplyr::select(-c(Status)) |>
#   pivot_longer(-PUID) |>
#   mutate(place = "all") |>
#   filter(name == "WoodlandForest_multi"|
#            name == "WoodlandForest_prod"|
#            name == "WoodlandForest_primary") |>
#   group_by(PUID) |>
#   summarise(landarea = ifelse(sum(value, na.rm = TRUE)>10000,10000,sum(value, na.rm = TRUE)))
# hist(lc_all_forests$landarea)
#
# forestry_bounds <- lc_all_forests |>
#   #left_join(lc_n2k_forests) |>
#   #left_join(PU_potential_lc, by = "PUID") |>
#   mutate(upper = landarea*1.1) |>
#   mutate(upper = ifelse(upper > 10000, 10000, upper)) |>
#   mutate(upper = ifelse(upper <0,0,upper)) |>
#   mutate(lower = 0) |> dplyr::select(-c(landarea)) |>
#   mutate(zone = "WoodlandForest_prod_production")


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
  mutate(upper = ifelse(upper >1,1, upper)) #|>
 # mutate(lower = ifelse(lower<0.01, lower, lower-0.01))
  # mutate(lower = round(lower, 2)) |>
  # mutate(upper = round(upper, 2))

manual_bounded_constraints |>
  ggplot(aes(x = upper)) + geom_histogram() +
  facet_wrap(~zone) + theme_classic()

unique(manual_bounded_constraints$zone)
if(apply_initialglobiom){
  write_csv(manual_bounded_constraints, "data/formatted-data/manual_bounded_constraints_production_globiom_initialGLOBIOM.csv")
} else {
  write_csv(manual_bounded_constraints, "data/formatted-data/manual_bounded_constraints_production_globiom.csv")
}


################## production flexibility ###############

nuts2 <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  st_transform(crs = st_crs(rast_template))

forestry <- st_read("data/production_targets/NUTS_TotalSumHa__FIT455_ver3.shp") |>
  as_tibble() |>
  dplyr::select(-c(geometry,NUTS2))|>
  rename(NUTS_ID = NURGCDL,
         WoodlandForest_prod = prdct__,
         WoodlandForest_multi = mltfn__,
         WoodlandForest_primary = stsd_f_) |>
  pivot_longer(-NUTS_ID) |> mutate(value = value/100/100) |>
  arrange(-value) |>
  dplyr::filter(name == "WoodlandForest_prod")

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
if(apply_initialglobiom){
  of <- "data/formatted-data/manual_bounded_constraints_production_globiom_initialGLOBIOM.csv"
} else {
  of <- "data/formatted-data/manual_bounded_constraints_production_globiom.csv"
}
manual_bounded_constraints <- read_csv(of) |>
  rename(pu = PUID) |>
  left_join(pu_in_EU) |>
  drop_na(EU_id) |>
  mutate(pu = EU_id) |>
  left_join(zones) |>
  #dplyr::select(-c(zone)) |>
  #rename(zone=name) |>
  dplyr::select(pu, zone, lower, upper, nuts2id)

n2k_pu <- manual_bounded_constraints |>
  separate(zone, c('maes_label', 'intensity', 'action'), sep = "_") |>
  filter(action == "conserve") |>
  group_by(pu) |> summarise(n2k = sum(lower))

manual_bounded_constraints <- manual_bounded_constraints |>
  left_join(n2k_pu)
#
nuts2_shp <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  mutate(country = substr(NUTS2, start = 1, stop = 2))
#
nuts2_all <- nuts2_shp |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nuts2id = seq(1:260)) |>
  dplyr::filter(country != "UK") |>
  dplyr::select(NUTS_ID, nuts2id)
#
#
forestry_conflicts <- manual_bounded_constraints |> as_tibble() |>
  filter(zone == "WoodlandForest_prod_production") |>
  group_by(nuts2id) |>
  summarise(upper = sum(upper-n2k),
            count = n()) |>
  left_join(nuts2_all) |>
  dplyr::select(-geometry) |>
  left_join(forestry) |> mutate(difference = upper -value) |>
  arrange(difference) |>
  dplyr::select(nuts2id, difference, NUTS_ID)

# crop conflicts
crop_intense <- read_csv("data/production_targets/EUCropland__primes_ref_2020REFERENCE_ver3.csv") |>
  group_by(Intensity_reclass,NUTS2) |>
  mutate(area_1000ha = replace_na(area_1000ha,0)) |>
  summarise(area_1000ha = sum(area_1000ha, na.rm = TRUE)) |> ungroup() |>
  pivot_wider(names_from = "Intensity_reclass", values_from = "area_1000ha") |>
  mutate(IntenseCropland = replace_na(IntenseCropland,0),
         LightCropland = replace_na(LightCropland,0),
         MinimalCropland = replace_na(MinimalCropland,0),
         PermanentCropland = replace_na(PermanentCropland,0)) |>
  mutate(IntenseCropland = IntenseCropland + PermanentCropland/3,
         LightCropland = LightCropland + PermanentCropland/3,
         MinimalCropland = MinimalCropland + PermanentCropland/3) |>
  dplyr::select(-PermanentCropland) |>
  pivot_longer(-NUTS2) |>
  rename(NUTS_ID = NUTS2) |>
  dplyr::filter(name == "IntenseCropland")

crop_conflicts <- manual_bounded_constraints |> as_tibble() |>
  filter(zone == "Cropland_high_production") |>
  group_by(nuts2id) |>
  summarise(upper = sum(upper-n2k),
            count = n()) |>
  left_join(nuts2_all) |>
  dplyr::select(-geometry) |>
  left_join(crop_intense) |> mutate(difference = upper -value/10) |>
  arrange(difference) |>
  dplyr::select(nuts2id, difference, NUTS_ID)
#
# pasture conflicts
pasture <- read_csv("data/production_targets/EUPasture__primes_MIX55_V2GHG_CO2_10_FIX_BLTrd_ver3.csv") |>
  dplyr::select(NUTS2, Intensity, area_1000ha) |>
  rename(NUTS_ID = NUTS2,
         name = Intensity,
         value = area_1000ha) |>
  mutate(name = ifelse(name == "LowIntensityPasture", "pasture_low", "pasture_high")) |>
  filter(name == "pasture_high")

pasture_conflicts <- manual_bounded_constraints |> as_tibble() |>
  filter(zone == "Pasture_high_production") |>
  group_by(nuts2id) |>
  summarise(upper = sum(upper-n2k),
            count = n()) |>
  left_join(nuts2_all) |>
  dplyr::select(-geometry) |>
  left_join(pasture) |> mutate(difference = upper -value/10) |>
  arrange(difference) |>
  dplyr::select(nuts2id, difference, NUTS_ID)


production_flex <- manual_bounded_constraints |>
  left_join(forestry_conflicts) |>
  mutate(upper = ifelse((difference < 0 & zone == "WoodlandForest_prod_production"),
                        1, upper)) |> dplyr::select(-difference) |>
  left_join(pasture_conflicts) |>
  mutate(upper = ifelse((difference < 0 & zone == "Pasture_high_production"),
                        1, upper)) |> dplyr::select(-c(difference, NUTS_ID)) |>
  left_join(crop_conflicts) |>
  mutate(upper = ifelse((difference < 0 & zone == "Cropland_high_production"),
                        1, upper)) |> dplyr::select(-c(difference, NUTS_ID, n2k))
#
if(apply_initialglobiom){
  write_csv(production_flex, "data/formatted-data/manual_bounded_constraints_production_globiom_flex_initialGLOBIOM.csv")
} else {
  write_csv(production_flex, "data/formatted-data/manual_bounded_constraints_production_globiom_flex.csv")
}


############ production flexibility ###########
############
################## forestry infeasibility ###############
#
# nuts2 <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
#   st_transform(crs = st_crs(rast_template))
#
# pu_in_EU <- read_csv("data/formatted-data/pu_in_EU.csv")
#
# pu <- read_fst("data/formatted-data/pu_data.fst") |>
#   left_join(pu_in_EU) |>
#   rename(id = EU_id) |>
#   dplyr::select(-c(pu, nuts2id)) |>
#   drop_na(id)
#
# zones <- read_csv("data/formatted-data/zone_id.csv") |>
#   mutate(name = paste0("z", id))
#
# manual_bounded_constraints <- read_csv("data/formatted-data/manual_bounded_constraints_production_globiom.csv") |>
#   rename(pu = PUID) |>
#   left_join(pu_in_EU) |>
#   drop_na(EU_id) |>
#   mutate(pu = EU_id) |>
#   left_join(zones) |>
#   #dplyr::select(-c(zone)) |>
#   #rename(zone=name) |>
#   dplyr::select(pu, zone, lower, upper, nuts2id) |>
#   drop_na()
#
# nuts2_shp <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
#   mutate(country = substr(NUTS2, start = 1, stop = 2))
#
# nuts2_all <- nuts2_shp |>
#   rename(NUTS_ID = NURGCDL2) |>
#   mutate(nuts2id = seq(1:260)) |> dplyr::select(NUTS_ID, nuts2id)
#
# production <- c("Cropland_high_production",
#                 "WoodlandForest_prod_production",
#                 "Cropland_med_production",
#                 "WoodlandForest_multi_production",
#                 "Pasture_high_production")
#
# production_flex <- manual_bounded_constraints |>
#   mutate(forests = ifelse(zone %in% c("WoodlandForest_primary_restore",
#                             "WoodlandForest_multi_restore",
#                             "WoodlandForest_multi_conserve", "WoodlandForest_multi_production"), 1, 0)) |>
#   group_by(pu) |>
#   mutate(forests = sum(forests)) |>
#   ungroup() |>
#   mutate(forests = ifelse(forests>0,1,0)) |>
#   mutate(upper = ifelse(zone %in% "WoodlandForest_prod_production" &
#                           forests > 0, 1, upper))
#
# unique(production_flex$zone)
# write_csv(production_flex, "data/formatted-data/manual_bounded_constraints_production_flex2.csv")
#
# read_csv("data/formatted-data/manual_bounded_constraints_production_flex2.csv")
########## Get names of all zones ################3
# write out for future use!
all_zones <- unique(manual_bounded_constraints$zone)
zone_ID <- as.data.frame(all_zones) |>
  rename(zone = all_zones) |>
  mutate(id = seq(1:length(all_zones)))
write_csv(zone_ID, "data/formatted-data/zone_id.csv")

