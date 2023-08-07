library(fst)
library(tidyverse)
library(raster)
library(terra)
library(sf)
library(fasterize)
library(exactextractr)
rm(list = ls())

# Apply initial globiom corrections
# Default should be FALSE
apply_initialglobiom <- TRUE

############## Planning unit costs per zone #############
# read in zone data
zone_id <- read_csv("data/formatted-data/zone_id.csv")
all_zones <- paste0("z", zone_id$id)

# Create a value of cost = 1 (area) for each PU
PU_template <- st_read("data/formatted-data/PU_template.shp")
raster <- raster("data/landcover/10km/Corine_2018_cropland.tif")

# make it into a raster
PU_raster <- fasterize(PU_template, raster, field = "PUID")
PU_raster_cost <- PU_raster/PU_raster
PU_zones_cost <- PU_raster_cost

# make cost raster for each zone (all = area for now)
for (i in seq(1:(length(all_zones)-1))) {
  PU_zones_cost <- stack(PU_zones_cost, PU_raster_cost)
}

# add names for each zone
names(PU_zones_cost) <- all_zones

# format and export for problem
PU_zones_cost <- as.data.frame(PU_zones_cost)

pu_costs_data_format <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(PU_zones_cost))) |>
  dplyr::select(-id) |>
  rename(pu = layer) |>
  drop_na(pu)

write_fst(pu_costs_data_format, "data/formatted-data/pu_data.fst")


##### intensity #######

pasture <- read_csv("data/production_targets/EUPasture__primes_MIX55_V2GHG_CO2_10_FIX_BLTrd_ver3.csv") |>
  #read_csv("data/prod_constraints/EUPasture__primes_MIX55_V2GHG_CO2_10_ver2.csv") |>
  dplyr::select(NUTS2, Intensity, area_1000ha) |>
  rename(NUTS_ID = NUTS2,
         name = Intensity,
         value = area_1000ha) |>
  mutate(name = ifelse(name == "LowIntensityPasture", "pasture_low", "pasture_high"))|>
  mutate(country = substr(NUTS_ID, start = 1, stop = 2)) |>
  filter(country != "UK") |> dplyr::select(-country)


# crpland_mgmt <- read_csv("data/2030_CrpLndMgmt_FitFor55Scenario_11Aug2022.csv") |>
#   group_by(NUTS2, ALLTECH) |>
#   summarise(area = sum(value)) |> ungroup()|>
#   rename(NUTS_ID = NUTS2)

crpland_mgmt <- read_csv("data/production_targets/EUCropland__primes_MIX55_V2GHG_CO2_10_FIX_BLTrd_ver3.csv") |>
  #read_csv("data/prod_constraints/EUCropland__primes_MIX55_V2GHG_CO2_10_ver2.csv") |>
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
  rename(NUTS_ID = NUTS2)|>
  mutate(country = substr(NUTS_ID, start = 1, stop = 2)) |>
  filter(country != "UK") |> dplyr::select(-country)

forestry <- st_read("data/production_targets/NUTS_TotalSumHa__FIT455_ver3.shp") |>
  #st_read("data/prod_constraints/NUTS_TotalSumHa_ver2.gpkg") |>
  as_tibble() |>
  dplyr::select(-c(geometry,NUTS2))|>
  rename(NUTS_ID = NURGCDL,
         WoodlandForest_prod = prdct__,
         WoodlandForest_multi = mltfn__,
         WoodlandForest_primary = stsd_f_) |>
  pivot_longer(-NUTS_ID) |> mutate(value = value/1000) |> # to 1000 ha, same as crop
  arrange(-value) |>
  mutate(country = substr(NUTS_ID, start = 1, stop = 2)) |>
  filter(country != "UK") |> dplyr::select(-country)

crpland_mgmt |>
  group_by(name) |>
  summarise(sum = sum(value)/10)

############## Production linear budgets per NUTS jurisdiction ###############

# Get NUTS2 jurisdiction data per planning unit and make a binary stack
nuts2 <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nutsIDnum = seq(1:260)) |>
  st_transform(crs = crs(PU_template)) |>
  mutate(country = substr(NUTS_ID, start = 1, stop = 2)) |>
  filter(country != "UK") |> dplyr::select(-country)
nuts2_raster <- fasterize(nuts2, PU_raster, field = "nutsIDnum", fun = "first")
nuts2_ID <- as_tibble(nuts2) |>
  dplyr::select(NUTS_ID, nutsIDnum) |>
  rename(nuts2id = nutsIDnum)
names(nuts2_raster) <- "nuts2id"
writeRaster(nuts2_raster, "data/formatted-data/nuts2_raster.tif", overwrite = TRUE)

# add a new id for problem formatting
pu_in_EU <- read_csv("data/formatted-data/pu_in_EU.csv")

############## bioregion constraints ##########
##############
#
# # Get NUTS2 jurisdiction data per planning unit and make a binary stack
# bioregion <- st_read("data/EU_BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") |>
#   rename(BR_ID = short_name) |>
#   mutate(BRIDnum = PK_UID) |>
#   st_transform(crs = crs(PU_template)) |>
#   st_simplify(dTolerance = 100, preserveTopology = TRUE) |>
#   st_make_valid()
# #glimpse(as_tibble(bioregion))
# BR_raster <- fasterize(bioregion, PU_raster, field = "BRIDnum", fun = "max")
# plot(BR_raster)
# names(BR_raster) <- "BR_ID"
# plot(BR_raster)
# writeRaster(BR_raster, "data/formatted-data/BR_raster.tif", overwrite = TRUE)
#
# # filter out PU outside of the EU
# pu_in_EU_BR <-
#   ### add in indices for planning units in raster to be organized
#   ### not totally necessary bc we will use PUID to reduce PU
#   tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
#   ### add in cost data
#   mutate(cost = 1) %>%
#   ### add in PUID
#   bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
#   ### add in SPP potential data
#   bind_cols(as_tibble(raster::as.data.frame(nuts2_raster))) |>
#   bind_cols(as_tibble(raster::as.data.frame(BR_raster))) |>
#   drop_na(nuts2id) |>
#   rename(pu = layer) |> dplyr::select(-id)
#
# bionames <- as_tibble(bioregion) |> dplyr::select(-geometry) |>
#   mutate(BR_ID = BRIDnum)
#
# pu_in_EU_BR |> group_by(BR_ID) |> count() |> full_join(bionames)
# # add a new id for problem formatting
#
# pu_in_EU_BR |>
#   mutate(EU_id = seq(1:nrow(pu_in_EU_BR))) |>
#   dplyr::select(-cost) |>
#   write_csv("data/formatted-data/pu_in_EU_BR.csv")

#### Pasture low #####

nuts_pasture_low <- pasture |>
  filter(name == "pasture_low")
write_csv(nuts_pasture_low, "data/formatted-data/linear_constraints/nuts_pasture_low_55.csv")

budget_pasture_low_data <- PU_zones_cost
# everything but the pasture_low zone = 0
nums <- zone_id %>%
  filter(zone != "Pasture_low_production" &
         zone != "Pasture_low_restore") |>
  dplyr::select(id)
for (i in nums$id) {
  budget_pasture_low_data[[i]] <- budget_pasture_low_data[[i]] * 0
}

pu_pasture_low_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_pasture_low_data))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(nuts2_raster))) |>
  left_join(nuts2_ID) |> dplyr::select(-nuts2id) |>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id))

write_csv(pu_pasture_low_budget_data,
          "data/formatted-data/linear_constraints/pu_pasture_low_budget_data_55.csv")

### Pasture high ###

nuts_pasture_high <- pasture|>
  filter(name == "pasture_high")
write_csv(nuts_pasture_high, "data/formatted-data/linear_constraints/nuts_pasture_high_55.csv")

budget_pasture_high_data <- PU_zones_cost
# everything but the pasture_low zone = 0
nums <- zone_id %>%
  filter(zone != "Pasture_high_production") |>
  dplyr::select(id)
for (i in nums$id) {
  budget_pasture_high_data[[i]] <- budget_pasture_high_data[[i]] * 0
}

pu_pasture_high_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_pasture_high_data))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(nuts2_raster))) |>
  left_join(nuts2_ID) |> dplyr::select(-nuts2id) |>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id))

write_csv(pu_pasture_high_budget_data,
          "data/formatted-data/linear_constraints/pu_pasture_high_budget_data_55.csv")

#### Crop high ###

nuts_crop_high <- crpland_mgmt |>
  filter(name == "IntenseCropland")
write_csv(nuts_crop_high, "data/formatted-data/linear_constraints/nuts_crop_high_55.csv")

budget_crop_high_data <- PU_zones_cost
# everything but the pasture_low zone = 0
nums <- zone_id %>%
  filter(zone != "Cropland_high_production") |>
  dplyr::select(id)
for (i in nums$id) {
  budget_crop_high_data[[i]] <- budget_crop_high_data[[i]] * 0
}

pu_crop_high_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_crop_high_data))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(nuts2_raster))) |>
  left_join(nuts2_ID) |> dplyr::select(-nuts2id) |>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id))

write_csv(pu_crop_high_budget_data,
          "data/formatted-data/linear_constraints/pu_crop_high_budget_data_55.csv")

##### crop medium ###

nuts_crop_med <- crpland_mgmt |>
  filter(name == "LightCropland")
write_csv(nuts_crop_med, "data/formatted-data/linear_constraints/nuts_crop_med_55.csv")

budget_crop_med_data <- PU_zones_cost
# everything but the pasture_low zone = 0
nums <- zone_id %>%
  filter(zone != "Cropland_med_production" &
           zone != "Cropland_med_restore" ) |>
  dplyr::select(id)
for (i in nums$id) {
  budget_crop_med_data[[i]] <- budget_crop_med_data[[i]] * 0
}

pu_crop_med_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_crop_med_data))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(nuts2_raster))) |>
  left_join(nuts2_ID) |> dplyr::select(-nuts2id) |>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id))

write_csv(pu_crop_med_budget_data,
          "data/formatted-data/linear_constraints/pu_crop_med_budget_data_55.csv")


# crop low

nuts_crop_low <- crpland_mgmt |>
  filter(name == "MinimalCropland")
write_csv(nuts_crop_low, "data/formatted-data/linear_constraints/nuts_crop_low_55.csv")

budget_crop_low_data <- PU_zones_cost
# everything but the pasture_low zone = 0
nums <- zone_id %>%
  filter(zone != "Cropland_low_production" &
           zone != "Cropland_low_restore") |>
  dplyr::select(id)
for (i in nums$id) {
  budget_crop_low_data[[i]] <- budget_crop_low_data[[i]] * 0
}

pu_crop_low_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_crop_low_data))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(nuts2_raster))) |>
  left_join(nuts2_ID) |> dplyr::select(-nuts2id) |>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id))

write_csv(pu_crop_low_budget_data,
          "data/formatted-data/linear_constraints/pu_crop_low_budget_data_55.csv")

# forestry - primary
nuts_forest_pri <- forestry |>
  filter(name == "WoodlandForest_primary")
write_csv(nuts_forest_pri, "data/formatted-data/linear_constraints/nuts_forest_pri_55.csv")

budget_forest_pri_data <- PU_zones_cost
# everything but the pasture_low zone = 0
nums <- zone_id %>%
  filter(zone != "WoodlandForest_primary_restore",
         zone != "WoodlandForest_primary_conserve") |>
  dplyr::select(id)
for (i in nums$id) {
  budget_forest_pri_data[[i]] <- budget_forest_pri_data[[i]] * 0
}

pu_forest_pri_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_forest_pri_data))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(nuts2_raster))) |>
  left_join(nuts2_ID) |> dplyr::select(-nuts2id) |>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id))

write_csv(pu_forest_pri_budget_data,
          "data/formatted-data/linear_constraints/pu_forest_pri_budget_data_55.csv")


# forestry multi
nuts_forest_multi <- forestry |>
  filter(name == "WoodlandForest_multi")
write_csv(nuts_forest_multi, "data/formatted-data/linear_constraints/nuts_forest_multi_55.csv")

budget_forest_multi_data <- PU_zones_cost
# everything but the pasture_low zone = 0
nums <- zone_id %>%
  filter(zone != "WoodlandForest_multi_production",
         zone != "WoodlandForest_multi_restore") |>
  dplyr::select(id)
for (i in nums$id) {
  budget_forest_multi_data[[i]] <- budget_forest_multi_data[[i]] * 0
}

pu_forest_multi_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_forest_multi_data))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(nuts2_raster))) |>
  left_join(nuts2_ID) |> dplyr::select(-nuts2id) |>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id))

write_csv(pu_forest_multi_budget_data,
          "data/formatted-data/linear_constraints/pu_forest_multi_budget_data_55.csv")

# forestry prod
nuts_forest_prod <- forestry |>
  filter(name == "WoodlandForest_prod")
write_csv(nuts_forest_prod, "data/formatted-data/linear_constraints/nuts_forest_prod_55.csv")

budget_forest_prod_data <- PU_zones_cost
# everything but the pasture_low zone = 0
nums <- zone_id %>%
  filter(zone != "WoodlandForest_prod_production") |>
  dplyr::select(id)
for (i in nums$id) {
  budget_forest_prod_data[[i]] <- budget_forest_prod_data[[i]] * 0
}

pu_forest_prod_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_forest_prod_data))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(nuts2_raster))) |>
  left_join(nuts2_ID) |> dplyr::select(-nuts2id) |>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id))

write_csv(pu_forest_prod_budget_data,
          "data/formatted-data/linear_constraints/pu_forest_prod_budget_data_55.csv")

# Add linear constraint for wetlands!
wetland_targets <- read_csv("data/wetland_restoration_targets.csv") |>
  drop_na() |>
  separate(name, into = c("country", "NUTS_ID", "LC", "type")) |>
  filter(type == "Area") |>
  mutate(name = "Wetlands_natural_restore") |>
  dplyr::select(c("NUTS_ID", "name", "value")) |>
  arrange(-value) |>
  mutate(country = substr(NUTS_ID, start = 1, stop = 2)) |>
  filter(country != "UK") |> dplyr::select(-country)

write_csv(wetland_targets, "data/formatted-data/linear_constraints/wetland_targets.csv")
sum(wetland_targets$value)
budget_wetland_rest_data <- PU_zones_cost
# everything but the pasture_low zone = 0
nums <- zone_id %>%
  filter(zone != "Wetlands_natural_restore") |>
  dplyr::select(id)
for (i in nums$id) {
  budget_wetland_rest_data[[i]] <- budget_wetland_rest_data[[i]] * 0
}

pu_wetland_rest_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_wetland_rest_data))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(nuts2_raster))) |>
  left_join(nuts2_ID) |> dplyr::select(-nuts2id) |>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id))

write_csv(pu_wetland_rest_budget_data,
          "data/formatted-data/linear_constraints/pu_wetland_rest_budget_data.csv")

#################### set restoration budget stacks ###############
# up to 20% restoration
budget_restoration_data <- PU_zones_cost
# everything but the restoration zones = 0
nums <- zone_id %>%
  dplyr::filter(!grepl('restore', zone)) |>
  dplyr::select(id)
# everything not restoration * 0
for (i in nums$id) {
  budget_restoration_data[[i]] <- budget_restoration_data[[i]] * 0
}

pu_restoration_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_restoration_data)))

write_csv(pu_restoration_budget_data,
          "data/formatted-data/linear_constraints/pu_restoration_budget_data.csv")


## Cropland restoration > 6%
budget_restore_crop_data <- PU_zones_cost
# everything but the restoration zones = 0
nums <- zone_id %>%
  dplyr::filter(zone != "Cropland_low_restore" &
                zone != "Cropland_med_restore" &
                zone != "Pasture_low_restore") |>
  dplyr::select(id)
# everything not restoration * 0
for (i in nums$id) {
  budget_restore_crop_data[[i]] <- budget_restore_crop_data[[i]] * 0
}

pu_restore_crop_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_restore_crop_data))) |>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id))


write_csv(pu_restore_crop_budget_data,
          "data/formatted-data/linear_constraints/pu_restore_crop_budget_data.csv")


## Forest restoration > 4%
budget_restore_forest_data <- PU_zones_cost
# everything but the restoration zones = 0
nums <- zone_id %>%
  dplyr::filter(zone != "WoodlandForest_multi_restore") |>
  dplyr::select(id)
# everything not restoration * 0
for (i in nums$id) {
  budget_restore_forest_data[[i]] <- budget_restore_forest_data[[i]] * 0
}

pu_restore_forest_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_restore_forest_data))) |>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id))


write_csv(pu_restore_forest_budget_data,
          "data/formatted-data/linear_constraints/pu_restore_forest_budget_data.csv")

## Natural restoration < 5%
budget_restore_natural_data <- PU_zones_cost
# everything but the restoration zones = 0
nums <- zone_id %>%
  dplyr::filter(zone != "Grassland_natural_restore" &
                  zone != "HeathlandShrub_natural_restore" &
                  zone != "MarineTransitional_natural_restore" &
                  zone != "SparseVeg_natural_restore" &
                  # zone != "Wetlands_natural_restore" &
                  zone != "WoodlandForest_primary_restore") |>
  dplyr::select(id)
# everything not restoration * 0
for (i in nums$id) {
  budget_restore_natural_data[[i]] <- budget_restore_natural_data[[i]] * 0
}

pu_restore_natural_budget_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(budget_restore_natural_data))) |>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id))


write_csv(pu_restore_natural_budget_data,
          "data/formatted-data/linear_constraints/pu_restore_natural_budget_data.csv")


# tests
#
cropland1 <- read_csv("data/formatted-data/linear_constraints/nuts_crop_high_55.csv")
cropland2 <- read_csv("data/formatted-data/linear_constraints/nuts_crop_med_55.csv")
cropland3 <- read_csv("data/formatted-data/linear_constraints/nuts_crop_low_55.csv")
pasture1 <- read_csv("data/formatted-data/linear_constraints/nuts_pasture_high_55.csv")
pasture2 <- read_csv("data/formatted-data/linear_constraints/nuts_pasture_low_55.csv")
forestry1 <- read_csv("data/formatted-data/linear_constraints/nuts_forest_multi_55.csv")
forestry2 <- read_csv("data/formatted-data/linear_constraints/nuts_forest_prod_55.csv")
forestry3 <- read_csv("data/formatted-data/linear_constraints/nuts_forest_pri_55.csv")

nuts2_shp <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  mutate(country = substr(NUTS2, start = 1, stop = 2))
nuts2_all <- nuts2_shp |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nuts2id = seq(1:260)) |> dplyr::select(NUTS_ID, nuts2id)

nuts_PU <- pu_in_EU |>
  left_join(as_tibble(nuts2_all))  |> # rename(nuts2id = nutsIDnum))) |>
  group_by(NUTS_ID) |>
  count()  |> arrange(n)

small_nuts2 <- nuts_PU |> filter(n < 41)
write_csv(small_nuts2, "small_jurisdictions.csv")

zones <- read_csv("data/formatted-data/zone_id.csv") |>
  mutate(name = paste0("z", id))

manual_bounded_constraints <- read_csv("data/formatted-data/manual_bounded_constraints_production_globiom_flex.csv") |>
  # rename(pu = PUID) |>
  # left_join(pu_in_EU) |>
  # drop_na(EU_id) |>
  # mutate(pu = EU_id) |>
  left_join(zones) |>
  #dplyr::select(-c(zone)) |>
  #rename(zone=name) |>
  dplyr::select(pu, zone, lower, upper, nuts2id)

n2k_nuts2 <- manual_bounded_constraints |>   left_join(as_tibble(nuts2_all)) |> dplyr::select(-geometry) |>
  separate(zone, c('maes_label', 'intensity', 'action'), sep = "_") |>
  filter(action == "conserve" | maes_label == "Urban") |>
  group_by(NUTS_ID) |> summarise(n2k = sum(lower))


conflict_nuts <- cropland1 |> dplyr::select(NUTS_ID, name, value) |>
  rbind(cropland2, cropland3,
                   pasture1, pasture2, forestry1, forestry2) |>
  group_by(NUTS_ID) |>
  summarise(value = sum(value)/10) |>
  left_join(nuts_PU) |>
  mutate(perc_prod = value/n*100) |>
  arrange(-perc_prod) |>
  left_join(n2k_nuts2) |>
  mutate(conflict = n-(value + n2k)) |>
  arrange(conflict) |>
  filter(conflict < 0) |>
  dplyr::select(c(NUTS_ID, conflict))

conflict_nuts |> arrange(NUTS_ID)


read_csv("data/formatted-data/linear_constraints/nuts_crop_high_55.csv") |>
  left_join(conflict_nuts) |> mutate(conflict = replace_na(conflict,0)) |>
  mutate(value = value+conflict) |>
  mutate(value = ifelse(value < 0.1, 0, value)) |>
  dplyr::select(c(NUTS_ID, name, value)) |>
  write_csv("data/formatted-data/linear_constraints/nuts_crop_high_55_adjusted.csv")


read_csv("data/formatted-data/linear_constraints/nuts_crop_med_55.csv") |>
  left_join(conflict_nuts) |> mutate(conflict = replace_na(conflict,0)) |>
  mutate(value = value+conflict) |>
  mutate(value = ifelse(value < 0.1, 0, value)) |>
  dplyr::select(c(NUTS_ID, name, value))|>
  write_csv("data/formatted-data/linear_constraints/nuts_crop_med_55_adjusted.csv")

read_csv("data/formatted-data/linear_constraints/nuts_crop_low_55.csv") |>
  left_join(conflict_nuts) |> mutate(conflict = replace_na(conflict,0)) |>
  mutate(value = value+conflict) |>
  mutate(value = ifelse(value < 0.1, 0, value)) |>
  dplyr::select(c(NUTS_ID, name, value))|>
  write_csv("data/formatted-data/linear_constraints/nuts_crop_low_55_adjusted.csv")

read_csv("data/formatted-data/linear_constraints/nuts_forest_prod_55.csv") |>
  left_join(conflict_nuts) |> mutate(conflict = replace_na(conflict,0)) |>
  mutate(value = value+conflict) |>
  mutate(value = ifelse(value < 0.1, 0, value)) |>
  dplyr::select(c(NUTS_ID, name, value)) |>
  write_csv("data/formatted-data/linear_constraints/nuts_forest_prod_55_adjusted.csv")


read_csv("data/formatted-data/linear_constraints/nuts_forest_multi_55.csv") |>
  left_join(conflict_nuts) |> mutate(conflict = replace_na(conflict,0)) |>
  mutate(value = value+conflict) |>
  mutate(value = ifelse(value < 0.1, 0, value)) |>
  dplyr::select(c(NUTS_ID, name, value))|>
  write_csv("data/formatted-data/linear_constraints/nuts_forest_multi_55_adjusted.csv")

read_csv("data/formatted-data/linear_constraints/nuts_pasture_low_55.csv") |>
  left_join(conflict_nuts) |> mutate(conflict = replace_na(conflict,0)) |>
  mutate(value = value+conflict) |>
  mutate(value = ifelse(value < 0.1, 0, value)) |>
  dplyr::select(c(NUTS_ID, name, value))|>
  write_csv("data/formatted-data/linear_constraints/nuts_pasture_low_55_adjusted.csv")

read_csv("data/formatted-data/linear_constraints/nuts_pasture_high_55.csv") |>
  left_join(conflict_nuts) |> mutate(conflict = replace_na(conflict,0)) |>
  mutate(value = value+conflict) |>
  mutate(value = ifelse(value < 0.1, 0, value)) |>
  dplyr::select(c(NUTS_ID, name, value))|>
  write_csv("data/formatted-data/linear_constraints/nuts_pasture_high_55_adjusted.csv")

# ###### refining to limit infeas - MJ loop ########
# ######
#
# ############ production flexibility - MJ Spain ###########
# # Plan:
# # Take the pasture, cropland and forestry production estimates and
# # cap them by the maximum amount of land available in the
#

#
# pu_in_EU <- read_csv("data/formatted-data/pu_in_EU.csv")
#
# nuts2_shp <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp", quiet = TRUE) |>
#   mutate(country = substr(NUTS2, start = 1, stop = 2)) |>
#   rename(NUTS_ID = NURGCDL2) |>
#   mutate(nuts2id = seq(1:260))
#
# nuts_PU <- pu_in_EU |>
#   left_join(as_tibble(nuts2_shp)) |>
#   group_by(NUTS_ID) |>
#   count()
#
# manbound <- read_csv("data/formatted-data/manual_bounded_constraints_production_globiom_flex.csv") |>
#   left_join(zones) |>
#   dplyr::select(pu, zone, nuts2id, lower, upper)
#
# # Manual bounded conservation constraints
# manual_bounded_constraints_conserve <- manbound |>
#   left_join(zones) |>
#   dplyr::select(-c(zone)) |>
#   rename(zone=name) |>
#   filter(zone %in% zones$name[grep("conserve", zones$zone)]) |>
#   group_by(pu) |>
#   summarise(lower_cons = sum(lower))
#
# #assertthat::assert_that(all(manual_bounded_constraints_conserve$lower_cons<1))
#
# # Write a loopie
# for(produce in grep("production", unique(manbound$zone), value = T) ){
#   print(produce)
#
#   input <- switch (produce,
#                    "WoodlandForest_multi_production" = "data/formatted-data/linear_constraints/nuts_forest_multi_55.csv",
#                    "WoodlandForest_prod_production" = "data/formatted-data/linear_constraints/nuts_forest_prod_55.csv",
#                    "Pasture_high_production" = "data/formatted-data/linear_constraints/nuts_pasture_high_55.csv",
#                    "Pasture_low_production" = "data/formatted-data/linear_constraints/nuts_pasture_low_55.csv",
#                    "Cropland_med_production" = "data/formatted-data/linear_constraints/nuts_crop_med_55.csv",
#                    "Cropland_low_production" = "data/formatted-data/linear_constraints/nuts_crop_low_55.csv",
#                    "Cropland_high_production" = "data/formatted-data/linear_constraints/nuts_crop_high_55.csv"
#   )
#   #input = "data/formatted-data/linear_constraints/nuts_crop_low_ref.csv"
#   f <- read.csv(input)
#
#   target <- manbound |> as_tibble() |>
#     filter(zone == produce) |>
#     left_join(manual_bounded_constraints_conserve) |>
#     mutate(lower_cons = replace_na(lower_cons, 0))|>
#     mutate(upper = replace_na(upper, 0)) |>
#     #filter(nuts2id %in% a$nuts2id) |>
#     group_by(nuts2id) |>
#     summarise(upper = ifelse(sum((upper-lower_cons))<0,0,sum((upper-lower_cons))) *10) |>
#     left_join(nuts2_shp) |>
#     dplyr::select(-geometry) |>
#     left_join(f) |>
#     mutate(difference = upper - value) |>
#     rename(new = value) |>
#     dplyr::select(NUTS_ID, upper, new, difference) |>
#     arrange(-difference) |>
#     mutate(name = unique(f$name),
#            value = ifelse(difference>0, new*0.95, upper*0.95)) |>
#     left_join(nuts_PU) |>
#     dplyr::select(NUTS_ID,name,value)
#
#   new <- target
#   write_csv(new, paste0(tools::file_path_sans_ext(input),"_l2.csv"))
# }
#
# sum(read_csv("data/formatted-data/linear_constraints/nuts_crop_high_55.csv")$value)/10
#
# ##### bit more refining for spain
#
# problematic_nutsregions <- c("ES22", "ES23", "ES41", "ES43")
# a <- nuts2_shp |> filter(NUTS_ID %in% problematic_nutsregions)
#
# ok <- manbound |> filter(nuts2id != a$nuts2id)
# check <- manbound |> filter(nuts2id == a$nuts2id)
#
# # Write a loopie
# for(produce in grep("production", unique(manbound$zone), value = T) ){
#   print(produce)
#
#   input <- switch (produce,
#                    "WoodlandForest_multi_production" = "data/formatted-data/linear_constraints/nuts_forest_multi_55_l2.csv",
#                    "WoodlandForest_prod_production" = "data/formatted-data/linear_constraints/nuts_forest_prod_55_l2.csv",
#                    "Pasture_high_production" = "data/formatted-data/linear_constraints/nuts_pasture_high_55_l2.csv",
#                    "Pasture_low_production" = "data/formatted-data/linear_constraints/nuts_pasture_low_55_l2.csv",
#                    "Cropland_med_production" = "data/formatted-data/linear_constraints/nuts_crop_med_55_l2.csv",
#                    "Cropland_low_production" = "data/formatted-data/linear_constraints/nuts_crop_low_55_l2.csv",
#                    "Cropland_high_production" = "data/formatted-data/linear_constraints/nuts_crop_high_55_l2.csv"
#   )
#   f <- read.csv(input)
#
#   target <- manbound |> as_tibble() |>
#     filter(zone == produce) |>
#     left_join(manual_bounded_constraints_conserve) |>
#     mutate(lower_cons = replace_na(lower_cons, 0))|>
#     filter(nuts2id %in% a$nuts2id) |>
#     group_by(nuts2id) |>
#     summarise(upper = ifelse(sum((upper-lower_cons))<0,0,sum((upper-lower_cons))) *10) |>
#     left_join(nuts2_shp) |>
#     dplyr::select(-geometry) |>
#     left_join(f) |> mutate(difference = upper - value) |>
#     rename(new = value) |>
#     dplyr::select(NUTS_ID, upper, new, difference) |>
#     arrange(-difference) |>
#     mutate(name = unique(f$name),
#            value = ifelse(difference>0, new*0.7, upper*0.7)) |>
#     dplyr::select(NUTS_ID,name,value) |> arrange(-value)
#
#   new <- bind_rows(target, f |> filter(!NUTS_ID %in% problematic_nutsregions))
#   write_csv(new, paste0(tools::file_path_sans_ext(input),"_ESTuned.csv"))
# }
#
#
# ##### bit more refining for forestry
#
# problematic_nutsregions <- c("ES22", "ES23", "ES41", "ES43")
# a <- nuts2_shp |> filter(NUTS_ID %in% problematic_nutsregions)
#
# ok <- manbound |> filter(nuts2id != a$nuts2id)
# check <- manbound |> filter(nuts2id == a$nuts2id)
#
# # Write a loopie
#
# f <- read.csv("data/formatted-data/linear_constraints/nuts_forest_prod_55_l2_ESTuned.csv")
#
# target <- manbound |> as_tibble() |>
#   filter(zone == produce) |>
#   left_join(manual_bounded_constraints_conserve) |>
#   mutate(lower_cons = replace_na(lower_cons, 0))|>
#   # filter(nuts2id %in% a$nuts2id) |>
#   group_by(nuts2id) |>
#   summarise(upper = ifelse(sum((upper-lower_cons))<0,0,sum((upper-lower_cons))) *10) |>
#   left_join(nuts2_shp) |>
#   dplyr::select(-geometry) |>
#   left_join(f) |> mutate(difference = upper - value) |>
#   rename(new = value) |>
#   dplyr::select(NUTS_ID, upper, new, difference) |>
#   arrange(-difference) |>
#   mutate(name = unique(f$name),
#          value = ifelse(difference>0, new*0.9, upper*0.9)) |>
#   dplyr::select(NUTS_ID,name,value) |> arrange(-value)
#
# new <- target#bind_rows(target, f |> filter(!NUTS_ID %in% problematic_nutsregions))
# write_csv(new,"data/formatted-data/linear_constraints/nuts_forest_prod_55_l_ESTuned3.csv")
#
#
#
# ########## refining constraints to limit infeasibilities  ##########
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
# manual_bounded_constraints <-
#   read_csv("data/formatted-data/manual_bounded_constraints_production_globiom_flex.csv") |>
#   #rename(pu = PUID) |>
#   # left_join(pu_in_EU) |>
#   # drop_na(EU_id) |>
#   # mutate(pu = EU_id) |>
#   # left_join(zones) |>
#   # #dplyr::select(-c(zone)) |>
#   # #rename(zone=name) |>
#   dplyr::select(pu, zone, lower, upper, nuts2id) #|>
#   # drop_na() #|>
#   # mutate(lower = ifelse(zone == "z8", 0, lower), # wetland restore
#   #        lower = ifelse(zone == "z5", 0, lower)) |> # marin transitional restore
#   # mutate(upper = ifelse(zone == "z8", 0, upper),
#   #        upper = ifelse(zone == "z5", 0, upper))
#

# manual_bounded_constraints_conserve <- manual_bounded_constraints |>
#   left_join(zones) |>
#   dplyr::select(-c(zone)) |>
#   rename(zone=name) |>
#   filter(zone %in% c("z19", "z20", "z21", "z22", "z23", "z24", "z25", "z25")) |>
#   group_by(pu) |>
#   summarise(lower_cons = sum(lower))
# hist(manual_bounded_constraints_conserve$lower_cons)
#
# ########## forestry production limited
# ##########
#
# forestry <- read_csv("data/formatted-data/linear_constraints/nuts_forest_prod_55.csv")
# #forestry <- read_csv("data/formatted-data/linear_constraints/nuts_forest_prod.csv")
#
# forestry_nuts <- manual_bounded_constraints |> as_tibble() |>
#   filter(zone == "WoodlandForest_prod_production") |>
#   left_join(manual_bounded_constraints_conserve) |>
#   mutate(lower_cons = replace_na(lower_cons,0))|>
#   group_by(nuts2id) |>
#   summarise(upper = ifelse(sum((upper-lower_cons))<0,0,sum((upper-lower_cons))) *10) |>
#   left_join(nuts2_all) |>
#   dplyr::select(-geometry) |>
#   left_join(forestry) |> mutate(difference = upper - value) |>
#   rename(G4M = value) |>
#   dplyr::select(NUTS_ID, upper, G4M, difference) |>
#   arrange(-difference) |>
#   mutate(name = "WoodlandForest_prod",
#          value = ifelse(difference>0, G4M*0.9,
#                         upper*0.9)) |>
#   #mutate(value = ifelse(upper*0.9 < value, upper*0.9, value)) |>
#   dplyr::select(NUTS_ID,name,value) |> arrange(-value)
#
# sum(forestry_nuts$value, na.rm = TRUE)/10
# sum(forestry$value)/10
#
# ########## cropland production limited ##########
#
# croplandhigh <- read_csv("data/formatted-data/linear_constraints/nuts_crop_high.csv")
#
# crophigh_nuts <- manual_bounded_constraints |> as_tibble() |>
#   filter(zone == "Cropland_high_production") |>
#   left_join(manual_bounded_constraints_conserve) |>
#   mutate(lower_cons = replace_na(lower_cons,0))|>
#   group_by(nuts2id) |>
#   summarise(upper = ifelse(sum(upper-lower_cons)<0,0,sum(upper-lower_cons)) *10) |>
#   left_join(nuts2_all) |>
#   dplyr::select(-geometry) |>
#   left_join(croplandhigh) |> mutate(difference = upper - value) |>
#   rename(G4M = value) |>
#   dplyr::select(NUTS_ID, upper, G4M, difference) |>
#   arrange(difference) |>
#   mutate(name = "IntenseCropland",
#          value = ifelse(difference>0, G4M*.98,
#                         upper*0.95)) |>
#   dplyr::select(NUTS_ID,name,value)
#
# sum(crophigh_nuts$value, na.rm=TRUE)/10
# sum(croplandhigh$value, na.rm=TRUE)/10
#
# write_csv(crophigh_nuts, "data/formatted-data/linear_constraints/nuts_crop_high_55_l.csv")
#
#
# # mid
# croplandmed <- read_csv("data/formatted-data/linear_constraints/nuts_crop_med.csv")
#
#
# cropmed_nuts <- manual_bounded_constraints |> as_tibble() |>
#   filter(zone == "Cropland_med_production" |
#            zone == "Cropland_med_restore") |>
#   left_join(manual_bounded_constraints_conserve) |>
#   mutate(lower_cons = replace_na(lower_cons,0))|>
#   group_by(nuts2id) |>
#   summarise(upper = sum(upper)*10) |>
#   left_join(nuts2_all) |>
#   dplyr::select(-geometry) |>
#   left_join(croplandmed) |> mutate(difference = upper - value) |>
#   rename(G4M = value) |>
#   dplyr::select(NUTS_ID, upper, G4M, difference) |>
#   arrange(difference) |>
#   mutate(name = "LightCropland",
#          value = ifelse(difference>0, G4M*0.98,
#                         upper*0.95)) |>
#   dplyr::select(NUTS_ID,name,value)
# sum(cropmed_nuts$value, na.rm=TRUE)/10
# sum(croplandmed$value, na.rm=TRUE)/10
#
# write_csv(cropmed_nuts, "data/formatted-data/linear_constraints/nuts_crop_med_55_l.csv")
#
# # low
#
# croplandlow <- read_csv("data/formatted-data/linear_constraints/nuts_crop_low.csv")
#
# croplow_nuts <- manual_bounded_constraints |> as_tibble() |>
#   filter(zone == "Cropland_low_production" |
#            zone == "Cropland_low_restore") |>
#   group_by(nuts2id) |>
#   summarise(upper = sum(upper)*10) |>
#   left_join(nuts2_all) |>
#   dplyr::select(-geometry) |>
#   left_join(croplandlow) |> mutate(difference = upper - value) |>
#   rename(G4M = value) |>
#   dplyr::select(NUTS_ID, upper, G4M, difference) |>
#   arrange(difference) |>
#   mutate(name = "MinimalCropland",
#          value = ifelse(difference>0, G4M*0.98,
#                         upper*0.95)) |>
#   dplyr::select(NUTS_ID,name,value)
#
# sum(croplow_nuts$value, na.rm=TRUE)/10
#
# write_csv(croplow_nuts, "data/formatted-data/linear_constraints/nuts_crop_low_55_l.csv")
#
#
# # pasture low
#
# pasturelow <- read_csv("data/formatted-data/linear_constraints/nuts_pasture_low.csv")
#
# pasturelow_nuts <- manual_bounded_constraints |> as_tibble() |>
#   filter(zone == "Pasture_low_production" |
#            zone == "Pasture_low_restore") |>
#   group_by(nuts2id) |>
#   summarise(upper = sum(upper, na.rm= TRUE)*10) |>
#   left_join(nuts2_all) |>
#   dplyr::select(-geometry) |>
#   left_join(pasturelow) |> mutate(difference = upper - value) |>
#   rename(G4M = value) |>
#   dplyr::select(NUTS_ID, upper, G4M, difference) |>
#   arrange(difference) |>
#   mutate(name = "pasture_low",
#          value = ifelse(difference>0, G4M*0.98,
#                         upper*0.95)) |>
#   dplyr::select(NUTS_ID,name,value)
# sum(pasturelow$value, na.rm=TRUE)/10
#
# write_csv(pasturelow_nuts, "data/formatted-data/linear_constraints/nuts_pasture_low_55_l.csv")
#
#
# # pasture high
#
# pasturehigh <- read_csv("data/formatted-data/linear_constraints/nuts_pasture_high.csv")
#
# pasturehigh_nuts <- manual_bounded_constraints |> as_tibble() |>
#   filter(zone == "Pasture_high_production") |>
#   left_join(manual_bounded_constraints_conserve) |>
#   mutate(lower_cons = replace_na(lower_cons,0),
#          upper = replace_na(upper,0))|>
#   group_by(nuts2id) |>
#   summarise(upper = ifelse(sum(upper-lower_cons)<0,0,sum(upper-lower_cons)) *10) |>
#   left_join(nuts2_all) |>
#   dplyr::select(-geometry) |>
#   left_join(pasturehigh) |> mutate(difference = upper - value) |>
#   rename(G4M = value) |>
#   dplyr::select(NUTS_ID, upper, G4M, difference) |>
#   arrange(difference) |>
#   mutate(name = "pasture_high",
#          value = ifelse(difference>0, G4M*0.98,
#                         upper*0.95)) |>
#   dplyr::select(NUTS_ID,name,value)
# sum(pasturehigh$value,na.rm=TRUE)/10
# write_csv(pasturehigh_nuts, "data/formatted-data/linear_constraints/nuts_pasture_high_55_l.csv")
#
# # forestry multi
#
# forestry_multi <- read_csv("data/formatted-data/linear_constraints/nuts_forest_multi.csv")
#
# forestrymulti_nuts <- manual_bounded_constraints |> as_tibble() |>
#   filter(zone == "WoodlandForest_multi_production" |
#            zone == "WoodlandForest_multi_restore" ) |>
#   group_by(nuts2id) |>
#   summarise(upper = sum(upper, na.rm= TRUE)*10) |>
#   left_join(nuts2_all) |>
#   dplyr::select(-geometry) |>
#   left_join(forestry_multi) |> mutate(difference = upper - value) |>
#   rename(G4M = value) |>
#   dplyr::select(NUTS_ID, upper, G4M, difference) |>
#   arrange(difference) |>
#   mutate(name = "WoodlandForest_multi",
#          value = ifelse(difference>0, G4M*0.98,
#                         upper*0.95)) |>
#   dplyr::select(NUTS_ID,name,value)
#
# write_csv(forestrymulti_nuts, "data/formatted-data/linear_constraints/nuts_forest_multi_55_l.csv")
#
#
# ##### bit more refining for spain
#
# problematic_nutsregions <- c("ES22", "ES23", "ES41", "ES43")
# a <- nuts2_shp |> filter(NUTS_ID %in% problematic_nutsregions)
#
# ok <- manbound |> filter(nuts2id != a$nuts2id)
# check <- manbound |> filter(nuts2id == a$nuts2id)
#
# # Write a loopie
# for(produce in grep("production", unique(manbound$zone), value = T) ){
#   print(produce)
#
#   input <- switch (produce,
#                    "WoodlandForest_multi_production" = "data/formatted-data/linear_constraints/nuts_forest_multi_ref_l.csv",
#                    "WoodlandForest_prod_production" = "data/formatted-data/linear_constraints/nuts_forest_prod_ref_l.csv",
#                    "Pasture_high_production" = "data/formatted-data/linear_constraints/nuts_pasture_high_ref_l.csv",
#                    "Pasture_low_production" = "data/formatted-data/linear_constraints/nuts_pasture_low_ref_l.csv",
#                    "Cropland_med_production" = "data/formatted-data/linear_constraints/nuts_crop_med_ref_l.csv",
#                    "Cropland_low_production" = "data/formatted-data/linear_constraints/nuts_crop_low_ref_l.csv",
#                    "Cropland_high_production" = "data/formatted-data/linear_constraints/nuts_crop_high_ref_l.csv"
#   )
#   f <- read.csv(input)
#
#   target <- manbound |> as_tibble() |>
#     filter(zone == produce) |>
#     left_join(manual_bounded_constraints_conserve) |>
#     mutate(lower_cons = replace_na(lower_cons, 0))|>
#     filter(nuts2id %in% a$nuts2id) |>
#     group_by(nuts2id) |>
#     summarise(upper = ifelse(sum((upper-lower_cons))<0,0,sum((upper-lower_cons))) *10) |>
#     left_join(nuts2_shp) |>
#     dplyr::select(-geometry) |>
#     left_join(f) |> mutate(difference = upper - value) |>
#     rename(new = value) |>
#     dplyr::select(NUTS_ID, upper, new, difference) |>
#     arrange(-difference) |>
#     mutate(name = unique(f$name),
#            value = ifelse(difference>0, new*0.95, upper*0.75)) |>
#     dplyr::select(NUTS_ID,name,value) |> arrange(-value)
#
#   new <- bind_rows(target, f |> filter(!NUTS_ID %in% problematic_nutsregions))
#   write_csv(new, paste0(tools::file_path_sans_ext(input),"_ESTuned.csv"))
# }
#
