library(fst)
library(tidyverse)
library(raster)
library(terra)
library(sf)
library(fasterize)
library(exactextractr)
rm(list = ls())


############### Load data #############

# map PU onto EU and nuts2
pu_in_EU <- read_csv("data/formatted-data/pu_in_EU.csv")
PU_template <- st_read("data/formatted-data/PU_template.shp") |>
  left_join(pu_in_EU |> rename(PUID= pu))
raster <- raster("data/landcover/10km/Corine_2018_cropland.tif")
PU_raster <- fasterize(PU_template, raster, field = "PUID")
#plot(PU_raster)

# grab nuts2 info
nuts2 <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nutsIDnum = seq(1:260)) |>
  #st_transform(crs = CRS(PU_template)) |>
  mutate(country = substr(NUTS_ID, start = 1, stop = 2)) |>
  filter(country != "UK") |> dplyr::select(-country)
nuts2_raster <- fasterize(nuts2, PU_raster, field = "nutsIDnum", fun = "first")
nuts2_ID <- as_tibble(nuts2) |>
  dplyr::select(NUTS_ID, nutsIDnum) |>
  rename(nuts2id = nutsIDnum)
names(nuts2_raster) <- "nuts2id"

nuts_shp <- st_read("data/GlobiomJune2024/NUTS_TotalSumHa__2020__BAU_FV02.shp") |>
  dplyr::select(NURGCDL, NUTS2, geometry)

# grab pu_lc data (initial conditions)
PU_lc <- read_csv("data/outputs/2-zones/PU_lc_intensity.csv")

###################### proportional constraints ###############

# Forest
forest_2020_ref <- as_tibble(st_read("data/GlobiomJune2024/NUTS_TotalSumHa__2020__BAU_FV02.shp")) |>
  dplyr::select(-geometry) |>
  pivot_longer(-c(NURGCDL, NUTS2)) |>
  rename(area_2020_ref = value)

forest_2020_f455 <- as_tibble(st_read("data/GlobiomJune2024/NUTS_TotalSumHa__2020__FIT455_FV02.shp")) |>
  dplyr::select(-geometry) |>
  pivot_longer(-c(NURGCDL, NUTS2)) |>
  rename(area_2020_f455 = value)

forest_2030_ref <- as_tibble(st_read("data/GlobiomJune2024/NUTS_TotalSumHa__2030__BAU_FV02.shp")) |>
  dplyr::select(-geometry) |>
  pivot_longer(-c(NURGCDL, NUTS2)) |>
  rename(area_2030_ref = value)

forest_2030_f455 <- as_tibble(st_read("data/GlobiomJune2024/NUTS_TotalSumHa__2030__FIT455_FV02.shp")) |>
  dplyr::select(-geometry) |>
  pivot_longer(-c(NURGCDL, NUTS2)) |>
  rename(area_2030_f455 = value)

forest_constraints <-  forest_2020_f455 |>
  left_join(forest_2020_ref) |>
  left_join(forest_2030_ref) |>
  left_join(forest_2030_f455) |>
  mutate_if(is.numeric, function(x) x / 1000) |>
  mutate(f455_change = (area_2030_f455 - area_2020_f455),
         ref_change = (area_2030_ref - area_2020_ref)) |>
  filter(name != "stsd_f_") |>
  arrange(-ref_change) |>
  dplyr::select(-NURGCDL)

forest_constraints |> group_by(name) |>
  summarise(area_2020_f455 = sum(area_2020_f455, na.rm = T),
            area_2020_ref = sum(area_2020_ref, na.rm = T))
# crop

crop_2020_ref <- read_csv("data/GlobiomJune2024/EUCropland__2020primes_ref_2020REFERENCE_FV02.csv") |>
  dplyr::select(NUTS2,area_1000ha, Intensity_reclass) |>
  group_by(NUTS2, Intensity_reclass) |>
  summarise(area_1000ha = sum(area_1000ha)) |>
  pivot_wider(names_from = Intensity_reclass, values_from = area_1000ha) |>
  mutate(IntenseCropland = IntenseCropland + PermanentCropland/3,
         LightCropland = LightCropland + PermanentCropland/3,
         MinimalCropland = MinimalCropland + PermanentCropland/3) |>
  dplyr::select(-PermanentCropland) |>
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) |>
  pivot_longer(-NUTS2) |>
  rename(area_2020_ref = value)

crop_2020_f455 <- read_csv("data/GlobiomJune2024/EUCropland__2020primes_MIX55_V2GHG_CO2_10_FIX_BLTrd_FV02.csv") |>
  dplyr::select(NUTS2,area_1000ha, Intensity_reclass) |>
  group_by(NUTS2, Intensity_reclass) |>
  summarise(area_1000ha = sum(area_1000ha)) |>
  pivot_wider(names_from = Intensity_reclass, values_from = area_1000ha) |>
  mutate(IntenseCropland = IntenseCropland + PermanentCropland/3,
         LightCropland = LightCropland + PermanentCropland/3,
         MinimalCropland = MinimalCropland + PermanentCropland/3) |>
  dplyr::select(-PermanentCropland) |>
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) |>
  pivot_longer(-NUTS2) |>
  rename(area_2020_f455 = value)

crop_2030_ref <- read_csv("data/GlobiomJune2024/EUCropland__2030primes_ref_2020REFERENCE_FV02.csv") |>
  dplyr::select(NUTS2,area_1000ha, Intensity_reclass) |>
  group_by(NUTS2, Intensity_reclass) |>
  summarise(area_1000ha = sum(area_1000ha)) |>
  pivot_wider(names_from = Intensity_reclass, values_from = area_1000ha) |>
  mutate(IntenseCropland = IntenseCropland + PermanentCropland/3,
         LightCropland = LightCropland + PermanentCropland/3,
         MinimalCropland = MinimalCropland + PermanentCropland/3) |>
  dplyr::select(-PermanentCropland) |>
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) |>
  pivot_longer(-NUTS2) |>
  rename(area_2030_ref = value)

crop_2030_f455 <- read_csv("data/GlobiomJune2024/EUCropland__2030primes_MIX55_V2GHG_CO2_10_FIX_BLTrd_FV02.csv") |>
  dplyr::select(NUTS2,area_1000ha, Intensity_reclass) |>
  group_by(NUTS2, Intensity_reclass) |>
  summarise(area_1000ha = sum(area_1000ha)) |>
  pivot_wider(names_from = Intensity_reclass, values_from = area_1000ha) |>
  mutate(IntenseCropland = IntenseCropland + PermanentCropland/3,
         LightCropland = LightCropland + PermanentCropland/3,
         MinimalCropland = MinimalCropland + PermanentCropland/3) |>
  dplyr::select(-PermanentCropland) |>
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) |>
  pivot_longer(-NUTS2) |>
  rename(area_2030_f455 = value)

crop_constraints <-  crop_2020_ref |>
  left_join(crop_2020_f455) |>
  left_join(crop_2030_f455) |>
  left_join(crop_2030_ref) |>
  mutate(f455_change = (area_2030_f455 - area_2020_f455),
         ref_change = (area_2030_ref - area_2020_ref))

# pasture

pasture_2020_ref <- read_csv("data/GlobiomJune2024/EUPasture__2020primes_ref_2020REFERENCE_FV02.csv") |>
  dplyr::select(NUTS2,area_1000ha, Intensity) |>
  rename(area_2020_ref = area_1000ha)

pasture_2020_f455 <- read_csv("data/GlobiomJune2024/EUPasture__2020primes_MIX55_V2GHG_CO2_10_FIX_BLTrd_FV02.csv") |>
  dplyr::select(NUTS2,area_1000ha, Intensity) |>
  rename(area_2020_f455 = area_1000ha)

pasture_2030_ref <- read_csv("data/GlobiomJune2024/EUPasture__2030primes_ref_2020REFERENCE_FV02.csv") |>
  dplyr::select(NUTS2,area_1000ha, Intensity) |>
  rename(area_2030_ref = area_1000ha)

pasture_2030_f455 <- read_csv("data/GlobiomJune2024/EUPasture__2030primes_MIX55_V2GHG_CO2_10_FIX_BLTrd_FV02.csv") |>
  dplyr::select(NUTS2,area_1000ha, Intensity) |>
  rename(area_2030_f455 = area_1000ha)

pasture_constraints <-  pasture_2020_ref |>
  left_join(pasture_2020_f455) |>
  left_join(pasture_2030_f455) |>
  left_join(pasture_2030_ref) |>
  mutate(f455_change = (area_2030_f455 - area_2020_f455),
         ref_change = (area_2030_ref - area_2020_ref)) |>
  rename(name = Intensity)

## All constraints
all_constraints <- forest_constraints |>
  bind_rows(crop_constraints, pasture_constraints) |>
  left_join(nuts_shp) |>
  mutate(area = as.numeric(st_area(geometry))*0.0001) |>
  mutate(
    name = case_when(
      name == "mltfn__" ~ "WoodlandForest_multi",
      name == "prdct__" ~ "WoodlandForest_prod",
      name == "IntenseCropland" ~ "Cropland_high",
      name == "LightCropland" ~ "Cropland_med",
      name == "MinimalCropland" ~ "Cropland_low",
      name == "HighIntensityPasture" ~ "Pasture_high",
      name == "LowIntensityPasture" ~ "Pasture_low",
      TRUE ~ name  # Keep other values as they are
    )
  ) |>
  mutate(area = area/1000) |>
  mutate(perc_f455_change = f455_change/area,
         perc_ref_change = ref_change/area)

all_constraints |> rename(zone = name) |> group_by(zone) |>
  summarise(area_2020_f455 = sum(area_2020_f455, na.rm = T),
            area_2020_ref = sum(area_2020_ref,  na.rm = T),
            area_2030_ref = sum(area_2030_ref,  na.rm = T),
            area_2030_f455 = sum(area_2030_f455,  na.rm = T)) |>
  pivot_longer(-c(zone))|>
  ggplot(aes(x = name, y = value, fill = zone))  + geom_col()

perc_production <- all_constraints |> group_by(NUTS2) |>
  summarise(a2020f = sum(area_2030_f455),
            a2020r = sum(area_2030_ref),
            area = mean(area)) |>
  mutate(perc_prod_2020_f455 = a2020f/area,
         perc_prod_2020_ref = a2020r/area) |>
  dplyr::select(perc_prod_2020_f455, perc_prod_2020_ref, NUTS2) |>
  arrange(-perc_prod_2020_f455)

diff_table <- all_constraints |>
  mutate(diff_2020 = area_2020_ref- area_2020_f455) |>
  arrange(-diff_2020) |>
  dplyr::select(-NURGCDL)
write_csv(as_tibble(diff_table)|>dplyr::select(-geometry), "data/diff_table_ic.csv")


a <- all_constraints |> mutate(area = area) |>
  mutate(perc_ref_change = ref_change/area) |>
  ggplot(aes(fill = perc_ref_change*100)) +
  geom_sf(aes(geometry = geometry)) +
  #theme_map() +
  facet_wrap(~name) +
  theme(legend.position = "bottom") +
  scale_fill_gradient2()

ggsave("figures/SIchange.png",a, height = 11, width = 5)

# NUTS_ID, name, value
percentage_constraints <- PU_lc |> left_join(pu_in_EU |>
                          rename(PUID= pu) |>
                          group_by(nuts2id) |> mutate(nuts_area = n())) |>
  left_join(nuts2_ID |> rename(NUTS2 = NUTS_ID)) |>
  filter(!is.na(EU_id)) |>
  pivot_longer(-c(PUID, nuts2id, EU_id, NUTS2, Status, nuts_area)) |>
  group_by(NUTS2) |> mutate(nuts_area = sum(value, na.rm = T)) |> ungroup() |>
  group_by(NUTS2, name, nuts_area) |>
  summarise(area_2020_IC = sum(value, na.rm = T)) |>
  mutate(percent_area_2020 = area_2020_IC/nuts_area) |>
  left_join(all_constraints |>
              dplyr::select(perc_ref_change, perc_f455_change, NUTS2, name)) |>
  filter(name %in% unique(all_constraints$name)) |>
  mutate(ref_min = (percent_area_2020 + perc_ref_change)*nuts_area,
         f455_min = (percent_area_2020 + perc_f455_change)*nuts_area) |>
  rename(NUTS_ID = NUTS2) |>
  mutate(ref_min = ifelse(ref_min <0, 0, ref_min),
         f455_min = ifelse(f455_min<0,0,f455_min))

mins <- percentage_constraints |> group_by(nuts_area, NUTS_ID) |>
  summarise(area_2020_IC = sum(area_2020_IC),
            ref_min = sum(ref_min, na.rm = T),
            f455_min = sum(f455_min, na.rm = T)) |>
  mutate(diff = nuts_area-f455_min) |>
  arrange(diff)

######################### Save formatted CSV #############

# crops

percentage_constraints |>
  filter(name == "Cropland_high") |>
  rename(value = ref_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_crop_high_ref_proportional.csv")

percentage_constraints |>
  filter(name == "Cropland_high") |>
  rename(value = f455_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_crop_high_55_proportional.csv")

percentage_constraints |>
  filter(name == "Cropland_med") |>
  rename(value = ref_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_crop_med_ref_proportional.csv")

percentage_constraints |>
  filter(name == "Cropland_med") |>
  rename(value = f455_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_crop_med_55_proportional.csv")

percentage_constraints |>
  filter(name == "Cropland_low") |>
  rename(value = ref_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_crop_low_ref_proportional.csv")

percentage_constraints |>
  filter(name == "Cropland_low") |>
  rename(value = f455_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_crop_low_55_proportional.csv")

## pastures

percentage_constraints |>
  filter(name == "Pasture_high") |>
  rename(value = ref_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_pasture_high_ref_proportional.csv")

percentage_constraints |>
  filter(name == "Pasture_high") |>
  rename(value = f455_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_pasture_high_55_proportional.csv")

percentage_constraints |>
  filter(name == "Pasture_low") |>
  rename(value = ref_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_pasture_low_ref_proportional.csv")

percentage_constraints |>
  filter(name == "Pasture_low") |>
  rename(value = f455_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_pasture_low_55_proportional.csv")

## forests

percentage_constraints |>
  filter(name == "WoodlandForest_multi") |>
  rename(value = ref_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_forest_multi_ref_proportional.csv")

percentage_constraints |>
  filter(name == "WoodlandForest_multi") |>
  rename(value = f455_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_forest_multi_55_proportional.csv")


percentage_constraints |>
  filter(name == "WoodlandForest_prod") |>
  rename(value = ref_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_forest_prod_ref_proportional.csv")

percentage_constraints |>
  filter(name == "WoodlandForest_prod") |>
  rename(value = f455_min) |> dplyr::select(NUTS_ID, name, value) |>
  write_csv("data/formatted-data/linear_constraints/nuts_forest_prod_55_proportional.csv")


## country constraints

cc <- percentage_constraints |> mutate(country = substr(NUTS_ID, 1, 2)) |>
  group_by(country, name) |>
  summarise(
    area_2020_IC_prod = sum(area_2020_IC),
    nuts_area = sum(nuts_area),
    value = sum(f455_min, na.rm = T),0) |>
  group_by(country) |>
  mutate(perc_production = sum(value)/nuts_area*100)|>
  dplyr::select(country, name, value)

write_csv(cc, "data/formatted-data/linear_constraints/country_constraints_55_proportional.csv")


cc <- percentage_constraints |> mutate(country = substr(NUTS_ID, 1, 2)) |>
  group_by(country, name) |>
  summarise(
    area_2020_IC_prod = sum(area_2020_IC),
    nuts_area = sum(nuts_area),
    value = sum(ref_min, na.rm = T),0) |>
  group_by(country) |>
  mutate(perc_production = sum(value)/nuts_area*100)|>
  dplyr::select(country, name, value)

write_csv(cc, "data/formatted-data/linear_constraints/country_constraints_ref_proportional.csv")


################# Create PU budget data #######################

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

#### Pasture low #####

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
          "data/formatted-data/linear_constraints/pu_pasture_low_budget_data.csv")

### Pasture high ###

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
          "data/formatted-data/linear_constraints/pu_pasture_high_budget_data.csv")

#### Crop high ###

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
          "data/formatted-data/linear_constraints/pu_crop_high_budget_data.csv")

##### crop medium ###

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
          "data/formatted-data/linear_constraints/pu_crop_med_budget_data.csv")

# crop low
budget_crop_low_data <- PU_zones_cost

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
          "data/formatted-data/linear_constraints/pu_crop_low_budget_data.csv")

# forestry multi

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
          "data/formatted-data/linear_constraints/pu_forest_multi_budget_data.csv")

# forestry prod

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
          "data/formatted-data/linear_constraints/pu_forest_prod_budget_data.csv")

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

