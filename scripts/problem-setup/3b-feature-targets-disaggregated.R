library(terra)
library(sf)
library(tidyverse)
library(raster)
library(fasterize)
library(fst)
rm(list = ls())

## PU template
filelist_temp <- list.files("data/SpeciesData/CurrentSDMs/", recursive = TRUE)
spp <- rast(paste0("data/SpeciesData/CurrentSDMs/", filelist_temp[1]))
# Grab the PU template (same used for LC analysis)
PU_template <- raster("data/landcover/10km/Corine_2018_cropland.tif") |>
  rasterToPolygons() |>
  st_as_sf() |>
  #st_transform(crs = st_crs(natura)) |>
  dplyr::mutate(PUID = seq(1:length(geometry))) |>
  dplyr::select(-Corine_2018_cropland) 
# make it into a raster
PU_raster <- fasterize(sf = PU_template, raster = raster(spp[[1]]), field = "PUID")

# bioregions
bioregion <- st_read("data/EU_BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") |>
  rename(BR_ID = short_name) |>
  mutate(BRIDnum = PK_UID) |>
  st_transform(crs = crs(PU_template))
BR_raster <- fasterize(bioregion, PU_raster, field = "BRIDnum", fun = "max")

#BR <- raster("data/formatted-data/BR_raster.tif", overwrite = TRUE)

# countries
nuts2_shp <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  mutate(country = substr(NUTS2, start = 1, stop = 2)) |>
  dplyr::select(country) |>
  st_transform(crs = crs(BR_raster)) |>
  filter(country != "UK")

country_code <- as_tibble(nuts2_shp) |>
  dplyr::select(-geometry) |> unique() |>
  mutate(country_code = seq(1:27))

country_raster <- nuts2_shp |> left_join(country_code) |>
  fasterize(BR_raster, field = "country_code") |>
  crop(BR_raster) 

# bioregion-country-PU

stack <- stack(PU_raster,BR_raster,  country_raster)
names(stack) <- c("pu","BR", "country")

bioregion_country_split <- as.data.frame(stack)

species_id <- readRDS("data/formatted-data/MinSpeciesToCover.rds") |>
  dplyr::select(speciesname, taxon_id) |>
  rename(feature = taxon_id)
species_id$speciesname <- gsub(" ", "_", species_id$speciesname, fixed=TRUE)

PU_lc <- read_csv("data/outputs/2-zones/PU_lc_intensity.csv") |>
  rename(pu = PUID) |>
  dplyr::select(-c(Status)) |>
  pivot_longer(-pu) |>
  rename(MAES= name, area = value)


zone_id <- read_csv("data/formatted-data/zone_id.csv") |>
  mutate(zone_sep = zone) |>
  separate(zone_sep, c('maes_label', 'level', 'action'), sep = "_") |>
  mutate(hab_pref = 1) |>
  filter(action %in% c("conserve", "production", "lockin")) |>
  mutate(MAES = paste0(maes_label, "_", level)) |>
  dplyr::select(id, MAES) |> rename(zone = id)


bioregion_country_split <- bioregion_country_split |> 
  drop_na() |>
  mutate(BR = sprintf("%02d",BR),
         country = sprintf("%02d",country)) |>
  mutate(BR_cou = sprintf("%04d",as.numeric(paste0(BR, country)))) |>
  dplyr::select(pu, BR_cou)

features_spp <- read_fst("data/formatted-data/features_data_all.fst") |>
  filter(feature != 999999) |>
  left_join(bioregion_country_split) 

features_split_spp <- features_spp |>
  drop_na() |>
  mutate(feature = paste0(feature, BR_cou)) |>
  mutate(feature = as.numeric(feature)) |>
  dplyr::select(-BR_cou)

carbon <- 
  read_fst("data/formatted-data/features_data_all.fst") |>
  filter(feature == 999999)

hist(carbon$amount)
features_split_spp |> bind_rows(carbon) |>
  write_fst("data/formatted-data/features_split.fst", compress = 75)

### recalculate with all zones and potential for setting targets

# landcover area per PU relevant to zone - conservation and production
# mutate to get total AOH for species
# mutate to get the AOH percentage in the country/biome
# target = AOH*AOH percentage 
features_split_spp_targets <- read_fst("data/formatted-data/features_data_spp_fortargets.fst") |>
  filter(feature != 999999) |>
  left_join(bioregion_country_split) |>
  drop_na() |>
  mutate(feature = paste0(feature, BR_cou)) |>
  mutate(feature = as.numeric(feature)) |>
  dplyr::select(-BR_cou)

features_split_spp_targets 


lc_targets <- PU_lc |>
  left_join(zone_id) |> dplyr::select(-MAES) |>
  mutate(zone = as.integer(zone))

hist(lc_targets$area)

targets_split_spp <- features_split_spp_targets |>
  left_join(lc_targets)

hist(targets_split_spp$area)

targets_split <- targets_split_spp |>
  mutate(spp = round(feature, -4)) |>
  drop_na(area) |>
  group_by(spp) |>
  mutate(AOH = sum(area*amount)) |> ungroup() |>
  group_by(feature) |>
  summarise(AOH_split = sum(area*amount),
            AOH = mean(AOH)) |>
  ungroup() |>
  mutate(perc_target = AOH_split/AOH) |>
  mutate(target = AOH) |>
  mutate(target = ifelse(target > 10^4, 10^4, target))  |> 
  mutate(target = ifelse(target <22, 22, target)) |>
  # mutate(target = target/100) |>
  mutate(target = target*perc_target) |>
  dplyr::select(target, feature)

hist(targets_split$target)
# ZONES
z <- read_csv("data/formatted-data/zone_id.csv") |>
  mutate(name = paste0("z", id)) |>
  relocate(id, 1) |>
  dplyr::select(-zone)

filelist_temp <- list.files("data/PotentialCarbonDistribution/")
carbon_potential <- rast(paste0("data/PotentialCarbonDistribution/", filelist_temp))
carbon_potential <- stack(carbon_potential)
carbon_target <- max(carbon_potential)
carbon_target <- carbon_target/100
carbon_target <- cellStats(carbon_target,stat="sum")
carbon_target <- tibble::tibble(id = 999999, name = "carbon", target = carbon_target)

carbon_target <- carbon_target |> mutate(feature = id) |>
  dplyr::select(-c(id, name))

targets_split |> 
  bind_rows(carbon_target) |>
  mutate(sense = ">=",
         type = "absolute",
         zone = list(z$name)) |>
  write_csv("data/formatted-data/targets_split.csv")

## formatted targets
## 


pu_in_EU <- read_csv("data/formatted-data/pu_in_EU.csv")

pu <- read_fst("data/formatted-data/pu_data.fst") |>
  left_join(pu_in_EU) |>
  rename(id = EU_id) |>
  dplyr::select(-c(pu, nuts2id)) |>
  drop_na(id)


rij <- read_fst("data/formatted-data/features_split.fst") |>
  rename(species = feature) |>
  #mutate(amount = round(amount)) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  dplyr::select(pu, species, zone, amount) |>
  drop_na(pu) |>
  mutate(amount = ifelse(amount<0.001, 0, amount))


# COST COLUMNS: names of cost cols from pu_data
cost_columns <- colnames(pu)[1:(ncol(pu)-1)]

# FEATURES
feat_rij <- data.frame(id = unique(rij$species),
                       prop = 1)

#carbon_id
carbon_id <- data.frame(
  name = "carbon",
  id = 999999
)

feat <-feat_rij |>
  #left_join(species_id) |> 
  mutate(name = as.factor(id)) |>
  dplyr::select(id,  name, prop) |>
  drop_na(id) |>
  drop_na(name)

# ZONES
z <- read_csv("data/formatted-data/zone_id.csv") |>
  mutate(name = paste0("z", id)) |>
  relocate(id, 1) |>
  dplyr::select(-zone)

# BOUNDED CONSTRAINTS
# # zone ids for this..

zones <- read_csv("data/formatted-data/zone_id.csv") |>
  mutate(name = paste0("z", id)) 

# FEATURE TARGETS
names <- readRDS("data/SpeciesData/SDMNameMatching.rds") |>
  mutate(speciesname= current_sname) |>
  rename(id = taxon_id)|>
  dplyr::select(id, speciesname)

names$speciesname <- sub(" ", "_", names$speciesname)


targs_existing <- read_csv("data/formatted-data/targets_split.csv") |>
  mutate(zone = list(z$name)) |>
  #dplyr::filter(feature %in% feat$name) |>
  mutate(target = ifelse(target > 100000, target/100, target),
         target = ifelse(target <0.001, 0, target))

setdiff(feat$id,targs_existing$feature)

# to line up with rij spp (2 spp that dont have targets bc ~0)
targs_default <- data.frame(
  feature = setdiff(feat$name, targs_existing$feature))|>
  mutate( target = 22,
          sense = ">=",
          type = "absolute",
          zone = list(targs_existing$zone[[1]]))

targs <- targs_existing |>
  mutate(feature = as.factor(feature)) |>
  bind_rows(targs_default) |>
  mutate(weight = 1) #|>
# dplyr::select(-id)

targs_rm <- data.frame(
  feature = setdiff(targs$feature, feat$id))|>
  mutate(remove = TRUE)

targs <- targs |> left_join(targs_rm) |>
  mutate(remove = replace_na(remove, FALSE)) |>
  filter(remove == FALSE) |>
  dplyr::select(-remove)

targs_join <- targs |>
  rename(species = feature) |>
  dplyr::select(species, target) |>
  mutate(species = as.numeric(species))
write_csv(targs_join, "data/formatted-data/targets_split_formatted.csv")

  
# pu, feature, amount, zone
# pu, feature, target, zones(all)
# take current target and percentage in a place

