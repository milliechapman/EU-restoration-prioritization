# Set zone contributions for each planning unit
# for each of the species features
library(terra)
library(sf)
library(tidyverse)
library(raster)
library(fasterize)
library(fst)
rm(list = ls())


MAES <- readRDS("data/EEA_habitatpreferences.rds") |>
  filter(typeasso %in% c("Preferred", "Suitable")) |>
  mutate(maes_label = recode(maes_label,
                             `Sparserly vegetated land` = "SparseVeg",
                             `Heathland and shrub` = "HeathlandShrub",
                             `Rivers and lakes` = "RiversLakes",
                             `Woodland and forest` = "WoodlandForest",
                             `Rivers and lakes` = "RiversLakes",
                             `Marine inlets and transitional waters` = "MarineTransitional"))
## fix naming of species to match above
MAES$speciesname <- gsub(" ", "_", MAES$speciesname, fixed=TRUE)

################3 threats data
################
zone_id <- read_csv("data/formatted-data/zone_id.csv") |>
  mutate(zone_sep = zone) |>
  separate(zone_sep, c('maes_label', 'level', 'action'), sep = "_") |>
  mutate(hab_pref = 1)

species_id <- readRDS("data/formatted-data/MinSpeciesToCover.rds") |>
  dplyr::select(speciesname, taxon_id)
species_id$speciesname <- gsub(" ", "_", species_id$speciesname, fixed=TRUE)
# spp level threats
threats_17 <- read_csv("data/SpeciesData/HabitatsDirectiveSpeciesThreats.csv") |>
  dplyr::select(taxon_id, pressure_code)
threats_12 <- read_csv("data/SpeciesData/BirdsDirectiveSpeciesThreats.csv") |>
  dplyr::select(taxon_id, pressure_code)
spp_threats <- rbind(threats_17, threats_12)

# threats by land production intensity
lc_threats <- read_csv("data/Pressures_Threats.csv") |>
  rename(Cropland_high = Croplandhigh,
         Cropland_low = Croplandlow,
         Cropland_med = Croplandmid,
         Pasture_low = PastureLow,
         Pasture_high = PastureHigh,
         WoodlandForest_primary = G4Mconservation,
         WoodlandForest_multi = G4Mmultipurpose,
         WoodlandForest_prod = G4Mproduction) |>
  pivot_longer(-Code) |>
  drop_na(value)

lc_threats_wide <- lc_threats |>
  pivot_wider(values_from = value,
              names_from = name) |>
  rename(pressure_code = Code)

threat_quantiles <- spp_threats |>
  left_join(lc_threats_wide) |>
  pivot_longer(-c(pressure_code, taxon_id)) |>
  mutate(value = replace_na(value, 0)) |>
  separate(name, into = c("type", "intensity"), "_") |>
  filter(intensity == "high" | intensity == "prod") |>
  group_by(type, taxon_id) |>
  summarise(value = sum(value)) |>
  filter(value >0) |>
  group_by(type) |>
  summarise(low_threat = quantile(value, 0.10),
            high_threat = quantile(value, 0.50))

################## Calculate zone contribution to targets by PU #################
################# Restoration zones #########################

filelist_temp <- list.files("data/SpeciesData/BIOCLIMA_PotentialSDMs_Art1217/", recursive = TRUE)
spp <- rast(paste0("data/SpeciesData/BIOCLIMA_PotentialSDMs_Art1217/", filelist_temp))
spp <- stack(spp)
# just species names
names(spp) <- substring(names(spp), 20)

#feature_targets <- read_rds("data/formatted-data/feature_targets.rds")

names_sdms <- names(spp)

names_sdms <- as_tibble(names_sdms)|>
  mutate(sdm = 1) |>
  rename(name = value)

# feature_targets |> left_join(names_sdms) |>
#   drop_na(sdm)


# Grab the PU template (same used for LC analysis)
PU_template <- raster("data/landcover/10km/Corine_2018_cropland.tif") |>
  rasterToPolygons() |>
  st_as_sf() |>
  #st_transform(crs = st_crs(natura)) |>
  dplyr::mutate(PUID = seq(1:length(geometry))) |>
  dplyr::select(-Corine_2018_cropland)

# make it into a raster
PU_raster <- fasterize(PU_template, spp[[1]], field = "PUID")

# Is SPP potential range within PU
pu_features_data <-
    ### add in indices for planning units in raster to be organized
    ### not totally necessary bc we will use PUID to reduce PU
    tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
    ### add in cost data
    mutate(cost = 1) %>%
    ### add in PUID
    bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
    ### add in SPP potential data
    bind_cols(as_tibble(raster::as.data.frame(spp)))

pu_features_data |>
    dplyr::select(-id) |>
    write_fst("data/intermediate-data/pu_features_data_spp.fst")

pu_features_data <- read_fst("data/intermediate-data/pu_features_data_spp.fst")

# tidy data table
feature_data_lookup <- MAES |>
  dplyr::select(speciesname, taxon_id) |>
  unique()

write_csv(feature_data_lookup, "data/formatted-data/feat_table.csv")

# make a list of zones and their land cover

pu_zone_features_data <- pu_features_data |>
  pivot_longer(-c(cost, layer)) |>
  rename(PUID = layer) |>
  drop_na(PUID) |>
  dplyr::select(-cost) |>
  rename(speciesname = name) |>
  left_join(species_id)

zone_id_rest <- zone_id |>
  filter(action == "restore")
zone_id_rest$id


## for each zone
for (i in zone_id_rest$id) {
  # make stack for zone
  zone_features <- pu_zone_features_data |>
    drop_na(value) |> filter(value >0)
  # get zone name
  zonename <- zone_id$zone[i]
  zoneid <- zone_id$id[i]

  # get the habitat of the zone
  habitat <- zone_id |>
    dplyr::select(maes_label, level, hab_pref, action)
  habitat <- habitat[i,] # change to i
  habitat_join <- habitat
  if(habitat$maes_label == "Pasture"){
    habitat_join$maes_label <- "Grassland"
  }

  # set habitat and action to 0 if not habitat
  spp_in_habitat <- MAES |>
    left_join(habitat_join) |> drop_na(hab_pref) |>
    group_by(speciesname) |>
    summarise(hab_pref = 1)

  # pu, zone, spp, value
  zone_features <- zone_features |>
    left_join(spp_in_habitat, by = "speciesname") |>
    mutate(value = replace_na(value, 0),
           hab_pref = replace_na(hab_pref,0)) |>
    mutate(value = value*hab_pref) |>
    dplyr::select(-c(speciesname, hab_pref)) |>
    mutate(zone = zoneid) |>
    mutate(value = as.integer(value),
           PUID = as.integer(PUID),
           zone = as.integer(zone)) |>
    rename(pu = PUID,
           species = taxon_id,
           amount = value) |>
    dplyr::select(pu, species, zone, amount) |>
    filter(amount > 0)

  # nuance with threats data
  zone_threats_i <- paste0(habitat$maes_label, "_", habitat$level)

  lc_threats_i <- lc_threats |> filter(name == zone_threats_i)

  quantiles <- threat_quantiles |> filter(type == habitat$maes_label)

  threat_impact <- spp_threats |> unique() |>
    rename(Code = pressure_code) |>
    left_join(lc_threats_i, by = "Code") |>
    filter(value >0) |>
    group_by(taxon_id) |>
    summarise(value = sum(value)) |>
    rename(species = taxon_id) |>
    mutate(value = replace_na(value, 0)) |>
    mutate(value = ifelse(value > quantiles$high_threat,0,
                          ifelse(value < quantiles$low_threat, 0.6, 0.3)))

  zone_features <- zone_features |> left_join(threat_impact) |>
    mutate(value = replace_na(value, 1)) |>
    mutate(amount = amount*value) |>
    dplyr::select(-value) |>
    filter(amount > 0)

  print(i)
  # write out for a given zone
  write_fst(zone_features,
            paste0("data/intermediate-data/zone-features-spp/", zonename, ".fst" ))
}

################# All zones - potential for targets #########################
## for each zone
for (i in zone_id$id) {
  # make stack for zone
  zone_features <- pu_zone_features_data |>
    drop_na(value) |> filter(value >0)
  # get zone name
  zonename <- zone_id$zone[i]
  zoneid <- zone_id$id[i]

  # get the habitat of the zone
  habitat <- zone_id |>
    dplyr::select(maes_label, level, hab_pref, action)
  habitat <- habitat[i,] # change to i
  habitat_join <- habitat
  if(habitat$maes_label == "Pasture"){
    habitat_join$maes_label <- "Grassland"
  }

  # set habitat and action to 0 if not habitat
  spp_in_habitat <- MAES |>
    left_join(habitat_join) |> drop_na(hab_pref) |>
    group_by(speciesname) |>
    summarise(hab_pref = 1)

  # pu, zone, spp, value
  zone_features <- zone_features |>
    left_join(spp_in_habitat, by = "speciesname") |>
    mutate(value = replace_na(value, 0),
           hab_pref = replace_na(hab_pref,0)) |>
    mutate(value = value*hab_pref) |>
    dplyr::select(-c(speciesname, hab_pref)) |>
    mutate(zone = zoneid) |>
    mutate(value = as.integer(value),
           PUID = as.integer(PUID),
           zone = as.integer(zone)) |>
    rename(pu = PUID,
           species = taxon_id,
           amount = value) |>
    dplyr::select(pu, species, zone, amount) |>
    filter(amount > 0)

  # nuance with threats data
  zone_threats_i <- paste0(habitat$maes_label, "_", habitat$level)

  lc_threats_i <- lc_threats |> filter(name == zone_threats_i)

  quantiles <- threat_quantiles |> filter(type == habitat$maes_label)

  threat_impact <- spp_threats |> unique() |>
    rename(Code = pressure_code) |>
    left_join(lc_threats_i, by = "Code") |>
    filter(value >0) |>
    group_by(taxon_id) |>
    summarise(value = sum(value)) |>
    rename(species = taxon_id) |>
    mutate(value = replace_na(value, 0)) |>
    mutate(value = ifelse(value > quantiles$high_threat,0,
                          ifelse(value < quantiles$low_threat, 0.6, 0.3)))

  zone_features <- zone_features |> left_join(threat_impact) |>
    mutate(value = replace_na(value, 1)) |>
    mutate(amount = amount*value) |>
    dplyr::select(-value) |>
    filter(amount > 0)

  print(i)
  # write out for a given zone
  write_fst(zone_features,
            paste0("data/intermediate-data/zone-features-spp-targets/", zonename, ".fst" ))
}

################### conservation and production zones #############

filelist_temp <- list.files("data/SpeciesData/CurrentSDMs//",
                            pattern = "*.tif$",
                            recursive = TRUE)
spp <- rast(paste0("data/SpeciesData/CurrentSDMs/", filelist_temp))
spp <- stack(spp)
# just species names

# Function to split the string and keep the part after "__"
split_and_keep_after <- function(x) {
  split_string <- strsplit(x, split = "__")[[1]]
  return(substr(split_string[2], 1, nchar(split_string[2]) - 4))
}

# Apply the function to each element in the list
names_sdms <- lapply(filelist_temp, split_and_keep_after)
names_sdms <- unlist(names_sdms)

names(spp) <- names_sdms

#feature_targets <- read_rds("data/formatted-data/feature_targets.rds")

names_sdms <- as_tibble(names_sdms)|>
  mutate(sdm = 1) |>
  rename(name = value)

# make it into a raster
PU_raster <- fasterize(PU_template, spp[[1]], field = "PUID")

# Is SPP current range within PU
pu_features_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in SPP potential data
  bind_cols(as_tibble(raster::as.data.frame(spp)))

pu_features_data |>
  dplyr::select(-id) |>
  write_fst("data/intermediate-data/pu_features_data_spp_current.fst")

pu_features_data <- read_fst("data/intermediate-data/pu_features_data_spp_current.fst")

# tidy data table
feature_data_lookup <- MAES |>
  dplyr::select(speciesname, taxon_id) |>
  unique()

write_csv(feature_data_lookup, "data/formatted-data/feat_table.csv")

# make a list of zones and their land cover

pu_zone_features_data <- pu_features_data |>
  pivot_longer(-c(cost, layer)) |>
  rename(PUID = layer) |>
  drop_na(PUID) |>
  dplyr::select(-cost) |>
  rename(speciesname = name) |>
  left_join(species_id)

zone_id_cons <- zone_id |>
  filter(action != "restore")
zone_id_cons$id

## for each zone
for (i in zone_id_cons$id) {
  # make stack for zone
  zone_features <- pu_zone_features_data |>
    drop_na(value) |> filter(value >0)
  # get zone name
  zonename <- zone_id$zone[i]
  zoneid <- zone_id$id[i]

  # get the habitat of the zone
  habitat <- zone_id |>
    dplyr::select(maes_label, level, hab_pref, action)
  habitat <- habitat[i,] # change to i
  habitat_join <- habitat
  if(habitat$maes_label == "Pasture"){
    habitat_join$maes_label <- "Grassland"
  }

  # set habitat and action to 0 if not habitat
  spp_in_habitat <- MAES |>
    left_join(habitat_join) |> drop_na(hab_pref) |>
    group_by(speciesname) |>
    summarise(hab_pref = 1)

  # pu, zone, spp, value
  zone_features <- zone_features |>
    left_join(spp_in_habitat, by = "speciesname") |>
    mutate(value = replace_na(value, 0),
           hab_pref = replace_na(hab_pref,0)) |>
    mutate(value = value*hab_pref) |>
    dplyr::select(-c(speciesname, hab_pref)) |>
    mutate(zone = zoneid) |>
    mutate(value = as.integer(value),
           PUID = as.integer(PUID),
           zone = as.integer(zone)) |>
    rename(pu = PUID,
           species = taxon_id,
           amount = value) |>
    dplyr::select(pu, species, zone, amount) |>
    filter(amount > 0)

  # nuance with threats data
  zone_threats_i <- paste0(habitat$maes_label, "_", habitat$level)

  lc_threats_i <- lc_threats |> filter(name == zone_threats_i)

  quantiles <- threat_quantiles |> filter(type == habitat$maes_label)

  threat_impact <- spp_threats |> unique() |>
    rename(Code = pressure_code) |>
    left_join(lc_threats_i, by = "Code") |>
    filter(value >0) |>
    group_by(taxon_id) |>
    summarise(value = sum(value)) |>
    rename(species = taxon_id) |>
    mutate(value = replace_na(value, 0)) |>
    mutate(value = ifelse(value > quantiles$high_threat,0,
                          ifelse(value < quantiles$low_threat, 0.6, 0.3)))

  zone_features <- zone_features |> left_join(threat_impact) |>
    mutate(value = replace_na(value, 1)) |>
    mutate(amount = amount*value) |>
    dplyr::select(-value) |>
    filter(amount > 0)
  print(i)
  # write out for a given zone
  write_fst(zone_features,
            paste0("data/intermediate-data/zone-features-spp/", zonename, ".fst" ))
}


### Zone features carbon #########
### pu, species, zone, amount
# current carbon

## production
## restoration (potential)
## conservation

filelist_temp <- list.files("data/NatureMap_currentCarbon/")
carbon_NatureMap <- rast(paste0("data/NatureMap_currentCarbon/", filelist_temp))
carbon_NatureMap <- stack(carbon_NatureMap)
# just carbon names
names(carbon_NatureMap) <- substring(names(carbon_NatureMap), 11)

# new carbon data
filelist_temp <- list.files("data/CurrentCarbonDensity_Corine_/")
carbon_current <- rast(paste0("data/CurrentCarbonDensity_Corine_/", filelist_temp))
carbon_current <- stack(carbon_current)

names(carbon_current) <- substring(filelist_temp, 40)


carbon_currentSOC <- resample(carbon_NatureMap[[3]], carbon_current[[1]])

names(carbon_currentSOC) <-"carbonSOC_10km"

carbon_current <- stack(carbon_current, carbon_currentSOC)
plot(carbon_current)
PU_raster <- fasterize(PU_template, carbon_current[[1]], field = "PUID")

pu_carbon_current_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in carbon potential data
  bind_cols(as_tibble(raster::as.data.frame(carbon_current))) |>
  drop_na(layer)


z <- read_csv("data/formatted-data/zone_id.csv") |>
  mutate(name = paste0("z", id)) |>
  relocate(id, 1)

glimpse(pu_carbon_current_data)

pu_carbon_conservation <- pu_carbon_current_data |>
  mutate(
    Grassland = NaturalGrassland_laea_tCha.tif + carbonSOC_10km,
    HeathlandShrub = HeathlandShrub_laea_tCha.tif + carbonSOC_10km,
    MarineTransitional = MarineTransitional_laea_tCha.tif + carbonSOC_10km,
    SparseVeg = SparseVeg_laea_tCha.tif + carbonSOC_10km,
    Wetlands = Wetlands_laea_tCha.tif + carbonSOC_10km,
    WoodlandForest = WoodlandForest_laea_tCha.tif + carbonSOC_10km
  ) |>
  dplyr::select(layer, Grassland:WoodlandForest) |>
  pivot_longer(-layer) |>
  mutate(carbon = replace_na(value, 0)) |>
  mutate(action = "conserve") |>
  rename(maes_label = name) |>
  left_join(zone_id, multiple = "all") |>
  dplyr::select(-c(maes_label, zone,  hab_pref, action)) |>
  rename(pu = layer, amount = carbon, zone = id) |>
  mutate(feature = 999999) |>
  dplyr::select(pu, feature, zone, amount) |>
  mutate(amount = amount) |>
  mutate(amount = ifelse(amount <0.001, 0,amount))

ggplot(pu_carbon_conservation) + geom_histogram(aes(x = amount)) +
  facet_wrap(~zone)


# potential carbon
filelist_temp <- list.files("data/PotentialCarbonDistribution/")
carbon_potential <- rast(paste0("data/PotentialCarbonDistribution/", filelist_temp))
carbon_potential <- stack(carbon_potential)
# just species names
names(carbon_potential) <- substring(names(carbon_potential), 37)
carbon_currentSOC <- crop(carbon_NatureMap[[3]], carbon_potential[[1]])
carbon_potential <- stack(carbon_potential, carbon_currentSOC, carbon_currentSOC)

PU_raster <- fasterize(PU_template, carbon_potential[[1]], field = "PUID")

pu_carbon_potential_data <-
  ### add in indices for planning units in raster to be organized
  ### not totally necessary bc we will use PUID to reduce PU
  tibble(id = as.list(seq_len(ncell(PU_raster)))) %>%
  ### add in cost data
  mutate(cost = 1) %>%
  ### add in PUID
  bind_cols(as_tibble(raster::as.data.frame(PU_raster))) |>
  ### add in carbon potential data
  bind_cols(as_tibble(raster::as.data.frame(carbon_potential))) |>
  drop_na(layer)

z <- read_csv("data/formatted-data/zone_id.csv") |>
  mutate(name = paste0("z", id)) |>
  relocate(id, 1)

pu_carbon_restoration <- pu_carbon_potential_data |>
  rename(
    Grassland = grassland__tCha,
    HeathlandShrub = heathland.and.shrub__tCha,
    MarineTransitional = marine.inlets.and.transitional.waters__tCha,
    SparseVeg = sparsely.vegetated.areas__tCha,
    Wetlands = wetlands__tCha,
    WoodlandForest = woodland.and.forest__tCha,
    Cropland = carbonSOC_10km.1,
    Pasture = carbonSOC_10km.2
  ) |>
  dplyr::select(layer, Grassland:Pasture) |>
  pivot_longer(-layer) |>
  mutate(carbon = replace_na(value, 0)) |>
  mutate(action = "restore") |>
  rename(maes_label = name) |>
  left_join(zone_id, multiple = "all") |>
  dplyr::select(-c(maes_label, zone,  hab_pref, action)) |>
  rename(pu = layer, amount = carbon, zone = id) |>
  mutate(feature = 999999) |>
  dplyr::select(pu, feature, zone, amount) |>
  mutate(amount = amount) |>
  mutate(amount = ifelse(amount <0.001, 0,amount)) |>
  mutate(amount = ifelse(zone == 2, amount*0.95,
                         ifelse(zone == 9, amount*0.7, amount)))

ggplot(pu_carbon_restoration) + geom_histogram(aes(x = amount)) +
  facet_wrap(~zone)

# get production carbon values
pu_carbon_current_production <- pu_carbon_current_data |>
  mutate(carbon_soc = carbonSOC_10km) |>
  dplyr::select(layer, carbon_soc, WoodlandForest_laea_tCha.tif) |>
  mutate(carbon_soc = replace_na(carbon_soc, 0)) |>
  mutate(action = "production")

### pu, species, zone, amount
pu_carbon_current_production <- pu_carbon_current_production |>
  left_join(zone_id) |>
  dplyr::select(-c(maes_label, zone, hab_pref, action)) |>
  mutate(carbon = ifelse(id == 18, carbon_soc*0.93, # crop high
                         ifelse(id == 17 |id == 15, carbon_soc, # low intensity crop and pasture
                                ifelse(id == 14, carbon_soc* 0.76, #high intensity pasture
                                       ifelse(id == 16, carbon_soc*0.95, #mid crop
                                            ifelse(id == 12, (carbon_soc + WoodlandForest_laea_tCha.tif)*0.7,# forest multi
                                                   ifelse(id == 13, (carbon_soc + WoodlandForest_laea_tCha.tif)*0.4, # forest production
                                                          0))))))) |>
  rename(pu = layer, amount = carbon, zone = id) |>
  mutate(feature = 999999) |>
  dplyr::select(pu, feature, zone, amount)

ggplot(pu_carbon_current_production) + geom_histogram(aes(x = amount)) +
  facet_wrap(~zone)


carbon_pu_data <-   bind_rows(pu_carbon_conservation,
                              pu_carbon_restoration,
                              pu_carbon_current_production) |>
  mutate(amount = amount/100)

carbon_pu_data |> ggplot() + geom_histogram(aes(x = amount)) +
  facet_wrap(~zone)

####################### Bringing data together#####################
filelist_temp <- list.files("data/intermediate-data/zone-features-spp/")
zone_features <- list()
for (i in 1:length(filelist_temp)) {
  zone_features[[i]] <- read_fst(paste0("data/intermediate-data/zone-features-spp/", filelist_temp[i]))
}

zone_features_spp <- zone_features |> bind_rows() |>
  rename(feature = species)


zone_features_spp |> bind_rows(carbon_pu_data) |>
  write_fst("data/formatted-data/features_data_all.fst", compress = 75)

####################### For calculating targets #####################
filelist_temp <- list.files("data/intermediate-data/zone-features-spp-targets/")
zone_features <- list()
for (i in 1:length(filelist_temp)) {
  zone_features[[i]] <- read_fst(paste0("data/intermediate-data/zone-features-spp-targets/", filelist_temp[i]))
}

zone_features_spp <- zone_features |> bind_rows() |>
  rename(feature = species)

zone_features_spp |>
  write_fst("data/formatted-data/features_data_spp_fortargets.fst", compress = 75)

