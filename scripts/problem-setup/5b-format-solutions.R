library(sf)
library(terra)
library(tidyverse)
library(fst)
library(raster)
library(fasterize)
#rm(list = ls())

targets <- read_csv("data/formatted-data/targets_split_formatted.csv")
pu_in_EU <- read_csv("data/formatted-data/pu_in_EU.csv")
raster <- raster("data/landcover/10km/Corine_2018_cropland.tif")

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
  mutate(amount = ifelse(amount<0.001, 0, amount)) |>
  mutate(amount = replace_na(amount, 0))

rij |> filter(species == 999999) |> ggplot(aes(x = amount)) + geom_histogram() + facet_wrap(~zone)
## Calculate representation and save

filelist_temp <- list.files("data/solutions/sol/") #, pattern = c("*_ref.csv|*_f455.csv"))
shortfalls <- list()
#  join targets and spp info
for (i in 1:length(filelist_temp)) {
  filenm <- unlist(strsplit(filelist_temp[[i]], "_"))
  carbon_weight <- filenm[[3]]
  country_contraint <- filenm[[9]]
  scenario <- filenm[[15]]
  future <- str_sub(filenm[[18]], 1,-5)
  
  solution <- read_csv(paste0("data/solutions/sol/", filelist_temp[[i]]))
  
  solution <- solution |>
    rowid_to_column("pu") |>
    dplyr::select(pu, solution_1_z1:solution_1_z26) |>
    pivot_longer(-pu) |>
    mutate(zone = str_sub(name, 13)) |>
    dplyr::select(-name) |>
    mutate(zone = as.numeric(zone))
  
  rij_sol <- rij |>
    left_join(solution, by = c("pu", "zone"))
  
  contribution <- rij_sol |>
    mutate(contribution = amount*value) |>
    group_by(species) |>
    summarise(rep = sum(contribution))|> 
    left_join(targets) |>
    mutate(shortfall = rep-target) |>
    mutate(shortfall_perc = shortfall/target) |>
    mutate(carbon_weight = carbon_weight,
           scenario = scenario,
           country_contraint = country_contraint,
           future = future) |>
    mutate(target_met = ifelse(shortfall<0,"no", "yes")) |>
    mutate(spp = ifelse(species != 999999, "biodiv", "carbon"))
  
  shortfalls[[i]] <- contribution
}


shortfalls <- bind_rows(shortfalls, .id = "column_label") 
write_csv(shortfalls, "data/solutions/representation_REF_f455_scenarios_version3.csv")

glimpse(shortfalls)

zones <- read_csv("data/formatted-data/zone_id.csv") |>
  filter(id >11) |>
  separate(zone, sep = "_", into = c("A", "B", "C")) |>
  mutate(zone = paste0(A, "_", B)) |>
  dplyr::select(-c(A,B,C))


LC <- read_csv("data/outputs/2-zones/PU_lc_intensity.csv") |> 
  dplyr::select(-Status) |>  
  pivot_longer(-PUID) |> rename(pu = PUID) |>
  left_join(pu_in_EU) |>
  rename(id = EU_id) |>
  dplyr::select(-c(pu, nuts2id)) |>
  drop_na(id) |> rename(zone = name) |>
  rename(pu = id) |>
  left_join(zones) |> 
  dplyr::select(-zone) |>
  rename(zone = id)

rij_sol <- rij |>
  left_join(LC, by = c("pu", "zone"))

contribution <- rij_sol |>
  mutate(contribution = amount*value) |>
  group_by(species) |>
  summarise(rep = sum(contribution, na.rm = TRUE))|> 
  left_join(targets) |>
  mutate(shortfall = rep-target) |>
  mutate(shortfall_perc = shortfall/target) |>
  mutate(carbon_weight = "IC",
         scenario = "IC",
         country_contraint = "IC") |>
  mutate(target_met = ifelse(shortfall<0,"no", "yes")) |>
  mutate(spp = ifelse(species != 999999, "biodiv", "carbon"))

write_csv(contribution, "data/solutions/representation_IC.csv")




## Save as raster stack
## 
## 
## 
#### plotting data

zone_id <- read_csv("data/formatted-data/zone_id.csv") 
zone_id$zone

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

PU_template <- raster("data/landcover/10km/Corine_2018_cropland.tif") |>
  rasterToPolygons() |>
  st_as_sf() |>
  #st_transform(crs = st_crs(natura)) |>
  dplyr::mutate(id = seq(1:length(geometry))) |>
  dplyr::select(-Corine_2018_cropland) 
raster <- raster("data/landcover/10km/Corine_2018_cropland.tif")

pu_in_EU <- read_csv("data/formatted-data/pu_in_EU.csv")

PU_plot <- pu_in_EU |> rename(id = pu) |> #plot_data |> 
  left_join(PU_template) |>
  dplyr::select(EU_id, geometry) |> st_as_sf()

PU_plot <- fasterize(st_as_sf(PU_plot), raster, field = "EU_id") 

bioregion <-read_csv("data/formatted-data/pu_in_EU_BR.csv") |>
  dplyr::select(pu, BR_ID) |>
  rename(id = pu)

country <- read_csv("data/formatted-data/linear_constraints/pu_crop_high_budget_data_55.csv") |>
  dplyr::select(NUTS_ID, pu) |>
  mutate(country = substr(NUTS_ID, start = 1, stop = 2))|>
  rename(id = pu) |> dplyr::select(-NUTS_ID)

PU_template_EU <- pu_in_EU |> 
  rename(id = pu) |> #plot_data |> 
  left_join(PU_template) |> 
  dplyr::select(-id) |>
  rename(id = EU_id)

filelist_temp <- list.files("data/solutions/sol/") #, pattern = c("*_ref.csv|*_f455.csv"))
solution_scenarios <- list()
zones <- read_csv("data/formatted-data/zone_id.csv") |>
  mutate(name = paste0("z", id)) 
#  join targets and spp info
#  
for (i in 1:length(filelist_temp)) {
  name <- substr(filelist_temp[[i]],1,nchar(filelist_temp[[i]])-4)
  solution <- read_csv(paste0("data/solutions/sol/", filelist_temp[[i]])) |>
    dplyr::select(id, solution_1_z1:solution_1_z26)
  
  colnames(solution) <- c("id", (zones$zone))
  
  solution_table <- pu_in_EU |> rename(id = pu) |> #plot_data |> 
    left_join(PU_template) |> 
    dplyr::select(-id) |>
    rename(id = EU_id) |>
    left_join(solution) |> 
    pivot_longer(-c(nuts2id:geometry)) 
  
  solution_table_plot <- solution_table |> 
    mutate(zone = name) |>
    separate(name, c('maes_label', 'action', 'place'), sep = "_")
  
  solution_raster <- fasterize(st_as_sf(solution_table_plot), raster, field = "value", by = "zone")
  solution_raster <- rast(solution_raster)
  
  writeRaster(solution_raster, paste0("data/solutions/rasters/", name, ".tif"), overwrite = TRUE)
}

plot(solution_raster)


# bioregions
bioregion <- st_read("data/EU_BiogeoRegions2016_shapefile/BiogeoRegions2016.shp") |>
  rename(BR_ID = short_name) |>
  mutate(BRIDnum = PK_UID)|>
  as_tibble() |> dplyr::select(-geometry)
write_csv(bioregion, "data/plotting/bioregion_code.csv")

#BR <- raster("data/formatted-data/BR_raster.tif", overwrite = TRUE)

# countries
nuts2_shp <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  mutate(country = substr(NUTS2, start = 1, stop = 2)) |>
  dplyr::select(country) 

country_code <- as_tibble(nuts2_shp) |>
  dplyr::select(-geometry) |> unique() |>
  mutate(country_code = seq(1:28))
write_csv(country_code, "data/plotting/country_code.csv")


species_id <- readRDS("data/formatted-data/MinSpeciesToCover.rds") |>
  dplyr::select(speciesname, taxon_id) |>
  rename(feature = taxon_id)
species_id$speciesname <- gsub(" ", "_", species_id$speciesname, fixed=TRUE)

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

write_csv(MAES, "data/plotting/spp_lookup.csv")



# Read  in solutions

plot_table <- list()
files <- list.files(path = "data/solutions/sol/")
for(i in 1:length(files)){
  nm <- files[i] #"sol_carbon_0.1_restoration_0.2_production_TRUE_country_TRUE_bioregion_TRUE_wetlands_TRUE_onlyrestoration_FALSE_TRIAL"
  carbon <- str_split(nm, "_")[[1]][3]
  country_TF <- str_split(nm, "_")[[1]][9]
  scenario <- str_split(nm, "_")[[1]][15]
  future <- str_split(nm, "_")[[1]][18]

  solution <- read_csv(paste0("data/solutions/sol/", nm)) |>
    dplyr::select(id, solution_1_z1:solution_1_z26)
  colnames(solution) <- c("id", (zone_id$zone))
  
  solution_table <- pu_in_EU |> 
    rename(id = pu) |> #plot_data |> 
    #left_join(PU_template) |> 
    #left_join(bioregion) |>
    dplyr::select(-id) |>
    rename(id = EU_id) |>
    left_join(solution) |> 
    pivot_longer(-c(nuts2id:id)) 
  
  solution_table_plot <- solution_table |> 
    mutate(zone = name) |>
    separate(name, c('maes_label', 'intensity', 'action'), sep = "_") |>
    mutate(carbon = carbon,
           country_TF = country_TF,
           scenario = scenario,
           future = future)
  
  plot_table[[i]] <- solution_table_plot
}

plotting_data <- bind_rows(plot_table)
write_csv(plotting_data, "data/plotting/plotting_data_ref_f455.csv")

