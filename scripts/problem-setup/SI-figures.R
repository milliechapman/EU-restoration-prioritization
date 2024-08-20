library(sf)
library(terra)
library(tidyverse)
library(fst)
library(raster)
library(fasterize)


### shortfall for SI solutions
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

## Calculate representation and save

filelist_temp <- list.files("data/solutions/sol/SI-sol/",
                           # pattern = c("*proportional_gurobi_f455.csv|*proportional_gurobi_ref.csv"),
                            full.names = T)
shortfalls <- list()

#  join targets and spp info
for (i in 1:length(filelist_temp)) {
  filenm <- unlist(strsplit(filelist_temp[[i]], "_"))
  carbon_weight <- filenm[[3]]
  country_contraint <- filenm[[9]]
  scenario <- filenm[[15]]
  future <- str_sub(filenm[[18]], 1,-5)
  name = filenm[[16]]

  solution <- read_csv(filelist_temp[[i]])

  solution <- solution |>
    rowid_to_column("pu") |>
    dplyr::select(pu, solution_1_z1:solution_1_z26) |>
    pivot_longer(-pu) |>
    mutate(zone = str_sub(name, 13)) |>
    dplyr::select(-name) |>
    mutate(zone = as.numeric(zone))

  if(name == "onlyRest"){
    rij_rel = rij
  }

  if(name == "proportional"){
    rij_rel = rij
  }

  if(name == "carbonSI"){
    rij_rel = rij |>
      mutate(amount = ifelse(zone %in% c(16, 17,18,1,2) & species == 999999,
                             1, amount))
  }
  if(name == "biodivLOW"){
    rij_rel = rij |>
      mutate(amount = ifelse(amount == 0.3 & species != 999999, 0.1,
                             ifelse(amount == 0.6 & species != 999999, 0.3, amount)))
  }
  if(name == "biodivHIGH"){
    rij_rel = rij |>
      mutate(amount = ifelse(amount == 0.3 & species != 999999, 0.5,
                             ifelse(amount == 0.6 & species != 999999, 0.8, amount)))
  }

  rij_sol <- rij_rel |>
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
           future = future,
           name = name) |>
    mutate(target_met = ifelse(shortfall<0,"no", "yes")) |>
    mutate(spp = ifelse(species != 999999, "biodiv", "carbon"))

  shortfalls[[i]] <- contribution
}


shortfalls <- bind_rows(shortfalls, .id = "column_label")

write_csv(shortfalls, "data/solutions/representation_SI_scenarios_proportional.csv")


### IC rep scenarios

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

name = "carbonSI"

if(name == "carbonSI"){
  rij_rel = rij |>
    mutate(amount = ifelse(zone %in% c(16, 17,18,1,2) & species == 999999,
                           1, amount))
  rij_sol <- rij_rel |>
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

  write_csv(contribution, "data/solutions/representation_IC_cSI.csv")
}

name = "biodivLOW"

if(name == "biodivLOW"){
  rij_rel = rij |>
    mutate(amount = ifelse(amount == 0.3, 0.1,
                           ifelse(amount == 0.6, 0.3, amount)))
  rij_sol <- rij_rel |>
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

  write_csv(contribution, "data/solutions/representation_IC_bLOWSI.csv")
}

name = "biodivHIGH"

if(name == "biodivHIGH"){
  rij_rel = rij |>
    mutate(amount = ifelse(amount == 0.3, 0.5,
                           ifelse(amount == 0.6, 0.8, amount)))

  rij_sol <- rij_rel |>
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

  write_csv(contribution, "data/solutions/representation_IC_HIGHSI.csv")
}


####################
### grab other relevant data
zone_id <- read_csv("data/formatted-data/zone_id.csv")
spp_lookup_count <- read_csv("data/plotting/spp_lookup.csv") |>
  filter(typeasso %in% c("Preferred", "Suitable"))|>
  group_by(speciesname) |> count() |>
  ungroup() |>
  mutate(n = ifelse(n>0,1,0))

spp_lookup <- read_csv("data/plotting/spp_lookup.csv") |>
  dplyr::select(speciesname, speciesgroup, taxon_id) |>
  left_join(spp_lookup_count) |> unique() #|>
#pivot_wider(names_from = maes_label, values_from = n)

BR_lookup <- read_csv("data/plotting/bioregion_code.csv") |>
  rename(bioregion = BRIDnum) |>
  dplyr::select(code, bioregion)
country_lookup <- read_csv("data/plotting/country_code.csv")
# plotting_data <- read_csv("data/plotting/plotting_data_ref_f455.csv")

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

country <- read_csv("data/formatted-data/linear_constraints/pu_pasture_low_budget_data.csv") |>
  dplyr::select(NUTS_ID, pu) |>
  mutate(country = substr(NUTS_ID, start = 1, stop = 2))|>
  rename(id = pu) |> dplyr::select(-NUTS_ID)

PU_template_EU <- pu_in_EU |>
  rename(id = pu) |> #plot_data |>
  left_join(PU_template) |>
  dplyr::select(-id) |>
  rename(id = EU_id)

colors <- c("grey", met.brewer(name="Isfahan1",type="continuous"))

#### make rep layers

rep_IC <- read_csv("data/solutions/representation_IC.csv")
rep_IC_bLOW <- read_csv("data/solutions/representation_IC_bLOWSI.csv")
rep_IC_bHIGH <- read_csv("data/solutions/representation_IC_cSI.csv")
rep_IC_carbon <- read_csv("data/solutions/representation_IC.csv")
rep <- read_csv("data/solutions/representation_SI_scenarios_proportional.csv")

ggplot(unique(rep_IC)) + geom_histogram(aes(x = rep))

rep_formatting <- function(x) {
  x |>
    #dplyr::select(-c(rep, target, shortfall)) |>
    #filter(species != 999999) |>
    mutate(taxon_id = ifelse(species > 10000000000,
                             substr(species, start = 1, stop = 7),
                             substr(species, start = 1, stop = 6)),
           bioregion = ifelse(species > 10000000000,
                              substr(species, start = 8, stop = 9),
                              substr(species, start = 7, stop = 8)),
           country_code = ifelse(species > 10000000000,
                                 substr(species, start = 10, stop = 11),
                                 substr(species, start = 9, stop = 10))) |>
    mutate(taxon_id = as.numeric(taxon_id),
           bioregion = as.numeric(bioregion),
           country_code = as.numeric(country_code)) |>
    left_join(BR_lookup, multiple = "all") |>
    left_join(spp_lookup, multiple = "all")
}

rep_IC_proportional <- rep_formatting(rep_IC) |>
  mutate(name = "proportional") |>
  rename(rep_ic = rep,
         shortfall_ic = shortfall,
         shortfall_perc_ic = shortfall_perc) |>
  dplyr::select(species, rep_ic, shortfall_ic, shortfall_perc_ic, name)

rep_IC_bLOW <- rep_formatting(rep_IC_bLOW) |>
  mutate(name = "biodivLOW") |>
  rename(rep_ic = rep,
         shortfall_ic = shortfall,
         shortfall_perc_ic = shortfall_perc) |>
  dplyr::select(species, rep_ic, shortfall_ic, shortfall_perc_ic, name)

rep_IC_bHIGH <- rep_formatting(rep_IC_bHIGH) |>
  mutate(name = "biodivHIGH")|>
  rename(rep_ic = rep,
         shortfall_ic = shortfall,
         shortfall_perc_ic = shortfall_perc) |>
  dplyr::select(species, rep_ic, shortfall_ic, shortfall_perc_ic, name)

rep_IC_carbon <- rep_formatting(rep_IC_carbon) |>
  mutate(name = "carbonSI")|>
  rename(rep_ic = rep,
         shortfall_ic = shortfall,
         shortfall_perc_ic = shortfall_perc) |>
  dplyr::select(species, rep_ic, shortfall_ic, shortfall_perc_ic, name)

rep_IC <- rep_IC_carbon |>
  bind_rows(rep_IC_bHIGH) |>
  bind_rows(rep_IC_bLOW) |>
  bind_rows(rep_IC_proportional)


rep_formatted <- rep |>
  dplyr::select(-c(column_label, rep, target, shortfall)) |>
  #filter(species != 999999) |>
  mutate(taxon_id = ifelse(species > 10000000000,
                           substr(species, start = 1, stop = 7),
                           substr(species, start = 1, stop = 6)),
         bioregion = ifelse(species > 10000000000,
                            substr(species, start = 8, stop = 9),
                            substr(species, start = 7, stop = 8)),
         country_code = ifelse(species > 10000000000,
                               substr(species, start = 10, stop = 11),
                               substr(species, start = 9, stop = 10))) |>
  mutate(taxon_id = as.numeric(taxon_id),
         bioregion = as.numeric(bioregion),
         country_code = as.numeric(country_code)) |>
  left_join(BR_lookup, multiple = "all") |>
  left_join(spp_lookup, multiple = "all")

##################### biodiversity sensitivity ##############################

# ic_join <- rep_formatted |>
#   left_join(rep_IC) |>
#   rename(#rep_ic = rep,
#          shortfall_ic = shortfall,
#          shortfall_perc_ic = shortfall_perc) |>
#   dplyr::select(species, rep_ic, shortfall_ic, shortfall_perc_ic)

## conservation status
b_status <- rep |>
  left_join(rep_IC) |>
  filter(species != 999999) |>
  mutate(target_met_ic = ifelse(shortfall_perc_ic < 0, 0, 1),
         target_met = ifelse(shortfall_perc < 0, 0, 1)) |>
  filter(target_met_ic<1) |>
  mutate(improved = (target_met-target_met_ic)) |>
  group_by(spp, scenario, carbon_weight, country_contraint, future, name) |>
  summarise(improved = sum(improved)/length(unique(species))) |>
  arrange(-improved)

b_status |> filter(carbon_weight == 0.5)

a <- b_status |> filter(carbon_weight>0) |>
  ggplot(aes(x = carbon_weight, y = improved*100,
             col = name)) +
  geom_line(aes(linetype = name)) +
  geom_point(aes()) +
  theme_classic() +
  #scale_color_manual(values = c("grey2", "#0072B2", "orange")) +
  scale_alpha_manual(values = c(1, 0.6)) +
  labs(x = "carbon weight \n (relative to biodiversity)",
       y = "% of species with improved \n conservation status") +
  facet_wrap(~future)

## carbon percentage change
c_perc <- rep |>
  filter(species == 999999) |>
  left_join(rep_IC) |>
  filter(country_contraint == "EVEN") |>
  filter(spp == "carbon") |>
  filter(rep_ic > 0) |>
  mutate(rep_diff = (rep - rep_ic)/rep_ic) |>
  #filter(name == "carbonSI") |>
  group_by(spp, scenario, carbon_weight, country_contraint, future, name) |>
  summarise(rep_diff = rep_diff*100) |>
  arrange(-rep_diff)

b <- c_perc |> filter(carbon_weight>0) |>
  ggplot(aes(x = carbon_weight, y = rep_diff,
             col = name)) +
  geom_line(aes(linetype = name)) +
  geom_point(aes()) +
  theme_classic() +
  #scale_color_manual(values = c("grey2", "#0072B2", "orange")) +
  labs(x = "carbon weight \n (relative to biodiversity)",
       y = "% increase in land carbon stock \n (compared to initial conditions)") +
  facet_wrap(~future)


################## Carbon Sensitivity ######################

rep_IC <- read_csv("data/solutions/representation_IC_cSI.csv")

rep_IC_formatted <- rep_IC |>
  dplyr::select(-c(rep, target, shortfall)) |>
  filter(species != 999999) |>
  mutate(taxon_id = ifelse(species > 10000000000,
                           substr(species, start = 1, stop = 7),
                           substr(species, start = 1, stop = 6)),
         bioregion = ifelse(species > 10000000000,
                            substr(species, start = 8, stop = 9),
                            substr(species, start = 7, stop = 8)),
         country_code = ifelse(species > 10000000000,
                               substr(species, start = 10, stop = 11),
                               substr(species, start = 9, stop = 10))) |>
  mutate(taxon_id = as.numeric(taxon_id),
         bioregion = as.numeric(bioregion),
         country_code = as.numeric(country_code)) |>
  left_join(BR_lookup, multiple = "all") |>
  left_join(spp_lookup, multiple = "all")

ic_join <- rep_IC |>
  rename(rep_ic = rep,
         shortfall_ic = shortfall,
         shortfall_perc_ic = shortfall_perc) |>
  dplyr::select(species, rep_ic, shortfall_ic, shortfall_perc_ic)

## carbon percentage change
c_perc <- rep |>
  left_join(ic_join) |>
  filter(spp == "carbon") |>
  filter(rep_ic > 0) |>
  mutate(rep_diff = (rep - rep_ic)/rep_ic) |>
  filter(name == "carbonSI") |>
  group_by(spp, scenario, carbon_weight, country_contraint, future, name) |>
  summarise(rep_diff = rep_diff*100) |>
  arrange(-rep_diff)

b <- c_perc |> filter(carbon_weight>0) |>
  ggplot(aes(x = carbon_weight, y = rep_diff,
             col = name)) +
  geom_line(aes(linetype = name)) +
  geom_point(aes()) +
  theme_classic() +
  #scale_color_manual(values = c("grey2", "#0072B2", "orange")) +
  labs(x = "carbon weight \n (relative to biodiversity)",
       y = "% increase in land carbon stock \n (compared to initial conditions)") +
  facet_wrap(~future)

