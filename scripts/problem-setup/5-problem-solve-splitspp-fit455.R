# formulate problem with all the pu and feature data
# work in progress!
library(prioritizr) # remotes::install_github("cboettig/prioritizr")
library(fst)
#library(cplexAPI)
library(rcbc)
#library(highs)
# library(gurobi)
library(tidyverse)
library(sf)
# rm(list = ls())

fs::dir_create("data/sol")
cores <- parallel::detectCores()
fst::threads_fst(nr_of_threads = cores, reset_after_fork = NULL)  # Set number for fst

########## function to set up, solve, and save solution #################

problem_setup <- function(
    cores = parallel::detectCores(),
    p, # conservation problem + min shortfall objective/ speeds things up to do this first...
    pu, # pu data (x from p)
    targs,
    manual_bounded_constraints = manual_bounded_constraints,
    carbon_weight = 0.5,
    restoration_constraint = 0.141,
    #conservation_constraint = 0.3,
    production_constraints = TRUE,
    country_constraints = "EVEN", # EVEN, FLEX, UNEVEN
    restoration_only = FALSE,
    wetlands = TRUE,
    nuts2_shp,
    pu_restoration_budget_data,
    pu_restoration_budget_data_country,
    restoration_scenario,
    name
)
{
  # FEATURE WEIGHTS matrix
  features_weighted <- targs |>
    mutate(weight = ifelse(feature == "999999", (nrow(targs)-1)*2*carbon_weight, weight)) |>
    # mutate(weight = ifelse(feature == "urban", 1000, weight)) |>
    dplyr::select(weight) |>
    as.matrix()

  targs_filtered <- targs |> dplyr::select(-weight)

  p_solve <- p |>
    add_proportion_decisions() |>
    add_manual_targets(targs_filtered) |>
#    add_cplex_solver(gap=0.2, threads = cores, time_limit = 60*60*12) |>
#   add_highs_solver(gap=0.2, threads = cores, time_limit = 60*60*12) |>
   add_cbc_solver(gap=0.2, threads = cores, time_limit = 60*60*12) |>
#    add_gurobi_solver(gap = 0.2, threads = cores, numeric_focus = FALSE,  verbose = TRUE) |>
    add_manual_bounded_constraints(manual_bounded_constraints) |>
    add_feature_weights(features_weighted)

  # restoration budgets
  restoration <-  pu_restoration_budget_data |>
    dplyr::select(-c(pu)) |>
    as.matrix() |>
    replace_na(0)


  p_solve <- p_solve |>
    add_linear_constraints(threshold = nrow(pu)*restoration_constraint*1.01,
                           sense = "<=",
                           data = restoration)

  p_solve <- p_solve |>
    add_linear_constraints(threshold = nrow(pu)*restoration_constraint*0.9,
                           sense = ">=",
                           data = restoration)

  if(restoration_scenario == "Baseline"){
    p_solve <- p_solve |>
      add_linear_constraints(threshold = nrow(pu)*0.07,
                             sense = "<=",
                             data = restoration_crop)
    p_solve <- p_solve |>
      add_linear_constraints(threshold = nrow(pu)*0.09,
                             sense = "<=",
                             data = restoration_forest)
    p_solve <-  p_solve |>
      add_linear_constraints(threshold = nrow(pu)*0.0013,
                             sense = "<=",
                             data = restoration_natural)
    p_solve <-  p_solve |>
      add_linear_constraints(threshold = nrow(pu)*0.0086,
                             sense = "<=",
                             data = restoration_wetland)
  }
  if(restoration_scenario == "HN"){
    p_solve <- p_solve |>
      add_linear_constraints(threshold = nrow(pu)*0.05,
                             sense = ">=",
                             data = restoration_crop)
    p_solve <- p_solve |>
      add_linear_constraints(threshold = nrow(pu)*0.05,
                             sense = ">=",
                             data = restoration_forest)
    p_solve <-  p_solve |>
      add_linear_constraints(threshold = nrow(pu)*0.05,
                             sense = "<=",
                             data = restoration_natural)
    p_solve <-  p_solve |>
      add_linear_constraints(threshold = nrow(pu)*0.0086,
                             sense = "<=",
                             data = restoration_wetland)
  }
  # country constraints for restoration
  if(country_constraints == "EVEN"){

    country2 <- nuts2_shp

    country_names <- unique(country2$country)

    for(i in 1:length(country_names)) {
      country2 <- country_names[[i]]

      #if(nrow(nuts_crop)>0) {

      country_budget_restore <- pu_restoration_budget_data_country |>
        mutate(m = ifelse(country == country2, 1, 0))
      country_budget_restore <- sum(country_budget_restore$m)

      country_budget <- pu_restoration_budget_data_country |>
        mutate(m = ifelse(country == country2, 1, 0)) |>
        mutate_at(vars(starts_with("z")), ~.*m) |>
        dplyr::select(-c(m, country, pu)) |>
        as.matrix() |>
        replace_na(0)

      p_solve <- p_solve |>
        add_linear_constraints(threshold = country_budget_restore*restoration_constraint*1.1,
                               sense = "<=",
                               data = country_budget)
      #}
      #else(print("no country constraint"))
    }
  }
    # country constraints for restoration
  if(country_constraints == "FLEX"){

      country2 <- nuts2_shp

      country_names <- unique(country2$country)

      for(i in 1:length(country_names)) {
        country2 <- country_names[[i]]

        #if(nrow(nuts_crop)>0) {

        country_budget_restore <- pu_restoration_budget_data_country |>
          mutate(m = ifelse(country == country2, 1, 0))
        country_budget_restore <- sum(country_budget_restore$m)

        country_budget <- pu_restoration_budget_data_country |>
          mutate(m = ifelse(country == country2, 1, 0)) |>
          mutate_at(vars(starts_with("z")), ~.*m) |>
          dplyr::select(-c(m, country, pu)) |>
          as.matrix() |>
          replace_na(0)

        p_solve <- p_solve |>
          add_linear_constraints(threshold = country_budget_restore*(restoration_constraint+0.1),
                                 sense = "<=",
                                 data = country_budget)
        #}
        #else(print("no country constraint"))
      }
      if(country_constraints == "UNCONSTRAINED"){
        p_solve <- p_solve
      }
  }
  # solve
  s <- solve(p_solve)

  #ss <- s |> dplyr::select(solution_1_z1:solution_1_z26)

  # save!!
  write_csv(s, paste0("data/solutions/sol/sol_carbon_",carbon_weight,
                      "_restoration_", restoration_constraint,
                      "_production_", production_constraints,
                      "_country_", country_constraints,
                      "_wetlands_", wetlands,
                      "_onlyrestoration_", restoration_only,
                      "_scenario_", restoration_scenario,
                      "_", name,
                      ".csv"))
}

######### load in and some minor formating of data ################

pu_in_EU <- read_csv("data-formatted/pu_in_EU.csv")

pu <- read_fst("data-formatted/pu_data.fst") |>
  left_join(pu_in_EU) |>
  rename(id = EU_id) |>
  dplyr::select(-c(pu, nuts2id)) |>
  drop_na(id)


rij <- read_fst("data-formatted/features_split.fst") |>
  rename(species = feature) |>
  #mutate(amount = round(amount)) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  dplyr::select(pu, species, zone, amount) |>
  drop_na(pu) |>
  mutate(amount = ifelse(amount<0.001, 0, amount)) |>
  mutate(amount = replace_na(amount,0))


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
z <- read_csv("data-formatted/zone_id.csv") |>
  mutate(name = paste0("z", id)) |>
  relocate(id, 1) |>
  dplyr::select(-zone)

# BOUNDED CONSTRAINTS
# # zone ids for this..

zones <- read_csv("data-formatted/zone_id.csv") |>
  mutate(name = paste0("z", id))

# FEATURE TARGETS
names <- readRDS("data-formatted/SDMNameMatching.rds") |>
  mutate(speciesname= current_sname) |>
  rename(id = taxon_id)|>
  dplyr::select(id, speciesname)

names$speciesname <- sub(" ", "_", names$speciesname)

#targs <- read_csv("data-formatted/targets_split_formatted.csv")
# target formatting to align with RIJ table -- this isn't perfect, but there are
# <4% mismatch between the features we have in the rij table and the features
# we have targets for (due to small outstanding mismatches in spp names)

targs_existing <- read_csv("data-formatted/targets_split.csv") |>
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

# all the constraints

manual_bounded_constraints <- read_csv("data-formatted/manual_bounded_constraints_production_globiom_flex.csv") |>
  left_join(zones) |>
  dplyr::select(-c(zone, nuts2id)) |>
  rename(zone=name) |>
  dplyr::select(pu, zone, lower, upper) |>
  drop_na() |>
  # allow production expansion to meet targets
  mutate(upper = ifelse(zone == "z13", 1, upper)) |>
  mutate(upper = ifelse(zone == "z14", 1, upper)) |>
  mutate(upper = ifelse(zone == "z18", 1, upper)) |>
  mutate(lower = ifelse(zone == "z5", 0, lower)) |>
  mutate(upper = ifelse(zone == "z5", 0, upper))

# all the constraints
#nuts2_crop_low <- read_csv("data-formatted/linear_constraints/nuts_crop_low_95.csv")
nuts2_crop_low <- read_csv("data-formatted/linear_constraints/nuts_crop_low_55_adjusted.csv")
pu_cropland_low_budget_data <- read_csv("data-formatted/linear_constraints/pu_crop_low_budget_data_55.csv")

#nuts2_crop_med <- read_csv("data-formatted/linear_constraints/nuts_crop_med_95.csv")
nuts2_crop_med <- read_csv("data-formatted/linear_constraints/nuts_crop_med_55_adjusted.csv")
pu_cropland_med_budget_data <- read_csv("data-formatted/linear_constraints/pu_crop_med_budget_data_55.csv")

#nuts2_crop_high <- read_csv("data-formatted/linear_constraints/nuts_crop_high_95.csv")
nuts2_crop_high <- read_csv("data-formatted/linear_constraints/nuts_crop_high_55_adjusted.csv")
pu_cropland_high_budget_data <- read_csv("data-formatted/linear_constraints/pu_crop_high_budget_data_55.csv")

#nuts2_pasture_high <- read_csv("data-formatted/linear_constraints/nuts_pasture_high_95.csv")
nuts2_pasture_high <- read_csv("data-formatted/linear_constraints/nuts_pasture_high_55_adjusted.csv")
pu_pasture_high_budget_data <- read_csv("data-formatted/linear_constraints/pu_pasture_high_budget_data_55.csv")

#nuts2_pasture_low <- read_csv("data-formatted/linear_constraints/nuts_pasture_low_95.csv")
nuts2_pasture_low <- read_csv("data-formatted/linear_constraints/nuts_pasture_low_55_adjusted.csv")
pu_pasture_low_budget_data <- read_csv("data-formatted/linear_constraints/pu_pasture_low_budget_data_55.csv")

nuts2_forest_multi <- read_csv("data-formatted/linear_constraints/nuts_forest_multi_55_adjusted.csv")
pu_forest_multi_budget_data <- read_csv("data-formatted/linear_constraints/pu_forest_multi_budget_data_55.csv")

nuts2_forest_prod <- read_csv("data-formatted/linear_constraints/nuts_forest_prod_55_adjusted.csv")
pu_forest_prod_budget_data <- read_csv("data-formatted/linear_constraints/pu_forest_prod_budget_data_55.csv")

nuts2_wetland_rest <- read_csv("data-formatted/linear_constraints/wetland_targets.csv")
pu_wetland_rest_budget_data <- read_csv("data-formatted/linear_constraints/pu_wetland_rest_budget_data.csv")

country_nuts <- pu_pasture_low_budget_data |> dplyr::select(pu, NUTS_ID)

# Restoration budget breakdown
pu_restore_natural_budget_data <- read_csv("data-formatted/linear_constraints/pu_restore_natural_budget_data.csv")

pu_restore_forest_budget_data <- read_csv("data-formatted/linear_constraints/pu_restore_forest_budget_data.csv")

pu_restore_crop_budget_data <- read_csv("data-formatted/linear_constraints/pu_restore_crop_budget_data.csv")

# restoration
pu_restoration_budget_data <- read_csv("data-formatted/linear_constraints/pu_restoration_budget_data.csv")|>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id)) |>
  mutate(z8 = 0)

country_nuts <- pu_pasture_low_budget_data |> dplyr::select(pu, NUTS_ID)

pu_restoration_budget_data_country <- read_csv("data-formatted/linear_constraints/pu_restoration_budget_data.csv")|>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |> left_join(country_nuts) |>
  mutate(country = substr(NUTS_ID, start = 1, stop = 2)) |>
  dplyr::select(-c(id, cost, nuts2id, NUTS_ID, EU_id))


restoration_crop <-  pu_restore_crop_budget_data |>
  dplyr::select(-c(pu)) |>
  as.matrix() |>
  replace_na(0)

restoration_forest <-  pu_restore_forest_budget_data |>
  dplyr::select(-c(pu)) |>
  as.matrix() |>
  replace_na(0)

restoration_natural <-  pu_restore_natural_budget_data |>
  dplyr::select(-c(pu)) |>
  as.matrix() |>
  replace_na(0)

restoration_wetland <-  pu_wetland_rest_budget_data |>
  dplyr::select(-c(pu, NUTS_ID)) |>
  as.matrix() |>
  replace_na(0)

# add production constraints
nuts2_shp <- st_read("data-formatted/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  mutate(country = substr(NUTS2, start = 1, stop = 2)) |>
  filter(country != "UK")
nuts2_all <- nuts2_shp |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nutsIDnum = seq(1:nrow(nuts2_shp)))

nuts2_names <- nuts2_all$NUTS_ID

# deal with infeasible nuts2

infeas_nuts2 <- read_csv("data-formatted/feasibility-tests/nuts2_feasible.csv") |>
  filter(TF == "TRUE")

# pre set up basic problem
p <- problem(x = pu,
             features = feat,
             zones = z,
             cost_column = cost_columns,
             rij = rij) |>
  add_min_shortfall_objective(nrow(pu))

# wetland restoration
for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
  nuts2 <-  nuts2_names[[i]]

  nuts_wetland <- nuts2_wetland_rest |>
    filter(NUTS_ID == nuts2) |>
    mutate(wetlandarea = value)

  if(nrow(nuts_wetland)>0) {

    wetland_budget <- pu_wetland_rest_budget_data |>
      mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
      mutate_at(vars(starts_with("z")), ~.*m) |>
      dplyr::select(-c(m, NUTS_ID, pu)) |>
      as.matrix() |>
      replace_na(0)


    p <- p |>
      add_linear_constraints(threshold = nuts_wetland$wetlandarea[[1]]/10,
                             sense = "<=",
                             data = wetland_budget)
    #print(nuts2)
    #threshold = nuts_wetland$wetlandarea[[1]]/10
#    print(threshold)
  }
  else(p <- p |>
         add_linear_constraints(threshold = 0,
                                sense = "<=",
                                data = wetland_budget))
}


# crop low

for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
  nuts2 <-  nuts2_names[[i]]
  adj <- ifelse(nuts2 %in% setdiff(infeas_nuts2$nuts, c("SE02", "NL23", "RO08")), 0.8,
                ifelse(nuts2 %in% c("SE02", "NL23"), 0.5,
                       ifelse(nuts2 == "RO08", 0.25, 0.98)))

  nuts_crop <- nuts2_crop_low |>
    mutate(croparea = replace_na(value,0)) |>
    filter(NUTS_ID == nuts2) |>
    mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here

  if(nrow(nuts_crop)>0) {

    crop_budget <- pu_cropland_low_budget_data |>
      mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
      mutate_at(vars(starts_with("z")), ~.*m) |>
      dplyr::select(-c(m, NUTS_ID, pu)) |>
      as.matrix() |>
      replace_na(0)

    p <- p |>
      add_linear_constraints(threshold = nuts_crop$croparea[[1]]/10*adj,
                             sense = ">=",
                             data = crop_budget)
    #print(nuts2)
  }
  else(print("no crop low constraint"))
}

# crop med

for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)

  nuts2 <- nuts2_names[[i]]
  adj <- ifelse(nuts2 %in% setdiff(infeas_nuts2$nuts, c("SE02", "NL23", "RO08")), 0.8,
                ifelse(nuts2 %in% c("SE02", "NL23"), 0.5,
                       ifelse(nuts2 == "RO08", 0.25, 0.98)))
  nuts_crop <- nuts2_crop_med |>
    mutate(croparea = replace_na(value,0)) |>
    filter(NUTS_ID == nuts2) |>
    mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here


  if(nrow(nuts_crop)>0) {

    crop_budget <- pu_cropland_med_budget_data |>
      mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
      mutate_at(vars(starts_with("z")), ~.*m) |>
      dplyr::select(-c(m, NUTS_ID, pu)) |>
      as.matrix() |>
      replace_na(0)

    p <- p |>
      add_linear_constraints(threshold = nuts_crop$croparea[[1]]/10*adj,
                             sense = ">=",
                             data = crop_budget)
    #print(nuts2)
  }
  else(print("no crop med constraint"))
}

#crop high

for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
  nuts2 <- nuts2_names[[i]]
  adj <- ifelse(nuts2 %in% setdiff(infeas_nuts2$nuts, c("SE02", "NL23", "RO08")), 0.8,
                ifelse(nuts2 %in% c("SE02", "NL23"), 0.5,
                       ifelse(nuts2 == "RO08", 0.25, 0.98)))
  nuts_crop <- nuts2_crop_high |>
    mutate(croparea = replace_na(value,0)) |>
    filter(NUTS_ID == nuts2) |>
    mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here

  if(nrow(nuts_crop)>0) {

    crop_budget <- pu_cropland_high_budget_data |>
      mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
      mutate_at(vars(starts_with("z")), ~.*m) |>
      dplyr::select(-c(m, NUTS_ID, pu)) |>
      as.matrix() |>
      replace_na(0)

    p <- p |>
      add_linear_constraints(threshold = nuts_crop$croparea[[1]]/10*adj,
                             sense = ">=",
                             data = crop_budget)
    #print(nuts2)
  }
  else(print("no crop high constraint"))
}

#   # Pasture low
#
for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
  nuts2 <- nuts2_names[[i]]
  adj <- ifelse(nuts2 %in% setdiff(infeas_nuts2$nuts, c("SE02", "NL23", "RO08")), 0.8,
                ifelse(nuts2 %in% c("SE02", "NL23"), 0.5,
                       ifelse(nuts2 == "RO08", 0.25, 0.98)))

  nuts_pasture <- nuts2_pasture_low |>
    mutate(pasturearea = replace_na(value,0)) |>
    filter(NUTS_ID == nuts2) |>
    mutate(value = replace_na(value, 0),
           pasturearea = replace_na(pasturearea,0))|>
    mutate(pasturearea = ifelse(value < 1, 0, value))

  if(nrow(nuts_pasture)>0) {

    pasture_budget <- pu_pasture_low_budget_data |>
      mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
      mutate_at(vars(starts_with("z")), ~.*m) |>
      dplyr::select(-c(m, NUTS_ID, pu)) |>
      as.matrix() |>
      replace_na(0)

    p <- p |>
      add_linear_constraints(threshold = nuts_pasture$pasturearea[[1]]/10*adj,
                             sense = ">=",
                             data = pasture_budget)
    #print(nuts2)
  }
  else(print("no pasture low constraint"))
}

#   # Pasture high
#
for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
  nuts2 <- nuts2_names[[i]]
  adj <- ifelse(nuts2 %in% setdiff(infeas_nuts2$nuts, c("SE02", "NL23", "RO08")), 0.8,
                ifelse(nuts2 %in% c("SE02", "NL23"), 0.5,
                       ifelse(nuts2 == "RO08", 0.25, 0.98)))
  nuts_pasture <- nuts2_pasture_high |>
    mutate(pasturearea = replace_na(value,0)) |>
    filter(NUTS_ID == nuts2) |>
    mutate(value = replace_na(value, 0),
           pasturearea = replace_na(pasturearea,0))|>
    mutate(pasturearea = ifelse(value < 1, 0, value))

  if(nrow(nuts_pasture)>0) {

    pasture_budget <- pu_pasture_high_budget_data |>
      mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
      mutate_at(vars(starts_with("z")), ~.*m) |>
      dplyr::select(-c(m, NUTS_ID, pu)) |>
      as.matrix() |>
      replace_na(0)

    p <- p |>
      add_linear_constraints(threshold = nuts_pasture$pasturearea[[1]]/10*adj,
                             sense = ">=",
                             data = pasture_budget)
    #print(nuts2)
  }
  else(print("no pasture high constraint"))
}

#
# # # Forest multi
for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
  nuts2 <- nuts2_names[[i]]
  adj <- ifelse(nuts2 %in% setdiff(infeas_nuts2$nuts, c("SE02", "NL23", "RO08")), 0.8,
                ifelse(nuts2 %in% c("SE02", "NL23"), 0.5,
                       ifelse(nuts2 == "RO08", 0.25, 0.98)))
  nuts_forest <- nuts2_forest_multi |>
    mutate(forestarea = replace_na(value,0)) |>
    filter(NUTS_ID == nuts2)

  if(nrow(nuts_forest)>0.01) {

    forest_budget <- pu_forest_multi_budget_data |>
      mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
      mutate_at(vars(starts_with("z")), ~.*m) |>
      dplyr::select(-c(m, NUTS_ID, pu)) |>
      as.matrix() |>
      replace_na(0)

    forest_threshold <- ifelse(nuts_forest$forestarea[[1]]/10<1, 0,
                               nuts_forest$forestarea[[1]]/10*adj)

    p <- p |>
      add_linear_constraints(threshold = forest_threshold,
                             sense = ">=",
                             data = forest_budget)
    #print(nuts2)
  }
  else(print("no forest multi constraint"))
}
#
#   # # Forest prod

for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
  nuts2 <- nuts2_names[[i]]
  adj <- ifelse(nuts2 %in% setdiff(infeas_nuts2$nuts, c("SE02", "NL23", "RO08")), 0.8,
                ifelse(nuts2 %in% c("SE02", "NL23"), 0,
                       ifelse(nuts2 == "RO08", 0, 0.98)))
  nuts_forest <- nuts2_forest_prod |>
    mutate(forestarea = replace_na(value,0)) |>
    filter(NUTS_ID == nuts2)

  if(nrow(nuts_forest)>0.01) {

    forest_budget <- pu_forest_prod_budget_data |>
      mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
      mutate_at(vars(starts_with("z")), ~.*m) |>
      dplyr::select(-c(m, NUTS_ID, pu)) |>
      as.matrix() |>
      replace_na(0)

    forest_threshold <- ifelse(nuts_forest$forestarea[[1]]/10< 1, 0,
                               nuts_forest$forestarea[[1]]/10*adj)
    p <- p |>
      add_linear_constraints(threshold = forest_threshold,
                             sense = ">=",
                             data = forest_budget)
    #print(nuts2)

  }
  else(print("no forest prod constraint"))
}

### solve ######
scenarios <-
  tidyr::crossing(# Do all combinations of:
    carbon_weight = seq(0,2, by = 0.5),
    country_constraints = c("EVEN", "FLEX", "UNCONSTRAINED"),
    restoration_constraint = 0.14,
    restoration_scenario = c("Baseline", "HN"),
  )

for (i in 2:nrow(scenarios)){ #nrow(scenarios)){
  problem_setup(
    cores = parallel::detectCores(),
    p, # conservation problem + min shortfall objective/ speeds things up to do this first...
    pu, # pu data (x from p)
    targs,
    manual_bounded_constraints = manual_bounded_constraints,
    carbon_weight = scenarios$carbon_weight[[i]],
    restoration_constraint = 0.141,
    #conservation_constraint = 0.3,
    production_constraints = TRUE,
    country_constraints = scenarios$country_constraints[[i]], # EVEN, FLEX, UNEVEN
    restoration_only = FALSE,
    wetlands = TRUE,
    nuts2_shp = nuts2_shp,
    pu_restoration_budget_data = pu_restoration_budget_data,
    pu_restoration_budget_data_country = pu_restoration_budget_data_country,
    restoration_scenario = scenarios$restoration_scenario[[i]],
    name = "globiom_IC_f455")
}

