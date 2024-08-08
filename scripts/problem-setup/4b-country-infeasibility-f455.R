### set production constraints ensuring feasibility
library(prioritizr)
library(fst)
library(highs)
library(gurobi)
library(tidyverse)
library(sf)
cores <- 7
fst::threads_fst(nr_of_threads = cores, reset_after_fork = NULL)  # Set number for fst
setwd(dir = "~/EU-restoration-prioritization/")

####### Data setup ##########
pu_in_EU <- read_csv("data/formatted-data/pu_in_EU.csv")

pu <- read_fst("data/formatted-data/pu_data.fst") |>
  left_join(pu_in_EU) |>
  rename(id = EU_id) |>
  dplyr::select(-c(pu, nuts2id)) |>
  drop_na(id)

# simplify features for quick solves
rij <- read_fst("data/formatted-data/features_split.fst") |>
  rename(species = feature) |>
  #mutate(amount = round(amount)) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  dplyr::select(pu, species, zone, amount) |>
  drop_na(pu) |>
  mutate(amount = ifelse(amount<0.001, 0, amount)) |>
  filter(species %in% c(10795880127, 999999)) |>
  mutate(amount = replace_na(amount, 0))

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

#targs <- read_csv("data/formatted-data/targets_split_formatted.csv")
# target formatting to align with RIJ table -- this isn't perfect, but there are
# <4% mismatch between the features we have in the rij table and the features
# we have targets for (due to small outstanding mismatches in spp names)

targs_existing <- read_csv("data/formatted-data/targets_split.csv") |>
  mutate(zone = list(z$name)) |>
  dplyr::filter(feature %in% feat$name) |>
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

manual_bounded_constraints <- read_csv("data/formatted-data/manual_bounded_constraints_CLC_adj.csv")

# all the constraints
#nuts2_crop_low <- read_csv("data/formatted-data/linear_constraints/nuts_crop_low_95.csv")
country_prod_constraints <- read_csv("data/formatted-data/linear_constraints/country_constraints_55_proportional.csv")

pu_cropland_low_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_crop_low_budget_data.csv")

pu_cropland_med_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_crop_med_budget_data.csv")

pu_cropland_high_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_crop_high_budget_data.csv")

pu_pasture_high_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_pasture_high_budget_data.csv")

pu_pasture_low_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_pasture_low_budget_data.csv")

pu_forest_multi_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_forest_multi_budget_data.csv")

pu_forest_prod_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_forest_prod_budget_data.csv")

country_nuts <- pu_pasture_low_budget_data |> dplyr::select(pu, NUTS_ID)

########### pre set up basic problem  ##########
problem <- problem(x = pu,
                   features = feat,
                   zones = z,
                   cost_column = cost_columns,
                   rij = rij) |>
  add_min_shortfall_objective(nrow(pu))


p1 <- problem



#remotes::install_github("cboettig/prioritizr")
library(prioritizr)
library(fst)
library(highs)
library(tidyverse)
library(sf)
library(furrr)
library(gurobi)

rm(list = ls())

cores <- 7

fst::threads_fst(nr_of_threads = cores, reset_after_fork = NULL)  # Set number for fst
options(future.globals.maxSize = 8000 * 1024^2)

setwd(dir = "~/EU-restoration-prioritization/") # to run in terminal...

########## function to set up, solve, and save solution #################

problem_setup <- function(
    cores = 7,
    p, # conservation problem + min shortfall objective/ speeds things up to do this first...
    pu, # pu data (x from p)
    targs,
    manual_bounded_constraints = manual_bounded_constraints,
    carbon_weight = 0.5,
    restoration_constraint = 0.141,
    #conservation_constraint = 0.3,
    production_constraints = "NUTS2adj", #NUTS2, NUTS2adj, COUNTRY
    country_constraints = "UNEVEN", # EVEN, FLEX, UNEVEN
    restoration_only = FALSE,
    wetlands = TRUE,
    nuts2_shp,
    pu_restoration_budget_data,
    pu_restoration_budget_data_country,
    restoration_scenario,
    future = "f455", #f455, ref
    name,
    solver = "gurobi",
    SI = FALSE
)
{
  # FEATURE WEIGHTS matrix
  features_weighted <- targs |>
    mutate(weight = ifelse(feature == "999999", (nrow(targs)-1)*2*carbon_weight, weight)) |>
    # mutate(weight = ifelse(feature == "urban", 1000, weight)) |>
    dplyr::select(weight) |>
    as.matrix()

  targs_filtered <- targs |> dplyr::select(-weight)

  if(restoration_only == TRUE) {
    manual_bounded_constraints <- manual_bounded_constraints |>
      mutate(upper = ifelse(zone %in% c("z12", "z13", "z14", "z15",
                                        "z16", "z17", "z18", "z19", "z20",
                                        "z21", "z22", "z23", "z23", "z24",
                                        "z25"), 0,  upper)) |>
      mutate(lower = ifelse(zone %in% c("z12", "z13", "z14", "z15",
                                        "z16", "z17", "z18", "z19", "z20",
                                        "z21", "z22", "z23", "z23", "z24",
                                        "z25"), 0,  lower))
  }

  # if no gurobi license, not recommended for lots of scenarios
  if(solver == "highs") {
    p_solve <- p |>
      add_proportion_decisions() |>
      add_manual_targets(targs_filtered) |>
      add_highs_solver(gap = 0.2, threads = cores,
                       verbose = TRUE, time_limit = 3600) |>
      add_manual_bounded_constraints(manual_bounded_constraints) |>
      add_feature_weights(features_weighted)
  }

  if(solver == "gurobi") {
    p_solve <- p |>
      add_proportion_decisions() |>
      add_manual_targets(targs_filtered) |>
      add_gurobi_solver(gap = 0.2, threads = cores,
                        verbose = TRUE) |>
      add_manual_bounded_constraints(manual_bounded_constraints) |>
      add_feature_weights(features_weighted)
  }

  # restoration budgets
  restoration <-  pu_restoration_budget_data |>
    dplyr::select(-c(pu)) |>
    as.matrix() |>
    replace_na(0)

  # no more than a bit over the restoration goal of ~14%
  p_solve <- p_solve |>
    add_linear_constraints(threshold = nrow(pu)*restoration_constraint*1.01,
                           sense = "<=",
                           data = restoration)

  # no less than restoration 25% restoration target
  p_solve <- p_solve |>
    add_linear_constraints(threshold = nrow(pu)*restoration_constraint*0,
                           sense = ">=",
                           data = restoration)


  # deintensify ag, little nature restoration
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

  # high nature scenario
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
  }


  if(country_constraints == "UNEVEN"){
    p_solve <- p_solve
  }

  # wetland restoration

  if(wetlands == "TRUE"){
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

        #threshold = nuts_wetland$wetlandarea[[1]]/10

        p_solve <- p_solve |>
          add_linear_constraints(threshold = nuts_wetland$wetlandarea[[1]]/10,
                                 sense = "<=",
                                 data = wetland_budget)
      }
      else(p_solve <- p_solve |>
             add_linear_constraints(threshold = 0,
                                    sense = "<=",
                                    data = wetland_budget))
    }

  }


  # allow restoration only option

  if(restoration_only == FALSE){
    # define production scenario
    if(future == "f455"){
      nuts2_crop_low <- read_csv("data/formatted-data/linear_constraints/nuts_crop_low_55_proportional.csv")

      nuts2_crop_med <- read_csv("data/formatted-data/linear_constraints/nuts_crop_med_55_proportional.csv")

      nuts2_crop_high <- read_csv("data/formatted-data/linear_constraints/nuts_crop_high_55_proportional.csv")

      nuts2_pasture_high <- read_csv("data/formatted-data/linear_constraints/nuts_pasture_high_55_proportional.csv")

      nuts2_pasture_low <- read_csv("data/formatted-data/linear_constraints/nuts_pasture_low_55_proportional.csv")

      nuts2_forest_multi <- read_csv("data/formatted-data/linear_constraints/nuts_forest_multi_55_proportional.csv")

      nuts2_forest_prod <- read_csv("data/formatted-data/linear_constraints/nuts_forest_prod_55_proportional.csv")

      country_prod_constraints <- read_csv("data/formatted-data/linear_constraints/country_constraints_55_proportional.csv")|>
        mutate(value = ifelse(country %in% c("FR", "ES", "IT"), 0, value))

      sum(country_prod_constraints$value)

      # deal with in-feasible nuts2

      infeas_nuts2 <- read_csv("data/feasibility-tests/nuts2_feasible_f455_adj.csv") |>
        mutate(adj = ifelse(TF =="FALSE", 1,0))

      if(production_constraints %in% c("NUTS2adj", "NUTS2")){

        # crop low nuts2
        for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
          nuts2 <-  nuts2_names[[i]]
          adj <- infeas_nuts2 |> filter(nuts == nuts2)
          adj <- adj$adj[1]

          nuts_crop <- nuts2_crop_low |>
            mutate(croparea = replace_na(value,0)) |>
            filter(NUTS_ID == nuts2) |> drop_na() |>
            mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here

          if(nrow(nuts_crop)>0) {

            crop_budget <- pu_cropland_low_budget_data |>
              mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
              mutate_at(vars(starts_with("z")), ~.*m) |>
              dplyr::select(-c(m, NUTS_ID, pu)) |>
              as.matrix() |>
              replace_na(0)

            p_solve <- p_solve |>
              add_linear_constraints(threshold = nuts_crop$croparea[[1]]*adj,#/10*adj,
                                     sense = ">=",
                                     data = crop_budget)
            print(nuts2)
          }
          else(print("no crop low constraint"))
        }

        # crop med nuts2
        for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)

          nuts2 <- nuts2_names[[i]]
          adj <- infeas_nuts2 |> filter(nuts == nuts2)
          adj <- adj$adj[1]

          nuts_crop <- nuts2_crop_med |>
            mutate(croparea = replace_na(value,0)) |>
            filter(NUTS_ID == nuts2) |>drop_na() |>
            mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here


          if(nrow(nuts_crop)>0) {

            crop_budget <- pu_cropland_med_budget_data |>
              mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
              mutate_at(vars(starts_with("z")), ~.*m) |>
              dplyr::select(-c(m, NUTS_ID, pu)) |>
              as.matrix() |>
              replace_na(0)

            p_solve <- p_solve |>
              add_linear_constraints(threshold = nuts_crop$croparea[[1]]*adj, #/10*adj,
                                     sense = ">=",
                                     data = crop_budget)
            print(nuts2)
          }
          else(print("no crop med constraint"))
        }

        #crop high nuts2
        for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
          nuts2 <- nuts2_names[[i]]
          adj <- infeas_nuts2 |> filter(nuts == nuts2)
          adj <- adj$adj[1]

          nuts_crop <- nuts2_crop_high |>
            mutate(croparea = replace_na(value,0)) |>
            filter(NUTS_ID == nuts2) |>drop_na() |>
            mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here

          if(nrow(nuts_crop)>0) {

            crop_budget <- pu_cropland_high_budget_data |>
              mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
              mutate_at(vars(starts_with("z")), ~.*m) |>
              dplyr::select(-c(m, NUTS_ID, pu)) |>
              as.matrix() |>
              replace_na(0)

            p_solve <- p_solve |>
              add_linear_constraints(threshold = nuts_crop$croparea[[1]]*adj,#/10*adj,
                                     sense = ">=",
                                     data = crop_budget)
            print(nuts2)
          }
          else(print("no crop high constraint"))
        }

        # Pasture low nuts2
        for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
          nuts2 <- nuts2_names[[i]]
          adj <- infeas_nuts2 |> filter(nuts == nuts2)
          adj <- adj$adj[1]

          nuts_pasture <- nuts2_pasture_low |>
            mutate(pasturearea = replace_na(value,0)) |>
            filter(NUTS_ID == nuts2) |>drop_na() |>
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

            p_solve <- p_solve |>
              add_linear_constraints(threshold = nuts_pasture$pasturearea[[1]]*adj,
                                     sense = ">=",
                                     data = pasture_budget)
            print(nuts2)
          }
          else(print("no pasture low constraint"))
        }

        # pasture high nuts2
        for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
          nuts2 <- nuts2_names[[i]]
          adj <- infeas_nuts2 |> filter(nuts == nuts2)
          adj <- adj$adj[1]

          nuts_pasture <- nuts2_pasture_high |>
            mutate(pasturearea = replace_na(value,0)) |>drop_na() |>
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

            p_solve <- p_solve |>
              add_linear_constraints(threshold = nuts_pasture$pasturearea[[1]]*adj,#/10*adj,
                                     sense = ">=",
                                     data = pasture_budget)
            print(nuts2)
          }
          else(print("no pasture high constraint"))
        }

        # forest multi nuts2
        for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
          nuts2 <- nuts2_names[[i]]
          adj <- infeas_nuts2 |> filter(nuts == nuts2)
          adj <- adj$adj[1]

          nuts_forest <- nuts2_forest_multi |>drop_na() |>
            mutate(forestarea = replace_na(value,0)) |>
            filter(NUTS_ID == nuts2)

          if(nrow(nuts_forest)>0.01) {

            forest_budget <- pu_forest_multi_budget_data |>
              mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
              mutate_at(vars(starts_with("z")), ~.*m) |>
              dplyr::select(-c(m, NUTS_ID, pu)) |>
              as.matrix() |>
              replace_na(0)

            forest_threshold <- nuts_forest$forestarea[[1]]*adj#/10*adj

            p_solve <- p_solve |>
              add_linear_constraints(threshold = forest_threshold,
                                     sense = ">=",
                                     data = forest_budget)
            print(nuts2)
          }
          else(print("no forest multi constraint"))
        }

        # forest production nuts2
        for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
          nuts2 <- nuts2_names[[i]]
          adj <- infeas_nuts2 |> filter(nuts == nuts2)
          adj <- adj$adj[1]

          nuts_forest <- nuts2_forest_prod |>
            mutate(forestarea = replace_na(value,0)) |>
            filter(NUTS_ID == nuts2)

          if(nrow(nuts_forest)>0.01) {

            forest_budget <- pu_forest_prod_budget_data |>drop_na() |>
              mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
              mutate_at(vars(starts_with("z")), ~.*m) |>
              dplyr::select(-c(m, NUTS_ID, pu)) |>
              as.matrix() |>
              replace_na(0)

            forest_threshold <- nuts_forest$forestarea[[1]]*adj #/10*adj)
            p_solve <- p_solve |>
              add_linear_constraints(threshold = forest_threshold,
                                     sense = ">=",
                                     data = forest_budget)
            print(nuts2)

          }
          else(print("no forest prod constraint"))
        }


      }

      if(production_constraints %in% c("NUTS2adj", "COUNTRY")) {

        country_names <- unique(nuts2_shp$country)

        # crop low country
        for(i in 1:length(country_names)) {
          country2 <- country_names[[i]]

          #if(nrow(nuts_crop)>0) {
          production_zone <- "Cropland_low"
          country_budget <- pu_cropland_low_budget_data |>
            mutate(country = substr(NUTS_ID, 1, 2)) |>
            mutate(m = ifelse(country == country2, 1, 0)) |>
            mutate_at(vars(starts_with("z")), ~.*m) |>
            dplyr::select(-c(m, country, pu, NUTS_ID)) |>
            as.matrix() |>
            replace_na(0)

          country_threshhold <- country_prod_constraints |>
            filter(country == country2) |>
            filter(name == production_zone)

          p_solve <- p_solve |>
            add_linear_constraints(threshold = country_threshhold$value[[1]],
                                   sense = ">=",
                                   data = country_budget)
          #}
          #else(print("no country constraint"))
        }

        # crop med country
        for(i in 1:length(country_names)) {
          country2 <- country_names[[i]]

          #if(nrow(nuts_crop)>0) {
          production_zone <- "Cropland_med"
          country_budget <- pu_cropland_med_budget_data |>
            mutate(country = substr(NUTS_ID, 1, 2)) |>
            mutate(m = ifelse(country == country2, 1, 0)) |>
            mutate_at(vars(starts_with("z")), ~.*m) |>
            dplyr::select(-c(m, country, pu, NUTS_ID)) |>
            as.matrix() |>
            replace_na(0)

          country_threshhold <- country_prod_constraints |>
            filter(country == country2) |>
            filter(name == production_zone)

          p_solve <- p_solve |>
            add_linear_constraints(threshold = country_threshhold$value[[1]],
                                   sense = ">=",
                                   data = country_budget)
          #}
          #else(print("no country constraint"))
        }

        # crop high country
        for(i in 1:length(country_names)) {
          country2 <- country_names[[i]]

          #if(nrow(nuts_crop)>0) {
          production_zone <- "Cropland_high"
          country_budget <- pu_cropland_high_budget_data |>
            mutate(country = substr(NUTS_ID, 1, 2)) |>
            mutate(m = ifelse(country == country2, 1, 0)) |>
            mutate_at(vars(starts_with("z")), ~.*m) |>
            dplyr::select(-c(m, country, pu, NUTS_ID)) |>
            as.matrix() |>
            replace_na(0)

          country_threshhold <- country_prod_constraints |>
            filter(country == country2) |>
            filter(name == production_zone)

          p_solve <- p_solve |>
            add_linear_constraints(threshold = country_threshhold$value[[1]],
                                   sense = ">=",
                                   data = country_budget)
          #}
          #else(print("no country constraint"))
        }

        # pasture low country
        for(i in 1:length(country_names)) {
          country2 <- country_names[[i]]

          #if(nrow(nuts_crop)>0) {
          production_zone <- "Pasture_low"
          country_budget <- pu_pasture_low_budget_data |>
            mutate(country = substr(NUTS_ID, 1, 2)) |>
            mutate(m = ifelse(country == country2, 1, 0)) |>
            mutate_at(vars(starts_with("z")), ~.*m) |>
            dplyr::select(-c(m, country, pu, NUTS_ID)) |>
            as.matrix() |>
            replace_na(0)

          country_threshhold <- country_prod_constraints |>
            filter(country == country2) |>
            filter(name == production_zone)

          p_solve <- p_solve |>
            add_linear_constraints(threshold = country_threshhold$value[[1]],
                                   sense = ">=",
                                   data = country_budget)
          #}
          #else(print("no country constraint"))
        }

        # pasture high country
        for(i in 1:length(country_names)) {
          country2 <- country_names[[i]]

          #if(nrow(nuts_crop)>0) {
          production_zone <- "Pasture_high"
          country_budget <- pu_pasture_high_budget_data |>
            mutate(country = substr(NUTS_ID, 1, 2)) |>
            mutate(m = ifelse(country == country2, 1, 0)) |>
            mutate_at(vars(starts_with("z")), ~.*m) |>
            dplyr::select(-c(m, country, pu, NUTS_ID)) |>
            as.matrix() |>
            replace_na(0)

          country_threshhold <- country_prod_constraints |>
            filter(country == country2) |>
            filter(name == production_zone)

          p_solve <- p_solve |>
            add_linear_constraints(threshold = country_threshhold$value[[1]],
                                   sense = ">=",
                                   data = country_budget)
          #}
          #else(print("no country constraint"))
        }

        # forest multi country
        for(i in 1:length(country_names)) {
          country2 <- country_names[[i]]

          #if(nrow(nuts_crop)>0) {
          production_zone <- "WoodlandForest_multi"
          country_budget <- pu_forest_multi_budget_data |>
            mutate(country = substr(NUTS_ID, 1, 2)) |>
            mutate(m = ifelse(country == country2, 1, 0)) |>
            mutate_at(vars(starts_with("z")), ~.*m) |>
            dplyr::select(-c(m, country, pu, NUTS_ID)) |>
            as.matrix() |>
            replace_na(0)

          country_threshhold <- country_prod_constraints |>
            filter(country == country2) |>
            filter(name == production_zone)

          p_solve <- p_solve |>
            add_linear_constraints(threshold = country_threshhold$value[[1]],
                                   sense = ">=",
                                   data = country_budget)
          #}
          #else(print("no country constraint"))
        }

        # forest production country
        for(i in 1:length(country_names)) {
          country2 <- country_names[[i]]

          #if(nrow(nuts_crop)>0) {
          production_zone <- "WoodlandForest_prod"
          country_budget <- pu_forest_prod_budget_data |>
            mutate(country = substr(NUTS_ID, 1, 2)) |>
            mutate(m = ifelse(country == country2, 1, 0)) |>
            mutate_at(vars(starts_with("z")), ~.*m) |>
            dplyr::select(-c(m, country, pu, NUTS_ID)) |>
            as.matrix() |>
            replace_na(0)

          country_threshhold <- country_prod_constraints |>
            filter(country == country2) |>
            filter(name == production_zone)

          p_solve <- p_solve |>
            add_linear_constraints(threshold = country_threshhold$value[[1]],
                                   sense = ">=",
                                   data = country_budget)
          #}
          #else(print("no country constraint"))
        }
      }

    }

  }
}

######### load in and some minor formatting of data ################
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
z <- read_csv("data/formatted-data/zone_id.csv") |>
  mutate(name = paste0("z", id)) |>
  relocate(id, 1) |>
  dplyr::select(-zone)

# add production constraints
nuts2_shp <- st_read("data/formatted-data/EU_GLOBIOM_NUTS2.shp") |>
  mutate(country = substr(NUTS2, start = 1, stop = 2)) |>
  filter(country != "UK")
nuts2_all <- nuts2_shp |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nutsIDnum = seq(1:nrow(nuts2_shp)))

nuts2_names <- nuts2_all$NUTS_ID

# BOUNDED CONSTRAINTS
# # zone ids for this..

zones <- read_csv("data/formatted-data/zone_id.csv") |>
  mutate(name = paste0("z", id))

# FEATURE TARGETS
names <- readRDS("data/formatted-data/SDMNameMatching.rds") |>
  mutate(speciesname= current_sname) |>
  rename(id = taxon_id)|>
  dplyr::select(id, speciesname)

names$speciesname <- sub(" ", "_", names$speciesname)

# targs <- read_csv("data-formatted/targets_split_formatted.csv")
# target formatting to align with RIJ table -- this isn't perfect, but there are
# <4% mismatch between the features we have in the rij table and the features
# we have targets for (due to small outstanding mismatches in spp names)

targs_existing <- read_csv("data/formatted-data/targets_split.csv") |>
  mutate(zone = list(z$name)) |>
  #dplyr::filter(feature %in% feat$name) |>
  mutate(target = ifelse(target > 100000, target/100, target),
         target = ifelse(target <0.001, 0, target))

# setdiff(feat$id,targs_existing$feature)

# to line up with rij spp
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
IC_mismatch <- read_fst("data/plotting/plotting_data_ic.fst") |> group_by(id) |>
  summarise(total = sum(value, na.rm = T)) |>
  rename(pu = id) |>
  mutate(diff = 1-total) |> # add to rivers and lakes LB to avoid benefits etc.
  dplyr::select(pu, diff) |>
  mutate(diff = diff - 0.02) |>
  mutate(diff = ifelse(diff < 0, 0, diff))

manual_bounded_constraints <- read_csv("data/formatted-data/manual_bounded_constraints_CLC_adj.csv")

# manual_bounded_constraints <- read_csv("data/formatted-data/manual_bounded_constraints_CLC.csv") |>
#   rename(pu = PUID) |>
#   left_join(pu_in_EU) |>
#   mutate(pu = EU_id) |>
#   left_join(zones) |>
#   mutate(zone = name) |>
#   dplyr::select(-c(name, EU_id, nuts2id)) |>
#   #rename(zone=name) |>
#   dplyr::select(pu, zone, lower, upper) |>
#   drop_na() |>
#   # allow production expansion to meet targets
#   mutate(upper = ifelse(zone == "z13", 1, upper)) |>
#   mutate(upper = ifelse(zone == "z14", 1, upper)) |>
#   mutate(upper = ifelse(zone == "z18", 1, upper)) |>
#   mutate(lower = ifelse(zone == "z5", 0, lower)) |>
#   mutate(upper = ifelse(zone == "z5", 0, upper)) |>
#   mutate(lower = ifelse(lower > upper, upper, lower)) |>
#   left_join(IC_mismatch) |>
#   mutate(upper = ifelse(zone == "z23", upper+diff, upper)) |>
#   mutate(lower = ifelse(zone %in% c("z23", "z24"), upper, lower)) |>
#   dplyr::select(-diff)

# all the constraints
pu_cropland_low_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_crop_low_budget_data.csv")

pu_cropland_med_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_crop_med_budget_data.csv")

pu_cropland_high_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_crop_high_budget_data.csv")

pu_pasture_high_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_pasture_high_budget_data.csv")

pu_pasture_low_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_pasture_low_budget_data.csv")

pu_forest_multi_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_forest_multi_budget_data.csv")

pu_forest_prod_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_forest_prod_budget_data.csv")

nuts2_wetland_rest <- read_csv("data/formatted-data/linear_constraints/wetland_targets.csv")
pu_wetland_rest_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_wetland_rest_budget_data.csv")

country_nuts <- pu_pasture_low_budget_data |> dplyr::select(pu, NUTS_ID)

# Restoration budget breakdown
pu_restore_natural_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_restore_natural_budget_data.csv")

pu_restore_forest_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_restore_forest_budget_data.csv")

pu_restore_crop_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_restore_crop_budget_data.csv")

# restoration
pu_restoration_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_restoration_budget_data.csv")|>
  rename(pu = layer) |>
  left_join(pu_in_EU) |>
  mutate(pu = EU_id) |>
  drop_na(pu) |>
  dplyr::select(-c(id, cost, nuts2id, EU_id)) |>
  mutate(z8 = 0)

country_nuts <- pu_pasture_low_budget_data |> dplyr::select(pu, NUTS_ID)

pu_restoration_budget_data_country <- read_csv("data/formatted-data/linear_constraints/pu_restoration_budget_data.csv")|>
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

# pre set up basic problem
p <- problem(x = pu,
             features = feat,
             zones = z,
             cost_column = cost_columns,
             rij = rij) |>
  add_min_shortfall_objective(nrow(pu))

# ### solve ######
scenarios <-
  tidyr::crossing(# Do all combinations of:
    carbon_weight = 2.5,#c(0.1,0.5,0.9,1.5, 2),#c(0.1, 0.3,seq(0.5,2, by = 0.5)),
    country_constraints = "EVEN", #c("EVEN", "FLEX", "UNEVEN"),
    restoration_constraint = 0.141,
    restoration_only = FALSE,
    production_constraints = "NUTS2adj",
    wetlands = TRUE,
    restoration_scenario = "Baseline",#c("Baseline", "HN"),
    future = "f455",
    name = "proportional",
    SI = FALSE
  )


i=1
print(i)
problem_setup(
  cores = parallel::detectCores(),
  p, # conservation problem + min shortfall objective/ speeds things up to do this first...
  pu, # pu data (x from p)
  targs,
  manual_bounded_constraints = manual_bounded_constraints,
  carbon_weight = scenarios$carbon_weight[[i]],
  restoration_constraint = scenarios$restoration_constraint[[i]],
  #conservation_constraint = 0.3,
  production_constraints = scenarios$production_constraints[[i]],
  country_constraints = scenarios$country_constraints[[i]], # EVEN, FLEX, UNEVEN
  restoration_only = scenarios$restoration_only[[i]],
  wetlands = scenarios$wetlands[[i]],
  nuts2_shp = nuts2_shp,
  pu_restoration_budget_data = pu_restoration_budget_data,
  pu_restoration_budget_data_country = pu_restoration_budget_data_country,
  restoration_scenario = scenarios$restoration_scenario[[i]],
  future = scenarios$future[[i]],
  name = scenarios$name[[i]],
  solver = "gurobi",
  SI = scenarios$SI[[i]] )




######### solve by nuts2 all constraints ##########

nuts2_shp <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  mutate(country = substr(NUTS2, start = 1, stop = 2))
nuts2_all <- nuts2_shp |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nutsIDnum = seq(1:260)) |>
  filter(country != "UK")
nuts2_names_all <- nuts2_all$NUTS_ID #setdiff(nuts2_all$NUTS_ID
feas <- read_csv("data/feasibility-tests/nuts2_feasible_f455.csv")

n = 1
adj <- 1
for(j in 1:length(nuts2_names_all)){
  # pre set up basic problem
  nuts <- nuts2_names_all[[j]]
  nuts2_all <- nuts2_shp |>
    rename(NUTS_ID = NURGCDL2) |>
    mutate(nutsIDnum = seq(1:260)) |>
    filter(NUTS_ID == nuts)
  nuts2_names <- nuts2_all$NUTS_ID

  # crop low

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <-  nuts2_names[[i]]
    nuts_crop <- nuts2_crop_low |>
      mutate(croparea = replace_na(value,0)) |>
      filter(NUTS_ID == nuts2) |> drop_na() |>
      mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here

    if(nrow(nuts_crop)>0) {

      crop_budget <- pu_cropland_low_budget_data |>
        mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
        mutate_at(vars(starts_with("z")), ~.*m) |>
        dplyr::select(-c(m, NUTS_ID, pu)) |>
        as.matrix() |>
        replace_na(0)

      p <- p1 |>
        add_linear_constraints(threshold = nuts_crop$croparea[[1]],#/10*adj,
                               sense = ">=",
                               data = crop_budget)
      #print(nuts2)
    }
    #else(print("no crop low constraint"))
  }

  # crop med

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)

    nuts2 <- nuts2_names[[i]]
    nuts_crop <- nuts2_crop_med |>
      mutate(croparea = replace_na(value,0)) |>
      filter(NUTS_ID == nuts2) |> drop_na() |>
      mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here


    if(nrow(nuts_crop)>0) {

      crop_budget <- pu_cropland_med_budget_data |>
        mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
        mutate_at(vars(starts_with("z")), ~.*m) |>
        dplyr::select(-c(m, NUTS_ID, pu)) |>
        as.matrix() |>
        replace_na(0)

      p <- p |>
        add_linear_constraints(threshold = nuts_crop$croparea[[1]],#/10*adj,
                               sense = ">=",
                               data = crop_budget)
      #print(nuts2)
    }
    # else(print("no crop med constraint"))
  }

  #crop high

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
    nuts_crop <- nuts2_crop_high |>
      mutate(croparea = replace_na(value,0)) |>
      filter(NUTS_ID == nuts2) |>drop_na() |>
      mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here

    if(nrow(nuts_crop)>0) {

      crop_budget <- pu_cropland_high_budget_data |>
        mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
        mutate_at(vars(starts_with("z")), ~.*m) |>
        dplyr::select(-c(m, NUTS_ID, pu)) |>
        as.matrix() |>
        replace_na(0)

      p <- p |>
        add_linear_constraints(threshold = nuts_crop$croparea[[1]],#/10*adj,
                               sense = ">=",
                               data = crop_budget)
      #print(nuts2)
    }
    # else(print("no crop high constraint"))
  }


  #   # Pasture low
  #
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
    nuts_pasture <- nuts2_pasture_low |>
      mutate(pasturearea = replace_na(value,0)) |>
      filter(NUTS_ID == nuts2) |>drop_na() |>
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
        add_linear_constraints(threshold = nuts_pasture$pasturearea[[1]],#/10*adj,
                               sense = ">=",
                               data = pasture_budget)
      #print(nuts2)
    }
    #else(print("no pasture low constraint"))
  }

  #   # Pasture high
  #
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
    nuts_pasture <- nuts2_pasture_high |>
      mutate(pasturearea = replace_na(value,0)) |>
      filter(NUTS_ID == nuts2) |>drop_na() |>
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
        add_linear_constraints(threshold = nuts_pasture$pasturearea[[1]],#/10*adj,
                               sense = ">=",
                               data = pasture_budget)
      #print(nuts2)
    }
    #else(print("no pasture high constraint"))
  }

  #
  # # # Forest multi
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
    nuts_forest <- nuts2_forest_multi |>drop_na() |>
      mutate(forestarea = replace_na(value,0)) |>
      filter(NUTS_ID == nuts2)

    if(nrow(nuts_forest)>0.01) {

      forest_budget <- pu_forest_multi_budget_data |>
        mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
        mutate_at(vars(starts_with("z")), ~.*m) |>
        dplyr::select(-c(m, NUTS_ID, pu)) |>
        as.matrix() |>
        replace_na(0)

      forest_threshold <- nuts_forest$forestarea[[1]]#/10*adj)

      p <- p |>
        add_linear_constraints(threshold = forest_threshold,
                               sense = ">=",
                               data = forest_budget)
      #print(nuts2)
    }
    #else(print("no forest multi constraint"))
  }
  #
  #   # # Forest prod

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
    nuts_forest <- nuts2_forest_prod |>drop_na() |>
      mutate(forestarea = replace_na(value,0)) |>
      filter(NUTS_ID == nuts2)

    if(nrow(nuts_forest)>0.01) {

      forest_budget <- pu_forest_prod_budget_data |>
        mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
        mutate_at(vars(starts_with("z")), ~.*m) |>
        dplyr::select(-c(m, NUTS_ID, pu)) |>
        as.matrix() |>
        replace_na(0)

      forest_threshold <- nuts_forest$forestarea[[1]]#/10*adj)
      p <- p |>
        add_linear_constraints(threshold = forest_threshold,
                               sense = ">=",
                               data = forest_budget)
      #print(nuts2)

    }
    # else(print("no forest prod constraint"))
  }


  features_weighted <- targs |>
    mutate(weight = ifelse(feature == "999999", (nrow(targs)-1)*2*0.5, weight)) |>
    # mutate(weight = ifelse(feature == "urban", 1000, weight)) |>
    dplyr::select(weight) |>
    as.matrix()

  targs_filtered <- targs |> dplyr::select(-weight)

  p_solve <- p |>
    add_proportion_decisions() |>
    add_manual_targets(targs_filtered) |>
    add_gurobi_solver(gap = 0.8, threads = cores,
                      verbose = TRUE, first_feasible = FALSE) |>
    add_manual_bounded_constraints(manual_bounded_constraints) |>
    add_feature_weights(features_weighted)
  print(nuts2_names)
  skip_to_next <- FALSE

  # Note that print(b) fails since b doesn't exist

  tryCatch(solve(p_solve), error = function(e) { skip_to_next <<- TRUE})
  feas$TF[n] <- skip_to_next
  feas$nuts[n] <- nuts2_names

  write_csv(feas, "data/feasibility-tests/nuts2_feasible_f455_adj.csv", append = FALSE)
  n = n+1
  print(feas)
  if(skip_to_next) { next }
}


# # ######### % for infeasible nuts2 ########
#
# p1 <- problem(x = pu,
#               features = feat,
#               zones = z,
#               cost_column = cost_columns,
#               rij = rij) |>
#   add_min_shortfall_objective(nrow(pu))
#
# # adj <- 0.9
# adj_list <- rev(seq(0,0.99, 0.01))
# l=1
#
#
# for(l in 76:length(adj_list)){
#   if(l>1){
#     infeas_nuts2 <- read_csv(paste0("data/feasibility-tests/nuts2_feasible_f455_", (l-1),".csv"))|>
#                                filter(TF == "TRUE") |>
#                                mutate(adj = adj_list[[l]])
#   }
#   if(l==1){
#   infeas_nuts2 <- read_csv("data/feasibility-tests/nuts2_feasible_f455.csv") |>
#     filter(TF == "TRUE") |>
#     mutate(adj = adj_list[[l]])
#   }
#   infeas2 <- infeas_nuts2
#   for(k in 1:length(infeas_nuts2$nuts)){
#     adj = adj_list[[l]]
#     nuts2_names <- infeas_nuts2$nuts[k] #setdiff(nuts2_all$NUTS_ID, small_jur$NUTS_ID)
#     # crop low
#     #nuts2_names <- nuts2_all$NUTS_ID #setdiff(nuts2_all$NUTS_ID
#     for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
#       nuts2 <- nuts2_names[[i]]
#
#       nuts_crop <- nuts2_crop_low |>
#         mutate(croparea = replace_na(value,0)) |>
#         filter(NUTS_ID == nuts2) |>
#         mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here
#
#       if(nrow(nuts_crop)>0) {
#
#         crop_budget <- pu_cropland_low_budget_data |>
#           mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
#           mutate_at(vars(starts_with("z")), ~.*m) |>
#           dplyr::select(-c(m, NUTS_ID, pu)) |>
#           as.matrix() |>
#           replace_na(0)
#
#         p <- p1 |>
#           add_linear_constraints(threshold = nuts_crop$croparea[[1]]*adj,#/10*adj,
#                                  sense = ">=",
#                                  data = crop_budget)
#         #print(nuts2)
#       }
#       #else(print("no crop low constraint"))
#     }
#
#     # crop med
#
#     for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
#       nuts2 <- nuts2_names[[i]]
#       nuts_crop <- nuts2_crop_med |>
#         mutate(croparea = replace_na(value,0)) |>
#         filter(NUTS_ID == nuts2) |>
#         mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here
#
#
#       if(nrow(nuts_crop)>0) {
#
#         crop_budget <- pu_cropland_med_budget_data |>
#           mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
#           mutate_at(vars(starts_with("z")), ~.*m) |>
#           dplyr::select(-c(m, NUTS_ID, pu)) |>
#           as.matrix() |>
#           replace_na(0)
#
#         p <- p |>
#           add_linear_constraints(threshold = nuts_crop$croparea[[1]]*adj, #/10*adj,
#                                  sense = ">=",
#                                  data = crop_budget)
#         #print(nuts2)
#       }
#       # else(print("no crop med constraint"))
#     }
#
#     #crop high
#
#     for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
#       nuts2 <- nuts2_names[[i]]
#       nuts_crop <- nuts2_crop_high |>
#         mutate(croparea = replace_na(value,0)) |>
#         filter(NUTS_ID == nuts2) |>
#         mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here
#
#       if(nrow(nuts_crop)>0) {
#
#         crop_budget <- pu_cropland_high_budget_data |>
#           mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
#           mutate_at(vars(starts_with("z")), ~.*m) |>
#           dplyr::select(-c(m, NUTS_ID, pu)) |>
#           as.matrix() |>
#           replace_na(0)
#
#         p <- p |>
#           add_linear_constraints(threshold = ifelse(nuts_crop$croparea[[1]]< 1, 0,
#                                                     nuts_crop$croparea[[1]]*adj),#/10*adj),
#                                  sense = ">=",
#                                  data = crop_budget)
#         #print(nuts2)
#       }
#       # else(print("no crop high constraint"))
#     }
#
#
#     #   # Pasture low
#     #
#     for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
#       nuts2 <- nuts2_names[[i]]
#       nuts_pasture <- nuts2_pasture_low |>
#         mutate(pasturearea = replace_na(value,0)) |>
#         filter(NUTS_ID == nuts2) |>
#         mutate(value = replace_na(value, 0),
#                pasturearea = replace_na(pasturearea,0))|>
#         mutate(pasturearea = ifelse(value < 1, 0, value))
#
#       if(nrow(nuts_pasture)>0) {
#
#         pasture_budget <- pu_pasture_low_budget_data |>
#           mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
#           mutate_at(vars(starts_with("z")), ~.*m) |>
#           dplyr::select(-c(m, NUTS_ID, pu)) |>
#           as.matrix() |>
#           replace_na(0)
#
#         p <- p |>
#           add_linear_constraints(threshold = nuts_pasture$pasturearea[[1]]*adj,#/10*adj,
#                                  sense = ">=",
#                                  data = pasture_budget)
#         #print(nuts2)
#       }
#       #else(print("no pasture low constraint"))
#     }
#
#     #   # Pasture high
#     #
#     for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
#       nuts2 <- nuts2_names[[i]]
#       nuts_pasture <- nuts2_pasture_high |>
#         mutate(pasturearea = replace_na(value,0)) |>
#         filter(NUTS_ID == nuts2) |>
#         mutate(value = replace_na(value, 0),
#                pasturearea = replace_na(pasturearea,0))|>
#         mutate(pasturearea = ifelse(value < 10, 0, value))
#
#       if(nrow(nuts_pasture)>0) {
#
#         pasture_budget <- pu_pasture_high_budget_data |>
#           mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
#           mutate_at(vars(starts_with("z")), ~.*m) |>
#           dplyr::select(-c(m, NUTS_ID, pu)) |>
#           as.matrix() |>
#           replace_na(0)
#
#         p <- p |>
#           add_linear_constraints(threshold = nuts_pasture$pasturearea[[1]]*adj,#/10*adj,
#                                  sense = ">=",
#                                  data = pasture_budget)
#         #print(nuts2)
#       }
#       #else(print("no pasture high constraint"))
#     }
#
#     #
#     # # # Forest multi
#     for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
#       nuts2 <- nuts2_names[[i]]
#       nuts_forest <- nuts2_forest_multi |>
#         mutate(forestarea = replace_na(value,0)) |>
#         filter(NUTS_ID == nuts2)
#
#       if(nrow(nuts_forest)>0.01) {
#
#         forest_budget <- pu_forest_multi_budget_data |>
#           mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
#           mutate_at(vars(starts_with("z")), ~.*m) |>
#           dplyr::select(-c(m, NUTS_ID, pu)) |>
#           as.matrix() |>
#           replace_na(0)
#
#         forest_threshold <- nuts_forest$forestarea[[1]]*adj
#
#         p <- p |>
#           add_linear_constraints(threshold = forest_threshold,
#                                  sense = ">=",
#                                  data = forest_budget)
#         #print(nuts2)
#       }
#       #else(print("no forest multi constraint"))
#     }
#     #
#     #   # # Forest prod
#
#     for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
#       nuts2 <- nuts2_names[[i]]
#       nuts_forest <- nuts2_forest_prod |>
#         mutate(forestarea = replace_na(value,0)) |>
#         filter(NUTS_ID == nuts2)
#
#       if(nrow(nuts_forest)>0.01) {
#
#         forest_budget <- pu_forest_prod_budget_data |>
#           mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
#           mutate_at(vars(starts_with("z")), ~.*m) |>
#           dplyr::select(-c(m, NUTS_ID, pu)) |>
#           as.matrix() |>
#           replace_na(0)
#
#         forest_threshold <- nuts_forest$forestarea[[1]]*adj
#         p <- p |>
#           add_linear_constraints(threshold = forest_threshold,
#                                  sense = ">=",
#                                  data = forest_budget)
#         #print(nuts2)
#
#       }
#       # else(print("no forest prod constraint"))
#     }
#     features_weighted <- targs |>
#       mutate(weight = ifelse(feature == "999999", (nrow(targs)-1)*2*0.5, weight)) |>
#       # mutate(weight = ifelse(feature == "urban", 1000, weight)) |>
#       dplyr::select(weight) |>
#       as.matrix()
#
#     targs_filtered <- targs |> dplyr::select(-weight)
#
#     p_solve <- p |>
#       add_proportion_decisions() |>
#       add_manual_targets(targs_filtered) |>
#       add_gurobi_solver(gap = 0.8, threads = cores,
#                         verbose = TRUE) |>
#       add_manual_bounded_constraints(manual_bounded_constraints) |>
#       add_feature_weights(features_weighted)
#     print(nuts2_names)
#     skip_to_next <- FALSE
#
#     # Note that print(b) fails since b doesn't exist
#
#     tryCatch(solve(p_solve), error = function(e) { skip_to_next <<- TRUE})
#     infeas2$TF[k] <- skip_to_next
#     infeas2$nuts[k] <- nuts2_names
#
#     write_csv(infeas2, paste0("data/feasibility-tests/nuts2_feasible_f455_", l, ".csv"), append = FALSE)
#     k = k+1
#     if(skip_to_next) { next }
#   }
# }
#
