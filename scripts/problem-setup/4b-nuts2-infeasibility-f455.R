### set production constraints ensuring feasibility
library(prioritizr)
library(fst)
library(gurobi)
library(tidyverse)
cores <- 7
fst::threads_fst(nr_of_threads = cores, reset_after_fork = NULL)  # Set number for fst

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

# updated manual constraints with globiom initial condition data
# z8 = wetlands (no restoration in this case)
# z5 = marine transitional (no restoration ever..)
# z13 = forestry production -- add extra flexibility for feasibility...
manual_bounded_constraints <- read_csv("data/formatted-data/manual_bounded_constraints_production_globiom_flex.csv") |>
  left_join(zones) |>
  dplyr::select(-c(zone, nuts2id)) |>
  rename(zone=name) |>
  dplyr::select(pu, zone, lower, upper) |>
  drop_na() |>
  mutate(upper = ifelse(zone == "z13", 1, upper)) |>
  mutate(upper = ifelse(zone == "z14", 1, upper)) |>
  mutate(upper = ifelse(zone == "z18", 1, upper)) |>
  mutate(lower = ifelse(zone == "z5", 0, lower)) |>
  mutate(upper = ifelse(zone == "z5", 0, upper))

# all the constraints
#nuts2_crop_low <- read_csv("data/formatted-data/linear_constraints/nuts_crop_low_95.csv")
nuts2_crop_low <- read_csv("data/formatted-data/linear_constraints/nuts_crop_low_55_adjusted.csv")
pu_cropland_low_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_crop_low_budget_data_55.csv")

#nuts2_crop_med <- read_csv("data/formatted-data/linear_constraints/nuts_crop_med_95.csv")
nuts2_crop_med <- read_csv("data/formatted-data/linear_constraints/nuts_crop_med_55_adjusted.csv")
pu_cropland_med_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_crop_med_budget_data_55.csv")

#nuts2_crop_high <- read_csv("data/formatted-data/linear_constraints/nuts_crop_high_95.csv")
nuts2_crop_high <- read_csv("data/formatted-data/linear_constraints/nuts_crop_high_55_adjusted.csv")
pu_cropland_high_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_crop_high_budget_data_55.csv")

#nuts2_pasture_high <- read_csv("data/formatted-data/linear_constraints/nuts_pasture_high_95.csv")
nuts2_pasture_high <- read_csv("data/formatted-data/linear_constraints/nuts_pasture_high_55_adjusted.csv")
pu_pasture_high_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_pasture_high_budget_data_55.csv")

#nuts2_pasture_low <- read_csv("data/formatted-data/linear_constraints/nuts_pasture_low_95.csv")
nuts2_pasture_low <- read_csv("data/formatted-data/linear_constraints/nuts_pasture_low_55_adjusted.csv")
pu_pasture_low_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_pasture_low_budget_data_55.csv")

nuts2_forest_multi <- read_csv("data/formatted-data/linear_constraints/nuts_forest_multi_55_adjusted.csv")
pu_forest_multi_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_forest_multi_budget_data_55.csv")

nuts2_forest_prod <- read_csv("data/formatted-data/linear_constraints/nuts_forest_prod_55_adjusted.csv")
pu_forest_prod_budget_data <- read_csv("data/formatted-data/linear_constraints/pu_forest_prod_budget_data_55.csv")


country_nuts <- pu_pasture_low_budget_data |> dplyr::select(pu, NUTS_ID)



########### pre set up basic problem  ##########
problem <- problem(x = pu,
             features = feat,
             zones = z,
             cost_column = cost_columns,
             rij = rij) |>
  add_min_shortfall_objective(nrow(pu))


p1 <- problem


######### solve by nuts2 all constraints ##########

nuts2_shp <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  mutate(country = substr(NUTS2, start = 1, stop = 2))
nuts2_all <- nuts2_shp |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nutsIDnum = seq(1:260)) |>
  filter(country != "UK")
nuts2_names_all <- nuts2_all$NUTS_ID #setdiff(nuts2_all$NUTS_ID
feas <- data.frame(
  TF = rep("TRUE",length(nuts2_names_all)),
  nuts = rep("A", length(nuts2_names_all)))

#feas <- read_csv("data/feasibility-tests/nuts2_feasible.csv")
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
      filter(NUTS_ID == nuts2) |>
      mutate(croparea = ifelse(value < 1, 0, value)) #ES63 is <1 PU and is driving an infeasibility here

    if(nrow(nuts_crop)>0) {

      crop_budget <- pu_cropland_low_budget_data |>
        mutate(m = ifelse(NUTS_ID == nuts2, 1, 0)) |>
        mutate_at(vars(starts_with("z")), ~.*m) |>
        dplyr::select(-c(m, NUTS_ID, pu)) |>
        as.matrix() |>
        replace_na(0)

      p <- p1 |>
        add_linear_constraints(threshold = nuts_crop$croparea[[1]]/10*adj,
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
    # else(print("no crop med constraint"))
  }

  #crop high

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    # else(print("no crop high constraint"))
  }


  #   # Pasture low
  #
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    #else(print("no pasture low constraint"))
  }

  #   # Pasture high
  #
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    #else(print("no pasture high constraint"))
  }

  #
  # # # Forest multi
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    #else(print("no forest multi constraint"))
  }
  #
  #   # # Forest prod

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    add_gurobi_solver(gap = 0.8, threads = cores, numeric_focus = FALSE,
                      verbose = TRUE) |>
    add_manual_bounded_constraints(manual_bounded_constraints) |>
    add_feature_weights(features_weighted)
  print(nuts2_names)
  skip_to_next <- FALSE

  # Note that print(b) fails since b doesn't exist

  tryCatch(solve(p_solve), error = function(e) { skip_to_next <<- TRUE})
  feas$TF[n] <- skip_to_next
  feas$nuts[n] <- nuts2_names

  write_csv(feas, "data/feasibility-tests/nuts2_feasible_f455.csv", append = FALSE)
  n = n+1
  print(feas)
  if(skip_to_next) { next }
}




######### 90% for infeasible nuts2 ########
#small_jur <- read_csv("small_jurisdictions.csv")
infeas_nuts2 <- read_csv("data/feasibility-tests/nuts2_feasible_f455.csv") |>
  filter(TF == "TRUE")

k=1
infeas2 <- data.frame(
  TF = rep("TRUE",length(infeas_nuts2$nuts)),
  nuts = rep("A", length(infeas_nuts2$nuts)))


p1 <- problem(x = pu,
              features = feat,
              zones = z,
              cost_column = cost_columns,
              rij = rij) |>
  add_min_shortfall_objective(nrow(pu))

adj <- 0.9

for(k in 1:length(infeas_nuts2$nuts)){
  nuts2_names <- infeas_nuts2$nuts[k] #setdiff(nuts2_all$NUTS_ID, small_jur$NUTS_ID)
  # crop low
  #nuts2_names <- nuts2_all$NUTS_ID #setdiff(nuts2_all$NUTS_ID
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]

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

      p <- p1 |>
        add_linear_constraints(threshold = nuts_crop$croparea[[1]]/10*adj,
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
    # else(print("no crop med constraint"))
  }

  #crop high

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
        add_linear_constraints(threshold = ifelse(nuts_crop$croparea[[1]]< 1, 0,
                                                  nuts_crop$croparea[[1]]/10*adj),
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
    #else(print("no pasture low constraint"))
  }

  #   # Pasture high
  #
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
    nuts_pasture <- nuts2_pasture_high |>
      mutate(pasturearea = replace_na(value,0)) |>
      filter(NUTS_ID == nuts2) |>
      mutate(value = replace_na(value, 0),
             pasturearea = replace_na(pasturearea,0))|>
      mutate(pasturearea = ifelse(value < 10, 0, value))

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
    #else(print("no pasture high constraint"))
  }

  #
  # # # Forest multi
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    #else(print("no forest multi constraint"))
  }
  #
  #   # # Forest prod

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    add_gurobi_solver(gap = 0.8, threads = cores, numeric_focus = FALSE,
                      verbose = TRUE) |>
    add_manual_bounded_constraints(manual_bounded_constraints) |>
    add_feature_weights(features_weighted)
  print(nuts2_names)
  skip_to_next <- FALSE

  # Note that print(b) fails since b doesn't exist

  tryCatch(solve(p_solve), error = function(e) { skip_to_next <<- TRUE})
  infeas2$TF[k] <- skip_to_next
  infeas2$nuts[k] <- nuts2_names

  write_csv(infeas2, "data/feasibility-tests/nuts2_feasible_f455_90.csv", append = FALSE)
  k = k+1
  if(skip_to_next) { next }
}

######### 80% for infeasible nuts2 ########
infeas_nuts2 <- read_csv("data/feasibility-tests/nuts2_feasible_f455_90.csv") |>
  filter(TF == "TRUE")

infeas2 <- data.frame(
  TF = rep("TRUE",length(infeas_nuts2$nuts)),
  nuts = rep("A", length(infeas_nuts2$nuts)))


p1 <- problem(x = pu,
              features = feat,
              zones = z,
              cost_column = cost_columns,
              rij = rij) |>
  add_min_shortfall_objective(nrow(pu))

adj <- 0.8
k=1

for(k in 1:length(infeas_nuts2$nuts)){
  nuts2_names <- infeas_nuts2$nuts[k] #setdiff(nuts2_all$NUTS_ID, small_jur$NUTS_ID)
  # crop low
  #nuts2_names <- nuts2_all$NUTS_ID #setdiff(nuts2_all$NUTS_ID
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]

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

      p <- p1 |>
        add_linear_constraints(threshold = nuts_crop$croparea[[1]]/10*adj,
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
    # else(print("no crop med constraint"))
  }

  #crop high

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
        add_linear_constraints(threshold = ifelse(nuts_crop$croparea[[1]]< 1, 0,
                                                  nuts_crop$croparea[[1]]/10*adj),
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
    #else(print("no pasture low constraint"))
  }

  #   # Pasture high
  #
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
    nuts_pasture <- nuts2_pasture_high |>
      mutate(pasturearea = replace_na(value,0)) |>
      filter(NUTS_ID == nuts2) |>
      mutate(value = replace_na(value, 0),
             pasturearea = replace_na(pasturearea,0))|>
      mutate(pasturearea = ifelse(value < 10, 0, value))

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
    #else(print("no pasture high constraint"))
  }

  #
  # # # Forest multi
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    #else(print("no forest multi constraint"))
  }
  #
  #   # # Forest prod

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    add_gurobi_solver(gap = 0.8, threads = cores, numeric_focus = FALSE,
                      verbose = TRUE) |>
    add_manual_bounded_constraints(manual_bounded_constraints) |>
    add_feature_weights(features_weighted)
  print(nuts2_names)
  skip_to_next <- FALSE

  # Note that print(b) fails since b doesn't exist

  tryCatch(solve(p_solve), error = function(e) { skip_to_next <<- TRUE})
  infeas2$TF[k] <- skip_to_next
  infeas2$nuts[k] <- nuts2_names

  write_csv(infeas2, "data/feasibility-tests/nuts2_feasible_f455_80.csv", append = FALSE)
  k = k+1
  if(skip_to_next) { next }
}

######### 70% for infeasible nuts2 ########
infeas_nuts2 <- read_csv("data/feasibility-tests/nuts2_feasible_f455_80.csv") |>
  filter(TF == "TRUE")

infeas2 <- data.frame(
  TF = rep("TRUE",length(infeas_nuts2$nuts)),
  nuts = rep("A", length(infeas_nuts2$nuts)))

p1 <- problem(x = pu,
              features = feat,
              zones = z,
              cost_column = cost_columns,
              rij = rij) |>
  add_min_shortfall_objective(nrow(pu))

adj <- 0.7
k=1

for(k in 1:length(infeas_nuts2$nuts)){
  nuts2_names <- infeas_nuts2$nuts[k] #setdiff(nuts2_all$NUTS_ID, small_jur$NUTS_ID)
  # crop low
  #nuts2_names <- nuts2_all$NUTS_ID #setdiff(nuts2_all$NUTS_ID
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]

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

      p <- p1 |>
        add_linear_constraints(threshold = nuts_crop$croparea[[1]]/10*adj,
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
    # else(print("no crop med constraint"))
  }

  #crop high

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
        add_linear_constraints(threshold = ifelse(nuts_crop$croparea[[1]]< 1, 0,
                                                  nuts_crop$croparea[[1]]/10*adj),
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
    #else(print("no pasture low constraint"))
  }

  #   # Pasture high
  #
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
    nuts_pasture <- nuts2_pasture_high |>
      mutate(pasturearea = replace_na(value,0)) |>
      filter(NUTS_ID == nuts2) |>
      mutate(value = replace_na(value, 0),
             pasturearea = replace_na(pasturearea,0))|>
      mutate(pasturearea = ifelse(value < 10, 0, value))

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
    #else(print("no pasture high constraint"))
  }

  #
  # # # Forest multi
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    #else(print("no forest multi constraint"))
  }
  #
  #   # # Forest prod

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    add_gurobi_solver(gap = 0.8, threads = cores, numeric_focus = FALSE,
                      verbose = TRUE) |>
    add_manual_bounded_constraints(manual_bounded_constraints) |>
    add_feature_weights(features_weighted)
  print(nuts2_names)
  skip_to_next <- FALSE

  # Note that print(b) fails since b doesn't exist

  tryCatch(solve(p_solve), error = function(e) { skip_to_next <<- TRUE})
  infeas2$TF[k] <- skip_to_next
  infeas2$nuts[k] <- nuts2_names

  write_csv(infeas2, "data/feasibility-tests/nuts2_feasible_f455_70.csv", append = FALSE)
  k = k+1
  if(skip_to_next) { next }
}


######### 80% for infeasible nuts2 ########
infeas_nuts2 <- read_csv("data/feasibility-tests/nuts2_feasible_f455_70.csv") |>
  filter(TF == "TRUE")

infeas2 <- data.frame(
  TF = rep("TRUE",length(infeas_nuts2$nuts)),
  nuts = rep("A", length(infeas_nuts2$nuts)))

p1 <- problem(x = pu,
              features = feat,
              zones = z,
              cost_column = cost_columns,
              rij = rij) |>
  add_min_shortfall_objective(nrow(pu))

adj <- 0.3
k=1

for(k in 1:length(infeas_nuts2$nuts)){
  nuts2_names <- infeas_nuts2$nuts[k] #setdiff(nuts2_all$NUTS_ID, small_jur$NUTS_ID)
  # crop low
  #nuts2_names <- nuts2_all$NUTS_ID #setdiff(nuts2_all$NUTS_ID
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]

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

      p <- p1 |>
        add_linear_constraints(threshold = nuts_crop$croparea[[1]]/10*adj,
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
    # else(print("no crop med constraint"))
  }

  #crop high

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
        add_linear_constraints(threshold = ifelse(nuts_crop$croparea[[1]]< 1, 0,
                                                  nuts_crop$croparea[[1]]/10*adj),
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
    #else(print("no pasture low constraint"))
  }

  #   # Pasture high
  #
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
    nuts_pasture <- nuts2_pasture_high |>
      mutate(pasturearea = replace_na(value,0)) |>
      filter(NUTS_ID == nuts2) |>
      mutate(value = replace_na(value, 0),
             pasturearea = replace_na(pasturearea,0))|>
      mutate(pasturearea = ifelse(value < 10, 0, value))

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
    #else(print("no pasture high constraint"))
  }

  #
  # # # Forest multi
  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    #else(print("no forest multi constraint"))
  }
  #
  #   # # Forest prod

  for(i in 1:length(nuts2_names)) { #length(nuts2_AT_names)
    nuts2 <- nuts2_names[[i]]
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
    add_gurobi_solver(gap = 0.8, threads = cores, numeric_focus = FALSE,
                      verbose = TRUE) |>
    add_manual_bounded_constraints(manual_bounded_constraints) |>
    add_feature_weights(features_weighted)
  print(nuts2_names)
  skip_to_next <- FALSE

  # Note that print(b) fails since b doesn't exist

  tryCatch(solve(p_solve), error = function(e) { skip_to_next <<- TRUE})
  infeas2$TF[k] <- skip_to_next
  infeas2$nuts[k] <- nuts2_names

  write_csv(infeas2, "data/feasibility-tests/nuts2_feasible_f455_30.csv", append = FALSE)
  k = k+1
  if(skip_to_next) { next }
}


#### combine into dataframe ############

f455_feas <- read_csv("data/feasibility-tests/nuts2_feasible_f455.csv") |>
  filter(TF == FALSE) |>
  mutate(adj = 1)
f455_90 <- read_csv("data/feasibility-tests/nuts2_feasible_f455_90.csv") |>
  filter(TF == FALSE) |>
  mutate(adj = 0.9)
f455_80 <- read_csv("data/feasibility-tests/nuts2_feasible_f455_80.csv") |>
  filter(TF == FALSE) |>
  mutate(adj = 0.8)
f455_70 <- read_csv("data/feasibility-tests/nuts2_feasible_f455_70.csv") |>
  filter(TF == FALSE) |>
  mutate(adj = 0.7)
f455_30 <- read_csv("data/feasibility-tests/nuts2_feasible_f455_30.csv") |>
  filter(TF == FALSE) |>
  mutate(adj = 0.3)
f455_feas_adj <- rbind(f455_feas, f455_90, f455_80, f455_70, f455_30)
write_csv(f455_feas_adj, "data/feasibility-tests/f455_feas_adj.csv")


########### Solve entire problem ############

nuts2_shp <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  mutate(country = substr(NUTS2, start = 1, stop = 2))

nuts2_all <- nuts2_shp |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nutsIDnum = seq(1:260)) |>
  filter(country != "UK")

nuts2_names <- nuts2_all$NUTS_ID #setdiff(nuts2_all$NUTS_ID

# deal with infeasible nuts2

infeas_nuts2 <- read_csv("data/feasibility-tests/nuts2_feasible.csv") |>
  filter(TF == "TRUE")

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

    p <- p1 |>
      add_linear_constraints(threshold = nuts_crop$croparea[[1]]/10*adj,
                             sense = ">=",
                             data = crop_budget)
    print(nuts2)
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
    print(nuts2)
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
    print(nuts2)
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
    print(nuts2)
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
    print(nuts2)
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
    print(nuts2)
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
    print(nuts2)

  }
  else(print("no forest prod constraint"))
}


# solve problem
features_weighted <- targs |>
  mutate(weight = ifelse(feature == "999999", (nrow(targs)-1)*2*0.5, weight)) |>
  # mutate(weight = ifelse(feature == "urban", 1000, weight)) |>
  dplyr::select(weight) |>
  as.matrix()

targs_filtered <- targs |> dplyr::select(-weight)

p_solve <- p |>
  add_proportion_decisions() |>
  add_manual_targets(targs_filtered) |>
  add_gurobi_solver(gap = 0.8, threads = cores, numeric_focus = FALSE,
                    verbose = TRUE) |>
  add_manual_bounded_constraints(manual_bounded_constraints) |>
  add_feature_weights(features_weighted)

test_feas_f455 <- solve(p_solve)



#### Test outcome differences ######
sum(test_feas_f455$solution_1_z18)/sum(nuts2_crop_high$value)*10
sum(test_feas_f455$solution_1_z13)/sum(nuts2_forest_prod$value)*10
sum(test_feas_f455$solution_1_z14)/sum(nuts2_pasture_high$value)*10

