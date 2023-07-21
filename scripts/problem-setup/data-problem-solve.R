### save problem solve files to one folder

pu_in_EU <- read_csv("data/formatted-data/pu_in_EU.csv")
write_csv(pu_in_EU, "data-formatted/pu_in_EU.csv")

pu <- read_fst("data/formatted-data/pu_data.fst")
write_fst(pu, "data-formatted/pu_data.fst")

rij <- read_fst("data/formatted-data/features_split.fst")
write_fst(rij, "data-formatted/features_split.fst", compress = 90)

z <- read_csv("data/formatted-data/zone_id.csv")
write_csv(z, "data-formatted/zone_id.csv")

names <- readRDS("data/SpeciesData/SDMNameMatching.rds")
writeRDS(names, "data-formatted/SpeciesData/SDMNameMatching.rds")

targs_existing <- read_csv("data/formatted-data/targets_split.csv")
write_csv(targs_existing, "data-formatted/targets_split.csv")

manual_bounded_constraints <- read_csv("data/formatted-data/manual_bounded_constraints_production_globiom_flex.csv") 
write_csv(manual_bounded_constraints, "data-formatted/manual_bounded_constraints_production_globiom_flex.csv")


R.utils::copyDirectory("data/formatted-data/linear_constraints/", "data-formatted/linear_constraints/")

R.utils::copyDirectory("data/feasibility-tests/", "data-formatted/feasibility-tests")

