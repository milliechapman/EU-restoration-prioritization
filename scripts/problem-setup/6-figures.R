library(tidyverse)
library(ggthemes)
library(sf)
library(raster)
library(fasterize)
library(terra)
require(viridisLite)
library(MetBrewer)
library(ggridges)
library(rasterVis)
library(fst)
#devtools::install_github("yutannihilation/ggsflabel")
#library(ggsflabel)
library(RStoolbox)
library(scico)
library(geomtextpath)
library(ggpubr)
library(ggpattern)
library(viridis)
library(ggsci)
library(formattable)
#rm(list = ls())

######################### set up files ##############
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
plotting_data <- read_fst("data/plotting/plotting_data_ref_f455.fst")

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

rep_IC <- read_csv("data/solutions/representation_IC.csv")
rep <- read_fst("data/solutions/representation_REF_f455_scenarios_proportional.fst")
rep_norest <- read_fst("data/solutions/representation_norest.fst")

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

rep_formatted <- rep |>
  dplyr::select(-c(column_label, rep, target, shortfall)) |>
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

rep_norest_formatted <- rep_norest |>
  dplyr::select(-c(column_label, rep, target, shortfall)) |>
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


biodiv_ic <- rep_IC_formatted |>
  group_by(target_met, scenario) |> count() |>
  pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met = yes/(yes+no)*100) |>
  mutate(scenario == "2020") |>
  mutate(scenario == ifelse(scenario == "HN", "High \n Nature", scenario))


biodiv_sol <- rep_formatted |>
  filter(country_contraint == "FLEX") |>
  group_by(target_met, scenario, future, carbon_weight) |>
  count() |>
  pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met = yes/(yes+no)*100) |>
  mutate(future = ifelse(future == "f455", "Fit for 55", "BAU")) |>
  mutate(scenario = ifelse(scenario == "HN", "High \n nature", scenario))

biodiv_noNRL <-rep_norest_formatted |>
  group_by(future, target_met) |>
  count() |>
  pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met = yes/(yes+no)*100) |>
  mutate(future = ifelse(future == "f455", "Fit for 55", "BAU")) |>
  mutate(label = paste0("2030 no NRL (", future, ")"))

biodiv_sol$scenario
nrl_biodiv_plot <- ggboxplot(biodiv_sol, x = "scenario", y = "perc_met",
                                color = "scenario",
                                add = "jitter") +
  geom_hline(data = biodiv_ic,
             aes(yintercept = perc_met), color = "#2c456b", linetype = "dashed", lwd = 1) +
  geom_text(data = biodiv_ic, aes(x = 1.5, y = perc_met, label = "2020"),
            color = "#2c456b", hjust = 0.55,vjust = -0.4, size = 3.5) +
  geom_hline(data = biodiv_noNRL,
             aes(yintercept = perc_met, color = "#4779c4"), linetype = "dotdash", lwd =1) +
  geom_text(data = biodiv_noNRL, aes(x = 1.5, y = perc_met, label = label, color = "#4779c4"),
             hjust = 0.5, vjust = 2, size = 3.5) +
  # stat_compare_means(aes(label = ..p.signif..))  +
  theme_bw() +
  scale_color_manual(values = c( "#4779c4", "grey", "black", "#4779c4", "black")) +
  ylim(0,62) +
  theme(legend.position = "none") +
  labs(x = "", y = "% of species targets met") +
  theme(axis.text.x = element_text(size = 11)) +
  facet_wrap(~future, nrow = 1)#+

ggsave("figures/updated/nrl_biodiv_plot.png", nrl_biodiv_plot, width = 3.7, height = 4, dpi = 300)

c_perc <- rep |>
  left_join(ic_join) |>
  filter(country_contraint == "FLEX") |>
  filter(spp == "carbon") |>
  filter(rep_ic > 0) |>
  mutate(rep_diff = (rep - rep_ic)/rep_ic) |>
  group_by(spp, scenario, carbon_weight, country_contraint, future) |>
  summarise(rep_diff = rep_diff*100) |>
  arrange(-rep_diff)

c_perc_nonrl <- rep_norest |>
  left_join(ic_join) |>
  filter(spp == "carbon") |>
  filter(rep_ic > 0) |>
  mutate(rep_diff = (rep - rep_ic)/rep_ic) |>
  group_by(spp, scenario, carbon_weight, country_contraint, future) |>
  summarise(rep_diff = rep_diff*100) |>
  arrange(-rep_diff)

nrl_carbon_plot <- ggboxplot(c_perc, x = "scenario", y = "rep_diff",
                             color = "scenario",
                             add = "jitter") +
  geom_hline(data = biodiv_ic,
             aes(yintercept = 0), color = "#2c456b", linetype = "dashed", lwd = 1) +
  geom_text(data = biodiv_ic, aes(x = 1.5, y = 0, label = scenario),
            color = "#2c456b", hjust = 0.55,vjust = -0.4, size = 3.5) +
  geom_hline(data = c_perc_nonrl,
             aes(yintercept = rep_diff, color = future), linetype = "dotdash", lwd =1) +
  geom_text(data = c_perc_nonrl, aes(x = 1.5, y = rep_diff, label = future, color = future),
            hjust = 0.35, vjust = 1.4, size = 3.5) +
  # stat_compare_means(aes(label = ..p.signif..))  +
  theme_bw() +
  scale_color_manual(values = c( "grey", "#4779c4", "black", "#4779c4")) +
  theme(legend.position = "none") +
  labs(x = "", y = "% increase in land carbon \n (relative to 2020 conditions)") +
  theme(axis.text.x = element_text(size = 12)) +
  facet_wrap(~future, nrow = 1)

ggsave("figures/updated/nrl_carbon_plot.png", nrl_carbon_plot, width = 3.7, height = 4, dpi = 300)

### for later figures, just f455
rep_norest_formatted <- rep_norest_formatted |> filter(future == "f455")

##################### abstract statistics ##############################

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
  group_by(spp, scenario, carbon_weight, country_contraint, future) |>
  summarise(rep_diff = rep_diff*100) |>
  arrange(-rep_diff)

## conservation status
b_status <-rep |>
  left_join(ic_join) |>
  filter(species != 999999) |>
  mutate(target_met_ic = ifelse(shortfall_perc_ic < 0, 0, 1),
         target_met = ifelse(shortfall_perc < 0, 0, 1)) |>
  filter(target_met_ic<1) |>
  mutate(improved = (target_met-target_met_ic)) |>
  group_by(spp, scenario, carbon_weight, country_contraint, future) |>
  summarise(improved = sum(improved)/length(unique(species))) |>
  arrange(-improved)



  # rep |>
  # left_join(ic_join) |>
  # filter(species != 999999) |>
  # mutate(spp = substr(species, start = 1, stop =nchar(species)-4)) |>
  # mutate(n_spp = n_distinct(spp)) |>
  # mutate(target_met_ic = ifelse(shortfall_perc_ic < 0, 0, 1),
  #        target_met = ifelse(shortfall_perc < 0, 0, 1)) |>
  # group_by(spp, n_spp, scenario, carbon_weight, country_contraint, future) |>
  # summarise(IC_met_n = sum(target_met_ic, na.rm = T),
  #           sol_met_n = sum(target_met, na.rm = T)) |>
  # mutate(improved = ifelse((sol_met_n - IC_met_n)<1,0,1)) |>
  # mutate(worse = ifelse((sol_met_n - IC_met_n)<0,1,0)) |>
  # #filter(target_met_ic<1) |>
  # #mutate(improved = (target_met-target_met_ic)) |>
  # group_by(scenario, carbon_weight, country_contraint, future, n_spp) |>
  # summarise(improved = sum(improved)/n_spp,
  #           worse = sum(worse)/n_spp) |>
  # arrange(-improved) |>
  # mutate(diff = IC_shortfall_n - improved) |> #/length(unique(species))) |>
  # arrange(-diff)

## burden sharing
rep |>
  left_join(ic_join) |>
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
  mutate(change = rep - rep_ic) |>
  mutate(change = ifelse(change >0, 1, 0)) |>
  group_by(scenario, carbon_weight, country_contraint, taxon_id, future) |>
  summarise(change = sum(change)) |>
  pivot_wider(names_from = country_contraint, values_from = change) |>
  mutate(diff = ifelse((`EVEN` - `UNEVEN`)<0,1,0),
         n_spp = length(unique(taxon_id))) |>
  group_by(scenario, carbon_weight, future) |>
  summarise(diff_spp = sum(diff, na.rm = TRUE)/mean(n_spp,  na.rm = TRUE)) |>
  arrange(-diff_spp)


rep |>
  filter(species == 999999) |>
  dplyr::select(scenario, carbon_weight, country_contraint, future, species, rep) |>
  pivot_wider(names_from = country_contraint, values_from = rep) |>
  left_join(ic_join) |>
  mutate(diff = (`EVEN` - `UNEVEN`)/rep_ic*-100) |>
  arrange(diff)


################### figure 4 - all scenarios #################

write_csv(b_status |> filter(carbon_weight>0), "data/biodiv_results.csv")
write_csv(c_perc |> filter(carbon_weight>0), "data/carbon_results.csv")

b_status |> filter(carbon_weight == 0.5)
a <- b_status |> filter(carbon_weight>0) |>
  ggplot(aes(x = carbon_weight, y = improved*100,
                            col = country_contraint)) +
  geom_line(aes(linetype = scenario)) +
  geom_point(aes()) +
  theme_classic() +
  scale_color_manual(values = c("grey2", "#0072B2", "orange")) +
  scale_alpha_manual(values = c(1, 0.6)) +
  labs(x = "carbon weight \n (relative to biodiversity)",
       y = "% of species with improved \n conservation status") +
  facet_wrap(~future)

b <- c_perc |> filter(carbon_weight>0) |>
  ggplot(aes(x = carbon_weight, y = rep_diff,
                          col = country_contraint)) +
  geom_line(aes(linetype = scenario)) +
  geom_point(aes()) +
  theme_classic() +
  scale_color_manual(values = c("grey2", "#0072B2", "orange")) +
  labs(x = "carbon weight \n (relative to biodiversity)",
       y = "% increase in land carbon stock \n (compared to initial conditions)") +
  facet_wrap(~future)

fig4 <- ggarrange(a,b, ncol=1, nrow=2, common.legend = TRUE, legend="right")

ggsave("figures/updated/combined_scenarios.png", fig4, width = 7, height = 3.5)

unique(c_perc$carbon_weight)
c_perc |> ungroup() |>
  dplyr::select(-spp) |>
  left_join(b_status) |>
  filter(carbon_weight >0) |>
  ggplot(aes(x = improved*100, y = rep_diff,
             col = country_contraint)) +
  geom_line(aes(linetype = future)) +
  geom_point(aes()) +
  theme_classic() +
  scale_color_manual(values = c("orange", "#0072B2", "grey")) +
  labs(x = "% of species with improved \n conservation status",
       y = "% increase in land carbon stock \n (compared to initial conditions)") +
  facet_wrap(~scenario, scales = "free")

c_perc |> ungroup() |>
  dplyr::select(-spp) |>
  left_join(b_status) |>
  filter(carbon_weight >0) |>
  ggplot(aes(x = improved*100, y = rep_diff,
             col = country_contraint, alpha = scenario)) +
  geom_line(aes(linetype = future)) +
  geom_point(aes()) +
  theme_classic() +
  scale_color_manual(values = c("orange", "#0072B2", "brown")) +
  scale_alpha_manual(values = c(1, 0.6)) +
  labs(x = "% of species with improved \n conservation status",
       y = "% increase in land carbon stock \n (compared to initial conditions)")

############ SI pareto front all scenarios (shortfall) ###############
pareto_all <- rep |> filter(target >0) |>
  mutate(shortfall_perc = ifelse(shortfall_perc >0,0, shortfall_perc)) |>
  drop_na(shortfall_perc) |>
  #filter(country_contraint == TRUE)|>
  mutate(type = ifelse(species == 99999, "carbon", "biodiversity")) |>
  group_by(carbon_weight, spp, scenario, country_contraint, future) |>
  summarise(shortfall = mean(shortfall_perc, na.rm = TRUE)*(-1))

pareto_all_plot  <- pareto_all |>
  pivot_wider(names_from = spp, values_from = shortfall) |>
  filter(carbon_weight!= 0) |>
  #filter(scenario == "Baseline") |>
  ggplot() +
  geom_line(aes(x = biodiv, y = carbon, col =country_contraint,
                linetype = future)) +
  geom_point(aes(x = biodiv, y = carbon,
                 color = country_contraint, #as.factor(carbon_weight),
                 alpha = future)) +
  scale_color_manual(values = c("#E69F00", "#0072B2", "black")) +
  scale_alpha_manual(values = c(1,0.8)) +
  #xlim(0.4,0.41)+ ylim(0.52,0.57)+
  theme_classic()   +
  #theme(legend.position = "none") +
  #geom_point(data = shortfall_IC_pareto, aes(x = biodiv, y = carbon))
  #scale_y_continuous(limits = c(0.45, 0.55), breaks = seq(0.45, 0.55, 0.02)) +
  #scale_x_continuous(limits = c(0.15, 0.20), breaks = seq(0.15, 0.20, 0.01)) +
  labs(x = "mean habitat shortfall", y = "carbon shortfall") +
  facet_wrap(~scenario, scales = "free")

ggsave(filename = 'figures/updated/pareto-all.png', plot = pareto_all_plot, width = 5, height = 4, dpi = 400)



################### figure 1 - conceptual ##########################
solution <- read_csv("data/solutions/sol/NUTSadj/sol_carbon_0.5_restoration_0.141_production_NUTS2adj_country_FLEX_wetlands_TRUE_onlyrestoration_FALSE_scenario_Baseline_proportional_gurobi_f455.csv") |>
  dplyr::select(id, solution_1_z1:solution_1_z26)
colnames(solution) <- c("id", (zone_id$zone))
pu_adjustment <- read_csv("data/formatted-data/pu_adjustment.csv") |>
  rename(id = pu)

solution_table <- pu_in_EU |> rename(id = pu) |> #plot_data |>
  left_join(PU_template) |>
  dplyr::select(-id) |>
  rename(id = EU_id) |>
  left_join(solution) |>
  pivot_longer(-c(nuts2id:geometry)) |>
  left_join(pu_adjustment) |>
  mutate(value = ifelse(name == "RiversLakes_natural_conserve", value-diff, value))

solution_table_plot <- solution_table |>
  mutate(zone = name) |>
  separate(name, c('maes_label', 'intensity', 'action'), sep = "_")

#nuts2_crop_low <- read_csv("data-formatted/linear_constraints/nuts_crop_low_95

fig1_solution_bar <- as_tibble(solution_table_plot) |>
  group_by(zone, maes_label, intensity, action) |>
  summarise(area = sum(value)) |>
  mutate(action = ifelse(action == "lockin", "production", action)) |>
  ungroup() |> filter(area >0) |>
  mutate(intensity = ifelse(intensity == "primary", "natural", intensity)) |>
  mutate(maes_label = ifelse(maes_label %in% c("RiversLakes", "MarineTransitional"), "RiversLakesMarine", maes_label)) |>
  mutate(name = paste0(maes_label, " \n (", intensity, ")")) |>
  ggplot(aes(x = reorder(name, -area), y = area/10, fill = action)) + geom_bar(stat="identity", width = 0.7) +
  theme_classic() +
  labs(x = element_blank(), y = "area (1000 km^2)") +
  theme(legend.position = "none",#c(0.8,0.8),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, hjust =1)) +
  scale_fill_manual(values = c(met.brewer(name="Kandinsky",n=4,type="continuous"))[c(1,3,2)]) #+ coord_flip()

ggsave("figures/updated/fig1_solution_bar.png", dpi = 300, width = 7, height = 3)
solution_raster <- fasterize(st_as_sf(solution_table_plot), raster(PU_plot), field = "value", by = "zone")
solution_raster <- rast(solution_raster)
names(solution_raster)

restore <- c(1,2,6,3,4,5,7,8,9,10)
produce <- c(12:18, 26)
conserve <- c(19:25)

colors <- c("grey", met.brewer(name="VanGogh3",n=20,type="continuous"))
myPal <- colors
myTheme <- rasterTheme(region = myPal)


rest_agg <- app(solution_raster[[restore]], sum)
png("figures/updated/rest_agg.png", width = 5, height = 5, units = "in", res = 400)
plot(rest_agg, col = colors)
dev.off()

conserve_agg <- app(solution_raster[[conserve]], sum)
png("figures/updated/conserve_agg.png", width = 5, height = 5, units = "in", res = 400)
plot(conserve_agg, col = colors)
dev.off()

prod_agg <- app(solution_raster[[produce]], sum)
png("figures/updated/prod_agg.png", width = 5, height = 5, units = "in", res = 400)
plot(prod_agg, col = colors)
dev.off()

wtland_rest <- app(solution_raster[[8]], sum)
png("figures/updated/wtland_rest.png", width = 5, height = 5, units = "in", res = 400)
plot(wtland_rest, col = colors)
dev.off()

spr_cons <- app(solution_raster[[21]], sum)
png("figures/updated/spr_cons.png", width = 5, height = 5, units = "in", res = 400)
plot(spr_cons, col = colors)
dev.off()


########## Table ###################

test <- plotting_data |> filter(carbon == 0.5) |>
  mutate(type = ifelse(intensity == "natural" & maes_label !=
                         "Wetlands", "nature",
                       ifelse(maes_label == "Wetlands", "rewetting",
                       "deintensification"))) |>
  group_by(action, country_TF, scenario, future, type) |>
  summarise(area = sum(value))


rest_table <- test |> ungroup() |>
  filter(action == "restore") |>
  dplyr::select(-action) |>
  group_by(country_TF, scenario, future) |>
  mutate(area_restored = sum(area),
         percent_restored = sum(area)/41046) |>
  ungroup() |>
  pivot_wider(names_from = type, values_from = area)

write_csv(rest_table, "data/rest_table_results.csv")

library(formattable)
library("htmltools")
library("webshot")
plain_formatter <- formatter("span")
plain_formatter(c(1, 2, 3))
width_formatter <- formatter("span",
                             style = x ~ style(width = suffix(x, "px")))
width_formatter(c(10, 11, 12))
sign_formatter <- formatter("span",
                            style = x ~ style(color = ifelse(x > 0, "green",
                                                             ifelse(x < 0, "red", "black"))))
sign_formatter(c(-1, 0, 1))

t <- rest_table |>
  mutate(`% total area restored` = round(percent_restored*100, 1),
         `% deintensification` = round((deintensification/area_restored)*percent_restored*100,1),
         `% nature` = round((nature/area_restored)*percent_restored*100,1),
         `% rewetting` = round((rewetting/area_restored)*percent_restored*100,1)) |>
  mutate(country_TF = ifelse(country_TF == "FLEX", "Flexible",
                             ifelse(country_TF == "EVEN", "Even ", "Unconstrained")) ) |>
  mutate(future = ifelse(future == "f455.csv", "Fit for 55", "BAU")) |>
  rename(`burden sharing` = country_TF,
         `restoration scenario` = scenario,
         `production constraints` = future) |>
  dplyr::select(-c(area_restored, nature, rewetting, deintensification, percent_restored)) |>
  #dplyr::select(-c(area_restored)) |>
  formattable(list(
    `% total area restored` =color_tile("grey","#654321"),
    `% nature` = color_bar("darkgreen"),
    `% rewetting` = color_bar("lightblue"),
    `% deintensification` = color_bar("#FFCC00")))

export_formattable <- function(f, file, width = "100%", height = NULL,
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}
export_formattable(t,"figures/updated/table.png")

#################### burden_sharing #######################
########## pareto
pareto_burden <- rep |> filter(target >0) |>
  mutate(shortfall_perc = ifelse(shortfall_perc >0,0, shortfall_perc)) |>
  #drop_na(shortfall_perc) |>
  #filter(country_contraint == TRUE)|>
  mutate(type = ifelse(species == 99999, "carbon", "biodiversity")) |>
  group_by(carbon_weight, spp, scenario, country_contraint) |>
  summarise(shortfall = mean(shortfall_perc, na.rm = TRUE)*(-1))

colors = c(met.brewer(name="Kandinsky",n=11,type="continuous"))

shortfall_IC_pareto <- rep_IC |>
  filter(target >0) |>
  mutate(shortfall_perc = ifelse(shortfall_perc >0,0, shortfall_perc)) |>
  drop_na(shortfall_perc) |>
  #filter(country_contraint == TRUE)|>
  mutate(type = ifelse(species == 99999, "carbon", "biodiversity")) |>
  group_by(carbon_weight, spp, scenario, country_contraint) |>
  summarise(shortfall = mean(shortfall_perc, na.rm = TRUE)*(-1)) |>
  pivot_wider(names_from = spp, values_from = shortfall)

pareto_burden_plot <- pareto_burden |>
  pivot_wider(names_from = spp, values_from = shortfall) |>
  filter(scenario == "Baseline") |>
  ggplot() +
  geom_line(aes(x = biodiv, y = carbon, linetype = country_contraint)) +
  geom_point(aes(x = biodiv, y = carbon,
                 color = country_contraint, #as.factor(carbon_weight),
                 alpha = country_contraint, size = 1.5)) +
  scale_color_manual(values = c("grey2", "grey", "black")) +
  scale_alpha_manual(values = c(1,0.8, 0.2)) +
  #xlim(0.4,0.41)+ ylim(0.52,0.57)+
  theme_classic()   +
  #theme(legend.position = "none") +
  #geom_point(data = shortfall_IC_pareto, aes(x = biodiv, y = carbon))
  #scale_y_continuous(limits = c(0.45, 0.55), breaks = seq(0.45, 0.55, 0.02)) +
  #scale_x_continuous(limits = c(0.15, 0.20), breaks = seq(0.15, 0.20, 0.01)) +
  labs(x = "mean habitat shortfall", y = "carbon shortfall")

ggsave(filename = 'figures/updated/burden-sharing-pareto.png', plot = pareto_burden_plot, width = 4, height = 4, dpi = 400)


## pareto new metrics
pareto_burden_plot_new <- b_status |>
  filter(carbon_weight %in% c("0.1","0.5", "0.9","1.5", "2")) |>
  ungroup() |>
  dplyr::select(-spp) |> filter(carbon_weight > 0) |>
  filter(scenario == "Baseline",
         future == "f455") |>
  #filter(country_contraint == "FALSE") |>
  left_join(c_perc |> ungroup() |> dplyr::select(-spp)) |>
  ggplot(aes(x = improved*100, y = rep_diff, col = country_contraint)) +
  geom_line(aes(linetype = future), lwd = 1) +
  geom_point() +
  #geom_point(aes(alpha =scenario)) +
  theme_classic() +
  xlim(24.5, 26.2) +
  #ylim(11, 40) +
  #xlim(52, 62) +
  #scale_color_manual(values = c("grey2", "grey")) +
  scale_color_manual(values = c( "grey", "grey2", "darkred")) +
  #scale_color_manual(values = c("orange", "blue", "black")) +
  #scale_alpha_manual(values = c(1, 0)) +
  labs(x = "% of species with improved \n conservation status",
       y = "% increase in land carbon stock \n (compared to current conditions)")
  #facet_wrap(~scenario, scale = "free")

ggsave(filename = 'figures/updated/burden-sharing-pareto-new-s1.png', plot = pareto_burden_plot_new, width = 6, height = 4, dpi = 400)

ggsave(filename = 'figures/updated/burden-sharing-pareto-new-s2.png', plot = pareto_burden_plot_new, width = 6, height = 4, dpi = 400)

ggsave(filename = 'figures/updated/burden-sharing-pareto-new.png', plot = pareto_burden_plot_new, width = 6, height = 4, dpi = 400)

## map of difference
country_true <- plotting_data |> filter(action == "restore",
                                        carbon == "0.5",
                                        country_TF == "EVEN",
                                        scenario == "Baseline",
                                        future == "f455.csv") |>
  group_by(id, country_TF) |>
  summarise(value = sum(value)) |>
  left_join(PU_template_EU)

country_true <- fasterize(st_as_sf(country_true), raster, field = "value")

country_false <- plotting_data |> filter(action == "restore",
                                         carbon == "0.5",
                                         country_TF == "FLEX",
                                         scenario == "Baseline",
                                         future == "f455.csv") |>
  group_by(id, country_TF) |>
  summarise(value = sum(value)) |>
  left_join(PU_template_EU)

country_false <- fasterize(st_as_sf(country_false), raster, field = "value")

burden_map <- country_false - country_true

library(RStoolbox)
fig3b_plot <- ggR(burden_map, layer = 1, maxpixels = 1e10,forceCat = FALSE,
                  geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_scico(direction = 1,palette = "vik",midpoint = 0, na.value = '#FFFFFF',
                   guide = guide_colorbar(title = "Difference",ticks = FALSE) ) +
  theme(legend.position = c(.05,.15),legend.background = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 23, hjust = .5),plot.subtitle = element_text(hjust = .5)) +
  labs(x = "", y = "")

ggsave(filename = 'figures/updated/burden-sharing-map.png', plot = fig3b_plot, width = 5, height = 4, dpi = 400)

### country summary

fig3d_data <- plotting_data |>
  filter(scenario == "Baseline",
         future == "f455.csv") |>
  left_join(country) |>
  group_by(country, carbon, country_TF) |>
  mutate(area_country = sum(value)) |> ungroup() |>
  group_by(country, carbon, country_TF, area_country, action) |>
  summarise(amount = sum(value)) |> ungroup() |>
  mutate(perc = amount/area_country) |> dplyr::select(-c(amount, area_country)) |>
  pivot_wider(names_from = country_TF, values_from = perc) |>
  filter(action == "restore") |>
  mutate(diff = `FLEX`-`EVEN`) |> drop_na() |>
  pivot_longer(-c(country:action, diff))

fig3d_box <- fig3d_data |>
  group_by(country) |> mutate(fill_value = mean(diff)) |>
  ggplot(aes(x = reorder(country, diff), y = diff*100, fill = fill_value)) +
  geom_boxplot()  +
  scale_fill_scico(direction = 1,palette = "vik",midpoint = 0, na.value = '#FFFFFF',
                   guide = guide_colorbar(title = "Difference",ticks = FALSE) ) +  coord_flip() + theme_classic() + #theme(legend.position = "none") +
  labs(x = element_blank()) + labs(y = "difference restoration priority \n (% land area; II-I all scenarios)") + geom_hline(yintercept = 0, linetype="dashed") +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11))

ggsave("figures/updated/burden-sharing-countries.png", fig3d_box, width = 4.5, height = 8, dpi = 300)

## land cover
colors = c(met.brewer(name="Kandinsky",n=4,type="continuous"))[c(4,2,1)]
fig3c_data <- plotting_data |>
  filter(carbon == 0.5) |>
  filter(action == "restore",
         scenario == "Baseline",
         future == "f455.csv") |>
  group_by(carbon, country_TF, maes_label, action) |>
  summarise(amount = sum(value)) |> ungroup() |>
  ggplot(aes(x = reorder(maes_label,amount), y = amount, fill = as.factor(carbon), alpha = country_TF)) + geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  coord_flip() +
  #facet_grid(rows = vars(carbon)) +
  scale_alpha_manual(values = c(0.5,0.9,1)) +
  scale_fill_manual(values = colors[4]) +
  labs(x = element_blank(), y = "area (100 km^2)") +
  theme_minimal() + theme(legend.position = "none")

ggsave("figures/updated/burden-sharing-LC.png",fig3c_data, width= 3, height = 3, dpi = 300)

## Biodiversity

lj_biodiv_ic <- rep_IC_formatted |> filter(scenario == "IC") |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Vascular plants", "Non-vascular plants"), "Plants", speciesgroup)) |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Mammals", "Birds", "Plants", "Amphibians", "Reptiles"),speciesgroup, "Other")) |>
  group_by(speciesgroup, target_met) |>
  count() |> pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met_IC = yes/(yes+no)) |> dplyr::select(-c(no, yes))

norest_burden <- rep_norest_formatted |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Vascular plants", "Non-vascular plants"), "Plants", speciesgroup)) |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Mammals", "Birds", "Plants", "Amphibians", "Reptiles"),speciesgroup, "Other")) |>
  group_by(speciesgroup, target_met) |>
  count() |> pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met_norest = yes/(yes+no)) |> dplyr::select(-c(no, yes))

spp <- rep_formatted |> filter(scenario == "Baseline",
                               future == "f455", carbon_weight == "0.5") |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Vascular plants", "Non-vascular plants"), "Plants", speciesgroup)) |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Mammals", "Birds", "Plants", "Amphibians", "Reptiles"),speciesgroup, "Other")) |>
  dplyr::select(speciesgroup, taxon_id) |> unique() |>
  group_by(speciesgroup) |> count() |> rename(nspp = n) |> ungroup()

spp_list <- unique(rep_formatted$speciesname)
write_csv(as.data.frame(spp_list), "figures/updated/spp_list.csv")

burden_biodiv <- rep_formatted |>
  filter(scenario == "Baseline",
         future == "f455") |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Vascular plants", "Non-vascular plants"), "Plants", speciesgroup)) |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Mammals", "Birds", "Plants", "Amphibians", "Reptiles"),speciesgroup, "Other")) |>
 # group_by(speciesgroup, country_contraint, carbon_weight) |># |> mutate(nspp = n()) |> ungroup()|>
  group_by(speciesgroup, country_contraint, target_met, carbon_weight) |>
  count() |> pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met = yes/(yes+no)) |>
  #left_join(lj_biodiv_ic) |>
  #mutate(perc_met_rel = (perc_met - perc_met_IC)*100) |>
   left_join(spp) |>
  mutate(country_contraint = ifelse(country_contraint == "EVEN", "I",
                                    ifelse(country_contraint == "FLEX", "II", "III"))) |>
  mutate(perc_met = perc_met*100) |>
  mutate(speciesgroup = paste0(speciesgroup, "\n n spp=", nspp))

lj_biodiv_ic$speciesgroup <- unique(burden_biodiv$speciesgroup)
lj_biodiv_ic <- lj_biodiv_ic |>
  mutate(label = ifelse(speciesgroup != "Plants\n n spp=140", "2020", ""))

norest_burden$speciesgroup <- unique(burden_biodiv$speciesgroup)
norest_burden <- norest_burden |>
  mutate(label = ifelse(speciesgroup != "Plants\n n spp=140", "2030 no NRL", ""))

burden_biodiv_plot <- ggboxplot(burden_biodiv, x = "country_contraint", y = "perc_met",
                 color = "country_contraint",
                 add = "jitter") +
  geom_hline(data = lj_biodiv_ic,
             aes(yintercept = perc_met_IC*100), color = "#2c456b", linetype = "dashed", lwd = 1) +
  geom_text(data = lj_biodiv_ic, aes(x = 1.5, y = perc_met_IC*100, label = label),
            color = "#2c456b", hjust = 0.55,vjust = -0.4, size = 3.5) +
  geom_hline(data = norest_burden,
             aes(yintercept = perc_met_norest*100), color = "#4779c4", linetype = "dotdash", lwd =1) +
  geom_text(data = norest_burden, aes(x = 1.5, y = perc_met_norest*100, label = label),
            color = "#4779c4", hjust = 0.35, vjust = 1.4, size = 3.5) +
   # stat_compare_means(aes(label = ..p.signif..))  +
  theme_bw() +
  scale_color_manual(values = c( "grey", "grey2", "darkred")) +
  ylim(0,62) +
  theme(legend.position = "none") +
  labs(x = "", y = "percent of species targets met") +
  theme(axis.text.x = element_text(size = 12)) +
  facet_wrap(~speciesgroup, nrow = 1)#+
  #scale_y_continuous(breaks=equal_breaks(n=10, s=0.1))

ggsave("figures/updated/burden-sharing-spp.png",burden_biodiv_plot, width=8.2 , height = 3.3, dpi = 300)

###################### restoration scenarios #############

## pareto new metrics
pareto_restoration <- b_status |> ungroup() |> dplyr::select(-spp) |> filter(carbon_weight >0) |>
  filter(country_contraint == "FLEX",
         future == "f455"
         ) |>
  left_join(c_perc |> ungroup() |> dplyr::select(-spp)) |>
  mutate(carbon_weight = as.numeric(carbon_weight)) |>
  #filter(carbon_weight > 0.1) |>
  ggplot(aes(x = improved*100, y = rep_diff)) +
  geom_line(aes(linetype = scenario, alpha = future), lwd = 1) +
  geom_point(aes(col = carbon_weight), alpha = 0.9, size = 4) +
  theme_classic() +
  scale_color_gradientn(colors = met.brewer("OKeeffe2")) +
  #scale_color_manual(values = c("grey2", "grey")) +
  scale_alpha_manual(values = c(1, 0.4)) +
  labs(x = "% of species with improved \n conservation status",
       y = "% increase in land carbon stock \n (compared to 2020 conditions)")


png("figures/updated/pareto_rest.png", width = 5.5, height = 4, units = "in", res = 400)
pareto_restoration
dev.off()

## pareto new metrics
pareto_restoration2 <- b_status |> ungroup() |> dplyr::select(-spp) |>
  filter(carbon_weight >0) |>
  filter(future == "f455") |>
  left_join(c_perc |> ungroup() |> dplyr::select(-spp)) |>
  ggplot(aes(x = improved*100, y = rep_diff)) +
  geom_line(aes(linetype = scenario, alpha = country_contraint), lwd = 1) +
  geom_point(aes(col = as.numeric(carbon_weight)), alpha = 0.9, size = 4) +
  theme_classic() +
  scale_color_gradientn(colors = met.brewer("OKeeffe2")) +
  #scale_color_manual(values = c("grey2", "grey")) +
  scale_alpha_manual(values = c(1, 0.6, 0.3)) +
  labs(x = "% of species with improved \n conservation status",
       y = "% increase in land carbon stock \n (compared to current conditions)")

ggsave("figures/pareto_restoration.png",pareto_restoration2, width = 5.5, height = 4, dpi = 400)

############ restoration scenario barchart ######

basline <- plotting_data |> filter(action == "restore",
                                   scenario == "Baseline",
                                   country_TF == "FLEX",
                                   future == "f455.csv") |>
  mutate(type = ifelse(maes_label== "Cropland", "Cropland \n de-intensification",
                       ifelse(maes_label %in% c("Pasture"), "Pasture \n de-intensification",
                              ifelse(maes_label %in% c("WoodlandForest"), "Forestry \n de-intensification",
                                     "Nature \n landscape")))) |>
  group_by(type, carbon, scenario, future) |>
  summarise(value = sum(value))

hn <-plotting_data |> filter(action == "restore",
                             scenario == "HN",
                             country_TF == "FLEX",
                             future == "f455.csv") |>
  mutate(type = ifelse(maes_label== "Cropland", "Cropland \n de-intensification",
                       ifelse(maes_label %in% c("Pasture"), "Pasture \n de-intensification",
                              ifelse(maes_label %in% c("WoodlandForest"), "Forestry \n de-intensification",
                                     "Nature \n landscape")))) |>
  group_by(type, carbon, scenario, future) |>
  summarise(value = sum(value))


hn |>
  filter(type == "Nature \n landscape") |>
  mutate(perc = value/40000*100) |>
  arrange(-value)

basline |> filter() |> arrange(-value)
hn_baseline_box <- basline |> bind_rows(hn) |> filter(future == "f455.csv",
                                                      carbon >0) |>
  ggplot() + geom_boxplot(aes(y = value/10, x = type, col = scenario, )) +
  theme_classic() +
  coord_flip() +
  theme(axis.title.y =element_blank(),
                          axis.text.y = element_text(angle = 90))+
  labs(y = "area (1000 km^2)") +
  scale_color_manual(values = c("grey2", "grey"))

hn_baseline_box <- basline |> bind_rows(hn) |> filter(future == "f455.csv",
                                                      carbon >0,
                                                      carbon == 0.5) |>
  ggplot() + geom_col(aes(y = value/10, x = type, fill = scenario), position =position_dodge(.7),
                      width = 0.6) +
  theme_classic() +
  coord_flip() +
  theme(axis.title.y =element_blank(),
        axis.text.y = element_text(angle = 90))+
  labs(y = "area (1000 km^2)") +
  scale_fill_manual(values = c( "grey", "grey2"))

ggsave("figures/updated/hn_baseline_box.png", hn_baseline_box,
       height = 4, width = 3.75)

########### restoration scenarios plots #########

basline <- plotting_data |> filter(action == "restore",
                                      scenario == "Baseline",
                                      country_TF == "FLEX",
                                      carbon == "0.5") |>
  mutate(type = ifelse(maes_label== "Cropland", "Cropland de-intensification",
                       ifelse(maes_label %in% c("Pasture"), "Pasture de-intensification",
                              ifelse(maes_label %in% c("WoodlandForest"), "Forestry de-intensification",
                                     "Nature landscape")))) |>
  group_by(id, type) |>
  summarise(value = sum(value)) |>
  left_join(PU_template_EU)

basline <- fasterize(st_as_sf(basline), raster, field = "value", by = "type")

plot(basline)

hn <-plotting_data |> filter(action == "restore",
                                     scenario == "HN",
                                     country_TF == "FLEX",
                                     carbon == 0.5) |>
  mutate(type = ifelse(maes_label== "Cropland", "Cropland de-intensification",
                       ifelse(maes_label %in% c("Pasture"), "Pasture de-intensification",
                              ifelse(maes_label %in% c("WoodlandForest"), "Forestry de-intensification",
                                     "Nature landscape")))) |>
  group_by(id, type) |>
  summarise(value = sum(value)) |>
  left_join(PU_template_EU)

hn <- fasterize(st_as_sf(hn), raster, field = "value", by = "type")

plot(hn)

crop_diff <- basline[[1]] - hn[[1]]

cd <- ggR(crop_diff, layer = 1, maxpixels = 1e10,forceCat = FALSE,
          geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_scico(direction = 1,palette = "vik",midpoint = 0, na.value = '#FFFFFF',
                   guide = guide_colorbar(title = "crop Difference",ticks = FALSE) ) +
  theme(legend.position = c(.05,.15),legend.background = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 23, hjust = .5),plot.subtitle = element_text(hjust = .5)) +
  labs(x = "", y = "")

ggsave(filename = 'figures/updated/cd.png', plot = cd, width = 16, height = 14, dpi = 400)

pasture_diff <- basline[[2]] - hn[[2]]

pd <- ggR(pasture_diff, layer = 1, maxpixels = 1e10,forceCat = FALSE,
          geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_scico(direction = 1,palette = "vik",midpoint = 0, na.value = '#FFFFFF',
                   guide = guide_colorbar(title = "pasture Difference",ticks = FALSE) ) +
  theme(legend.position = c(.05,.15),legend.background = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 23, hjust = .5),plot.subtitle = element_text(hjust = .5)) +
  labs(x = "", y = "")

ggsave(filename = 'figures/updated/pd.png', plot = pd, width = 16, height = 14, dpi = 400)

nature_diff <- basline[[3]] - hn[[3]]

nd <- ggR(nature_diff, layer = 1, maxpixels = 1e10,forceCat = FALSE,
          geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_scico(direction = 1,palette = "vik",midpoint = 0, na.value = '#FFFFFF',
                   guide = guide_colorbar(title = "nature Difference",ticks = FALSE) ) +
  theme(legend.position = c(.05,.15),legend.background = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 23, hjust = .5),plot.subtitle = element_text(hjust = .5)) +
  labs(x = "", y = "")

ggsave(filename = 'figures/updated/nd.png', plot = nd, width = 16, height = 14, dpi = 400)

forest_diff <- basline[[4]] - hn[[4]]

fd <- ggR(forest_diff, layer = 1, maxpixels = 1e10,forceCat = FALSE,
          geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_scico(direction = 1,palette = "vik",midpoint = 0, na.value = '#FFFFFF',
                   guide = guide_colorbar(title = "forest Difference",ticks = FALSE) ) +
  theme(legend.position = c(.05,.15),legend.background = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 23, hjust = .5),plot.subtitle = element_text(hjust = .5)) +
  labs(x = "", y = "")

ggsave(filename = 'figures/updated/fd.png', plot = fd, width = 16, height = 14, dpi = 400)


######## rest biodiv ###########

lj_biodiv_ic <- rep_IC_formatted |> filter(scenario == "IC") |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Vascular plants", "Non-vascular plants"), "Plants", speciesgroup)) |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Mammals", "Birds", "Plants", "Amphibians", "Reptiles"),speciesgroup, "Other")) |>
  group_by(speciesgroup, target_met) |>
  count() |> pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met_IC = yes/(yes+no)) |> dplyr::select(-c(no, yes))

rest_biodiv <- rep_formatted |>
  filter(country_contraint == "FLEX",
         future == "f455") |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Vascular plants", "Non-vascular plants"), "Plants", speciesgroup)) |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Mammals", "Birds", "Plants", "Amphibians", "Reptiles"),speciesgroup, "Other")) |>
  # group_by(speciesgroup, country_contraint, carbon_weight) |># |> mutate(nspp = n()) |> ungroup()|>
  group_by(speciesgroup, scenario, target_met, carbon_weight) |>
  count() |> pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met = yes/(yes+no)) |>
  left_join(lj_biodiv_ic) |>
  mutate(perc_met_rel = (perc_met - perc_met_IC)*100) |>
  left_join(spp) |>
  #mutate(country_contraint = ifelse(country_contraint == "EVEN", "I",
                                  #  ifelse(country_contraint == "FLEX", "II", "III"))) |>
  mutate(perc_met = round(perc_met*100,0)) |>
  mutate(speciesgroup = paste0(speciesgroup, "\n n spp=", nspp)) |>
  mutate(scenario = ifelse(scenario == "Baseline", "Ba", scenario))


norest_burden <- rep_norest_formatted |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Vascular plants", "Non-vascular plants"), "Plants", speciesgroup)) |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Mammals", "Birds", "Plants", "Amphibians", "Reptiles"),speciesgroup, "Other")) |>
  group_by(speciesgroup, target_met) |>
  count() |> pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met_norest = yes/(yes+no)) |> dplyr::select(-c(no, yes))
lj_biodiv_ic$speciesgroup <- unique(burden_biodiv$speciesgroup)
lj_biodiv_ic <- lj_biodiv_ic |>
  mutate(label = ifelse(speciesgroup != "Plants\n n spp=140", "2020", ""))

norest_burden$speciesgroup <- unique(burden_biodiv$speciesgroup)
norest_burden <- norest_burden |>
  mutate(label = ifelse(speciesgroup != "Plants\n n spp=140", "2030 no NRL", ""))

# Box plot facetted
rest_biodiv_plot <- ggboxplot(rest_biodiv, x = "scenario", y = "perc_met",
                                color = "scenario",
                                add = "jitter",
                                facet.by = "speciesgroup", nrow = 1,
                                short.panel.labs = TRUE) +
  # stat_compare_means(aes(label = ..p.signif..))  +
  theme_bw() +
  geom_hline(data = lj_biodiv_ic,
             aes(yintercept = perc_met_IC*100), color = "#2c456b", linetype = "dashed", lwd = 1) +
  geom_text(data = lj_biodiv_ic, aes(x = 1.5, y = perc_met_IC*100, label = label),
            color = "#2c456b", hjust = 1,vjust = -0.4, size = 3.5) +
  geom_hline(data = norest_burden,
             aes(yintercept = perc_met_norest*100), color = "#4779c4", linetype = "dotdash", lwd =1) +
  geom_text(data = norest_burden, aes(x = 1.5, y = perc_met_norest*100, label = label),
            color = "#4779c4", hjust = 0.4, vjust = 1.4, size = 3.5) +
  scale_color_manual(values = c( "grey", "grey2", "darkred")) +
  ylim(0,62) +
  theme(legend.position = "none") +
  labs(x = "", y = "percent of species targets met") +
  theme(axis.text.x = element_text(size = 12)) #+
#scale_y_continuous(breaks=equal_breaks(n=10, s=0.1))

ggsave("figures/updated/rest-spp.png",rest_biodiv_plot, width=8 , height = 3.2, dpi = 300)

########## pareto ######
pareto_scenario <- rep |> filter(target >0) |>
  mutate(shortfall_perc = ifelse(shortfall_perc >0,0, shortfall_perc)) |>
  drop_na(shortfall_perc) |>
  #filter(country_contraint == TRUE)|>
  mutate(type = ifelse(species == 99999, "carbon", "biodiversity")) |>
  group_by(carbon_weight, spp, scenario, country_contraint) |>
  summarise(shortfall = mean(shortfall_perc, na.rm = TRUE)*(-1))

colors = c(met.brewer(name="Kandinsky",n=11,type="continuous"))

pareto_scenario_plot <- pareto_scenario |>
  pivot_wider(names_from = spp, values_from = shortfall) |>
  filter(country_contraint == "FLEX") |>
  ggplot() +
  geom_line(aes(x = biodiv, y = carbon, linetype = scenario)) +
  geom_point(aes(x = biodiv, y = carbon,
                 color = scenario, #as.factor(carbon_weight),
                 alpha = scenario, size = 1.5)) +
  scale_color_manual(values = c("grey2", "grey")) +
  scale_alpha_manual(values = c(1,0.8)) +
  #xlim(0.4,0.41)+ ylim(0.52,0.57)+
  theme_classic()   +
  #theme(legend.position = "none") +
  #geom_point(data = shortfall_IC_pareto, aes(x = biodiv, y = carbon))
  #scale_y_continuous(limits = c(0.45, 0.55), breaks = seq(0.45, 0.55, 0.02)) +
  #scale_x_continuous(limits = c(0.15, 0.20), breaks = seq(0.15, 0.20, 0.01)) +
  labs(x = "mean habitat shortfall", y = "carbon shortfall")

ggsave(filename = 'figures/updated/rest-scenario-pareto.png', plot = pareto_scenario_plot, width = 4, height = 4, dpi = 400)


##################### all scenarios simple outcomes ##################
full_ic <- rep_IC_formatted |> filter(scenario == "IC") |>
  group_by(spp, target_met) |>
  count() |> pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met_IC = ifelse(spp == "biodiv", yes/(yes+no), shortfall_perc)) |>
           dplyr::select(-c(no, yes))

## carbon percentage change
c_perc <- rep |>
  left_join(ic_join) |>
  filter(spp == "carbon") |>
  filter(rep_ic > 0) |>
  mutate(rep_diff = (rep - rep_ic)/rep_ic) |>
  group_by(spp, scenario, carbon_weight, country_contraint, future) |>
  summarise(rep_diff = rep_diff*100) |>
  arrange(-rep_diff)


full_benefits <- rep_formatted |>
  filter(country_contraint == "FLEX") |>
  # group_by(speciesgroup, country_contraint, carbon_weight) |># |> mutate(nspp = n()) |> ungroup()|>
  group_by(scenario, future, spp, target_met, carbon_weight) |>
  count() |> pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met = yes/(yes+no)) |>
  left_join(lj_biodiv_ic) |>
  mutate(perc_met_rel = (perc_met - perc_met_IC)*100) |>
  #mutate(country_contraint = ifelse(country_contraint == "EVEN", "I",
  #  ifelse(country_contraint == "FLEX", "II", "III"))) |>
  mutate(perc_met = round(perc_met*100,0))


norest_burden <- rep_norest_formatted |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Vascular plants", "Non-vascular plants"), "Plants", speciesgroup)) |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Mammals", "Birds", "Plants", "Amphibians", "Reptiles"),speciesgroup, "Other")) |>
  group_by(speciesgroup, target_met) |>
  count() |> pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met_norest = yes/(yes+no)) |> dplyr::select(-c(no, yes))
lj_biodiv_ic$speciesgroup <- unique(burden_biodiv$speciesgroup)
lj_biodiv_ic <- lj_biodiv_ic |>
  mutate(label = ifelse(speciesgroup != "Plants\n n spp=140", "2020", ""))

norest_burden$speciesgroup <- unique(burden_biodiv$speciesgroup)
norest_burden <- norest_burden |>
  mutate(label = ifelse(speciesgroup != "Plants\n n spp=140", "2030 no NRL", ""))

# Box plot facetted
rest_biodiv_plot <- ggboxplot(rest_biodiv, x = "scenario", y = "perc_met",
                              color = "scenario",
                              add = "jitter",
                              facet.by = "speciesgroup", nrow = 1,
                              short.panel.labs = TRUE) +
  # stat_compare_means(aes(label = ..p.signif..))  +
  theme_bw() +
  geom_hline(data = lj_biodiv_ic,
             aes(yintercept = perc_met_IC*100), color = "#2c456b", linetype = "dashed", lwd = 1) +
  geom_text(data = lj_biodiv_ic, aes(x = 1.5, y = perc_met_IC*100, label = label),
            color = "#2c456b", hjust = 1,vjust = -0.4, size = 3.5) +
  geom_hline(data = norest_burden,
             aes(yintercept = perc_met_norest*100), color = "#4779c4", linetype = "dotdash", lwd =1) +
  geom_text(data = norest_burden, aes(x = 1.5, y = perc_met_norest*100, label = label),
            color = "#4779c4", hjust = 0.4, vjust = 1.4, size = 3.5) +
  scale_color_manual(values = c( "grey", "grey2", "darkred")) +
  ylim(0,62) +
  theme(legend.position = "none") +
  labs(x = "", y = "percent of species targets met \n relative to initial conditions") +
  theme(axis.text.x = element_text(size = 12)) #+
#scale_y_continuous(breaks=equal_breaks(n=10, s=0.1))

ggsave("figures/updated/rest-spp.png",rest_biodiv_plot, width=8 , height = 3.2, dpi = 300)


############## bioregion #######################
# fig4E_legend <- plotting_data |> filter(carbon == "0.3",
#                                         country_TF == "FLEX") |>
#   #group_by(BR_ID) |>
#   #mutate(area_br = sum(value)) |>
#   left_join(BR_lookup |> rename(BR_ID = bioregion)) |>
#   filter(action == "restore") |> drop_na() |>
#   mutate(type = ifelse(maes_label== "Cropland", "Cropland \n de-intensification",
#                        ifelse(maes_label %in% c("Pasture"), "Pasture \n de-intensification",
#                               ifelse(maes_label %in% c("WoodlandForest"), "Forestry \n de-intensification",
#                                      "Natural landscape")))) |>
#   mutate(type = factor(type, levels = c("Natural landscape",
#                                         "Forestry \n de-intensification",
#                                         "Cropland \n de-intensification",
#                                         "Pasture \n de-intensification"))) |>
#   group_by(code, scenario) |>
#   mutate(area = sum(value)) |>
#   group_by(code, type, area, scenario) |>
#   summarise(value = sum(value)) |>
#   mutate(code = fct_rev(code)) |>
#   ggplot(aes(x = scenario, y = value/area, fill = type, alpha = scenario)) + facet_grid(rows = vars(code), scales = "free") + coord_flip() +
#   geom_bar(stat = "identity", width = 0.8) + #facet_wrap(~scenario) +
#   scale_fill_manual(values = c(met.brewer(name="Moreau",n=5,type="continuous"))) + #c("darkolivegreen3", "chocolate4", "snow2", "khaki3"))+
#   #scale_pattern_manual(values = c("stripe","none")) +
#   scale_alpha_manual(values = c(1,0.5)) +
#   theme_classic() + labs(y = "% restoration priority", x = element_blank())
# ggsave("figures/updated/fig4E_legend.png",fig4E_legend,  height = 4, width = 4, dpi = 300)
#
#
# fig4E <- fig4E_legend +
#   theme(legend.position = "none", legend.title = element_blank()) + coord_flip() +
#   labs(x = element_blank(), y = "area restored (100 km^2)")
#
# ggsave("figures/updated/fig4E.png",fig4E,  height = 3, width = 3, dpi = 300)
#

## Figure 2D:
##
# bioregion_plot <- plotting_data |> filter(action == "restore",
#                                           carbon == "0.3",
#                                           country_TF == "TRUE") |>
#   group_by(id, BR_ID) |> mutate(value = sum(value))
#
# fig2D <- bioregion_plot_sum |> left_join(bioregion) |> drop_na() |>
#   ggplot() + geom_sf(aes(geometry = geometry, fill = code), lwd = 0) +
#   scale_fill_manual(values = c(met.brewer(name="Navajo",n=8,type="continuous")))+
#   geom_sf_label_repel(aes(label = code),
#                       force = 100, nudge_x = -2, seed = 10, size = 2) +
#   theme_map() + theme(legend.position = "none")
#
# fig4D <- bioregion_plot_sum |> left_join(bioregion) |> drop_na() |>
#   ggplot() + geom_sf(aes(geometry = geometry, fill = code), lwd = 0) +
#   scale_fill_manual(values = c(met.brewer(name="Navajo",n=8,type="continuous")))+
#   #geom_sf_label_repel(aes(label = code),
#   #  force = 100, nudge_x = -2, seed = 10, size = 2) +
#   theme_map() + theme(legend.position = "none")
#
# ggsave("figures/updated/fig4D.png",fig4D,  height = 3, width = 3, dpi = 300)
# ## Figure 2E:
# ##
#
#
#
# data_pie <- plotting_data |>
#   filter(carbon == "0.3",
#          country_TF == TRUE) |>
#   drop_na() |>
#   group_by(BR_ID) |> mutate(area = sum(value)) |>
#   left_join(bioregion) |> #mutate(code = pre_2012) |>
#   drop_na() |>
#   filter(action == "restore") |>
#   mutate(type = ifelse(maes_label== "Cropland", "Cropland \n de-intensification",
#                        ifelse(maes_label %in% c("Pasture"), "Pasture \n de-intensification",
#                               ifelse(maes_label %in% c("WoodlandForest"), "Forestry \n de-intensification",
#                                      "Natural landscape")))) |>
#   mutate(type = factor(type, levels = c("Natural landscape",
#                                         "Forestry \n de-intensification",
#                                         "Cropland \n de-intensification",
#                                         "Pasture \n de-intensification"))) |>
#   group_by(scenario) |> mutate(area = sum(value)) |>
#   group_by(type, area, scenario) |>
#   summarise(value = sum(value)) |>
#   mutate(perc = value/area)
#
# # Compute the cumulative percentages (top of each rectangle)
# data_pieI <- data_pie |> filter(scenario == "I")
# data_pieI$ymax = cumsum(data_pieI$perc)
# # Compute the bottom of each rectangle
# data_pieI$ymin = c(0, head(data_pieI$ymax, n=-1))
#
# # Make the plot
# pieI <- ggplot(data_pieI, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
#   geom_rect() +
#   coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
#   xlim(c(2, 4)) +   theme_void() +
#   scale_fill_manual(values = c(met.brewer(name="Moreau",n=5,type="continuous")))
#
# ggsave("figures/updated/pie_s1.png", pieI, width = 3, height = 2, dpi = 300)
#
#
# # Compute the cumulative percentages (top of each rectangle)
# data_pieII <- data_pie |> filter(scenario == "II")
# data_pieII$ymax = cumsum(data_pieII$perc)
# # Compute the bottom of each rectangle
# data_pieII$ymin = c(0, head(data_pieII$ymax, n=-1))
#
# # Make the plot
# pieII <- ggplot(data_pieII, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type, alpha =0.6)) +
#   geom_rect() +
#   coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
#   xlim(c(2, 4)) +   theme_void() +
#   scale_fill_manual(values = c(met.brewer(name="Moreau",n=5,type="continuous")))
#
# ggsave("figures/updated/pie_s2.png", pieII, width = 3, height = 2, dpi = 300)
#
# fig4F <- as_tibble(plot_table[[2]]) |>
#   group_by(BR_ID) |>
#   mutate(n = sum(value)) |> ungroup() |>
#   group_by(action, BR_ID, n) |>
#   summarise(value = sum(value)) |>
#   mutate(perc_value = value/n) |>
#   filter(action == "restore") |> drop_na() |>
#   left_join(bioregion) |>
#   ggplot(aes(x = fct_rev(code), y = value, fill = code)) +
#   geom_bar(stat = "identity", width = 0.7) +
#   #geom_hline(yintercept = 0.2, linetype = "dashed") +
#   scale_fill_manual(values = c(met.brewer(name="Navajo",n=8,type="continuous"))) +
#   theme_classic() +
#   theme(legend.position = "none") + coord_flip() +
#   labs(x = element_blank(), y = "area restored (100 km^2)")
# ggsave("figures/updated/fig4F.png",fig4F,  height = 4, width = 3, dpi = 300)
#
#
#
#
#
#
# ############### make pareto fronts #################
#
# carbon_low <- plotting_data |> filter(action == "restore",
#                                       carbon == "0",
#                                       country_TF == "TRUE") |>
#   group_by(id, country_TF) |>
#   summarise(value = sum(value)) |>
#   left_join(PU_template_EU)
#
# carbon_low <- fasterize(st_as_sf(carbon_low), raster, field = "value")
#
# carbon_high <- plotting_data |> filter(action == "restore",
#                                        carbon == "1",
#                                        country_TF == "TRUE") |>
#   group_by(id, country_TF) |>
#   summarise(value = sum(value)) |>
#   left_join(PU_template_EU)
#
# carbon_high <- fasterize(st_as_sf(carbon_high), raster, field = "value")
#
# fig2b_data <- carbon_high - carbon_low
#
# fig2b_plot <- ggR(fig2b_data, layer = 1, maxpixels = 1e10,forceCat = FALSE,
#                   geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
#   ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
#   scale_fill_scico(direction = 1,palette = "vik",midpoint = 0, na.value = '#FFFFFF',
#                    guide = guide_colorbar(title = "Difference",ticks = FALSE) ) +
#   theme(legend.position = c(.05,.15),legend.background = element_rect(fill = 'transparent'),
#         plot.title = element_text(size = 23, hjust = .5),plot.subtitle = element_text(hjust = .5)) +
#   labs(x = "", y = "")
#
# ggsave(filename = 'figures/updated/fig2B.png', plot = fig2b_plot, width = 16, height = 14, dpi = 400)
#
#
#
# fig2c <- plotting_data |> filter(country_TF == TRUE) |>
#   filter(action == "restore") |>
#   mutate(type = ifelse(maes_label== "Cropland", "Cropland de-intensification",
#                        ifelse(maes_label %in% c("Pasture"), "Pasture de-intensification",
#                               ifelse(maes_label %in% c("WoodlandForest"), "Forestry de-intensification",
#                                      "Nature landscape")))) |>
#   group_by(maes_label, type, carbon) |>
#   summarise(value = sum(value)) |>
#   ggplot(aes(x = reorder(maes_label, value), y = value, fill = as.factor(carbon))) +
#   geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
#   coord_flip() +
#   scale_fill_manual(values = c(met.brewer(name="Kandinsky",n=11,type="continuous"))) +
#   labs(x = element_blank(), y = "area (100 km^2)") +
#   theme_minimal() + theme(legend.position = "none")
#
# ggsave("figures/updated/fig2c.png", fig2c, width = 3, height = 3, dpi = 300)


##################### restoration only ########################
colors <- c("grey", met.brewer(name="VanGogh3",n=20,type="continuous"))

#### A - overestimating potential
#### B - even not overestimating potential
solution <- read_csv("data/solutions/sol/NUTSadj/sol_carbon_0.5_restoration_0.141_production_NUTS2adj_country_FLEX_wetlands_TRUE_onlyrestoration_FALSE_scenario_Baseline_proportional_gurobi_f455.csv") |>
  dplyr::select(id, solution_1_z1:solution_1_z26)
colnames(solution) <- c("id", (zone_id$zone))

solution_table <- pu_in_EU |> rename(id = pu) |> #plot_data |>
  left_join(PU_template) |>
  dplyr::select(-id) |>
  rename(id = EU_id) |>
  left_join(solution) |>
  pivot_longer(-c(nuts2id:geometry))

solution_table_plot <- solution_table |>
  mutate(zone = name) |>
  separate(name, c('maes_label', 'intensity', 'action'), sep = "_")

solution_raster <- fasterize(st_as_sf(solution_table_plot), raster(PU_plot), field = "value", by = "zone")
solution_raster <- rast(solution_raster)
names(solution_raster)

restore <- c(1,2,6,3,4,5,7,8,9,10)

rest_agg <- app(solution_raster[[restore]], sum)

solution_niave <- read_csv("data/solutions/sol/SI-sol/sol_carbon_0.5_restoration_0.141_production_NUTS2adj_country_FLEX_wetlands_TRUE_onlyrestoration_TRUE_scenario_Baseline_onlyrest_highs_f455.csv")|>
  dplyr::select(id, solution_1_z1:solution_1_z26)
colnames(solution_niave) <- c("id", (zone_id$zone))

solution_table_niave <- pu_in_EU |> rename(id = pu) |> #plot_data |>
  left_join(PU_template) |>
  dplyr::select(-id) |>
  rename(id = EU_id) |>
  left_join(solution_niave) |>
  pivot_longer(-c(nuts2id:geometry))

solution_table_plot_n <- solution_table_niave |>
  mutate(zone = name) |>
  separate(name, c('maes_label', 'intensity', 'action'), sep = "_")

rest_bar_niave <- as_tibble(solution_table_plot_n) |>
  group_by(zone, maes_label, intensity, action) |>
  summarise(area_naive = sum(value, na.rm = TRUE)) |>
  mutate(action = ifelse(action == "lockin", "production", action)) |>
  ungroup() |> #filter(area >0) |>
  mutate(intensity = ifelse(intensity == "primary", "natural", intensity)) |>
  mutate(name = paste0(maes_label, " \n (", intensity, ")")) |>
  filter(action == "restore")

rest_bar_coord <- as_tibble(solution_table_plot) |>
  group_by(zone, maes_label, intensity, action) |>
  summarise(area = sum(value)) |>
  mutate(action = ifelse(action == "lockin", "production", action)) |>
  ungroup() |> filter(area >0) |>
  mutate(intensity = ifelse(intensity == "primary", "natural", intensity)) |>
  mutate(name = paste0(maes_label, " \n (", intensity, ")")) |>
  filter(action == "restore") |>
  left_join(rest_bar_niave)

bar_rest_SI <- rest_bar_coord |> dplyr::select(name, area, area_naive) |>
  rename(area_joint = area) |>
  pivot_longer(cols = c("area_joint", "area_naive"), names_to = "scenario") |>
  ggplot(aes(x = reorder(name, -value), y = value/10, fill = scenario)) +
  geom_bar(stat="identity", width = 0.5, position = position_dodge(width = 0.6)) +
  theme_classic() +
  ylim(0,350) +
  labs(x = element_blank(), y =  expression("area (1000 km"^2*")")) +
  theme(legend.position = c(0.8,0.8),
        legend.title = element_blank(),
        axis.text = element_text(size = 9, hjust =1)) + coord_flip() +
  scale_fill_manual(values = c(met.brewer(name="Kandinsky",n=4,type="continuous"))[c(1,3,2)])

ggsave("figures/updated/bar_rest_SI_A.png", bar_rest_SI, width = 3.5, height = 3.5, dpi = 300)

solution_raster_n <- fasterize(st_as_sf(solution_table_plot_n), raster(PU_plot), field = "value", by = "zone")
solution_raster_n <- rast(solution_raster_n)

rest_agg_naive <- app(solution_raster_n[[restore]], sum)
rst_diff = rest_agg_naive - rest_agg
plot(rst_diff)

cellStats(raster(abs(rst_diff)), "sum")/cellStats(raster(abs(rest_agg)), "sum")

rst_diff_plot <- ggR(rst_diff*100, layer = 1, maxpixels = 1e10, forceCat = FALSE,
                     geom_raster = TRUE, coord_equal = TRUE, stretch = "none", ggObj = TRUE) +
  ggthemes::theme_map(base_size = 11, base_family = 'Arial') +
  scale_fill_scico(direction = 1, palette = "vik",
                   midpoint = 0, limits = c(-100, 100),  # Set limits from -1 to 1
                   na.value = '#FFFFFF',
                   guide = guide_colorbar(title = "Difference", ticks = FALSE)) +
  theme(
    legend.text = element_text(size = 8),
    legend.position = "bottom",  # Move legend to the bottom
    legend.direction = "horizontal",  # Arrange legend items horizontally
    legend.title = element_blank(),  # Center the legend title
    legend.background = element_rect(fill = 'transparent'),
    plot.title = element_text(size = 23, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(x = "", y = "")

ggsave(filename = 'figures/updated/rst_diff_plot_A.png', plot = rst_diff_plot, width = 6, height = 4, dpi = 400)



#### B - even not overestimating potential
solution <- read_csv("data/solutions/sol/NUTSadj/sol_carbon_0.5_restoration_0.141_production_NUTS2adj_country_UNEVEN_wetlands_TRUE_onlyrestoration_FALSE_scenario_HN_proportional_gurobi_f455.csv") |>
  dplyr::select(id, solution_1_z1:solution_1_z26)
colnames(solution) <- c("id", (zone_id$zone))

solution_table <- pu_in_EU |> rename(id = pu) |> #plot_data |>
  left_join(PU_template) |>
  dplyr::select(-id) |>
  rename(id = EU_id) |>
  left_join(solution) |>
  pivot_longer(-c(nuts2id:geometry))

solution_table_plot <- solution_table |>
  mutate(zone = name) |>
  separate(name, c('maes_label', 'intensity', 'action'), sep = "_")

solution_raster <- fasterize(st_as_sf(solution_table_plot), raster(PU_plot), field = "value", by = "zone")
solution_raster <- rast(solution_raster)
names(solution_raster)

restore <- c(1,2,6,3,4,5,7,8,9,10)

rest_agg <- app(solution_raster[[restore]], sum)

solution_niave <- read_csv("data/solutions/sol/SI-sol/sol_carbon_0.5_restoration_0.141_production_NUTS2adj_country_UNEVEN_wetlands_TRUE_onlyrestoration_TRUE_scenario_HN_onlyrest_highs_f455.csv")|>
  dplyr::select(id, solution_1_z1:solution_1_z26)
colnames(solution_niave) <- c("id", (zone_id$zone))

solution_table_niave <- pu_in_EU |> rename(id = pu) |> #plot_data |>
  left_join(PU_template) |>
  dplyr::select(-id) |>
  rename(id = EU_id) |>
  left_join(solution_niave) |>
  pivot_longer(-c(nuts2id:geometry))

solution_table_plot_n <- solution_table_niave |>
  mutate(zone = name) |>
  separate(name, c('maes_label', 'intensity', 'action'), sep = "_")

rest_bar_niave <- as_tibble(solution_table_plot_n) |>
  group_by(zone, maes_label, intensity, action) |>
  summarise(area_naive = sum(value, na.rm = TRUE)) |>
  mutate(action = ifelse(action == "lockin", "production", action)) |>
  ungroup() |> #filter(area >0) |>
  mutate(intensity = ifelse(intensity == "primary", "natural", intensity)) |>
  mutate(name = paste0(maes_label, " \n (", intensity, ")")) |>
  filter(action == "restore")

rest_bar_coord <- as_tibble(solution_table_plot) |>
  group_by(zone, maes_label, intensity, action) |>
  summarise(area = sum(value)) |>
  mutate(action = ifelse(action == "lockin", "production", action)) |>
  ungroup() |> filter(area >0) |>
  mutate(intensity = ifelse(intensity == "primary", "natural", intensity)) |>
  mutate(name = paste0(maes_label, " \n (", intensity, ")")) |>
  filter(action == "restore") |>
  left_join(rest_bar_niave)

rest_bar_coord$name <- factor(rest_bar_coord$name,
                              levels = rev(c("Grassland \n (natural)",
                                             "HeathlandShrub \n (natural)",
                                             "WoodlandForest \n (natural)",
                                             "SparseVeg \n (natural)",
                                             "Wetlands \n (natural)",
                                             "Pasture \n (low)",
                                             "Cropland \n (med)",
                                             "Cropland \n (low)",
                                             "WoodlandForest \n (multi)")))


bar_rest_SI <- rest_bar_coord |> dplyr::select(name, area, area_naive) |>
  rename(area_joint = area) |>
  pivot_longer(cols = c("area_joint", "area_naive"), names_to = "scenario") |>
  ggplot(aes(x = name, y = value/10, fill = scenario)) +
  geom_bar(stat="identity", width = 0.5, position = position_dodge(width = 0.6)) +
  theme_classic() +
  ylim(0,350) +
  labs(x = element_blank(), y =  expression("area (1000 km"^2*")")) +
  theme(legend.position = c(0.8,0.8),
        legend.title = element_blank(),
        axis.text = element_text(size = 9, hjust =1)) + coord_flip() +
  scale_fill_manual(values = c(met.brewer(name="Kandinsky",n=4,type="continuous"))[c(1,3,2)])

ggsave("figures/updated/bar_rest_SI_B.png", bar_rest_SI, width = 3.5, height = 3.5, dpi = 300)

solution_raster_n <- fasterize(st_as_sf(solution_table_plot_n), raster(PU_plot), field = "value", by = "zone")
solution_raster_n <- rast(solution_raster_n)

rest_agg_naive <- app(solution_raster_n[[restore]], sum)
rst_diff = rest_agg_naive - rest_agg
plot(rst_diff)

cellStats(raster(abs(rst_diff)), "sum")/cellStats(raster(abs(rest_agg)), "sum")

rst_diff_plot <- ggR(rst_diff*100, layer = 1, maxpixels = 1e10, forceCat = FALSE,
                     geom_raster = TRUE, coord_equal = TRUE, stretch = "none", ggObj = TRUE) +
  ggthemes::theme_map(base_size = 11, base_family = 'Arial') +
  scale_fill_scico(direction = 1, palette = "vik",
                   midpoint = 0, limits = c(-100, 100),  # Set limits from -1 to 1
                   na.value = '#FFFFFF',
                   guide = guide_colorbar(title = "Difference", ticks = FALSE)) +
  theme(
    legend.text = element_text(size = 8),
    legend.position = "bottom",  # Move legend to the bottom
    legend.direction = "horizontal",  # Arrange legend items horizontally
    legend.title = element_blank(),  # Center the legend title
    legend.background = element_rect(fill = 'transparent'),
    plot.title = element_text(size = 23, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(x = "", y = "")

ggsave(filename = 'figures/updated/rst_diff_plot_B.png', plot = rst_diff_plot, width = 6, height = 4, dpi = 400)


############### land cover figure #############################
nuts2 <- st_read("data/EU_NUTS2_GLOBIOM/EU_GLOBIOM_NUTS2.shp") |>
  rename(NUTS_ID = NURGCDL2) |>
  mutate(nuts2id = seq(1:260)) |>
  mutate(country = str_sub(NUTS_ID,1,2)) |>
  filter(country != "UK") |>
  dplyr::select(nuts2id, NUTS_ID)

plotting_data_ic <- read_csv("data/outputs/2-zones/PU_lc_intensity.csv") |>
  dplyr::select(-Status) |>
  pivot_longer(-PUID) |> rename(pu = PUID) |>
  left_join(pu_in_EU) |>
  rename(id = EU_id) |>
  dplyr::select(-c(pu)) |>
  drop_na(id) |> rename(zone = name) |>
  separate(zone, c('maes_label', 'intensity'), sep = "_") |>
  left_join(as_tibble(nuts2) |> dplyr::select(-geometry))

write_fst(plotting_data_ic, "data/plotting/plotting_data_ic.fst")

##### land use change figure
plotting_data_ic <- read_fst("data/plotting/plotting_data_ic.fst")

solution_table_plot <- plotting_data |>
  filter(carbon == "0.5",
         country_TF == "FLEX")

LC_change_fig_sol <- solution_table_plot |>
  mutate(intensity = ifelse(intensity %in% c("low", "med", "multi"),"low",
                            ifelse(intensity %in% c("high", "prod"), "high",
                                   ifelse(intensity %in% c("primary", "natural"),
                                          "natural",intensity)))) |>
  group_by(nuts2id, intensity, future, scenario) |> summarise(value_sol = sum(value,na.rm = T))

LC_change_fig_ic <- plotting_data_ic |>
  mutate(intensity = ifelse(intensity %in% c("low", "med", "multi"),"low",
                            ifelse(intensity %in% c("high", "prod"), "high",
                                   ifelse(intensity %in% c("primary", "natural"),
                                          "natural",intensity)))) |>
  group_by(nuts2id, intensity) |> summarise(value_ic = sum(value,na.rm = T))

LC_change_fig <- LC_change_fig_sol |>
  left_join(LC_change_fig_ic) |>
  mutate(delta = value_sol-value_ic) |>
  left_join(nuts2_area) |>
  mutate(delta = delta/n)

sum(LC_change_fig$value_sol)/41046

library(scico)  # Load the scico package for the "vik" palette

LC_change_plot_f455_bn <- LC_change_fig |>
  left_join(nuts2) |>
  filter(
    scenario == "HN",
    intensity != "urban"
  ) |>
  ggplot() +
  geom_sf(aes(fill = delta, geometry = geometry)) +
  theme_classic() +
  scale_fill_scico(
    direction = 1,
    palette = "vik",
    midpoint = 0,
    na.value = '#FFFFFF',
    name = "solution - initial conditions \n (% of total land area)"  # Label the legend
  ) +
  facet_grid(
    rows = vars(intensity),
    cols = vars(future),
    labeller = labeller(future = c("f455.csv" = "F455", "ref.csv" = "BAU"))  # Customize column labels
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "white", color = "black", size = 0.5),
    legend.background = element_rect(fill = "white", color = "black", size = 0.5, linetype = "solid"),  # Box around the legend
    legend.title = element_text(size = 10),  # Bold the legend title
    legend.position = "bottom"  # Position the legend on the right
  )

LC_change_plot_f455_bn

ggsave("figures/updated/lc_change.png", LC_change_plot_f455_bn, width = 4, height = 6, dpi = 300)

LC_change_bar <- solution_table_plot |>
  group_by(maes_label, intensity, future, scenario) |>
  summarise(value_sol = sum(value,na.rm = T))

LC_change_fig_ic <- plotting_data_ic |>
  group_by(maes_label, intensity) |> summarise(value_ic = sum(value,na.rm = T))

unique(LC_change_bar_plot$lc)
LC_change_bar_plot <- LC_change_bar |>
  left_join(LC_change_fig_ic) |>
  mutate(diff = value_sol - value_ic) |>
  filter(intensity != "natural") |>
  filter(intensity != "primary") |>
  filter(intensity != "urban") |>
  filter(scenario == "Baseline") |>
  mutate(lc = paste0(maes_label, "\n", intensity)) |>
  mutate(lc = case_when(
    lc ==  "WoodlandForest\nprod" ~ "Forest \n production",
    lc == "WoodlandForest\nmulti" ~ "Forest \n multifunctional",
    lc == "Pasture\nlow" ~ "Pasture \n low-intensity",
    lc == "Pasture\nhigh" ~ "Pasture \n high-intensity",
    lc == "Cropland\nmed" ~ "Cropland \n mid-intensity",
    lc == "Cropland\nlow" ~ "Cropland \n low-intensity",
    lc == "Cropland\nhigh" ~ "Cropland \n high-intensity",
    TRUE ~ "F"
  )) |>
  arrange(lc) |>
  mutate(diff_perc = diff / 41046 * 100) |>
  mutate(future = ifelse(future == "f455.csv","F455", "BAU")) |>  # Customize column labels
  filter(intensity != "urban",
         maes_label != "RiversLakes",
         maes_label != "MarineTransitional") |>
  ggplot(aes(x = diff_perc, y = lc, fill = future)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  scale_alpha_manual(values = c(0.6, 1)) +  # Set alpha values manually
  labs(
    x =  "% change from 2020 conditions \n  (% of total land area)",
    y = "Land Cover Type",
    fill = "Production \n Scenario",
    alpha = "NRL \n Scenario",
    title = "Comparison of Land Cover Changes"
  ) +
  theme_bw(base_size = 10) +  # Adjust base font size for readability
  theme(
    axis.title.y = element_blank(),  # Remove y-axis title
    legend.position = c(0.95, 0.4),  # Position the legend inside the plot area
    legend.justification = c(1, 1),  # Justify the legend in the top-right corner
    legend.background = element_rect(fill = "white", color = "black", size = 0.5, linetype = "solid"),  # Add a white background with a black border to the legend
    #legend.title = element_blank(),  # Remove the legend title
    plot.title = element_blank()  # Center and bold the title
  ) + geom_vline(aes(xintercept = 0), linetype = "dashed") +
  #facet_grid(rows = vars(intensity), scales = "free", space = "free")  +
  theme(
    axis.text = element_text(size = 9.5),
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "white", color = "black", size = 0.5),
  ) + xlim(-3, 10)
# Allow facets to have different sizes

ggsave("figures/updated/lc_change_bar.png", LC_change_bar_plot, width = 3, height = 4, dpi = 300)


LC_change_bar_plot <- LC_change_bar |>
  filter(scenario == "Baseline") |>
  left_join(LC_change_fig_ic) |>
  mutate(diff = value_sol - value_ic) |>
  mutate(lc = paste0(maes_label, intensity)) |>
  arrange(-diff) |>
  mutate(diff_perc = diff / 41046 * 100) |>
  mutate(future = ifelse(future == "f455.csv","F455", "BAU")) |>  # Customize column labels
  filter(intensity != "urban",
         maes_label != "RiversLakes",
         maes_label != "MarineTransitional") |>
  ggplot(aes(x = diff_perc, y = fct_reorder(maes_label, diff), fill = future)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  scale_alpha_manual(values = c(0.6, 1)) +  # Set alpha values manually
  labs(
    x =  "change from 2020 initial conditions \n  (% of total land area)",
    y = "Land Cover Type",
    fill = "Production \n Scenario",
    alpha = "NRL \n Scenario",
    title = "Comparison of Land Cover Changes"
  ) +
  theme_classic(base_size = 10) +  # Adjust base font size for readability
  theme(
    axis.title.y = element_blank(),  # Remove y-axis title
    legend.position = c(0.95, 0.4),  # Position the legend inside the plot area
    legend.justification = c(1, 1),  # Justify the legend in the top-right corner
    legend.background = element_rect(fill = "white", color = "black", size = 0.5, linetype = "solid"),  # Add a white background with a black border to the legend
    #legend.title = element_blank(),  # Remove the legend title
    plot.title = element_blank()  # Center and bold the title
  )

ggsave("figures/updated/lc_change_bar_main.png", LC_change_bar_plot, width = 4, height = 6, dpi = 300)

