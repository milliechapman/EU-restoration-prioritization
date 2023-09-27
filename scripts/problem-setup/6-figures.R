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
#devtools::install_github("yutannihilation/ggsflabel")
library(ggsflabel)
library(RStoolbox)
library(scico)
library(geomtextpath)
library(ggpubr)
library(ggpattern)
library(viridis)
library(ggsci)
library(formattable)
rm(list = ls())

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
plotting_data <- read_csv("data/plotting/plotting_data_ref_f455.csv")

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

country <- read_csv("data/formatted-data/linear_constraints/pu_pasture_low_budget_data_55.csv") |>
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
rep <- read_csv("data/solutions/representation_REF_f455_scenarios_version3.csv")

rep |> filter(species == 999999)

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
b_status <- rep |>
  left_join(ic_join) |>
  filter(species != 999999) |>
  mutate(target_met_ic = ifelse(shortfall_perc_ic < 0, 0, 1),
         target_met = ifelse(shortfall_perc < 0, 0, 1)) |>
  filter(target_met_ic<1) |>
  mutate(improved = (target_met-target_met_ic)) |>
  group_by(spp, scenario, carbon_weight, country_contraint, future) |>
  summarise(improved = sum(improved)/length(unique(species))) |>
  arrange(-improved)

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
  mutate(diff = ifelse((`EVEN` - `UNCONSTRAINED`)<0,1,0),
         n_spp = length(unique(taxon_id))) |>
  group_by(scenario, carbon_weight, future) |>
  summarise(diff_spp = sum(diff, na.rm = TRUE)/mean(n_spp,  na.rm = TRUE)) |>
  arrange(-diff_spp)


rep |>
  filter(species == 999999) |>
  dplyr::select(scenario, carbon_weight, country_contraint, future, species, rep) |>
  pivot_wider(names_from = country_contraint, values_from = rep) |>
  left_join(ic_join) |>
  mutate(diff = (`EVEN` - `UNCONSTRAINED`)/rep_ic*-100) |>
  arrange(diff)


################### figure 4 - all scenarios #################

write_csv(b_status |> filter(carbon_weight>0), "data/biodiv_results.csv")
write_csv(c_perc |> filter(carbon_weight>0), "data/carbon_results.csv")

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
solution <- read_csv("data-formatted/sol/sol_carbon_0.1_restoration_0.141_production_TRUE_country_FLEX_wetlands_TRUE_onlyrestoration_FALSE_scenario_Baseline_globiomICflat_gurobi_f455.csv") |>
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

#nuts2_crop_low <- read_csv("data-formatted/linear_constraints/nuts_crop_low_95

fig1_solution_bar <- as_tibble(solution_table_plot) |>
  group_by(zone, maes_label, intensity, action) |>
  summarise(area = sum(value)) |>
  mutate(action = ifelse(action == "lockin", "production", action)) |>
  ungroup() |> filter(area >0) |>
  mutate(intensity = ifelse(intensity == "primary", "natural", intensity)) |>
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
produce <- c(13:18, 26)
conserve <- c(19:25,12)

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

#################### burden_sharing #######################
########## pareto
pareto_burden <- rep |> filter(target >0) |>
  mutate(shortfall_perc = ifelse(shortfall_perc >0,0, shortfall_perc)) |>
  drop_na(shortfall_perc) |>
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
pareto_burden_plot_new <- b_status |> ungroup() |>
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
  #ylim(11, 40) +
  #xlim(52, 62) +
  #scale_color_manual(values = c("grey2", "grey")) +
  scale_color_manual(values = c( "grey", "grey2", "red")) +
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
                                        carbon == "0.3",
                                        country_TF == "EVEN",
                                        scenario == "Baseline",
                                        future == "f455.csv") |>
  group_by(id, country_TF) |>
  summarise(value = sum(value)) |>
  left_join(PU_template_EU)

country_true <- fasterize(st_as_sf(country_true), raster, field = "value")

country_false <- plotting_data |> filter(action == "restore",
                                         carbon == "0.3",
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
  labs(x = element_blank()) + labs(y = "restoration priority area \n (% land area; II-I all scenarios)") + geom_hline(yintercept = 0, linetype="dashed") +
  theme(axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11))

ggsave("figures/updated/burden-sharing-countries.png", fig3d_box, width = 4.5, height = 8, dpi = 300)

## land cover
colors = c(met.brewer(name="Kandinsky",n=4,type="continuous"))[c(4,2,1)]
fig3c_data <- plotting_data |>
  filter(carbon == 0.3) |>
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
  group_by(speciesgroup, target_met) |>
  count() |> pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met_IC = yes/(yes+no)) |> dplyr::select(-c(no, yes))

library(ggpubr)

unique(rep_formatted$speciesgroup)

spp <- rep_formatted |> filter(scenario == "Baseline",
                               future == "f455", carbon_weight == "0.3") |>
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


  # Box plot facetted by "dose"
burden_biodiv_plot <- ggboxplot(burden_biodiv, x = "country_contraint", y = "perc_met",
                 color = "country_contraint",
                 add = "jitter",
                 facet.by = "speciesgroup", nrow = 1,
                 short.panel.labs = TRUE, scales = "free") +
   # stat_compare_means(aes(label = ..p.signif..))  +
  theme_classic() +
  scale_color_manual(values = c( "grey", "grey2", "red")) +
  theme(legend.position = "none") +
  labs(x = "", y = "percent of species targets met") +
  theme(axis.text.x = element_text(size = 12)) #+
  #scale_y_continuous(breaks=equal_breaks(n=10, s=0.1))

ggsave("figures/updated/burden-sharing-spp.png",burden_biodiv_plot, width=8 , height = 3.2, dpi = 300)

###################### restoration scenarios #############

## pareto new metrics
pareto_restoration <- b_status |> ungroup() |> dplyr::select(-spp) |> filter(carbon_weight >0) |>
  filter(country_contraint == "FLEX") |>
  left_join(c_perc |> ungroup() |> dplyr::select(-spp)) |>
  ggplot(aes(x = improved*100, y = rep_diff)) +
  geom_line(aes(linetype = scenario, alpha = future), lwd = 1) +
  geom_point(aes(col = carbon_weight), alpha = 0.9, size = 4) +
  theme_classic() +
  scale_color_gradientn(colors = met.brewer("OKeeffe2")) +
  #scale_color_manual(values = c("grey2", "grey")) +
  scale_alpha_manual(values = c(1, 0.4)) +
  labs(x = "% of species with improved \n conservation status",
       y = "% increase in land carbon stock \n (compared to current conditions)")

## pareto new metrics
pareto_restoration2 <- b_status |> ungroup() |> dplyr::select(-spp) |> filter(carbon_weight >0) |>
  filter(future == "ref") |>
  left_join(c_perc |> ungroup() |> dplyr::select(-spp)) |>
  ggplot(aes(x = improved*100, y = rep_diff)) +
  geom_line(aes(linetype = scenario, alpha = country_contraint), lwd = 1) +
  geom_point(aes(col = carbon_weight), alpha = 0.9, size = 4) +
  theme_classic() +
  scale_color_gradientn(colors = met.brewer("OKeeffe2")) +
  #scale_color_manual(values = c("grey2", "grey")) +
  scale_alpha_manual(values = c(1, 0.6, 0.3)) +
  labs(x = "% of species with improved \n conservation status",
       y = "% increase in land carbon stock \n (compared to current conditions)")

ggsave(filename = 'figures/updated/pareto_restoration.png', plot = pareto_restoration, width = 5.5, height = 4, dpi = 400)

############ restoration scenario barchart ######

basline <- plotting_data |> filter(action == "restore",
                                   scenario == "Baseline",
                                   country_TF == "FLEX") |>
  mutate(type = ifelse(maes_label== "Cropland", "Cropland \n de-intensification",
                       ifelse(maes_label %in% c("Pasture"), "Pasture \n de-intensification",
                              ifelse(maes_label %in% c("WoodlandForest"), "Forestry \n de-intensification",
                                     "Nature \n landscape")))) |>
  group_by(type, carbon, scenario, future) |>
  summarise(value = sum(value))

hn <-plotting_data |> filter(action == "restore",
                             scenario == "HN",
                             country_TF == "FLEX") |>
  mutate(type = ifelse(maes_label== "Cropland", "Cropland \n de-intensification",
                       ifelse(maes_label %in% c("Pasture"), "Pasture \n de-intensification",
                              ifelse(maes_label %in% c("WoodlandForest"), "Forestry \n de-intensification",
                                     "Nature \n landscape")))) |>
  group_by(type, carbon, scenario, future) |>
  summarise(value = sum(value))


hn |>
  filter(type == "Nature \n landscape") |>
  mutate(perc = value/41046*100) |>
  arrange(-value)

hn_baseline_box <- basline |> bind_rows(hn) |> filter(future == "f455.csv") |>
  ggplot() + geom_boxplot(aes(y = value/10, x = type, col = scenario, )) +
  theme_classic() +
  coord_flip() +
  theme(axis.title.y =element_blank(),
                          axis.text.y = element_text(angle = 90))+
  labs(y = "area (1000 km^2)") +
  scale_color_manual(values = c("grey2", "grey"))

ggsave("figures/updated/hn_baseline_box.png", hn_baseline_box,
       height = 4, width = 3.75)

########### restoration scenarios plots #########

basline <- plotting_data |> filter(action == "restore",
                                      scenario == "Baseline",
                                      country_TF == "FLEX",
                                      carbon == "0.3") |>
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
                                     carbon == 0.3) |>
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
                   guide = guide_colorbar(title = "pasture Difference",ticks = FALSE) ) +
  theme(legend.position = c(.05,.15),legend.background = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 23, hjust = .5),plot.subtitle = element_text(hjust = .5)) +
  labs(x = "", y = "")

ggsave(filename = 'figures/updated/nd.png', plot = nd, width = 16, height = 14, dpi = 400)

forest_diff <- basline[[4]] - hn[[4]]

fd <- ggR(forest_diff, layer = 1, maxpixels = 1e10,forceCat = FALSE,
          geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_scico(direction = 1,palette = "vik",midpoint = 0, na.value = '#FFFFFF',
                   guide = guide_colorbar(title = "pasture Difference",ticks = FALSE) ) +
  theme(legend.position = c(.05,.15),legend.background = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 23, hjust = .5),plot.subtitle = element_text(hjust = .5)) +
  labs(x = "", y = "")

ggsave(filename = 'figures/updated/fd.png', plot = fd, width = 16, height = 14, dpi = 400)


######## rest biodiv ###########

rest_biodiv <- rep_formatted |>
  filter(country_contraint == "FLEX",
         future == "f455") |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Vascular plants", "Non-vascular plants"), "Plants", speciesgroup)) |>
  mutate(speciesgroup = ifelse(speciesgroup %in% c("Mammals", "Birds", "Plants", "Amphibians", "Reptiles"),speciesgroup, "Other")) |>
  # group_by(speciesgroup, country_contraint, carbon_weight) |># |> mutate(nspp = n()) |> ungroup()|>
  group_by(speciesgroup, scenario, target_met, carbon_weight) |>
  count() |> pivot_wider(names_from = target_met, values_from = n) |>
  mutate(perc_met = yes/(yes+no)) |>
  #left_join(lj_biodiv_ic) |>
  #mutate(perc_met_rel = (perc_met - perc_met_IC)*100) |>
  left_join(spp) |>
  #mutate(country_contraint = ifelse(country_contraint == "EVEN", "I",
                                  #  ifelse(country_contraint == "FLEX", "II", "III"))) |>
  mutate(perc_met = round(perc_met*100,0)) |>
  mutate(speciesgroup = paste0(speciesgroup, "\n n spp=", nspp))


# Box plot facetted by "dose"
rest_biodiv_plot <- ggboxplot(rest_biodiv, x = "scenario", y = "perc_met",
                                color = "scenario",
                                add = "jitter",
                                facet.by = "speciesgroup", nrow = 1,
                                short.panel.labs = TRUE, scales = "free") +
  # stat_compare_means(aes(label = ..p.signif..))  +
  theme_classic() +
  scale_color_manual(values = c( "grey", "grey2", "red")) +
  theme(legend.position = "none") +
  labs(x = "", y = "percent of species targets met") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust =1)) #+
#scale_y_continuous(breaks=equal_breaks(n=10, s=0.1))

ggsave("figures/updated/rest-spp.png",burden_biodiv_plot, width=8 , height = 3.2, dpi = 300)

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


##################### SI restoration only ########################
colors <- c("grey", met.brewer(name="VanGogh3",n=20,type="continuous"))

solution <- read_csv("data/solutions/sol/sol_carbon_0.5_restoration_0.2_production_TRUE_country_TRUE_bioregion_FALSE_wetlands_TRUE_onlyrestoration_FALSE_scenario_Baseline_globiom_IC_nm.csv") |>
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

#+ coord_flip()

solution_raster <- fasterize(st_as_sf(solution_table_plot), raster(PU_plot), field = "value", by = "zone")
solution_raster <- rast(solution_raster)
names(solution_raster)

restore <- c(1,2,6,3,4,5,7,8,9,10)

rest_agg <- app(solution_raster[[restore]], sum)

png("figures/updated/rest_agg2.png", width = 5, height = 5, units = "in", res = 400)
plot(rest_agg, col = colors)
dev.off()

solution_niave <- read_csv("data/solutions/sol/rest_only/sol_carbon_0.5_restoration_0.2_production_FALSE_country_TRUE_bioregion_FALSE_wetlands_TRUE_onlyrestoration_TRUE_scenario_Baseline.csv")|>
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

bar_rest_SI <- rest_bar_coord |> select(name, area, area_naive) |>
  rename(area_joint = area) |>
  pivot_longer(cols = c("area_joint", "area_naive"), names_to = "scenario") |>
  ggplot(aes(x = reorder(name, -value), y = value/10, fill = scenario)) + geom_bar(position = "dodge", stat="identity", width = 0.6) +
  theme_classic() +
  labs(x = element_blank(), y = "area (1000 km^2)") +
  theme(legend.position = c(0.8,0.8),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, hjust =1)) +
  scale_fill_manual(values = c(met.brewer(name="Kandinsky",n=4,type="continuous"))[c(1,3,2)])

ggsave("figures/updated/bar_rest_SI.png", bar_rest_SI, width = 4, height = 4, dpi = 300)

solution_raster_n <- fasterize(st_as_sf(solution_table_plot_n), raster(PU_plot), field = "value", by = "zone")
solution_raster_n <- rast(solution_raster_n)


rest_agg_naive <- app(solution_raster_n[[restore]], sum)
png("figures/updated/rest_agg_naive.png", width = 5, height = 5, units = "in", res = 400)
plot(rest_agg_naive, col = colors)
dev.off()


rst_diff = rest_agg_naive - rest_agg
plot(rst_diff)

cellStats(raster(abs(rst_diff)), "sum")/cellStats(raster(abs(rest_agg)), "sum")

rst_diff_plot <- ggR(rst_diff, layer = 1, maxpixels = 1e10,forceCat = FALSE,
                  geom_raster = TRUE, coord_equal = TRUE,stretch = "none", ggObj = TRUE) +
  ggthemes::theme_map(base_size = 11,base_family = 'Arial') +
  scale_fill_scico(direction = 1,palette = "vik",midpoint = 0, na.value = '#FFFFFF',
                   guide = guide_colorbar(title = "Difference",ticks = FALSE) ) +
  theme(legend.position = c(.05,.15),legend.background = element_rect(fill = 'transparent'),
        plot.title = element_text(size = 23, hjust = .5),plot.subtitle = element_text(hjust = .5)) +
  labs(x = "", y = "")

ggsave(filename = 'figures/updated/rst_diff_plot.png', plot = rst_diff_plot, width = 5, height = 4, dpi = 400)

colors <- c("grey", met.brewer(name="VanGogh3",n=20,type="continuous"))
myPal <- colors
myTheme <- rasterTheme(region = myPal)

png("figures/updated/rest_agg.png", width = 5, height = 5, units = "in", res = 400)
plot(rest_agg, col = colors)
dev.off()





