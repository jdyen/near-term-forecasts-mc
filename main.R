# Analysis of flow scenarios for Murray cod (Maccullochella peelii)
# 
# Flow scenarios consider observed daily flows with and without
#   environmental water allocations from 2009-2023 and near-term (to 2025) 
#   forecasts of flows under different plausible climates and
#   flow management strategies
#
# Author: Jian Yen (jdl.yen [at] gmail.com)
# 
# Date created: 5 July 2023
# Date modified: 26 November 2025

# load some packages
library(qs)
library(dplyr)
library(tidyr)
library(lubridate)
library(aae.db)
library(aae.hydro)
library(aae.pop.templates)
library(sf)
library(ggplot2)
library(ggspatial)
library(ragg)
library(rstanarm)
library(bayesplot)
library(patchwork)

# and some load helpers
source("R/utils.R")
source("R/fish.R")
source("R/flow.R")
source("R/validation.R")

# settings
set.seed(2023-10-17)
nsim <- 1000
nburnin <- 10
simulate_again <- FALSE

# load data
cpue <- readRDS("data/cpue-compiled.rds")
cpue_recruits <- readRDS("data/cpue-recruits-compiled.rds")

# load stocking info
stocking <- read.csv("data/stocking.csv")

# fetch hydrological data
flow <- readRDS("data/discharge.rds")

# specify futures (individual years and events; can chop and change these to 
#    create specific scenarios)
flow_futures <- specify_flow_future(flow)

# and add some extra information on hypoxia risk and K by species (expands
#   over species)
flow_futures <- flow_futures |>
  left_join(
    fetch_hypoxia_risk(),
    by = "waterbody"
  ) |>
  left_join(
    fetch_carrying_capacity(), 
    by = "waterbody",
    relationship = "many-to-many"
  )

# calculate flow metrics
metrics <- calculate_metrics(flow_futures, recompile = FALSE)

# then left_join this to the pre-2024 metrics
metrics_observed <- metrics |>
  filter(future == "ave", scenario == "none", water_year < 2024) |>
  select(-future, -scenario)

# list all possible future scenarios 
#   (81 combinations per species and waterbody)
scenario_options <- flow_futures |>
  distinct(species, waterbody, future, scenario) |>
  mutate(future_next = future, scenario_next = scenario) |>
  complete(nesting(species, waterbody), future, future_next, scenario, scenario_next)

# for each future, grab the correct 2024 metrics and then grab the
#   2025 metrics based on the *_next settings, but replace antecedent
#   for that scenario with the appropriate 2024 value
metrics_2024 <- metrics |>
  filter(water_year == max(water_year))
metrics_2025 <- metrics_2024 |> mutate(water_year = water_year + 1)
metrics_future <- scenario_options |>
  left_join(
    metrics_2024,
    by = c("species", "waterbody", "future", "scenario")
  ) |>
  left_join(
    metrics_2025,
    by = c("species", "waterbody", "future_next" = "future", "scenario_next" = "scenario"),
    suffix = c("_2024", "_2025")
  ) |>
  mutate(
    proportional_antecedent_flow_2025 = proportional_annual_flow_2024,
    proportional_max_antecedent_2025 = proportional_max_annual_2024
  )
metrics_2024 <- metrics_future |>
  select(species, waterbody, future, future_next, scenario, scenario_next, contains("2024")) |>
  rename(water_year = water_year_2024) |>
  rename_with(\(x) gsub("_2024", "", x), contains("2024"))
metrics_2025 <- metrics_future |>
  select(species, waterbody, future, future_next, scenario, scenario_next, contains("2025")) |>
  rename_with(\(x) gsub("_2025", "", x), contains("2025"))
metrics_future <- bind_rows(metrics_2024, metrics_2025)

# simulate for each species in turn (if required)
if (simulate_again) {
  
  species_list <- metrics_observed |> pull(species) |> unique()
  for (i in seq_along(species_list)) {
    
    # pull out flow/covariate metrics for species
    metrics_observed_sp <- metrics_observed |>
      filter(species == species_list[i]) |>
      select(waterbody, water_year, all_of(get_metric_names(species_list[i])))

    # rename a few metrics
    if (species_list[i] == "maccullochella_peelii") {
      metrics_observed_sp <- metrics_observed_sp |>
        rename(blackwater_risk = hypoxia_risk)
    }

    # simulate for each waterbody in turn
    waterbodies <- metrics_observed_sp |> pull(waterbody) |> unique()  
    for (j in seq_along(waterbodies)) {
      
      # extract carrying capacity for species and reach
      k <- flow_futures |> 
        filter(species == species_list[i], waterbody == waterbodies[j]) |>
        pull(carrying_capacity) |>
        unique()
      
      # filter to each waterbody in turn
      metrics_observed_wb <- metrics_observed_sp |>
        filter(waterbody == waterbodies[j])

      # specify initial conditions
      initial <- specify_initial_conditions(
        species = species_list[i],
        waterbody = waterbodies[j],
        cpue = cpue,
        start = min(metrics_observed_wb$water_year),
        nsim = nsim,
        k = k
      )
      
      # grab stocking info if required (default to zero, otherwise)
      n_stocked <- rep(0, nrow(metrics_observed_wb))
      system_lu <- c(
        "broken_creek_r4" = "Broken Creek",
        "broken_river_r3" = "Broken River",
        "campaspe_river_r4" = "Campaspe River",
        "goulburn_river_r4" = "Goulburn River",
        "loddon_river_r4" = "Loddon River",
        "ovens_river_r5" = "Ovens River"
      )
      if (species_list[i] == "maccullochella_peelii") {
        stocking_rates <- stocking |>
          filter(
            Species == "Murray Cod",
            System == system_lu[waterbodies[j]]
          ) |>
          select(Year, Number) |>
          mutate(Year = Year + 1) |>
          rename(water_year = Year, number_stocked = Number)
        n_stocked <- metrics_observed_wb |>
          select(water_year) |>
          left_join(stocking_rates, by = c("water_year")) |>
          mutate(number_stocked = ifelse(is.na(number_stocked), 0, number_stocked)) |>
          pull(number_stocked)
      }
      
      # initialise population model for this species and waterbody
      pop <- specify_pop_model(
        species = species_list[i],
        waterbody = waterbodies[j],
        ntime = nrow(metrics_observed_wb), 
        nstocked = n_stocked,
        k = k
      )
      
      # and simulate from this model
      initial <- simulate_scenario(
        species = species_list[i],
        x = pop, 
        nsim = nsim, 
        init = initial,
        metrics = metrics_observed_wb[1, ],
        coefs = get_coefs(species_list[i], waterbodies[j]),
        nburnin = nburnin - 1
      )
      sims_observed <- simulate_scenario(
        species = species_list[i],
        x = pop, 
        nsim = nsim, 
        init = initial[, , dim(initial)[3]],
        metrics = metrics_observed_wb,
        coefs = get_coefs(species_list[i], waterbodies[j]),
        nburnin = 0
      )
      
      # save output
      qsave(
        sims_observed, 
        file = paste0("outputs/simulated/observed-", species_list[i], "-", waterbodies[j], ".qs")
      )

      # extract initial conditions for forecasts from sims_observed
      initial_future <- sims_observed[, , dim(sims_observed)[3]]
      
      # simulate futures
      future_sub <- metrics_future |>
        filter(
          species == species_list[i],
          waterbody == waterbodies[j]
        ) |>
        distinct(future, future_next, scenario, scenario_next)
      for (ff in seq_len(nrow(future_sub))) {
        
        # pull out metrics for a given scenario
        metrics_future_sub <- metrics_future |>
          filter(
            species == species_list[i],
            waterbody == waterbodies[j],
            future == future_sub$future[ff],
            future_next == future_sub$future_next[ff],
            scenario == future_sub$scenario[ff],
            scenario_next == future_sub$scenario_next[ff]
          )
        
        # rename a few metrics
        if (species_list[i] == "maccullochella_peelii") {
          metrics_future_sub <- metrics_future_sub |>
            rename(blackwater_risk = hypoxia_risk)
        }

        # simulate under the specific scenario
        sims_future <- simulate_scenario(
          species = species_list[i],
          x = pop, 
          nsim = nsim, 
          init = initial_future,
          metrics = metrics_future_sub,
          coefs = get_coefs(species_list[i], waterbodies[j]),
          nburnin = 0
        )
        
        # and save output
        future_name <- paste(
          future_sub$future[ff],
          future_sub$future_next[ff],
          future_sub$scenario[ff],
          future_sub$scenario_next[ff],
          sep = "_"
        )
        qsave(
          sims_future, 
          file = paste0("outputs/simulated/future-", species_list[i], "-", waterbodies[j], "-", future_name, ".qs")
        )
        
      }
      
    }
    
  }
  
}

# load all simulated models
mc_sim_obs <- load_simulated(type = "observed", species = "maccullochella")
mc_sim_future <- load_simulated(type = "future", species = "maccullochella")

# model CPUE using an AR1 model to estimate values (simple AR1 model with
#   random terms to soak up variation)
iter <- 4000
warmup <- 2000
chains <- 4
cores <- 4
use_cached <- TRUE
cpue_mc <- estimate_cpue(
  x = cpue, 
  use_cached = use_cached,
  species = "Maccullochella peelii",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)                       
cpue_recruit_mc <- estimate_cpue(
  x = cpue_recruits, 
  recruit = TRUE,
  use_cached = use_cached,
  species = "Maccullochella peelii",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores
)

# posterior checks (saved to figures)
pp_mc <- pp_check(cpue_mc) + scale_x_log10() + xlab("Catch (total)") + ylab("Density") + theme(legend.position = "none")
pp_mc_recruit <- pp_check(cpue_recruit_mc) + scale_x_log10() + xlab("Catch (young of year)") + ylab("Density") + theme(legend.position = "none")
pp_all <- (pp_mc | pp_mc_recruit) +
  plot_annotation(tag_levels = "a")
ggsave(
  filename = "outputs/figures/pp-checks.png",
  plot = pp_all,
  device = ragg::agg_png,
  width = 6,
  height = 3,
  units = "in",
  dpi = 600,
  bg = "white"
)

# model validation (using observed minus outputs)
mc_sim_metrics <- calculate_val_metrics(
  x = mc_sim_obs,
  cpue_mod = cpue_mc,
  subset = 1:50, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
mc_sim_metrics_recruit <- calculate_val_metrics(
  x = mc_sim_obs,
  cpue_mod = cpue_recruit_mc,
  recruit = TRUE,
  subset = 1, 
  sim_years = (min(metrics_observed$water_year) - 1L):max(metrics_observed$water_year)
)
sim_metrics <- bind_rows(
  mc_sim_metrics |> mutate(species = "Murray Cod"),
  mc_sim_metrics_recruit |> mutate(species = "Murray Cod (young of year)")
)
metrics_plot_mc <- plot_metric(sim_metrics)
ggsave(
  filename = "outputs/figures/metrics-mc.png",
  plot = metrics_plot_mc,
  device = ragg::agg_png,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

# plot all hindcast combinations
mc_hindcast <- plot_hindcasts(
  x = mc_sim_obs,
  cpue = cpue_mc,
  subset = 1:50, 
  sim_years = min(metrics_observed$water_year):max(metrics_observed$water_year)
)
mc_hindcast_recruit <- plot_hindcasts(
  x = mc_sim_obs,
  cpue = cpue_recruit_mc,
  recruit = TRUE,
  subset = 1, 
  sim_years = (min(metrics_observed$water_year) - 1L):max(metrics_observed$water_year)
)
ggsave(
  filename = "outputs/figures/hindcast-mc.png",
  plot = mc_hindcast,
  device = ragg::agg_png,
  width = 7.5,
  height = 5,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/hindcast-mc-recruits.png",
  plot = mc_hindcast_recruit,
  device = ragg::agg_png,
  width = 7.5,
  height = 5,
  units = "in",
  dpi = 600
)

# plot futures for each species and waterbody
#    - one step ahead (recruits)
#    - two steps ahead under all combos of climate and flow
mdb_systems <- c(
  "goulburn_river_r4", "campaspe_river_r4", "broken_creek_r4",
  "broken_river_r3", "ovens_river_r5", "loddon_river_r4"
)
mc_recruit_futures <- plot_forecasts(
  x = mc_sim_future,
  subset = 1,
  probs = c(0.1, 0.9),
  target = 2024
)
for (i in seq_along(mdb_systems)) {
  
  # plot future for all MC over all climates and flows, save to file
  mc_all_futures <- plot_forecasts(
    x = mc_sim_future,
    subset = 1:50,
    probs = c(0.1, 0.9),
    system = mdb_systems[i]
  )
  ggsave(
    filename = paste0("outputs/figures/futures-mc-", mdb_systems[i], ".png"),
    plot = mc_all_futures,
    device = ragg::agg_png,
    width = 7,
    height = 7,
    units = "in",
    dpi = 600
  )
  
}

# save recruit plots to file
ggsave(
  filename = "outputs/figures/futures-mc-recruits.png",
  plot = mc_recruit_futures,
  device = ragg::agg_png,
  width = 7.5,
  height = 5,
  units = "in",
  dpi = 600
)

# download 2024 survey data and validate one-step forecasts for recruitment
#   and pop growth
cpue_recruits24 <- readRDS("data/cpue-recruits24-compiled.rds")
cpue_adults24 <- readRDS("data/cpue-adults24-compiled.rds")

# refit the statistical models with these updated data
use_cached <- FALSE
iter <- 4000
warmup <- 2000
chains <- 4
cores <- 4
cpue_mc24 <- estimate_cpue(
  x = cpue_adults24, 
  use_cached = use_cached,
  species = "Maccullochella peelii",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores,
  cache_model = FALSE
)                       
cpue_mc_recruits24 <- estimate_cpue(
  x = cpue_recruits24, 
  recruit = TRUE,
  use_cached = use_cached,
  species = "Maccullochella peelii",
  iter = iter,
  warmup = warmup,
  chains = chains,
  cores = cores,
  cache_model = FALSE
)                       

# plot these future forecasts
p_adults <- plot_forecasts_update(
  x = list(mc_sim_future),
  cpue = list(cpue_mc24),
  sciname = c("Maccullochella peelii"),
  recruit = FALSE,
  subset = list(5:50),
  sim_years = 2023:2025,
  survey_max = 2024,
  scenario_set = "baseflow",
  future_set = "ave"
)
p_recruit <- plot_forecasts_update(
  x = list(mc_sim_future),
  cpue = list(cpue_mc_recruits24),
  sciname = c("Maccullochella peelii"),
  recruit = TRUE,
  subset = list(1),
  sim_years = 2023:2025,
  survey_max = 2024,
  scenario_set = "baseflow",
  future_set = "ave"
)

ggsave(
  filename = "outputs/figures/forecast-validation.png",
  plot = p_adults + theme(axis.title.y = element_text(margin = margin(l = 18), vjust = 9)),
  device = ragg::agg_png,
  width = 6,
  height = 3,
  units = "in",
  dpi = 600
)
ggsave(
  filename = "outputs/figures/forecast-validation-recruit.png",
  plot = p_recruit + theme(axis.title.y = element_text(margin = margin(l = 18), vjust = 9)),
  device = ragg::agg_png,
  width = 6,
  height = 3,
  units = "in",
  dpi = 600
)
