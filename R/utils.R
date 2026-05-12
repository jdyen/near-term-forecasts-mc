# function to tidy up lists by adding names as a column and then binding all rows
add_col_bind <- function(x, name = "scenario") {
  id <- names(x)
  out <- mapply(\(.x, .y) .x |> mutate(category = .y), .x = x, .y = id, SIMPLIFY = FALSE)
  bind_rows(out) |> rename(!!name := category)
}

# function to rescale values between [base, top]
rescale_fn <- function(x, base, top = 1, ...) {
  
  # shift to min = 0
  x <- x - min(x, ...)
  
  # scale to [0, 1]
  x <- x / max(x, ...)
  
  # shrink to [0, top - base]
  x <- (top - base) * x
  
  # and add base
  x <- base + x
  
  # return
  x
  
}

# function to set initial conditions based on assumed constant proportional
#   survival over all age classes, but modified based on survey data
#   to account for shifts to older/younger age classes
set_initial <- function(
    species,
    cpue,
    cpue_max,
    effort_h,
    n = 10000,
    nsim = 1,
    rescale = NULL
) {
  
  # define a population matrix and work out steady state age frequency
  if (species == "maccullochella_peelii") {
    adults <- 5:50
    mat <- murray_cod(k = n)$dynamics$matrix
    initial_age_frequency <- Re(eigen(mat)$vectors[, 1])
  }
  
  # work out proportion of fish relative to max catch for a reach
  n_scaled <- 0
  if (cpue_max > 0)
    n_scaled <- n * (cpue / cpue_max)
  
  # rescale if required
  if (!is.null(rescale))
    n_scaled <- 0.6 * cpue * rescale  # current assumption is 60% of catch is adults
  
  # set a min value if cpue_adult = 0
  if (n_scaled == 0) {
    if (cpue_max > 0) {
      cpue_min <- min(c(cpue_max, 1)) / effort_h
      n_scaled <- n * (cpue_min / cpue_max)
    } else {
      n_scaled <- 1 / effort_h
    }
  }
  
  # standardise this for adult abundances only
  initial_age_frequency <- initial_age_frequency / sum(initial_age_frequency[adults])
  
  # return
  matrix(
    rpois(
      nsim * length(initial_age_frequency),
      lambda = n_scaled * initial_age_frequency
    ),
    nrow = nsim,
    byrow = TRUE
  )
  
}

# function to initialise populations
specify_initial_conditions <- function(species, waterbody, cpue, start, nsim, k, ...) {
  
  # stop if species isn't one of the three targets
  stopifnot(
    species %in% c("maccullochella_peelii")
  )
  
  # rename waterbody to avoid conflicts below
  wb <- waterbody
  
  # pull out CPUE for the target species
  cpue_sub <- cpue |> 
    mutate(
      sciname = tolower(gsub(" ", "_", scientific_name)),
      waterbody = paste(
        tolower(gsub(" ", "_", waterbody)),
        reach_no,
        sep = "_r"
      )
    ) |>
    filter(
      sciname == species,
      waterbody == wb
    ) |>
    mutate(cpue_max = max(catch / effort_h))
  
  # work out the minimum starting year that's actually in the data
  if (start %in% cpue_sub$survey_year) {
    cpue_sub <- cpue_sub |>
      filter(survey_year == start)
  } else {
    start <- abs(start - unique(cpue_sub$survey_year))
    start <- unique(cpue_sub$survey_year)[which.min(start)]
    cpue_sub <- cpue_sub |>
      filter(survey_year == start)
  }
  
  # calculate total CPUE per reach
  cpue_sub <- cpue_sub |>
    group_by(survey_year) |>
    summarise(
      cpue_max = unique(cpue_max),
      catch = sum(catch),
      effort_h = sum(effort_h)
    ) |>
    mutate(
      cpue = catch / effort_h
    )
  
  # pull out reach length for target waterbody
  reach_length <- .vewh_reach_lengths |> 
    filter(waterbody == wb) |>
    pull(reach_length) |>
    as.numeric()
  
  # work out the rescaling rate for catch to abundance based on a capture
  #    rate of 40%
  adult_rescale <- (1 / 0.4) * (reach_length / 100)
  init <- set_initial(
    species = species,
    cpue = cpue_sub$cpue,
    cpue_max = cpue_sub$cpue_max,
    effort_h = cpue_sub$effort_h,
    n = k,
    nsim = nsim,
    rescale = adult_rescale
  )
  
  # return
  init
  
}

# function to set up all three population models
specify_pop_model <- function(species, waterbody, ntime, nstocked, k, ...) {
  
  # stop if species isn't one of the three targets
  stopifnot(
    species %in% c("maccullochella_peelii")
  )
  
  # murray cod model
  if (species == "maccullochella_peelii") {
    
    # use a lookup to define system from waterbody
    system <- switch(
      waterbody,
      "broken_creek_r4" = "broken_creek",
      "broken_river_r3" = "broken_river",
      "campaspe_river_r4" = "campaspe_river",
      "goulburn_river_r4" = "goulburn_river",
      "loddon_river_r4" = "campaspe_river",
      "ovens_river_r5" = "ovens_river",
      "murray_river"
    )
    
    mod <- murray_cod(
      k = k,
      system = system,
      n = list(
        # number stocked, accounting for fingerling mortality and 50:50 sex ratio
        nstocked[seq_len(ntime)],
        rep(0, ntime),
        rep(0, ntime)
      ),
      ntime = ntime, 
      start = rep(1, 3), 
      end = rep(ntime, 3), 
      add = c(TRUE, TRUE, TRUE),
      p_capture = 0.1,   # 10% capture probability for any fish in slot
      slot = c(550, 750) # slot in mm
    )
  }
  
  # return
  mod
  
}

# function to return coefficients for each species
get_coefs <- function(species, waterbody) {
  coefs <- list(
    "maccullochella_peelii" = list(
      "broken_creek_r4" = c(-10, 30, 6, -30, 100, 100),
      "broken_river_r3" = c(-15, 35, 45, -30, 80, 25),
      "campaspe_river_r4" = c(-60, 10, 20, -15, 45, 10),
      "goulburn_river_r4" = c(-5, 20, 6, -10, 30, 10),
      "loddon_river_r4" = c(-30, 20, 10, -10, 30, 10),
      "ovens_river_r5" = c(-10, 60, 20, -30, 60, 25)
    )
  )
  coefs[[species]][[waterbody]]
}

# function to pull out names of all metrics for relevant species
get_metric_names <- function(species) {
  metric_list <- list(
    "maccullochella_peelii" = c(
      "spawning_flow_variability",
      "proportional_spring_flow",
      "proportional_max_antecedent",
      "proportional_summer_flow",
      "proportional_winter_flow",
      "spawning_temperature",
      "hypoxia_risk",
      "minimum_daily_flow",
      "carrying_capacity"
    )
  )
  metric_list[[species]]
}

# function to simulate a pop model for a single species based on input
#    metrics and coefficients
simulate_scenario <- function(species, x, nsim, init, metrics, coefs, nburnin = 0, ...) {
  
  # stop if species isn't one of the three targets
  stopifnot(
    species %in% c("maccullochella_peelii")
  )
  
  # murray cod model
  if (species == "maccullochella_peelii") {
    
    # include a short loop to stabilise the initial conditions if 
    #   required
    if (nburnin > 0) {
      
      # simulate repeatedly to burn-in the inits
      sims <- simulate(
        x,
        nsim = nsim,
        init = init,
        args = list(
          covariates = c(
            format_covariates(metrics[rep(1, nburnin), ]),
            list(threshold = 0.05),
            list(coefs = coefs)
          ),
          density_dependence = list(
            kdyn = lapply(
              rep(1, nburnin), 
              function(.x) metrics$carrying_capacity[.x]
            )
          )
        ),
        options = list(
          update = update_binomial_leslie,
          tidy_abundances = floor
        )
      )
      
      # pull out final iteration as the starting point for main sims
      init <- sims[, , nburnin] 
      
    }
    
    # main simulation
    sims <- simulate(
      x,
      nsim = nsim,
      init = init,
      args = list(
        covariates = c(
          format_covariates(metrics),
          list(threshold = 0.05),
          list(coefs = coefs)
        ),
        density_dependence = list(
          kdyn = lapply(
            seq_len(nrow(metrics)), 
            function(.x) metrics$carrying_capacity[.x]
          )
        )
      ),
      options = list(
        update = update_binomial_leslie,
        tidy_abundances = floor
      )
    )
    
  }
  
  # return
  sims
  
}

# function to load simulated scenarios and parse scenarios from filenames
load_simulated <- function(type, species) {
  
  # list all saved simulations
  sim_list <- dir("outputs/simulated/")
  
  # filter to target type and species
  sim_list_sub <- sim_list[
    grepl(paste(type, species, sep = "-"), sim_list, ignore.case = TRUE)
  ]
  out <- lapply(sim_list_sub, \(x) qread(paste0("outputs/simulated/", x)))
  
  # parse names
  waterbody <- sapply(strsplit(sim_list_sub, split = "-"), \(x) x[3])
  if (type == "future") {
    
    # work out the scenario as well as the system
    scenario <- sapply(strsplit(sim_list_sub, split = "-"), \(x) x[4])
    scenario <- gsub(".qs", "", scenario)
    
    # split it up by climate and flows, and add waterbody info
    naming_fn <- function(x, wb) {
      tibble(
        waterbody = wb,
        future = x[1],
        future_next = x[2],
        scenario = x[3],
        scenario_next = x[4]
      )
    }
    scenario <- strsplit(scenario, "_") |>
      mapply(FUN = naming_fn, wb = waterbody, SIMPLIFY = FALSE) |>
      bind_rows()
    
  } else {
    
    # just strip the file suffix from hte system name
    waterbody <- gsub(".qs", "", waterbody)
    
    # reformat
    scenario <- tibble(waterbody = waterbody)
    
  }
  
  # and return
  list(
    scenario = scenario,
    sims = out
  )
  
}
