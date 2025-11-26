# function to return hypoxia risk as a neat tibble
fetch_hypoxia_risk <- function() {
  
  # set carrying capacity
  .sys_list <- c(
    "broken_creek_r4",
    "broken_river_r3",
    "campaspe_river_r4",
    "goulburn_river_r4",
    "loddon_river_r4",
    "ovens_river_r5"
  )
  
  # collate and return
  tibble(waterbody = .sys_list) |>
    mutate(
      hypoxia_threshold_discharge = NA,
      hypoxia_threshold_discharge = ifelse(
        waterbody == "broken_creek_r4", 
        350,
        hypoxia_threshold_discharge
      ),
      hypoxia_threshold_temp = NA,
      hypoxia_threshold_temp = ifelse(
        waterbody == "broken_creek_r4", 
        25,
        hypoxia_threshold_temp
      )
    )
  
}

# function calculate flow metrics from daily discharge and temperature data
calculate_metrics <- function(data, reference = NULL, recompile = FALSE) {
  
  # check if metrics exist
  metrics_exist <- "metrics-compiled.rds" %in% dir("data/")
  
  # load metrics from file if they exist and !recompile
  if (!recompile & metrics_exist) {
    
    # load from file
    out <- readRDS("data/metrics-compiled.rds")
    
  } else {
    
    # which year range do we want?  
    year <- min(year(data$date_formatted)):max(year(data$date_formatted))
    
    # and grab some reference values for rescaling
    if (is.null(reference)) {
      reference <- data |>
        group_by(waterbody) |>
        summarise(
          reference_discharge = median(stream_discharge_mld, na.rm = TRUE),
          reference_temperature = median(water_temperature_c, na.rm = TRUE),
          reference_max = max(stream_discharge_mld, na.rm = TRUE)
        )
    } else {
      reference <- reference |>
        group_by(waterbody) |>
        summarise(
          reference_discharge = median(stream_discharge_mld, na.rm = TRUE),
          reference_temperature = median(water_temperature_c, na.rm = TRUE),
          reference_max = max(stream_discharge_mld, na.rm = TRUE)
        )
    }
    
    # pull out carrying capacity and hypoxia levels by group
    carrying_capacity <- data |> 
      group_by(species, waterbody) |>
      summarise(carrying_capacity = unique(carrying_capacity))
    hypoxia_thresholds <- data |>
      group_by(waterbody) |>
      summarise(
        hypoxia_threshold_discharge = unique(hypoxia_threshold_discharge),
        hypoxia_threshold_temp = unique(hypoxia_threshold_temp)
      )
    
    # pull out relevant variables
    data <- data |> 
      select(
        species,
        waterbody,
        future,
        scenario,
        stream_discharge_mld,
        date_formatted,
        water_temperature_c
      ) |>
      rename(
        flow = stream_discharge_mld,
        date = date_formatted,
        temp = water_temperature_c
      )
    
    # calculate metrics
    out <- data |> 
      left_join(hypoxia_thresholds, by = "waterbody") |>
      group_by(species, waterbody, future, scenario) |>
      reframe(
        water_year = calculate(
          flow,
          date,
          resolution = survey(season = 9:11, subset = year), 
          na.rm = TRUE
        )$date,
        proportional_spring_flow = calculate(
          flow,
          date,
          resolution = survey(season = 9:11, subset = year), 
          na.rm = TRUE
        )$metric,
        relative_max_annual = calculate(
          flow, 
          date,
          resolution = survey(season = 7:18, subset = year),
          fun = max,
          na.rm = TRUE
        )$metric,
        proportional_annual_flow = calculate(
          flow,
          date,
          resolution = survey(season = c(7:18), subset = year),
          na.rm = TRUE
        )$metric,
        relative_max_antecedent = calculate(
          flow, 
          date,
          resolution = survey(season = 7:18, lag = 1, subset = year + 1),
          fun = max,
          na.rm = TRUE
        )$metric,
        proportional_antecedent_flow = calculate(
          flow,
          date,
          resolution = survey(season = c(7:18), subset = year + 1, lag = 1),
          na.rm = TRUE
        )$metric,
        proportional_summer_flow = calculate(
          flow,
          date,
          na.rm = TRUE,
          resolution = survey(season = c(12:15), subset = year)
        )$metric,
        proportional_winter_flow = calculate(
          flow, 
          date, 
          na.rm = TRUE, 
          resolution = survey(season = 6:8, subset = year)
        )$metric,
        spawning_temperature = calculate(
          temp, date, resolution = survey(season = 10:12, subset = year)
        )$metric,
        spawning_flow_variability = calculate(
          flow,
          date,
          resolution = survey(season = 10:12, subset = year),
          fun = rolling_range,
          lag = 3,
          na.rm = TRUE
        )$metric,
        minimum_daily_flow = calculate(
          flow, 
          date, 
          na.rm = TRUE, 
          resolution = survey(season = 7:18, subset = year),
          fun = min
        )$metric,
        annual_median_flow = calculate(
          flow,
          date,
          na.rm = TRUE,
          resolution = survey(season = 7:18, subset = year)
        )$metric,
        annual_median_flow_tm1 = calculate(
          flow,
          date,
          na.rm = TRUE,
          resolution = survey(season = 7:18, subset = year, lag = 1)
        )$metric,
        annual_median_flow_tm2 = calculate(
          flow,
          date,
          na.rm = TRUE,
          resolution = survey(season = 7:18, subset = year, lag = 2)
        )$metric,
        annual_median_flow_tm3 = calculate(
          flow,
          date,
          na.rm = TRUE,
          resolution = survey(season = 7:18, subset = year, lag = 3)
        )$metric,
        hypoxia_risk_discharge = calculate(
          flow,
          date,
          resolution = survey(season = 12:15, subset = year), 
          fun = blackwater_fn,
          trigger = hypoxia_threshold_discharge, 
          na.rm = TRUE
        )$metric,
        hypoxia_risk_temp = calculate(
          temp,
          date,
          resolution = survey(season = 12:15, subset = year), 
          fun = blackwater_fn,
          trigger = hypoxia_threshold_temp, 
          type = days_above,
          na.rm = TRUE
        )$metric
      ) |>
      left_join(carrying_capacity, by = c("species", "waterbody")) |>
      left_join(reference, by = "waterbody") |>
      mutate(
        lagged_median_flow = median(
          c(annual_median_flow_tm1, annual_median_flow_tm2, annual_median_flow_tm3),
          na.rm = TRUE
        ) / reference_discharge,
        kdyn = lagged_median_flow * carrying_capacity,
        kdyn = ifelse(kdyn > 3 * carrying_capacity, 3 * carrying_capacity, kdyn),
        hypoxia_risk = ifelse(hypoxia_risk_discharge & hypoxia_risk_temp, 1, 0)
      )
    
    # rescale based on reference values
    out <- out |>
      mutate(
        across(contains("proportional"), .fns = \(x) x / reference_discharge),
        spawning_temperature = spawning_temperature / reference_temperature,
        proportional_max_antecedent = relative_max_antecedent / reference_max,
        proportional_max_annual = relative_max_annual / reference_max
      ) |>
      select(!contains("relative"))
    
    # convert to tibble
    out <- as_tibble(out)
    
    # and save output
    saveRDS(out, "data/metrics-compiled.rds")
    
  }
  
  # and return
  out
  
}

# function to calculate a rolling median over specified lag
rolling_median <- function(x, lag, ...) {
  n <- length(x)
  idx <- sapply(rev(seq_len(lag)) - 1, function(x) rep(x, n))
  idx <- sweep(idx, 1, seq_len(n), "+")
  idx <- ifelse(idx > n, NA, idx)
  df <- matrix(x[idx], nrow = n)
  apply(df, 1, median, ...)
}

# function to calculate blackwater risk from a set of flow and
#    temperature triggers
blackwater_fn <- function(x, trigger, type = days_below, ...) {
  out <- 0
  if (length(trigger) > 1) 
    trigger <- unique(trigger)
  if (!is.na(trigger) & (type(x, trigger, ...) > 0))
    out <- 1
  out
}

# helper functions to define a flow scenario based on freshes or baseflows
get_season <- function(type, season) {
  season <- season
  if (type %in% c("none", "baseflow"))
    season <- c(1:12)
  season
}
get_duration <- function(type, duration) {
  duration <- duration
  if (type %in% c("none", "baseflow"))
    duration <- rep(Inf, 12)
  duration
}
get_magnitude <- function(type, magnitude) {
  magnitude <- magnitude
  if (type %in% c("none", "baseflow") & length(magnitude) == 1)
    magnitude <- rep(magnitude, 12)
  magnitude
}

# function to add scenarios to the of an observed flow sequence based
#   on years to resample and events to include
add_scenario <- function(x, scn, template_years, cutoff = NULL) {
  
  # set a default cutoff
  if (is.null(cutoff))
    cutoff <- max(x$date_formatted[!is.na(x$stream_discharge_mld)])
  
  # drop all values after the cutoff if specified
  x <- x |> filter(date_formatted <= cutoff)
  
  # apply each flow event for a given starting year
  apply_flow_event(
    x$stream_discharge_mld,
    x$date_formatted,
    scn,
    template_years,
    x$water_temperature_c,
    cutoff = cutoff
  )
  
}

# internal function to apply a series of flow events based on specified
#   season, duration, and magnitude
apply_flow_event <- function(x, y, scn, template_years, w, cutoff) {
  
  # initialise outer list (for different years)
  out <- vector("list", length = length(template_years))
  
  # add some identifiers to this list
  names(out) <- names(template_years)
  
  # work through each added year
  for (i in seq_along(template_years)) {
    
    # initialise inner list (for different e-flows)
    out[[i]] <- vector("list", length = length(scn))
    
    # and add some identifiers to this
    names(out[[i]]) <- names(scn)
    
    # pull out data for the year to be replicated
    cutoff_start <- cutoff + days(1)
    year(cutoff_start) <- template_years[i]
    cutoff_end <- cutoff_start + years(1) - days(1)
    idy <- y %within% interval(cutoff_start, cutoff_end)
    xbase <- x[idy]
    ybase <- y[idy]
    wbase <- w[idy]
    
    # work through each different scenario
    for (j in seq_along(scn)) {
      
      # set an initial data set for this scenario
      xsub <- xbase
      ysub <- ybase
      wsub <- wbase
      
      # reformat scenario definitions so they can be applied generically        
      season_set <- get_season(names(scn)[j], season = scn[[j]]$season)
      duration_set <- get_duration(names(scn)[j], duration = scn[[j]]$duration)
      magnitude_set <- get_magnitude(names(scn)[j], magnitude = scn[[j]]$magnitude)
      
      # work through each flow "event" and modify flows accordingly
      for (k in seq_along(season_set)) {
        
        # filter to target month
        idx <- month(ysub) == season_set[k]
        
        # and subsample the days if event has a duration shorter than the month
        if (duration_set[k] < sum(idx)) {
          start <- sample(which(idx), size = 1)
          if ((start + duration_set[k]) > max(which(idx)))
            start <- max(which(idx)) - duration_set[k] + 1L
          idx <- idx & (seq_along(idx) %in% start:(start + duration_set[k] - 1L))
        }
        
        xsub[idx] <- ifelse(
          xsub[idx] < magnitude_set[k],
          magnitude_set[k],
          xsub[idx]
        )
        
      }
      
      # combine xsub into a clean data.frame, adding 6 months of NAs at the end
      #   to keep things complete (sequence ends in June)
      out[[i]][[j]] <- data.frame(
        date_formatted = ysub,
        stream_discharge_mld = xsub,
        water_temperature_c = wsub
      )
      
      # update the date to run from the day after the cutoff
      #  to one year later
      cutoff_start <- cutoff + days(1)
      cutoff_end <- cutoff_start + years(1) - days(1)
      new_date <- seq(cutoff_start, cutoff_end, by = 1)
      
      # but drop Feb 29 from the resampled values or the imputed dates
      #  if the number of days doesn't match
      if (length(new_date) > nrow(out[[i]][[j]]))
        new_date <- new_date[!(month(new_date) == 2 & day(new_date) == 29)]
      if (length(new_date) < nrow(out[[i]][[j]])) {
        out[[i]][[j]] <- out[[i]][[j]] |>
          filter(!(month(date_formatted) == 2 & day(date_formatted) == 29))
      }
      
      # and add this new date to the output
      out[[i]][[j]]$date_formatted <- new_date
      
      # pull out any extra days to the end of the calendar year if required
      #   and fill with NAs (to keep aae.hydro::calculate happy)
      if (!(month(cutoff) == 12 & day(cutoff) == 31)) {
        fill_start <- cutoff + days(1)
        year(fill_start) <- max(template_years[i])
        fill_end <- fill_start
        month(fill_end) <- 12
        day(fill_end) <- 31
        idz <- y %within% interval(fill_start, fill_end)
        yextra <- y[idz]
        year(yextra) <- max(year(out[[i]][[j]]$date_formatted))
        na_extra <- rep(NA, sum(idz))
        out[[i]][[j]] <- rbind(
          out[[i]][[j]],
          data.frame(
            date_formatted = yextra,
            stream_discharge_mld = na_extra,
            water_temperature_c = na_extra
          )
        )
      }
      
      # and append earlier years of data
      out[[i]][[j]] <- rbind(
        data.frame(
          date_formatted = y,
          stream_discharge_mld = x,
          water_temperature_c = w
        ),
        out[[i]][[j]]
      )
      
    }
    
  }
  
  # return
  out
  
}

# define flow futures (scenarios to append to the initial data set to
#   demonstrate possible outcomes)
specify_flow_future <- function(x) {
  
  # specify flow events for each reach  
  system_details <- list(
    "broken_creek_r4" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = c(rep(250, 5), rep(40, 3), rep(250, 5))),
      "fresh" = list(
        magnitude = c(300, 450, 400), 
        duration = c(10, 8, 14),
        season = c(7, 10, 11)
      )
    ),
    "broken_river_r3" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = 30),
      "fresh" = list(
        magnitude = 500, 
        duration = 5,
        season = 3
      )
    ),
    "campaspe_river_r4" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = c(rep(50, 5), rep(200, 6), 50)),
      "fresh" = list(
        magnitude = c(200, 1800, 1800, 200), 
        duration = c(3, 5, 7, 3),
        season = c(1, 7, 10, 12)
      )
    ),
    "goulburn_river_r4" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = 800),
      "fresh" = list(
        magnitude = c(5700, 7300, 10500, 8000),
        duration = c(2, 2, 8, 3),
        season = c(5, 7, 10, 11)
      )
    ),
    "loddon_river_r4" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = c(rep(50, 5), rep(100, 6), 50)),
      "fresh" = list(
        magnitude = c(450),
        duration = c(9),
        season = c(10)
      )
    ),
    "ovens_river_r5" = list(
      "none" = list(magnitude = 0),
      "baseflow" = list(magnitude = 0),
      "fresh" = list(
        magnitude = c(80, 430, 260), 
        duration = c(2, 3, 3),
        season = c(2, 3, 4)
      )
    )  
  )
  
  # define a annual flow scenarios for each future type (can chop and
  #   change these to create multi-year sequences)
  template_years <- c("ave" = 2019, "dry" = 2009, "wet" = 2011)
  out <- vector("list", length = length(system_details))
  for (i in seq_along(system_details)) {
    
    # pull out flow for a specific system
    xsub <- x[[names(system_details)[i]]]
    
    # and append the relevant scenarios
    out[[i]] <- add_scenario(
      xsub, 
      system_details[[i]],
      template_years = template_years, 
      cutoff = dmy("30-06-2023")
    )
    
  }
  
  # add some identifiers to the output
  names(out) <- names(system_details)
  
  # flatten to a big tibble
  out <- lapply(out, \(x) lapply(x, add_col_bind))
  out <- lapply(out, add_col_bind, name = "future")
  out <- add_col_bind(out, name = "waterbody")
  
  # return
  out |> as_tibble()
  
}
