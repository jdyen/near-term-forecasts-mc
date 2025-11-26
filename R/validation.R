# helper functions for validation summaries
# function to estimate CPUE from catch data
estimate_cpue <- function(
    x,
    species, 
    iter = 2000,
    warmup = 1000, 
    chains = 2, 
    cores = 1, 
    recruit = FALSE,
    adult = FALSE,
    use_cached = TRUE,
    cache_model = TRUE
) {
  
  # check species are in correct list
  if (!species %in% unique(x$scientific_name))
    stop("species must be represented in data set", call. = FALSE)
  
  # use saved version if available, otherwise refit the model
  species_clean <- gsub(" ", "_", tolower(species))
  if (recruit) {
    species_clean <- paste0("recruit-", species_clean)
  }
  if (use_cached & any(grepl(species_clean, dir("outputs/fitted/")))) {
    
    # if exists, load fitted version
    mod <- qread(paste0("outputs/fitted/cpue-mod-", species_clean, ".qs"))
    
  } else {
    
    # subset data to target species and collate over all gear types and surveys
    #   within a year
    x <- x |>
      filter(scientific_name == species) |>
      group_by(id_site, waterbody, reach_no, survey_year) |>
      summarise(
        catch = sum(catch),
        effort_h = sum(effort_h)
      ) |>
      ungroup()
    
    # fit model
    if (recruit) {
      
      mod <- stan_glmer(
        catch ~ (1 | waterbody / reach_no) +
          (1 | id_site) +
          (1 | survey_year) +
          (1 | waterbody:survey_year) +
          offset(effort_h),
        family = poisson,
        data = x,
        iter = iter,
        warmup = warmup,
        chains = chains,
        cores = cores
      )
      
      if (cache_model)
        qsave(mod, file = paste0("outputs/fitted/cpue-mod-", species_clean, ".qs"))
      
    } else {
      
      mod <- stan_glmer(
        catch ~ log_cpue_ym1 +
          (1 | waterbody / reach_no) +
          (1 | id_site) +
          (1 | survey_year) +
          (1 | waterbody:survey_year) +
          offset(effort_h),
        family = poisson,
        data = x |> 
          left_join(
            x |>
              mutate(
                survey_year = survey_year + 1,
                log_cpue_ym1 = log(catch + 1) - log(effort_h)
              ) |>
              select(id_site, survey_year, log_cpue_ym1),
            by = c("id_site", "survey_year")
          ) |>
          filter(!is.na(log_cpue_ym1)),
        iter = iter,
        warmup = warmup,
        chains = chains,
        cores = cores
      )
      
      if (cache_model)
        qsave(mod, file = paste0("outputs/fitted/cpue-mod-", species_clean, ".qs"))
      
    }
    
  }
  
  # return
  mod
  
}  

# function to calculate mid, lower, and upper bounds from simulated
#   trajectories
summarise_sim <- function(x, y, subset, probs, growth_rate = TRUE, zscale = TRUE) {
  
  # pull out abundances for the subset of ages/stages
  abund <- apply(x[, subset, , drop = FALSE], c(1, 3), sum)
  
  # return growth rates if required (raw abundances returned otherwise)
  if (growth_rate) {
    
    # skip if all values are zero
    if (!all(abund == 0)) {
      
      # calculate the pop growth rate first, then summarise this if
      #   z-scaling
      # fill zeros with half of the min observed value
      abund[abund == 0] <- min(abund[abund > 0], na.rm = TRUE) / 2.0
      
      # calculate pop growth rate
      abund <- abund / abund[, c(1L, seq_len(ncol(abund) - 1L))]
      
      # calculate zscores
      if (zscale) {
        abund <- abund - mean(abund)
        abund <- abund / sd(abund)
      }
      
    }
    
    # remove first column (duplicated in previous to avoid errors when
    #    length(abund) == 1, suspect this will still error but not
    #    relevant to this analysis)
    abund <- abund[, -1]
    
  } else {
    
    # just calculate z-scores
    if (zscale) {
      abund <- abund - mean(abund)
      abund <- abund / sd(abund)
    }
    
  }
  
  # collate raw predicted values, dropping first column
  out <- tibble(
    y,
    mid = apply(abund, 2, median),
    lower = apply(abund, 2, quantile, probs = probs[1]),
    upper = apply(abund, 2, quantile, probs = probs[2])
  )
  
  # return
  out
  
}

# function to calculate pop growth rates from observed data and compare
#   to simulated values
add_cpue <- function(
    sim,
    cpue_mod, 
    sim_years = 2010:2023,
    probs = c(0.1, 0.9)
) {
  
  #    generate new samples from the fitted posterior for each year/waterbody,
  #    setting previous cpue to 0 to estimate growth rate directly
  #    (no need to divide by catch_ym1)
  
  newdata <- cpue_mod$data |> 
    distinct(waterbody, reach_no, survey_year) |>
    filter(!is.na(reach_no)) |>
    mutate(
      log_cpue_ym1 = 0,
      effort_h = 1,
      id_site = "abc"
    )
  cpue_pred <- posterior_epred(
    cpue_mod, 
    newdata = newdata,
    re.form = ~ (1 | waterbody / reach_no) +
      (1 | survey_year)
  )
  cpue_ar1 <- tibble(
    newdata,
    cpue = apply(cpue_pred, 2, median),
    lower = apply(cpue_pred, 2, quantile, probs = probs[1]),
    upper = apply(cpue_pred, 2, quantile, probs = probs[2])
  )
  
  # add reach info and rename cpue field
  cpue_ar1 <- cpue_ar1 |>
    mutate(
      waterbody = paste0(
        tolower(gsub(" ", "_", waterbody)),
        "_r",
        reach_no
      )
    ) |>
    select(-reach_no) |>
    rename(growth_rate = cpue)
  
  # add in survey year info to simulated values
  nwaterbody <- sim |> pull(waterbody) |> unique() |> length()
  survey_year_tmp <- rep(sim_years, nwaterbody)
  
  # but correct this if it's a future scenario which has many many layers to it
  if (length(survey_year_tmp) != nrow(sim))
    survey_year_tmp <- rep(sim_years, nrow(sim) / length(sim_years))
  
  # and add this value in  
  sim <- sim |> mutate(survey_year = survey_year_tmp)
  
  # z-scale it all
  cpue_std <- cpue_ar1 |>
    group_by(waterbody) |>
    summarise(
      center = mean(growth_rate, na.rm = TRUE),
      width = sd(growth_rate, na.rm = TRUE)
    )
  cpue_ar1 <- cpue_ar1 |>
    left_join(cpue_std, by = "waterbody") |>
    mutate(
      growth_rate_z = (growth_rate - center) / width,
      lower_z = (lower - center) / width,
      upper_z = (upper - center) / width
    ) |>
    select(waterbody, survey_year, growth_rate_z, lower_z, upper_z)
  
  # return this value joined to simulated pop growth rates
  sim |> 
    left_join(cpue_ar1, by = c("waterbody", "survey_year")) |>
    pivot_longer(
      cols = c(mid, growth_rate_z, lower, lower_z, upper, upper_z),
      values_to = "value",
      names_to = "type"
    ) |>
    mutate(
      category = ifelse(grepl("_z", type), "Observed", "Simulated"),
      type = gsub("_z", "", type),
      type = gsub("growth_rate", "mid", type)
    ) |>
    pivot_wider(
      id_cols = c(waterbody, survey_year, category),
      names_from = type,
      values_from = value
    )
  
}

# calculate summary metrics
calculate_val_metrics <- function(
    x, cpue_mod, 
    subset, 
    sim_years,
    probs = c(0.1, 0.9), 
    recruit = FALSE
) {
  
  # use functions above to summarise the simulated population trajectories
  x <- mapply(
    summarise_sim, 
    x = x$sims,
    y = lapply(
      seq_len(nrow(x$scenario)),
      \(i) x$scenario[i, ]
    ),
    MoreArgs = list(
      subset = subset, probs = probs, growth_rate = !recruit
    ),
    SIMPLIFY = FALSE
  )
  x <- bind_rows(x)
  
  # add estimated CPUE
  x <- add_cpue(
    sim = x,
    cpue_mod = cpue_mod,
    sim_years = sim_years,
    probs = probs
  )
  
  # split out the modelled and observed values and calculate residuals
  x <- x |>
    select(waterbody, survey_year, category, mid) |>
    pivot_wider(
      id_cols = c(waterbody, survey_year),
      names_from = category,
      values_from = mid
    ) |>
    mutate(eps = Simulated - Observed)
  
  # calculate all the metrics
  x |>
    group_by(waterbody) |>
    summarise(
      r = ifelse(
        !all(is.na(Observed)), 
        cor(Simulated, Observed, use = "complete"),
        NA
      ),
      md = mean(eps, na.rm = TRUE),
      rmse = sqrt(mean(eps ^ 2, na.rm = TRUE)),
      sign = sum(sign(Observed) == sign(Simulated), na.rm = TRUE) / length(Observed)
    )
  
}

# helpers to tidy names in a plot
tidy_names <- function(x) {
  x <- gsub("_", " ", x)
  x <- strsplit(x, split = " ")
  y <- sapply(x, \(.x) gsub("r", "Reach ", .x[3]))
  x <- lapply(x, \(.x) paste0(toupper(substr(.x[1:2], 1, 1)), substr(.x[1:2], 2, nchar(.x))))
  paste(sapply(x, paste, collapse = " "), y, sep = ": ")
}
metric_names <- c(
  "r" = "r",
  "md" = "MD",
  "rmse" = "RMSE",
  "sign" = "Sign"
)

# function to plot validation metrics for one or more rivers at a time
plot_metric <- function(x) {
  
  # prepare data
  x <- x |>
    pivot_longer(
      cols = c(r, md, rmse, sign),
      names_to = "name",
      values_to = "value"
    ) |>
    mutate(
      waterbody = tidy_names(waterbody),
      metric = metric_names[name],
      metric = factor(metric, levels = c("r", "Sign", "RMSE", "MD")),
      species = factor(
        species,
        levels = c(
          "Murray Cod", "Murray Cod (young of year)"
        )
      )
    )
  
  # work out level and labels for text
  x <- x |>
    left_join(
      x |> 
        group_by(metric) |> 
        summarise(level = median(value, na.rm = TRUE)),
      by = "metric"
    ) |>
    mutate(label = ifelse(is.na(value), "*", ""))
  
  # set a width based on species
  width_set <- 0.45

  # create a dummy data set that adds lines at 1/0 for the different facets
  dummy <- tibble(
    metric = factor(c("r", "MD", "RMSE", "Sign"), levels = c("r", "Sign", "RMSE", "MD")),
    height = c(1, 0, 0, 1)
  )
  
  # plot 
  p <- x |>
    ggplot(aes(y = value, x = waterbody, fill = species)) + 
    geom_bar(position = position_dodge(width = 0.9, preserve = "single"), stat = "identity") +
    geom_text(
      aes(y = level, label = label), 
      position = position_dodge(width = width_set, preserve = "single")
    ) +
    geom_hline(data = dummy, aes(yintercept = height), col = "gray30", linewidth = 1.25) +
    ylab("Value") +
    xlab("Waterbody") +
    facet_wrap( ~ metric, scales = "free") +
    scale_fill_brewer(palette = "Set2", name = "") +
    ggthemes::theme_hc() +
    theme(
      legend.text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 60, hjust = 1),
      panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
      strip.background = element_rect(fill = "white")
    )
  
  # remove legend if just one species
  if (length(unique(x$species)) == 1)
    p <- p + theme(legend.position = "none")
  
  # return
  p
  
}

# river names lookup
.river_lookup <- c(
  "broken_creek_r4" = "Broken Creek (Reach 4)",
  "broken_river_r3" = "Broken River (Reach 3)",
  "campaspe_river_r4" = "Campaspe River (Reach 4)",
  "goulburn_river_r4" = "Goulburn River (Reach 4)",
  "loddon_river_r4" = "Loddon River (Reach 4)",
  "ovens_river_r5" = "Ovens River (Reach 5)"
)

# function to create abundance hindcast plots from simulated and observed data
plot_hindcasts <- function(
    x, 
    cpue, 
    subset, 
    sim_years, 
    probs = c(0.1, 0.9), 
    recruit = FALSE,
    rb = FALSE
) {
  
  # use functions above to summarise the simulated population trajectories
  x <- mapply(
    summarise_sim, 
    x = x$sims,
    y = lapply(
      seq_len(nrow(x$scenario)),
      \(i) x$scenario[i, ]
    ),
    MoreArgs = list(
      subset = subset, probs = probs, growth_rate = !recruit
    ),
    SIMPLIFY = FALSE
  )
  x <- bind_rows(x)
  
  # add estimated CPUE
  x <- add_cpue(
    sim = x,
    cpue_mod = cpue,
    sim_years = sim_years,
    probs = probs
  )
  
  # set up base plot
  p <- x |>
    mutate(waterbody = .river_lookup[waterbody]) |>
    ggplot(aes(x = survey_year, y = mid, col = category, group = category)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_line(position = position_dodge(width = 0.2)) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper), 
      width = 0.2, 
      position = position_dodge(width = 0.2)
    ) +
    scale_color_brewer(
      name = "",
      palette = "Set2"
    ) +
    xlab("Water year") +
    ggthemes::theme_hc() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 8),
      panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
      strip.background = element_rect(fill = "white")
    ) + 
    facet_wrap( ~ waterbody, scales = "free")
  
  if (recruit) {
    p <- p + ylab("Scaled recruitment")
  } else {
    p <- p + ylab("Scaled population growth rate")
  }
  
  # and return
  p
  
}

# function to plot near-term forecasts from start to final observed year
plot_forecasts <- function(
    x, subset, probs, system = NULL, target = NULL, climate = NULL, marker = NULL, rb = FALSE
) {
  
  # use functions above to summarise the simulated population trajectories
  nscn <- nrow(x$scenario)
  x <- mapply(
    summarise_sim, 
    x = x$sims,
    y = lapply(
      seq_len(nscn),
      \(i) x$scenario[i, ]
    ),
    MoreArgs = list(
      subset = subset, probs = probs, growth_rate = FALSE, zscale = FALSE
    ),
    SIMPLIFY = FALSE
  )
  x <- bind_rows(x)
  
  # add year information
  x <- x |>
    mutate(survey_year = rep(c(2023:2025), times = nscn))
  
  # clean up variable values
  x <- x |>
    mutate(
      future = factor(
        future,
        levels = c("dry", "ave", "wet"), 
        labels = c("Dry (2023/2024)", "Ave. (2023/2024)", "Wet (2023/2024)")
      ),
      future_next = factor(
        future_next,
        levels = c("dry", "ave", "wet"),
        labels = c("Dry (2024/2025)", "Ave. (2024/2025)", "Wet (2024/2025)")
      ),
      scenario = factor(
        scenario,
        levels = c("none", "baseflow", "fresh"),
        labels = c("None", "Baseflows", "Freshes")
      ),
      scenario_next = factor(
        scenario_next,
        levels = c("none", "baseflow", "fresh"),
        labels = c("None", "Baseflows", "Freshes")
      )
    )
  
  # set target to latest year if not specified
  if (is.null(target))
    target <- max(x$survey_year)
  
  # and set a flag to work out plots below
  one_step_ahead <- target != max(x$survey_year)
  
  # filter to target water year
  x <- x |> filter(survey_year == target)
  
  # filter to a single system if required
  if (!is.null(system))
    x <- x |> filter(waterbody == system)
  
  # add labels to mark "poor" performers if required
  x <- x |>
    mutate(waterbody = .river_lookup[waterbody])
  if (!is.null(marker)) {
    new_names <- sort(unique(x$waterbody))
    names(new_names) <- new_names
    new_names[marker] <- paste0(new_names[marker], "*")
    x <- x |>
      mutate(waterbody = new_names[waterbody])
  }
  
  # if showing all climates, need a plot that expands out
  if (is.null(climate)) {
    
    # two options: simpler plot if only showing one step ahead
    if (one_step_ahead) {
      
      p <- x |>
        filter(
          future_next == "Ave. (2024/2025)",
          scenario_next == "None"
        ) |>
        mutate(
          future = factor(future, labels = c("Dry", "Ave.", "Wet"))
        ) |>
        ggplot(aes(y = mid, x = future, fill = scenario)) +
        geom_bar(position = position_dodge(0.9), stat = "identity") +
        geom_errorbar(
          aes(ymin = lower, ymax = upper),
          position = position_dodge(0.9),
          col = "black",
          width = 0.2
        ) +
        xlab("Climate state") +
        ylab("Abundance") +
        scale_fill_brewer(name = "Flow priority", palette = "Set2") +
        ggthemes::theme_hc() +
        theme(
          legend.position = "bottom",
          axis.text = element_text(size = 8),
          panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
          strip.background = element_rect(fill = "white")
        )
      
      # and add a facet wrap if system is not specified
      if (is.null(system))
        p <- p + facet_wrap( ~ waterbody, scales = "free")
      
    } else {
      
      # must specify system or this plot gets unwieldy
      if (is.null(system))
        stop(
          "system must be specified when plotting multiple ",
          "time steps ahead",
          call. = FALSE
        )
      
      #  plot it
      p <- x |>
        mutate(survey_year = factor(survey_year)) |>
        ggplot(aes(y = mid, x = scenario_next, fill = scenario)) +
        geom_bar(position = position_dodge(0.9), stat = "identity") +
        geom_errorbar(
          aes(ymin = lower, ymax = upper),
          position = position_dodge(0.9),
          col = "black",
          width = 0.2
        ) +
        xlab("Flow priority (2024/2025)") +
        ylab("Abundance") +
        scale_fill_brewer(name = "Flow priority (2023/2024)", palette = "Set2") +
        facet_grid(future_next ~ future) +
        ggthemes::theme_hc() +
        theme(
          legend.position = "bottom",
          axis.text = element_text(size = 8),
          axis.text.x = element_text(angle = 60, hjust = 1),
          panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
          strip.background = element_rect(fill = "white")
        )
      
    }
    
  } else {
    
    # show just a single climate for 2023/24 then all for 2024/25
    #  plot it
    p <- x |>
      filter(future == climate) |>
      mutate(survey_year = factor(survey_year)) |>
      ggplot(aes(y = mid, x = scenario_next, fill = scenario)) +
      geom_bar(position = position_dodge(0.9), stat = "identity") +
      geom_errorbar(
        aes(ymin = lower, ymax = upper),
        position = position_dodge(0.9),
        col = "black",
        width = 0.2
      ) +
      xlab("Flow priority (2024/2025)") +
      ylab("Abundance") +
      scale_fill_brewer(name = "Flow priority (2023/2024)", palette = "Set2") +
      facet_grid( ~ future_next) +
      ggthemes::theme_hc() +
      theme(
        legend.position = "bottom",
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
        strip.background = element_rect(fill = "white")
      )
    
  }
  
  # and return
  p
  
}

# function to extract near-term forecasts one year ahead
extract_forecasts <- function(
    x, subset, probs, climate = NULL
) {
  
  # use functions above to summarise the simulated population trajectories
  nscn <- nrow(x$scenario)
  x <- mapply(
    summarise_sim, 
    x = x$sims,
    y = lapply(
      seq_len(nscn),
      \(i) x$scenario[i, ]
    ),
    MoreArgs = list(
      subset = subset, probs = probs, growth_rate = FALSE, zscale = FALSE
    ),
    SIMPLIFY = FALSE
  )
  x <- bind_rows(x)
  
  # add year information
  x <- x |>
    mutate(survey_year = rep(c(2023:2025), times = nscn))
  
  # clean up variable values
  x <- x |>
    mutate(
      future = factor(
        future,
        levels = c("dry", "ave", "wet"), 
        labels = c("Dry (2023/2024)", "Ave. (2023/2024)", "Wet (2023/2024)")
      ),
      future_next = factor(
        future_next,
        levels = c("dry", "ave", "wet"),
        labels = c("Dry (2024/2025)", "Ave. (2024/2025)", "Wet (2024/2025)")
      ),
      scenario = factor(
        scenario,
        levels = c("none", "baseflow", "fresh"),
        labels = c("None", "Baseflows", "Freshes")
      ),
      scenario_next = factor(
        scenario_next,
        levels = c("none", "baseflow", "fresh"),
        labels = c("None", "Baseflows", "Freshes")
      )
    )
  
  # filter to target water year and scenarios
  x <- x |> 
    filter(
      survey_year %in% c(2023, 2024),
      future_next == "Ave. (2024/2025)",
      scenario_next == "None",
      scenario == "Baseflows",
      future == "Ave. (2023/2024)"
    )
  
  # calculate change from one year to the next
  x <- x |>
    mutate(
      lower = (lower - mean(mid)) / sd(mid),
      upper = (upper - mean(mid)) / sd(mid),
      mid = (mid - mean(mid)) / sd(mid)
    ) |>
    pivot_wider(
      id_cols = c(waterbody),
      names_from = survey_year,
      values_from = c(lower, upper, mid)
    ) |>
    mutate(
      abund_change_lower = lower_2024 - upper_2023,
      abund_change_upper = upper_2024 - lower_2023,
      abund_change = mid_2024 - mid_2023
    )
  
  # and return
  x
  
}

# plot forecasts against validation data
plot_onestep <- function(mod, obs, species, ...) {
  
  # do the work to line things up
  x <- mod |>
    select(waterbody, scientific_name, abund_change_lower, abund_change_upper, abund_change) |>
    left_join(
      obs |>
        ungroup() |>
        mutate(
          cpue_change_lower = (cpue_change_lower - mean(cpue_change)) / sd(cpue_change),
          cpue_change_upper = (cpue_change_upper - mean(cpue_change)) / sd(cpue_change),
          cpue_change = (cpue_change - mean(cpue_change)) / sd(cpue_change),
          waterbody = paste0(tolower(gsub(" ", "_", waterbody)), "_r", reach_no)
        ) |>
        select(waterbody, scientific_name, cpue_change_lower, cpue_change_upper, cpue_change),
      by = c("waterbody", "scientific_name")
    ) |>
    pivot_longer(
      cols = contains("_change"),
      values_to = "change",
      names_to = "type"
    ) |>
    mutate(
      level = ifelse(grepl("lower", type), "lower", ifelse(grepl("upper", type), "upper", "mid")),
      type = ifelse(grepl("abund_", type), "Modelled", "Observed")
    ) |>
    pivot_wider(
      id_cols = c(waterbody, scientific_name, type),
      values_from = change,
      names_from = level
    )
  
  # plot it
  x |>
    mutate(
      waterbody = .river_lookup[waterbody],
      scientific_name = factor(
        scientific_name,
        levels = c("Maccullochella peelii"),
        labels = c("Murray Cod")
      )
    ) |>
    ggplot(aes(y = mid, x = waterbody, fill = type, ymin = lower, ymax = upper)) +
    geom_hline(yintercept = 0, col = "gray60", linetype = "dashed") +
    geom_bar(position = position_dodge(0.9), stat = "identity") +
    geom_errorbar(position = position_dodge(0.9), width = 0.25) +
    xlab("") +
    ylab("Relative change (2023 to 2024)") +
    scale_fill_brewer(palette = "Set2", name = "") +
    facet_wrap( ~ scientific_name, scales = "free", nrow = 3) +
    theme(
      legend.position = "bottom", 
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  
}

# plot forecasts against validation data
plot_onestep_trend <- function(mod, obs, ytitle = "Relative recruitment", ...) {
  
  # do the work to line things up
  x <- mod |>
    select(waterbody, scientific_name, lower_2023, lower_2024, upper_2023, upper_2024, mid_2023, mid_2024) |>
    left_join(
      obs |>
        ungroup() |>
        mutate(
          cpue_2023_lower = (cpue_2023_lower - mean(cpue_2023)) / sd(cpue_2023),
          cpue_2023_upper = (cpue_2023_upper - mean(cpue_2023)) / sd(cpue_2023),
          cpue_2024_lower = (cpue_2024_lower - mean(cpue_2023)) / sd(cpue_2023),
          cpue_2024_upper = (cpue_2024_upper - mean(cpue_2023)) / sd(cpue_2023),
          cpue_2024 = (cpue_2024 - mean(cpue_2023)) / sd(cpue_2023),
          cpue_2023 = (cpue_2023 - mean(cpue_2023)) / sd(cpue_2023),
          waterbody = paste0(tolower(gsub(" ", "_", waterbody)), "_r", reach_no)
        ) |>
        select(
          waterbody, scientific_name, cpue_2023_lower,
          cpue_2023_upper, cpue_2024_lower,
          cpue_2024_upper, cpue_2024, cpue_2023
        ),
      by = c("waterbody", "scientific_name")
    ) |>
    pivot_longer(
      cols = c(contains("_2023"), contains("_2024")),
      values_to = "abund",
      names_to = "type"
    ) |>
    mutate(
      mod = ifelse(grepl("cpue_", type), "Observed", "Modelled"),
      type = gsub("cpue_", "", type),
      year = ifelse(
        mod == "Modelled", 
        sapply(strsplit(type, "_"), \(x) x[2]),
        sapply(strsplit(type, "_"), \(x) x[1])
      ),
      level = ifelse(
        mod == "Modelled", 
        sapply(strsplit(type, "_"), \(x) x[1]),
        sapply(strsplit(type, "_"), \(x) x[2])
      ),
      level = ifelse(is.na(level), "mid", level)
    ) |>
    pivot_wider(
      id_cols = c(waterbody, scientific_name, year, mod),
      values_from = abund,
      names_from = level
    )
  
  # plot it
  x |>
    filter(year == 2024) |>
    mutate(
      waterbody = .river_lookup[waterbody],
      scientific_name = factor(
        scientific_name,
        levels = c("Maccullochella peelii"),
        labels = c("Murray Cod")
      )
    ) |>
    ggplot(aes(y = mid, x = waterbody, fill = mod, ymin = lower, ymax = upper)) +
    geom_bar(stat = "identity", position = position_dodge(0.9)) +
    geom_errorbar(width = 0.25, position = position_dodge(0.9)) +
    xlab("") +
    ylab(ytitle) +
    scale_fill_brewer(palette = "Set2", name = "") +
    facet_wrap( ~ scientific_name, scales = "free", nrow = 1) +
    theme(
      legend.position = "bottom", 
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  
}


# function to create abundance forecast plots from simulated and observed data
#   once new observed data are available
plot_forecasts_update <- function(
    x, 
    cpue, 
    subset, 
    sciname,
    sim_years, 
    survey_max,
    probs = c(0.1, 0.9), 
    recruit = FALSE,
    scenario_set = "baseflow",
    future_set = "ave"
) {
  
  # use functions above to summarise the simulated population trajectories
  for (w in seq_along(x)) {
    
    x[[w]] <- mapply(
      summarise_sim, 
      x = x[[w]]$sims,
      y = lapply(
        seq_len(nrow(x[[w]]$scenario)),
        \(i) x[[w]]$scenario[i, ]
      ),
      MoreArgs = list(
        subset = subset[[w]], probs = probs, growth_rate = !recruit
      ),
      SIMPLIFY = FALSE
    )
    x[[w]] <- bind_rows(x[[w]])
    
    # add estimated CPUE
    x[[w]] <- add_cpue_update(
      sim = x[[w]],
      cpue_mod = cpue[[w]],
      sim_years = sim_years,
      probs = c(0.4, 0.6)
    )
    
    # add scientific names
    x[[w]] <- x[[w]] |> mutate(scientific_name = sciname[w])
    
  }
  
  # flatten
  x <- bind_rows(x)
  
  # remove any years wtihout new data
  x <- x |> filter(survey_year <= survey_max)
  
  # only need one set of future_next because we're not plotting it
  x <- x |> filter(scenario_next == "baseflow", future_next == "ave")
  
  # set up base plot
  p <- x |>
    filter(
      scenario == !!scenario_set, 
      future == !!future_set,
      # !(category == "Observed" & !future %in% !!future_set),
      survey_year == survey_max
    ) |>
    mutate(
      waterbody = .river_lookup[waterbody],
      scientific_name = factor(
        scientific_name,
        levels = c("Maccullochella peelii"),
        labels = c("Murray Cod")
      ),
      category = factor(category, levels = c("Simulated", "Observed"))
    ) |>
    ggplot(aes(x = waterbody, y = mid, fill = category)) +
    geom_bar(stat = "identity", position = position_dodge(0.9)) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper), 
      width = 0.2, 
      position = position_dodge(width = 0.9)
    ) +
    scale_fill_brewer(
      name = "",
      palette = "Set2"
    ) +
    xlab("Waterbody") +
    ggthemes::theme_hc() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 8, angle = 35, vjust = 0, hjust = 1),
      panel.border = element_rect(fill = NA, colour = "gray30", linetype = 1),
      strip.background = element_rect(fill = "white")
    ) +
    facet_wrap( ~ scientific_name, scales = "free", nrow = length(sciname))
  
  if (recruit) {
    p <- p + ylab("Scaled recruitment")
  } else {
    p <- p + ylab("Scaled population growth rate")
  }
  
  # and return
  p
  
}


# function to calculate pop growth rates from obsered data and compare
#   to simulated values
add_cpue_update <- function(
    sim,
    cpue_mod, 
    sim_years = 2010:2023,
    probs = c(0.1, 0.9)
) {
  
  #    generate new samples from the fitted posterior for each year/waterbody,
  #    setting previous cpue to 0 to estimate growth rate directly
  #    (no need to divide by catch_ym1)
  
  newdata <- cpue_mod$data |> 
    distinct(waterbody, reach_no, survey_year) |>
    filter(!is.na(reach_no)) |>
    mutate(
      log_cpue_ym1 = 0,
      effort_h = 1,
      id_site = "abc"
    )
  cpue_pred <- posterior_epred(
    cpue_mod, 
    newdata = newdata,
    re.form = ~ (1 | waterbody / reach_no) +
      (1 | survey_year)
  )
  cpue_ar1 <- tibble(
    newdata,
    cpue = apply(cpue_pred, 2, median),
    lower = apply(cpue_pred, 2, quantile, probs = probs[1]),
    upper = apply(cpue_pred, 2, quantile, probs = probs[2])
  )
  
  # add reach info and rename cpue field
  cpue_ar1 <- cpue_ar1 |>
    mutate(
      waterbody = paste0(
        tolower(gsub(" ", "_", waterbody)),
        "_r",
        reach_no
      )
    ) |>
    select(-reach_no) |>
    rename(growth_rate = cpue)
  
  # add in survey year info to simulated values
  nwaterbody <- sim |> pull(waterbody) |> unique() |> length()
  survey_year_tmp <- rep(sim_years, nwaterbody)
  
  # but correct this if it's a future scenario which has many many layers to it
  if (length(survey_year_tmp) != nrow(sim))
    survey_year_tmp <- rep(sim_years, nrow(sim) / length(sim_years))
  
  # and add this value in  
  sim <- sim |> mutate(survey_year = survey_year_tmp)
  
  # z-scale it all
  cpue_std <- cpue_ar1 |>
    group_by(waterbody) |>
    summarise(
      center = mean(growth_rate, na.rm = TRUE),
      width = sd(growth_rate, na.rm = TRUE)
    )
  cpue_ar1 <- cpue_ar1 |>
    left_join(cpue_std, by = "waterbody") |>
    mutate(
      growth_rate_z = (growth_rate - center) / width,
      lower_z = (lower - center) / width,
      upper_z = (upper - center) / width
    ) |>
    select(waterbody, survey_year, growth_rate_z, lower_z, upper_z)
  
  # return this value joined to simulated pop growth rates but remove any
  #   years not surveyed yet
  sim |> 
    left_join(cpue_ar1, by = c("waterbody", "survey_year")) |>
    pivot_longer(
      cols = c(mid, growth_rate_z, lower, lower_z, upper, upper_z),
      values_to = "value",
      names_to = "type"
    ) |>
    mutate(
      category = ifelse(grepl("_z", type), "Observed", "Simulated"),
      type = gsub("_z", "", type),
      type = gsub("growth_rate", "mid", type)
    ) |>
    pivot_wider(
      id_cols = c(waterbody, future, future_next, scenario, scenario_next, survey_year, category),
      names_from = type,
      values_from = value
    )
  
}
