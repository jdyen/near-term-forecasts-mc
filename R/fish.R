# function to return carrying capacity as a neat tibble
fetch_carrying_capacity <- function() {
  
  # set carrying capacity
  .carrying_capacity <- list(
    maccullochella_peelii = c(
      "broken_creek_r4" = 50000,
      "broken_river_r3" = 200000,
      "campaspe_river_r4" = 50000,
      "goulburn_river_r4" = 100000,
      "loddon_river_r4" = 200000,
      "ovens_river_r5" = 100000
    )
  )
  
  # collate and return
  tibble(
    species = rep(names(.carrying_capacity), times = sapply(.carrying_capacity, length)),
    waterbody = unlist(lapply(.carrying_capacity, names)),
    carrying_capacity = unlist(.carrying_capacity)
  )
  
}
