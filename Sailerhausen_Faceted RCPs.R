# Multi-Species Tree Ring Projection Analysis
# Sites: Sailershausen (B1-B2)

# Load required libraries
source("External/Vs.R")
library(VSLiteR)
library(dplR)
library(ggplot2)
library(tidyverse)
library(DEoptim)

# Load tree ring data
b1ac <- read.rwl("Data/b1ac.rwl")
b1apl <- read.rwl("Data/b1apl.rwl")
b1aps <- read.rwl("Data/b1aps.rwl")
b1fs <- read.rwl("Data/b1fs.rwl")
b1pa <- read.rwl("Data/b1pa.rwl")
b1qp <- read.rwl("Data/b1qp.rwl")
b2ac <- read.rwl("Data/b2ac.rwl")
b2pa <- read.rwl("Data/b2pa.rwl")
b2qp <- read.rwl("Data/b2qp.rwl")

# Load climate data
climate_b1 <- read_csv2("Data/B1_Sailerhausen.csv")
climate_b2 <- read_csv2("Data/B2_Sailerhausen.csv")
cordex_b1 <- read_csv2("Data/lfu_cordex_ensemble_monthly_b1.csv")
cordex_b2 <- read_csv2("Data/lfu_cordex_ensemble_monthly_b2.csv")

# Detrend tree ring data
b1ac_d <- detrend(b1ac, method = "Spline", nyrs = 32)
b1apl_d <- detrend(b1apl, method = "Spline", nyrs = 32)
b1aps_d <- detrend(b1aps, method = "Spline", nyrs = 32)
b1fs_d <- detrend(b1fs, method = "Spline", nyrs = 32)
b1pa_d <- detrend(b1pa, method = "Spline", nyrs = 32)
b1qp_d <- detrend(b1qp, method = "Spline", nyrs = 32)
b2ac_d <- detrend(b2ac, method = "Spline", nyrs = 32)
b2pa_d <- detrend(b2pa, method = "Spline", nyrs = 32)
b2qp_d <- detrend(b2qp, method = "Spline", nyrs = 32)

# Create chronologies
b1ac_c <- chron(b1ac_d)
b1apl_c <- chron(b1apl_d)
b1aps_c <- chron(b1aps_d)
b1fs_c <- chron(b1fs_d)
b1pa_c <- chron(b1pa_d)
b1qp_c <- chron(b1qp_d)
b2ac_c <- chron(b2ac_d)
b2pa_c <- chron(b2pa_d)
b2qp_c <- chron(b2qp_d)

# 2. PARAMETER ESTIMATION FUNCTION
estimate_params <- function(chronology, climate_data, phi = 50, .iter = 200) {
  input_historic <- make_vsinput_historic(chronology, climate_data)
  
  params <- vs_params(
    input_historic$trw,
    input_historic$tmean,
    input_historic$prec,
    input_historic$syear,
    input_historic$eyear,
    .phi = phi,
    iter = .iter
  )
  
  return(params)
}

# 3. PROJECTION FUNCTION
generate_projection <- function(params, cordex_data, phi = 50) {
  
  complete_years <- function(x) {
    cy <- x |> 
      group_by(year) |> 
      summarise(n_months = length(month)) |> 
      filter(n_months == 12) |> 
      _$year
    
    x |> filter(year %in% cy)
  }
  
  do_project <- function(x) {
    vs_run_forward(
      params,
      x$temp,
      x$prec,
      x$syear,
      x$eyear,
      .phi = phi
    )
  }
  
  do_smooth <- function(x) {
    x$trw_smooth <- lowess(x$trw)$y
    return(x)
  }
  
  projection <- cordex_data %>% 
    filter(year > 2021) %>% 
    group_by(rcp, gcm, rcm) |> 
    nest() |> 
    mutate(
      data = purrr::map(data, complete_years),
      input = purrr::map(data, make_vsinput_transient),
      projection = purrr::map(input, do_project),
      projection = purrr::map(projection, do_smooth)
      ) |> 
    unnest("projection")
  
  return(projection)
}

# 4. ESTIMATE PARAMETERS FOR ALL SPECIES

# B1 species
b1ac_params <- estimate_params(b1ac_c, climate_b1, .iter = 200)
b1apl_params <- estimate_params(b1apl_c, climate_b1, .iter = 200)
b1aps_params <- estimate_params(b1aps_c, climate_b1, .iter = 200)
b1fs_params <- estimate_params(b1fs_c, climate_b1, .iter = 200)
b1pa_params <- estimate_params(b1pa_c, climate_b1, .iter = 200)
b1qp_params <- estimate_params(b1qp_c, climate_b1, .iter = 200)

# B2 species
b2ac_params <- estimate_params(b2ac_c, climate_b2, .iter = 200)
b2pa_params <- estimate_params(b2pa_c, climate_b2, .iter = 200)
b2qp_params <- estimate_params(b2qp_c, climate_b2, .iter = 200)


# 5. GENERATE ALL PROJECTIONS
generate_all_projections <- function(species_params_list, cordex_data, species_names) {
  tibble(params = species_params_list, species = species_names) |> 
    mutate(projection = purrr::map(params, ~generate_projection(.x, cordex_data))) |> 
    unnest(projection) |> 
    select(-params)
}

# 6a. PLOTTING ALL IN ONE

plot_all_in_one <- function(cordex_data, 
                            species_params_list, 
                            species_names) {
  
  projections <- generate_all_projections(
    species_params_list,
    cordex_data,
    species_names
  )
  
  projections |> 
    group_by(species, rcp, year) |> 
    summarise(trw_mean = mean(trw_smooth),
              trw_sd = sd(trw_smooth),
              trw_ci_lower = trw_mean - trw_sd,
              trw_ci_upper = trw_mean + trw_sd) |> 
    ggplot(aes(x = year, y = trw_mean)) +
    geom_ribbon(aes(ymin = trw_ci_lower,
                    ymax = trw_ci_upper,
                    fill = species),
                alpha = 0.1) +
    geom_line(aes(colour = species)) +
    facet_wrap(. ~ rcp)
}


# 7. GENERATE PLOTS FOR BOTH SITES
rcp_scenarios <- c("RCP_26", "RCP_45", "RCP_85")

# Site B1
b1_species_params <- list(b1ac_params, b1apl_params, b1aps_params, 
                          b1fs_params, b1pa_params, b1qp_params)
b1_species_names <- c("Acer campestre", "Acer platanoides", "Acer pseudoplatanus",
                      "Fagus sylvatica", "Prunus avium", "Quercus petraea")

# Site B2
b2_species_params <- list(b2ac_params, b2pa_params, b2qp_params)
b2_species_names <- c("Acer campestre", "Prunus avium", "Quercus petraea")

# 8. CALL plot_all_in_one
plot_all_in_one(cordex_b1, b1_species_params, b1_species_names)
plot_all_in_one(cordex_b2, b2_species_params, b2_species_names)







