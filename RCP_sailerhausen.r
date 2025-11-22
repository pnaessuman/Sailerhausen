source("External/Vs.R")
library(VSLiteR)
library(dplR) # for reading in dendro data
library(ggplot2) # for plotting later
library(tidyverse)
library(DEoptim)

b1ac <- read.rwl("Data/b1ac.rwl")
b1apl  <- read.rwl("Data/b1apl.rwl")
b1aps <- read.rwl("Data/b1aps.rwl")
b1fs <- read.rwl("Data/b1fs.rwl")
b1pa <- read.rwl("Data/b1pa.rwl")
b1qp <- read.rwl("Data/b1qp.rwl")
climate_b1 <- read_csv2("Data/B1_Sailerhausen.csv")
cordex_b1 <- read_csv2("Data/cordex_b1_padded.csv")

b2ac <- read.rwl("Data/b2ac.rwl")
b2pa <- read.rwl("Data/b2pa.rwl")
b2qp <- read.rwl("Data/b2qp.rwl")
climate_b2 <- read_csv2("Data/B2_Sailerhausen.csv")
cordex_b2 <- read_csv2("Data/cordex_b2_padded.csv")

b1ac_d <- detrend(b1ac, method = "Spline", nyrs = 32)
b1apl_d  <- detrend(b1apl, method = "Spline", nyrs = 32)
b1aps_d <- detrend(b1aps, method = "Spline", nyrs = 32)
b1fs_d <- detrend(b1fs, method = "Spline", nyrs = 32)
b1pa_d <- detrend(b1pa, method = "Spline", nyrs = 32)
b1qp_d <- detrend(b1qp, method = "Spline", nyrs = 32)
b2ac_d <- detrend(b2ac, method = "Spline", nyrs = 32)
b2pa_d <- detrend(b2pa, method = "Spline", nyrs = 32)
b2qp_d <- detrend(b2qp, method = "Spline", nyrs = 32)

b1ac_c <- chron(b1ac_d)
b1apl_c <- chron(b1apl_d)
b1aps_c <- chron(b1aps_d)
b1fs_c <- chron(b1fs_d)
b1pa_c <- chron(b1pa_d)
b1qp_c <- chron(b1qp_d)
b2ac_c <- chron(b2ac_d)
b2pa_c <- chron(b2pa_d)
b2qp_c <- chron(b2qp_d)

# B1ac - Acer campestre
input_historic_b1ac <- make_vsinput_historic(b1ac_c, climate_b1)

b1ac_params <- vs_params(input_historic_b1ac$trw,
                         input_historic_b1ac$tmean,
                         input_historic_b1ac$prec,
                         input_historic_b1ac$syear,
                         input_historic_b1ac$eyear,
                         .phi = 50)

input_transient_b1ac <- make_vsinput_transient(climate_b1)

b1ac_forward <- vs_run_forward(b1ac_params,
                                input_transient_b1ac$tmean,
                                input_transient_b1ac$prec,
                                input_transient_b1ac$syear,
                                input_transient_b1ac$eyear,
                                .phi = 50)



# Method 1: Simple and reliable approach
all_rcp_projections <- list()

for(rcp_scenario in c("RCP_85", "RCP_45", "RCP_26")) {
  
  # Get data for this specific RCP
  rcp_data <- cordex_b1 %>% 
    filter(year > 2021, rcp == rcp_scenario) %>% 
    group_by(year, month) %>% 
    summarise(
      tmean = mean(tmean),
      prec = mean(prec)
    ) %>% 
    ungroup() %>% 
    make_vsinput_transient()
  
  # Run projection
  projection <- vs_run_forward(b1ac_params,
                               rcp_data$tmean,
                               rcp_data$prec,
                               rcp_data$syear,
                               rcp_data$eyear,
                               .phi = 50)
  
  names(projection) <- c("year", "projection")
  projection$rcp <- rcp_scenario
  
  all_rcp_projections[[rcp_scenario]] <- projection
}

# Combine into one dataframe
combined_data <- bind_rows(all_rcp_projections)

# Now plot
ggplot(combined_data, aes(year, projection)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_smooth(color = "red", se = FALSE) +
  labs(
    title = "B1ac Future Projection",
    subtitle = "All RCP Scenarios",
    x = "Year",
    y = "Ring Width Index"
  ) +
  theme_minimal() +
  facet_wrap(~ rcp, ncol = 3)
