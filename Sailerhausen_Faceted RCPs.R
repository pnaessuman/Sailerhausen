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
cordex_b1 <- read_csv2("Data/cordex_b1_padded.csv")
cordex_b2 <- read_csv2("Data/cordex_b2_padded.csv")

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
estimate_params <- function(chronology, climate_data, phi = 50) {
  input_historic <- make_vsinput_historic(chronology, climate_data)
  
  params <- vs_params(
    input_historic$trw,
    input_historic$tmean,
    input_historic$prec,
    input_historic$syear,
    input_historic$eyear,
    .phi = phi
  )
  
  return(params)
}

# 3. PROJECTION FUNCTION
generate_projection <- function(params, cordex_data, rcp_scenario, phi = 50) {
  rcp_input <- cordex_data %>% 
    filter(year > 2021, rcp == rcp_scenario) %>% 
    group_by(year, month) %>% 
    summarise(
      tmean = mean(tmean),
      prec = mean(prec),
      .groups = "drop"
    ) %>% 
    ungroup() %>% 
    make_vsinput_transient()
  
  projection <- vs_run_forward(
    params,
    rcp_input$tmean,
    rcp_input$prec,
    rcp_input$syear,
    rcp_input$eyear,
    .phi = phi
  )
  
  result <- data.frame(
    year = projection[, 1],
    projection = projection[, 2]
  )
  
  return(result)
}

# 4. ESTIMATE PARAMETERS FOR ALL SPECIES

cat("Estimating parameters for all species...\n")

# B1 species
b1ac_params <- estimate_params(b1ac_c, climate_b1)
b1apl_params <- estimate_params(b1apl_c, climate_b1)
b1aps_params <- estimate_params(b1aps_c, climate_b1)
b1fs_params <- estimate_params(b1fs_c, climate_b1)
b1pa_params <- estimate_params(b1pa_c, climate_b1)
b1qp_params <- estimate_params(b1qp_c, climate_b1)

# B2 species
b2ac_params <- estimate_params(b2ac_c, climate_b2)
b2pa_params <- estimate_params(b2pa_c, climate_b2)
b2qp_params <- estimate_params(b2qp_c, climate_b2)


# 5. GENERATE ALL PROJECTIONS
generate_all_projections <- function(species_params_list, cordex_data, 
                                     species_names, rcp_scenario) {
  all_projections <- list()
  
  for (i in seq_along(species_params_list)) {
    projection <- generate_projection(
      species_params_list[[i]], 
      cordex_data, 
      rcp_scenario
    )
    projection$species <- species_names[i]
    all_projections[[i]] <- projection
  }
  
  combined <- bind_rows(all_projections)
  return(combined)
}

# 6. PLOTTING FUNCTIONS
# Function to plot each RCP separately (original version)
plot_site_rcp_separate <- function(site_name, cordex_data, 
                                   species_params_list, 
                                   species_names, rcp_scenarios) {
  
  plot_list <- list()
  
  for (rcp in rcp_scenarios) {
    cat(paste("Generating projection for", site_name, "-", rcp, "\n"))
    
    projections <- generate_all_projections(
      species_params_list,
      cordex_data,
      species_names,
      rcp
    )
    
    p <- ggplot(projections, aes(x = year, y = projection, color = species)) +
      geom_line(linewidth = 0.8, alpha = 0.7) +
      geom_smooth(se = FALSE, linewidth = 1.2, method = "loess") +
      labs(
        title = paste(site_name, "- Future Projections"),
        subtitle = paste("RCP Scenario:", rcp),
        x = "Year",
        y = "Ring Width Index",
        color = "Species"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")
      ) +
      scale_color_brewer(palette = "Set2")
    
    plot_list[[rcp]] <- p
  }
  
  return(plot_list)
}

# Function to plot all RCPs as facets with species in each panel
plot_site_all_rcp_faceted <- function(site_name, cordex_data, species_params_list, 
                                      species_names, rcp_scenarios) {
  
  cat(paste("Generating faceted projection for", site_name, "\n"))
  
  # Generate projections for all RCP scenarios
  all_data <- list()
  
  for (rcp in rcp_scenarios) {
    projections <- generate_all_projections(
      species_params_list,
      cordex_data,
      species_names,
      rcp
    )
    projections$rcp <- rcp
    all_data[[rcp]] <- projections
  }
  
  # Combine all data
  combined_data <- bind_rows(all_data)
  
  # Create proper factor levels for RCP scenarios
  combined_data$rcp <- factor(combined_data$rcp, 
                              levels = c("RCP_26", "RCP_45", "RCP_85"),
                              labels = c("RCP_26", "RCP_45", "RCP_85"))
  
  # Create the faceted plot
  p <- ggplot(combined_data, aes(x = year, y = projection, color = species)) +
    geom_line(linewidth = 1.2, alpha = 0.85) +
    geom_smooth(se = TRUE, alpha = 0.15, linewidth = 1, method = "loess", span = 0.3) +
    facet_wrap(~ rcp, ncol = 3, scales = "fixed") +
    labs(
      title = paste(site_name, "- Future Projections Under All RCP Scenarios"),
      subtitle = "Species responses across different climate scenarios (2020-2100)",
      x = "Year",
      y = "Ring Width Index",
      color = "Species"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, color = "gray40", hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = "gray90", color = "gray70"),
      strip.text = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 9)
    ) +
    scale_color_brewer(palette = "Set2") +
    guides(color = guide_legend(nrow = 2))
  
  return(p)
}

# 7. GENERATE PLOTS FOR BOTH SITES
rcp_scenarios <- c("RCP_26", "RCP_45", "RCP_85")

# Site B1
cat("\n=== Processing Site B1 ===\n")
b1_species_params <- list(b1ac_params, b1apl_params, b1aps_params, 
                          b1fs_params, b1pa_params, b1qp_params)
b1_species_names <- c("Acer campestre", "Acer platanoides", "Acer pseudoplatanus",
                      "Fagus sylvatica", "Prunus avium", "Quercus petraea")

# Site B2
cat("\n=== Processing Site B2 ===\n")
b2_species_params <- list(b2ac_params, b2pa_params, b2qp_params)
b2_species_names <- c("Acer campestre", "Prunus avium", "Quercus petraea")

# 8. CREATE FACETED PLOTS (ALL RCPs AS PANELS IN ONE GRAPH PER SITE)
cat("\n=== Creating Faceted Plots ===\n")

# Create faceted plot for B1 (RCPs as panels, species as colored lines in each panel)
b1_faceted_plot <- plot_site_all_rcp_faceted(
  "Site B1 (Sailershausen)",
  cordex_b1,
  b1_species_params,
  b1_species_names,
  rcp_scenarios
)

# Create faceted plot for B2 (RCPs as panels, species as colored lines in each panel)
b2_faceted_plot <- plot_site_all_rcp_faceted(
  "Site B2 (Sailershausen)",
  cordex_b2,
  b2_species_params,
  b2_species_names,
  rcp_scenarios
)

# 9. DISPLAY FACETED PLOTS
cat("\n=== Displaying Faceted Plots ===\n")

# Display the two main plots
print(b1_faceted_plot)
print(b2_faceted_plot)

# 10. OPTIONAL: SAVE FACETED PLOTS TO FILES
save_faceted_plots <- function(output_dir = "plots") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save B1 faceted plot
  ggsave(
    file.path(output_dir, "B1_All_RCP_Faceted.png"),
    b1_faceted_plot,
    width = 16,
    height = 8,
    dpi = 300
  )
  cat("Saved: B1_All_RCP_Faceted.png\n")
  
  # Save B2 faceted plot
  ggsave(
    file.path(output_dir, "B2_All_RCP_Faceted.png"),
    b2_faceted_plot,
    width = 16,
    height = 8,
    dpi = 300
  )
  cat("Saved: B2_All_RCP_Faceted.png\n")
}

# Uncomment to save the faceted plots
# save_faceted_plots()

cat("\n=== Analysis Complete ===\n")
cat("You now have 2 faceted plots:\n")
cat("  1. Site B1 - 3 RCP panels, each showing all 6 species\n")
cat("  2. Site B2 - 3 RCP panels, each showing all 3 species\n")
cat("\nEach RCP panel contains colored lines for each species,\n")
cat("making it easy to compare species within and across scenarios!\n")

