estimate_parameters <- function(chrono, historic_climate) {
  input_historic <- make_vsinput_historic(chrono, historic_climate)
  
  vs_params_estimate <- vs_params(input_historic$trw,
                                  input_historic$tmean,
                                  input_historic$prec,
                                  input_historic$syear,
                                  input_historic$eyear,
                                  .phi = 50, iter = 20)
  vs_params_estimate
}

forward_historic <- function(vs_params_estimate, historic_climate) {
  input_transient_historic <- make_vsinput_transient(historic_climate)
  vs_run_forward(vs_params_estimate,
                 input_transient_historic$tmean,
                 input_transient_historic$prec,
                 input_transient_historic$syear,
                 input_transient_historic$eyear,
                 .phi = 50)
}

forward_projection <- function(vs_params_estimate, projected_climate) {
  
  make_forward_model <- function(x) {
    input <- make_vsinput_transient(x)
    vs_run_forward(vs_params_estimate,
                   input$temp,
                   input$prec,
                   input$syear,
                   input$eyear,
                   .phi = 50)
  }
  
  projected_climate %>% 
    filter(year > 2020) %>% 
    group_by(rcp, gcm, rcm) %>% 
    filter(length(year) == 80 * 12) %>% 
    nest() %>% 
    mutate(vs_forward = purrr::map(data, make_forward_model)) %>% 
    unnest(vs_forward)
}

smooth_trw <- function(x) {
  trw <- predict(loess(trw ~ year, data = x))
  data.frame(
    year = x$year,
    trw = trw
  )
}

forward_plot <- function(forward_model) {
  forward_model %>% 
    group_by(rcp, gcm, rcm) %>% 
    nest() %>% 
    mutate(trw = map(data, smooth_trw)) %>%
    unnest(trw) %>% 
    group_by(rcp, year) %>% 
    summarise(
      trw_sd = sd(trw),
      trw = mean(trw)
    ) %>% 
    ggplot(aes(year, trw)) +
    geom_ribbon(aes(ymin = trw - trw_sd, ymax = trw + trw_sd)) +
    geom_line(colour = "red") +
    facet_wrap(~ rcp)
}