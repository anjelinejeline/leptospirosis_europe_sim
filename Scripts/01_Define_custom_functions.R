#####################
### Author: Angela Fanelli
### Role: Epidemiologist
### Institution: JRC
### Project: APES

# Generate all  combinations of covars  
create_all_combinations<-function(covars){
  
  # Generate all combinations
  combinations <- list()
  for (i in 1:length(covars)) {
    comb <- combn(covars, i)
    for (j in 1:ncol(comb)) {
      combinations[[length(combinations) + 1]] <- comb[, j]
    }
  }
  
  # Print the list of combinations
  return(combinations)
  
}

# inla model function: include formula and set defaults

fit_INLA <- function(formula, data, family="nbinomial", config = FALSE, temp_dir,return_marginals=FALSE, verbose=TRUE)
  
{
  model <- inla(
    # formula, data and model family 
    formula = formula, data = data, family = family, 
    
    # offset 
    offset=log(population),
    
    # configure inla approx
    control.inla = list(strategy = "adaptive"), 
    
    # items to compute
    control.compute = list(dic = TRUE, waic = TRUE, config = config, 
                           cpo = TRUE,return.marginals = return_marginals), #  do not return marginals unless specified(saves memory)
    
    # fixed effects calibration
    control.fixed = list(correlation.matrix = TRUE, 
                         prec.intercept = 1, # precision 1
                         prec = 1), # weakly regularising on fixed effects (sd of 1)
    
    # save predicted values on response scale
    control.predictor = list(link = 1, compute = TRUE), 
    
    # Verbose 
    verbose = verbose,
    
    # Working directory 
    working.directory = temp_dir)
  
  model <- inla.rerun(model)
  return(model)
}

# Extract metrics from the model 
get_metrics_INLA <- function(model, label) {
  
  # Create a data frame with the metrics
  model_summary <- data.frame(
    # Name of the model 
    model = label,
    # Fixed effects 
    fixed = paste(model$names.fixed, collapse = " + "),
    # Random effects
    random= paste(names(model$summary.random), collapse = " + "),
    # DIC: Deviance Information Criterion, a measure of model fit and complexity
    dic = model$dic$dic,
    # WAIC: Watanabe-Akaike information criterion, a measure of model fit and complexity
    waic = model$waic$waic,
    # p.eff: Effective number of parameters, a measure of model complexity
    p_eff = model$waic$p.eff,
    # logscore: Log score, a measure of model predictive performance
    # The CPO is the probability of y_i knowing all other values of y.
    # The higher it is, the best. It a tool to compare models. 
    logscore = ifelse(all(model$cpo$cpo == 0), NA, -mean(log(model$cpo$cpo), na.rm = TRUE)),
    # cpo_fail: Number of observations for which CPO calculation failed
    cpo_fail = sum(model$cpo$failure == 1 & !is.na(model$cpo$failure)),
    stringsAsFactors = FALSE
  )
  
  # Return the data frame with the metrics
  return(model_summary)
  
}


# Function to extract the random effects 
get_random_INLA <- function(model, covars, transform=FALSE){
  
  summary_random <- model$summary.random
  
  # Select covariates 
  summary_random <- summary_random[names(summary_random)%in% covars]
  
  # Retrieve the estimates 
  summary_random <- map2(summary_random, names(summary_random), function(x, y) {
    x <- x |>
      rename("value" = 1, "lower" = 4, "median" = 5, "upper" = 6)
    x$covar <- y
    x <- x |> relocate(covar,.before=value)
    return(x)
  }) |> list_rbind()
  
  if (transform==FALSE){
    
    return(summary_random)
    
  } else {
    
    summary_random_transformed <- summary_random |> 
      mutate(across(c(mean, sd, lower, median, upper, mode), exp))
    
    return(summary_random_transformed)
    
  }
  
}


# Create a custom theme similar to matplotlib's style
theme_matplotlib <- function(base_size = 14) {
  theme(
    # Set background color
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    
    # Set text and title colors and font
    axis.title = element_text(color = "black", size = base_size * 1.2,  face = "plain"),
    
    
    axis.text = element_text(color = "black", size = base_size,  face = "plain"),
    plot.title = element_text(color = "black", size = base_size * 1.4, hjust = 0.5,  face = "plain"),
    
    # Set axis line colors and styles
    axis.line = element_line(size = 0.5, color = "black", linetype = "solid"),
    
    # Set grid line colors and styles
    panel.grid.major = element_line(size = 0.25, color = "grey90", linetype = "solid"),
    panel.grid.minor = element_line(size = 0.25, color = "grey95", linetype = "solid"),
    
    # Remove legend background and border
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    
    # Set legend text color and font
    legend.text = element_text(color = "black", size = base_size * 0.8,  face = "plain"),
    
    # Set plot margins
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )
}

# Function to build a year-month grid as a data.frame
make_time_grid <- function(start_year, end_year) {
  expand.grid(
    year = start_year:end_year,
    month = 1:12
  ) |>
    arrange(year, month)|>
    as.data.frame()
}

# Function to create padding for a time grid
make_padding <- function(time_grid, meta, by = NULL) {
  merge(
    as.data.frame(meta),
    as.data.frame(time_grid),
    by = by
  ) |>
    transform(
      cases = NA_real_,
      t2m_lag = NA_real_,
      spei_lag = NA_real_,
      population = NA_real_,
      fhn = NA_real_,
      mammal_richness = NA_real_
    )
}


# function to process a group of files (all from the same SSP)
process_group <- function(group) {
  group |>
    map(\(f) {
      model <- str_match(basename(f), "pred_(.*)_ssp[0-9]+_future.csv")[, 2]
      
      df <- read.csv(f)
      df$model <- model
      df
    }) |>
    list_rbind()
}
