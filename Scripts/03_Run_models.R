#####################
### Author: Angela Fanelli
### Role: Epidemiologist
### Institution: JRC
#####################

# ---------------------------
# Source required scripts
# ---------------------------
source("Scripts/00_Load_packages.R")
source("Scripts/01_Define_custom_functions.R")
source("Scripts/02_Prepare_and_load_data_model.R")

# ---------------------------
# Set output paths and create directories if they do not exist
# ---------------------------
output_path <- "Output"

dirs_to_create <- c(
  glue("{output_path}/Models"),
  glue("{output_path}/Predictions"),
  glue("{output_path}/Metrics")
)

walk(dirs_to_create, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# ---------------------------
# Define precision prior
# ---------------------------
alpha <- 0.01   # 1% chance that σ > U, weakly informative prior
u <- 1
precision_prior <- list(prec = list(prior = "pc.prec", param = c(u, alpha)))

# ---------------------------
# Define base formulas
# ---------------------------
base_formula <- Y ~ 1 +
  f(month, model = "rw1", cyclic = TRUE, constr = TRUE, scale.model = TRUE, hyper = precision_prior) +
  f(year, model = "iid", hyper = precision_prior) +
  f(nuts_index, model = "bym2", graph = g, scale.model = TRUE, hyper = precision_prior)

only_seasonal_formula <- Y ~ 1 +
  f(month, model = "rw1", cyclic = TRUE, constr = TRUE, scale.model = TRUE, hyper = precision_prior)

only_inter_annual_formula <- Y ~ 1 +
  f(year, model = "iid", hyper = precision_prior)

# ---------------------------
# Generate all combinations of covariates dynamically
# ---------------------------
covariates <- c("t2m_lag_g", "spei_lag_g", "fhn_g", "mammal_richness_g")
covar_combinations <- create_all_combinations(covariates)

# ---------------------------
# Create smooth term formulae for each combination
# ---------------------------
formulae <- map(covar_combinations, function(covars) {
  if (length(covars) == 0) return(NULL)
  
  smooth_terms <- map(covars, function(y) {
    paste0("f(", y, ", model = 'rw2', hyper = precision_prior, scale.model = TRUE, constr = TRUE)")
  }) |> paste(collapse = " + ")
  
  update.formula(base_formula, as.formula(paste("~ . +", smooth_terms)))
})

# Remove NULL (empty combinations) if any
formulae <- compact(formulae)

# Add simple models
formulae <- c(formulae, list(base_formula, only_inter_annual_formula, only_seasonal_formula))

# ---------------------------
# Automatically generate model labels
# ---------------------------
num_smooth_models <- length(formulae) - 3
labels <- c(
  paste0("model_smooth_", 1:num_smooth_models),
  "base_model",
  "model_only_inter_annual",
  "model_only_seasonal"
)

# ---------------------------
# Set INLA CPU threads
# ---------------------------
inla.setOption(num.threads = 70)

# ---------------------------
# Fit models iteratively
# ---------------------------
walk2(formulae, labels, function(f, l) {
  
  # Create temporary directory for INLA
  temp_dir <- paste0("/scratch/panelan/INLA/inla_tmp_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  dir.create(temp_dir, recursive = TRUE)
  
  # Log start time
  start_msg <- glue("{l}: processing started at {Sys.time()}\n")
  print(start_msg)
  write(start_msg, file = file.path(output_path, "Models/log_file.txt"), append = TRUE)
  
  # Fit INLA model
  model <- fit_INLA(
    formula = f, 
    data = data_model, 
    family = "nbinomial", 
    config = FALSE, 
    temp_dir = temp_dir, 
    return_marginals = TRUE,
    verbose = TRUE
  )
  
  # Log end time
  end_msg <- glue("{l} processing ended at {Sys.time()}\n")
  print(end_msg)
  write(end_msg, file = file.path(output_path, "Models/log_file.txt"), append = TRUE)
  
  # Delete temporary directory
  unlink(temp_dir, recursive = TRUE)
  
  # Save model
  saveRDS(model, file = file.path(output_path, "Models", paste0(l, ".RDS")))
  
  # Save predictions
  pred <- data_model |>
    mutate(
      mean = model$summary.fitted.values$mean,
      sd = model$summary.fitted.values$sd,
      upper = model$summary.fitted.values$`0.975quant`,
      lower = model$summary.fitted.values$`0.025quant`
    )
  fwrite(pred, file = file.path(output_path, "Predictions", paste0("pred_", l, ".csv")), row.names = FALSE)
  
  # Extract and save model metrics
  model_metrics <- get_metrics_INLA(model = model, label = l)
  write.csv(model_metrics, file = file.path(output_path, "Metrics", paste0("metrics_", l, ".csv")), row.names = FALSE)
  
}, .progress = TRUE)

