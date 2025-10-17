#####################
### Author: Angela Fanelli
### Role: Epidemiologist
### Institution: JRC
#####################

# ---------------------------
# Load required scripts
# ---------------------------
source("Scripts/00_Load_packages.R")
source("Scripts/01_Define_custom_functions.R")

# ---------------------------
# Load adjacency matrix
# ---------------------------
g <- inla.read.graph(filename = "Input/Data_model/nuts3_graph")

# ---------------------------
# Load input data
# ---------------------------
data <- fread("Input/Data_model/data_2010_2023.csv", header = TRUE)

# ---------------------------
# Define paths and models
# ---------------------------
output_path <- "Output"
nexgdpp_models <- c("FGOALS-g3","IPSL-CM6A-LR","CanESM5","NorESM2-MM","CMCC-ESM2")
scenario <- "reference"
timeframe <- "200912-202312"

# Create necessary output directories
dirs_to_create <- c(
  glue("{output_path}/NEX-GDDP/Refererence/Models"),
  glue("{output_path}/NEX-GDDP/Refererence/Predictions")
)
walk(dirs_to_create, ~ if(!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# ---------------------------
# Loop over climate models
# ---------------------------
walk(nexgdpp_models, function(model_name) {
  
  # ---------------------------
  # Load climate covariates
  # ---------------------------
  spei3 <- read.csv(glue("Input/NEX-GDPP/SPEI3/{model_name}_{scenario}_{timeframe}.csv"))
  t2m <- read.csv(glue("Input/NEX-GDPP/TMEAN/{model_name}_{scenario}_{timeframe}.csv"))
  
  label <- glue("{model_name}_{scenario}_{timeframe}")
  
  # ---------------------------
  # Create lagged covariates
  # ---------------------------
  t2m <- t2m |>
    arrange(id, year, month) |>
    group_by(id) |>
    mutate(t2m_lag = lag(t2m, 1)) |>
    ungroup() |>
    select(!t2m)
  
  spei3 <- spei3 |>
    arrange(id, year, month) |>
    group_by(id) |>
    mutate(spei_lag = lag(spei, 1)) |>
    ungroup() |>
    select(!spei)
  
  # ---------------------------
  # Merge covariates with main data
  # ---------------------------
  data_model <- data |>
    select(!c(t2m_lag, spei_lag)) |>
    left_join(t2m, by = c("nuts_id" = "id", "month", "year")) |>
    left_join(spei3, by = c("nuts_id" = "id", "month", "year")) |>
    
  # ---------------------------
  # Prepare data for INLA
  # ---------------------------
  mutate(
    t2m_lag_g = inla.group(t2m_lag, method = "cut", n = 20),
    spei_lag_g = inla.group(spei_lag, method = "cut", n = 20),
    fhn_g = inla.group(fhn, method = "cut", n = 20),
    mammal_richness_g = inla.group(mammal_richness, method = "cut", n = 20),
    population = population / 1e5  # rescale population
  ) |>
    rename(Y = cases) |>
    select(
      country_code, nuts_index, nuts_id, nuts_name,
      disease, month, year, Y,
      t2m_lag_g, spei_lag_g, fhn_g, mammal_richness_g,
      population
    )
  
  # ---------------------------
  # Define INLA precision prior
  # ---------------------------
  alpha <- 0.01  # 1% chance that sigma > u
  u <- 1
  precision_prior <- list(prec = list(prior = "pc.prec", param = c(u, alpha)))
  
  # ---------------------------
  # INLA formula
  # ---------------------------
  formula <- Y ~ 
    f(month, model = "rw1", cyclic = TRUE, constr = TRUE, scale.model = TRUE, hyper = precision_prior) +
    f(year, model = "iid", hyper = precision_prior) +
    f(nuts_index, model = "bym2", graph = g, scale.model = TRUE, hyper = precision_prior) +
    f(t2m_lag_g, model = "rw2", hyper = precision_prior, scale.model = TRUE, constr = TRUE) +
    f(spei_lag_g, model = "rw2", hyper = precision_prior, scale.model = TRUE, constr = TRUE) +
    f(fhn_g, model = "rw2", hyper = precision_prior, scale.model = TRUE, constr = TRUE) +
    f(mammal_richness_g, model = "rw2", hyper = precision_prior, scale.model = TRUE, constr = TRUE)
  
  # ---------------------------
  # Temporary directory for INLA
  # ---------------------------
  temp_dir <- paste0("/scratch/panelan/INLA/inla_tmp_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  dir.create(temp_dir, recursive = TRUE)
  
  # ---------------------------
  # Start modeling
  # ---------------------------
  start_msg <- glue("{label}: processing started at {Sys.time()}")
  print(start_msg)
  write(start_msg, glue('{output_path}/NEX-GDDP/Refererence/Models/log_file.txt'), append = TRUE)
  
  model <- fit_INLA(
    formula = formula,
    data = data_model,
    family = "nbinomial",
    config = FALSE,
    temp_dir = temp_dir,
    return_marginals = TRUE,
    verbose = TRUE
  )
  
  end_msg <- glue("{label}: processing ended at {Sys.time()}")
  print(end_msg)
  write(end_msg, glue("{output_path}/NEX-GDDP/Refererence/Models/log_file.txt"), append = TRUE)
  
  # ---------------------------
  # Clean up temporary directory
  # ---------------------------
  unlink(temp_dir, recursive = TRUE)
  
  # ---------------------------
  # Save model and predictions
  # ---------------------------
  saveRDS(model, glue("{output_path}/NEX-GDDP/Refererence/Models/model_{label}.RDS"))
  
  pred <- data_model |>
    mutate(
      mean = model$summary.fitted.values$mean,
      sd = model$summary.fitted.values$sd,
      upper = model$summary.fitted.values$`0.975quant`,
      lower = model$summary.fitted.values$`0.025quant`
    )
  
  fwrite(pred, glue("{output_path}/NEX-GDDP/Refererence/Predictions/pred_{label}.csv"), row.names = FALSE)
  
}, .progress = TRUE)

