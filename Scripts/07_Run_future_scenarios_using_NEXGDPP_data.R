#####################
### Author: Angela Fanelli
### Role: Epidemiologist
### Institution: JRC
### Description:
### This code runs the best-performing model while keeping environmental variables constant.
#####################

# ---------------------------
# Source required scripts
# ---------------------------
source("Scripts/00_Load_packages.R")
source("Scripts/01_Define_custom_functions.R")

# ---------------------------
# Define output path
# ---------------------------
output_path <- "Output"

# ---------------------------
# Create necessary output directories
# ---------------------------
dirs_to_create <- c(
  glue("{output_path}/NEX-GDDP/Future_scenarios/Models"),
  glue("{output_path}/NEX-GDDP/Future_scenarios/Predictions")
)

walk(dirs_to_create, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# ---------------------------
# Load data
# ---------------------------
data <- fread("Input/Data_model/data_2010_2023.csv", header = TRUE)

spei3_h <- read.csv("Input/Historical_climate/spei_3_v2.10_h.csv")
t2m_h <- read.csv("Input/Historical_climate/t2m_h.csv")

g <- inla.read.graph(filename = "Input/Data_model/nuts3_graph")

best_model <- readRDS(glue("{output_path}/Models/model_smooth_15.RDS"))

data_model <- fread("Input/Data_model/data_model.csv")

# ---------------------------
# Extract random effects for constant variables
# ---------------------------
fhn_g_re <- best_model$summary.random$fhn_g |> 
  select(ID, mean) |> 
  rename(fhn_g = ID)

mammal_richness_g_re <- best_model$summary.random$mammal_richness_g |> 
  select(ID, mean) |> 
  rename(mammal_richness_g = ID)

constant_variables <- data_model |>
  select(nuts_id, year, fhn_g, mammal_richness_g) |>
  left_join(fhn_g_re, by = "fhn_g") |>
  select(!fhn_g) |> rename(fhn_g = mean) |>
  left_join(mammal_richness_g_re, by = "mammal_richness_g") |>
  select(!mammal_richness_g) |> rename(mammal_richness_g = mean)

constant_variables <- constant_variables |>
  group_by(nuts_id) |>
  summarise(
    fhn_g = mean(fhn_g),
    mammal_richness_g = mean(mammal_richness_g)
  ) |>
  ungroup()

# ---------------------------
# Calculate lagged climate covariates for historical period
# ---------------------------
t2m_h <- t2m_h |>
  arrange(id, year, month) |>
  group_by(id) |>
  mutate(t2m_lag = lag(t2m, 1)) |>
  ungroup() |> select(!t2m) |> filter(year >= 2004 & year <= 2009)

spei3_h <- spei3_h |>
  arrange(id, year, month) |>
  group_by(id) |>
  mutate(spei_lag = lag(spei, 1)) |>
  ungroup() |> select(!spei) |> filter(year >= 2004 & year <= 2009)

# Only NUTS included in reference period
t2m_h <- t2m_h |> filter(id %in% unique(data$nuts_id))
spei3_h <- spei3_h |> filter(id %in% unique(data$nuts_id))

# ---------------------------
# Extend dataset with future and historical periods
# ---------------------------
meta <- data |> distinct(country_code, nuts_index, nuts_id, nuts_name, disease)

extra_1 <- make_time_grid(2004, 2009)
extra_2 <- make_time_grid(2041, 2060)
extra_3 <- make_time_grid(2081, 2100)

padding <- rbind(
  make_padding(extra_1, meta, by = NULL),
  make_padding(extra_2, meta, by = NULL),
  make_padding(extra_3, meta, by = NULL)
)

data_longer <- rbind(data, padding) |> arrange(nuts_index, year, month)

data_longer <- data_longer |>
  mutate(period = case_when(
    year < 2010 ~ "historical",
    year >= 2010 & year <= 2023 ~ "reference_period",
    year >= 2041 & year <= 2060 ~ "future_short_term",
    year >= 2081 & year <= 2100 ~ "future_long_term",
    TRUE ~ NA_character_
  ))

# Fix population at mean over reference period
mean_pop <- data |>
  mutate(population = population / 10^5) |>
  group_by(nuts_id) |>
  summarise(population = mean(population)) |> ungroup()

# ---------------------------
# Create new data model
# ---------------------------
newdata_model <- data_longer |>
  rename(Y = cases) |>
  select(country_code, nuts_index, nuts_id, nuts_name, disease, month, year, Y, t2m_lag, spei_lag, period) |>
  left_join(mean_pop, by = "nuts_id") |>
  left_join(constant_variables)

# Substitute NAs with historical values
newdata_model <- newdata_model |>
  left_join(spei3_h, by = c("nuts_id" = "id", "year" = "year", "month" = "month"), suffix = c("", "_new")) |>
  left_join(t2m_h,  by = c("nuts_id" = "id", "year" = "year", "month" = "month"), suffix = c("", "_new")) |>
  mutate(
    spei_lag = coalesce(spei_lag, spei_lag_new),
    t2m_lag  = coalesce(t2m_lag,  t2m_lag_new)
  ) |> select(-spei_lag_new, -t2m_lag_new)

# ---------------------------
# Define scenarios and models
# ---------------------------
nexgdpp_models <- c("FGOALS-g3","IPSL-CM6A-LR","CanESM5","NorESM2-MM","CMCC-ESM2")
scenarios <- c("ssp245","ssp585")
timeframe <- "future"

# ---------------------------
# Loop over scenarios and models
# ---------------------------
walk(scenarios, function(scenario) {
  walk(nexgdpp_models, function(model_name) {
    
    # Import NEX-GDPP models
    spei3_f <- read.csv(glue("Input/NEX-GDPP/SPEI3/{model_name}_{scenario}_{timeframe}.csv"))
    t2m_f <- read.csv(glue("Input/NEX-GDPP/TMEAN/{model_name}_{scenario}_{timeframe}.csv"))
    
    l <- glue("{model_name}_{scenario}_{timeframe}")
    
    # Create lagged covariates
    t2m_f <- t2m_f |>
      arrange(id, year, month) |>
      group_by(id) |>
      mutate(t2m_lag = lag(t2m,1)) |>
      ungroup() |> select(!t2m) |> filter(!year < 2041)
    
    spei3_f <- spei3_f |>
      arrange(id, year, month) |>
      group_by(id) |>
      mutate(spei_lag = lag(spei,1)) |>
      ungroup() |> select(!spei) |> filter(!year < 2041)
    
    # Substitute NAs with future scenario values
    newdata_model <- newdata_model |>
      left_join(spei3_f, by = c("nuts_id"="id","year"="year","month"="month"), suffix = c("", "_new")) |>
      left_join(t2m_f,  by = c("nuts_id"="id","year"="year","month"="month"), suffix = c("", "_new")) |>
      mutate(
        spei_lag = coalesce(spei_lag, spei_lag_new),
        t2m_lag  = coalesce(t2m_lag,  t2m_lag_new)
      ) |> select(-spei_lag_new, -t2m_lag_new)
    
    # Reclassify climate variables
    newdata_model <- newdata_model |>
      mutate(
        t2m_lag_g = inla.group(newdata_model$t2m_lag, method = "cut", n = 20),
        spei_lag_g = inla.group(newdata_model$spei_lag, method = "cut", n = 20)
      ) |> select(-c(t2m_lag, spei_lag))
    
    # ---------------------------
    # Define precision prior and formula
    # ---------------------------
    alpha <- 0.01
    u <- 1
    precision_prior <- list(prec = list(prior = "pc.prec", param = c(u, alpha)))
    
    formula <- Y ~ offset(fhn_g) + offset(mammal_richness_g) +
      f(month, model = "rw1", cyclic = TRUE, constr = TRUE, scale.model = TRUE, hyper = precision_prior) +
      f(year, model = "iid", hyper = precision_prior) +
      f(nuts_index, model = "bym2", graph = g, scale.model = TRUE, hyper = precision_prior) +
      f(t2m_lag_g, model = "rw2", hyper = precision_prior, scale.model = TRUE, constr = TRUE) +
      f(spei_lag_g, model = "rw2", hyper = precision_prior, scale.model = TRUE, constr = TRUE)
    
    # ---------------------------
    # Run the model
    # ---------------------------
    temp_dir <- paste0("/scratch/panelan/INLA/inla_tmp_", format(Sys.time(),"%Y%m%d_%H%M%S"))
    dir.create(temp_dir, recursive = TRUE)
    
    start_msg <- glue("{l}: processing started at {Sys.time()}\n")
    print(start_msg)
    write(start_msg, glue('{output_path}/NEX-GDDP/Future_scenarios/Models/log_file.txt'), append = TRUE)
    
    inla.setOption(num.threads = 7)
    
    model <- fit_INLA(
      formula = formula, data = newdata_model, family = "nbinomial", config = FALSE,
      temp_dir = temp_dir, return_marginals = TRUE, verbose = TRUE
    )
    
    end_msg <- glue("{l} processing ended at {Sys.time()}\n")
    print(end_msg)
    write(end_msg, glue("{output_path}/NEX-GDDP/Future_scenarios/Models/log_file.txt"), append = TRUE)
    
    unlink(temp_dir, recursive = TRUE)
    
    # Save model and predictions
    saveRDS(model, glue("{output_path}/NEX-GDDP/Future_scenarios/Models/model_{l}.RDS"))
    
    pred <- newdata_model |>
      mutate(
        mean = model$summary.fitted.values$mean,
        sd = model$summary.fitted.values$sd,
        upper = model$summary.fitted.values$`0.975quant`,
        lower = model$summary.fitted.values$`0.025quant`
      )
    
    fwrite(pred, glue("{output_path}/NEX-GDDP/Future_scenarios/Predictions/pred_{l}.csv"), row.names = FALSE)
    
  }, .progress = TRUE)
})
