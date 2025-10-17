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

# ---------------------------
# Define output path
# ---------------------------
output_path <- "Output"

# ---------------------------
# Create necessary output directories
# ---------------------------
dirs_to_create <- c(
  # glue("{output_path}/NEX-GDDP/Future_scenarios/Predictions"),
  glue("{output_path}/NEX-GDDP/Future_scenarios/Linear_models")
)

walk(dirs_to_create, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# ---------------------------
# Define the files of interest
# ---------------------------

files <- list.files(
  glue("{output_path}/NEX-GDDP/Future_scenarios/Predictions"),
  full.names = TRUE,
  pattern = ".csv"
)

# Extract scenario (ssp245 or ssp585) using regex
file_names <- sub(".*_(ssp[0-9]+)_.*", "\\1", files)

# Split into groups by scenario
# Process each SSP, return a list with two elements: $ssp245 and $ssp585
file_groups <- split(files, file_names)
file_groups

# ---------------------------
# Calculate prediction for each GCM
# ---------------------------

iwalk(file_groups, \(group, scenario) {
  message("Processing scenario: ", scenario)
  
  results <- process_group(group) |>
    group_by(period, model, month, nuts_id) |>
    summarise(
      n = n(),
      monthly_incidence = mean(mean, na.rm = TRUE),
      monthly_sd = sd(mean, na.rm = TRUE),
      monthly_cv = monthly_sd / monthly_incidence,
      monthly_se = monthly_sd / sqrt(n),
      monthly_upper = monthly_incidence + (1.96 * monthly_se),
      monthly_lower = monthly_incidence - (1.96 * monthly_se),
      .groups = "drop"
    ) |>
    as.data.frame()
  
  # Save to scenario-specific file
  saveRDS(
    results,
    glue("{output_path}/NEX-GDDP/Future_scenarios/Predictions/monthly_incidence_by_model_{scenario}.RDS")
  )
})

# ---------------------------
# Calculate average across models and run linear models
# ---------------------------

iwalk(file_groups, function(group, scenario) {
  
  message("Processing scenario: ", scenario)
  
  walk(c("future_short_term", "future_long_term"), function(timeframe) {
    
    # Use process_group function to read & combine all models for this SSP
    data <- process_group(group) |>
      filter(period %in% c(timeframe, "historical", "reference_period")) |>
      mutate(period_group = ifelse(period %in% c("historical", "reference_period"), 1, 2))
    
    # Compute the monthly mean incidence across models for all years
    pred_by_period_all_years <- data |>
      group_by(period_group, year, month, nuts_id) |>
      summarise(
        n = n(),
        monthly_sd = sd(mean, na.rm = TRUE),
        monthly_incidence = mean(mean, na.rm = TRUE),
        monthly_cv = monthly_sd / monthly_incidence,
        monthly_se = monthly_sd / sqrt(n),
        monthly_upper = monthly_incidence + (1.96 * monthly_se),
        monthly_lower = monthly_incidence - (1.96 * monthly_se),
        .groups = "drop"
      )
    
    # Save the monthly data for all years
    saveRDS(
      pred_by_period_all_years,
      glue("{output_path}/NEX-GDDP/Future_scenarios/Predictions/monthly_incidence_by_period_{scenario}_{timeframe}_all_years.RDS")
    )
    
    # Compute monthly mean incidence across all years
    pred_by_period_across_years <- pred_by_period_all_years |>
      select(period_group, month, nuts_id, monthly_incidence) |>
      group_by(period_group, month, nuts_id) |>
      summarise(
        n = n(),
        monthly_sd = sd(monthly_incidence, na.rm = TRUE),
        monthly_incidence = mean(monthly_incidence, na.rm = TRUE),
        monthly_cv = monthly_sd / monthly_incidence,
        monthly_se = monthly_sd / sqrt(n),
        monthly_upper = monthly_incidence + (1.96 * monthly_se),
        monthly_lower = monthly_incidence - (1.96 * monthly_se),
        .groups = "drop"
      )
    
    # Save the monthly data by period
    saveRDS(
      pred_by_period_across_years,
      glue("{output_path}/NEX-GDDP/Future_scenarios/Predictions/monthly_incidence_by_period_{scenario}_{timeframe}_across_years.RDS")
    )
    
    # Calculate % difference between future and historical/reference
    perc_diff <- pred_by_period_across_years |>
      select(period_group, month, nuts_id, monthly_incidence) |>
      pivot_wider(
        names_from = period_group,
        values_from = monthly_incidence,
        names_prefix = "pg"
      ) |>
      mutate(perc_diff = 100 * (pg2 - pg1) / pg1)
    
    # Linear model comparing period groups
    monthly_lm <- pred_by_period_all_years |>
      group_by(nuts_id, month) |>
      summarise(
        slope = coef(summary(lm(monthly_incidence ~ period_group)))[2, 1],
        se = coef(summary(lm(monthly_incidence ~ period_group)))[2, 2],
        pvalue = coef(summary(lm(monthly_incidence ~ period_group)))[2, 4],
        .groups = "drop"
      ) |>
      mutate(comparison = glue("historical_reference_vs_{timeframe}"))
    
    monthly_lm <- monthly_lm |>
      left_join(perc_diff, by = c("nuts_id", "month"))
    
    # Summarize results
    monthly_lm_summary <- monthly_lm |>
      group_by(comparison) |>
      summarise(
        n_all = n(),
        n_significant = sum(pvalue < 0.05),
        n_slope_pos_significant = sum(pvalue < 0.05 & slope > 0),
        n_slope_neg_significant = sum(pvalue < 0.05 & slope < 0),
        .groups = "drop"
      ) |>
      mutate(
        prop_pos_significant_out_of_all = round(n_slope_pos_significant / n_all, 2),
        prop_neg_significant_out_of_all = round(n_slope_neg_significant / n_all, 2),
        prop_pos_significant_out_of_significant = round(n_slope_pos_significant / n_significant, 2),
        prop_neg_significant_out_of_significant = round(n_slope_neg_significant / n_significant, 2)
      )
    
    # Return a list for this scenario and timeframe
    results <- list(
      monthly_lm = monthly_lm,
      monthly_lm_summary = monthly_lm_summary
    )
    
    results |>
      saveRDS(
        glue("{output_path}/NEX-GDDP/Future_scenarios/Linear_models/lm_results_{scenario}_{timeframe}.RDS")
      )
    
  }, .progress = TRUE)
})

# ---------------------------
# Create a unique df to describe results
# ---------------------------

lm_future_scenarios <- list.files(
  glue("{output_path}/NEX-GDDP/Future_scenarios/Linear_models"),
  full.names = TRUE,
  pattern = "\\.RDS$"
)

lm_summary <- map(lm_future_scenarios, function(x) {
  filename <- basename(x)
  parts <- str_match(filename, "^lm_results_(ssp[0-9]+)_(.*)\\.RDS$")
  ssp <- parts[, 2]
  incidence_change <- readRDS(x)
  
  df <- incidence_change$monthly_lm_summary |>
    mutate(scenario = ssp)
}) |> list_rbind()

# Save summary
lm_summary |>
  write.csv(
    glue("{output_path}/NEX-GDDP/Future_scenarios/Linear_models/lm_summary_all.csv"),
    row.names = FALSE
  )

# ---------------------------
# Extract estimates for discussion
# ---------------------------

lm_future_scenarios <- list.files(
  glue("{output_path}/NEX-GDDP/Future_scenarios/Linear_models"),
  full.names = TRUE,
  pattern = "\\.RDS$"
)

future_estimates <- map(lm_future_scenarios, function(x) {
  filename <- basename(x)
  parts <- str_match(filename, "^lm_results_(ssp[0-9]+)_(.*)\\.RDS$")
  ssp <- parts[, 2]
  incidence_change <- readRDS(x)
  
  df <- incidence_change$monthly_lm |>
    mutate(scenario = ssp)
}) |> list_rbind()

head(future_estimates)

future_estimates_wide <- future_estimates |>
  mutate(period = gsub("^historical_reference_vs_", "", comparison)) |>
  select(-c(comparison, pg1, pg2, slope, se, pvalue)) |>
  pivot_wider(
    names_from = period,
    values_from = perc_diff
  )

head(future_estimates_wide)

# ---------------------------
# Compare within-scenario differences (long vs short term)
# ---------------------------

future_estimates_wide |>
  group_by(scenario) |>
  summarise(
    n_total = n(),
    n_long_gt_short = sum(future_long_term > future_short_term, na.rm = TRUE),
    prop_long_gt_short = n_long_gt_short / n_total
  )

# ---------------------------
# Compare scenarios by period (ssp585 > ssp245)
# ---------------------------

future_estimates_wide |>
  select(nuts_id, month, scenario, future_long_term, future_short_term) |>
  pivot_longer(
    cols = starts_with("future_"),
    names_to = "period",
    values_to = "value"
  ) |>
  pivot_wider(
    names_from = scenario,
    values_from = value
  ) |>
  summarise(
    n_total = n(),
    n_ssp585_gt_ssp245 = sum(ssp585 > ssp245, na.rm = TRUE),
    prop_ssp585_gt_ssp245 = n_ssp585_gt_ssp245 / n_total,
    .by = period
  )

# ---------------------------
# Count significant changes by month and scenario
# ---------------------------

lm_future_scenarios <- list.files(
  glue("{output_path}/NEX-GDDP/Future_scenarios/Linear_models"),
  full.names = TRUE,
  pattern = "\\.RDS$"
)

nuts_change_significant_summary_grouped_by_month <- map(lm_future_scenarios, function(x) {
  incidence_change <- readRDS(x)
  
  filename <- basename(x)
  parts <- str_match(filename, "^lm_results_(ssp[0-9]+)_(.*)\\.RDS$")
  ssp <- parts[, 2]
  timeframe <- parts[, 3]
  
  change_significant_grouped <- incidence_change$monthly_lm |>
    mutate(
      ssp = ssp,
      period = timeframe
    ) |>
    group_by(ssp, period, month) |>
    summarise(
      total_nuts = n(),
      sig_pos = sum(pvalue < 0.05 & slope > 0, na.rm = TRUE),
      sig_neg = sum(pvalue < 0.05 & slope < 0, na.rm = TRUE),
      prop_pos = sig_pos / total_nuts,
      prop_neg = sig_neg / total_nuts,
      .groups = "drop"
    )
}, .progress = TRUE) |> list_rbind()

# Save it
nuts_change_significant_summary_grouped_by_month |>
  write.csv(
    glue("{output_path}/NEX-GDDP/Future_scenarios/Linear_models/nuts_change_significant_summary_grouped_by_month.csv"),
    row.names = FALSE
  )

