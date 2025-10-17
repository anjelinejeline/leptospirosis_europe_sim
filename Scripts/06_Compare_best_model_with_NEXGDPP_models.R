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

set.seed(123)  # Ensure reproducibility

# ---------------------------
# Define output path
# ---------------------------
output_path <- "Output"

# ---------------------------
# Create necessary output directories
# ---------------------------
#dirs_to_create <- c(
#  glue("{output_path}/Figures"))

#walk(dirs_to_create, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# ---------------------------
# Load prediction files
# ---------------------------
files <- list.files(
  path = glue("{output_path}/NEX-GDDP/Refererence/Predictions"),
  pattern = "\\.csv$",
  full.names = TRUE
)

# Extract model names from filenames
model_names <- str_match(basename(files), "pred_(.*?)_reference")[,2]

# Load each CSV into global environment as pred_{model_name}
walk2(files, model_names, function(file, model) {
  assign(glue("pred_{model}"), fread(file), envir = .GlobalEnv)
})

# ---------------------------
# Load observed predictions
# ---------------------------
pred_obs <- fread(glue("{output_path}/Predictions/pred_model_smooth_15.csv"), header = TRUE)

# ---------------------------
# Correlation with observed predictions
# ---------------------------
cor(pred_obs$mean, pred_CanESM5$mean)
cor(pred_obs$mean, `pred_CMCC-ESM2`$mean)
cor(pred_obs$mean, `pred_FGOALS-g3`$mean)
cor(pred_obs$mean, `pred_IPSL-CM6A-LR`$mean)
cor(pred_obs$mean, `pred_NorESM2-MM`$mean)

# ---------------------------
# Mean difference (model bias)
# ---------------------------
mean(pred_obs$mean - pred_CanESM5$mean)        
mean(pred_obs$mean - `pred_CMCC-ESM2`$mean)   
mean(pred_obs$mean - `pred_FGOALS-g3`$mean)   
mean(pred_obs$mean - `pred_IPSL-CM6A-LR`$mean)
mean(pred_obs$mean - `pred_NorESM2-MM`$mean)  

# ---------------------------
# Combine observed and model predictions
# ---------------------------
df_pred <- data.frame(
  observed = pred_obs$mean,
  CanESM5 = pred_CanESM5$mean,
  CMCC_ESM2 = `pred_CMCC-ESM2`$mean,
  FGOALS_g3 = `pred_FGOALS-g3`$mean,
  IPSL_CM6A_LR = `pred_IPSL-CM6A-LR`$mean,
  NorESM2_MM = `pred_NorESM2-MM`$mean
)

# ---------------------------
# Scatter plot: observed vs predicted (sample of 100)
# ---------------------------
scatter_plot<-df_pred |>
  slice_sample(n = 100) |>
  pivot_longer(
    cols = -observed,
    names_to = "model",
    values_to = "predicted"
  ) |>
  ggplot(aes(x = observed, y = predicted)) +
  geom_point(color = "navy", alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkred") +
  facet_wrap(~model, scales = "free_y") +
  labs(
    x = "Incidence (observed values)",
    y = "Incidence (NEX-GDPP)",
    title = "Observed vs NEX-GDPP (sample of 100)"
  ) +
  theme_matplotlib(base_size = 12) +
  theme(
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "white", size = 0.5)
  ) 

print(scatter_plot)

ggsave(
  glue("{output_path}/Figures/comparison_observed_nexgdpp.tiff"),
  plot = scatter_plot,
  width = 10, height = 8, dpi = 300
)

# ---------------------------
# Time series plot for sample NUTS3 regions
# ---------------------------
nuts_sample <- sample(unique(pred_obs$nuts_id), 4)

# Create a named list of datasets
datasets <- list(
  Observed = pred_obs,
  CanESM5 = pred_CanESM5,
  `CMCC-ESM2` = `pred_CMCC-ESM2`,
  `FGOALS-g3` = `pred_FGOALS-g3`,
  `IPSL-CM6A-LR` = `pred_IPSL-CM6A-LR`,
  `NorESM2-MM` = `pred_NorESM2-MM`
)

# Filter for sampled NUTS IDs and add model column
df_all <- datasets |>
  imap(function(df, name) {
    df |>
      filter(nuts_id %in% nuts_sample) |>
      mutate(model = name)
  }) |>
  bind_rows() |>
  mutate(date = as.Date(paste(year, month, "01", sep = "-")))

# Plot time series
timeseries_plot<-ggplot(df_all, aes(x = date, y = mean, color = model)) +
  geom_line(size = 1) +
  facet_wrap(~nuts_id, scales = "free_y") +
  scale_color_viridis_d(option = "viridis") +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  labs(
    x = "Time",
    y = "Incidence",
    title = "Time series of predicted incidence (sample of 4 NUTS3)"
  ) +
  theme_matplotlib(base_size = 12) +
  theme(
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "white", size = 0.5),
    legend.title = element_blank(),
    legend.position = "bottom"
  ) 

print(timeseries_plot)

ggsave(
  glue("{output_path}/Figures/comparison_observed_nexgdpp_timeseries.tiff"),
  plot = timeseries_plot,
  width = 10, height = 8, dpi = 300
)

