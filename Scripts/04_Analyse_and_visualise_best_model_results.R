#####################
### Author: Angela Fanelli
### Role: Epidemiologist
### Institution: JRC

# ---------------------------
# Load required scripts
# ---------------------------
source("Scripts/00_Load_packages.R")
source("Scripts/01_Define_custom_functions.R")
source("Scripts/02_Prepare_and_load_data_model.R")

# ---------------------------
# Define output paths and create directories if they do not exist
# ---------------------------
output_path <- "Output"

dirs_to_create <- c(
  output_path,
  glue("{output_path}/Figures")
)

walk(dirs_to_create, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# ---------------------------
# Load NUTS shapefiles
# ---------------------------
nuts3 <- gisco_get_nuts(nuts_level = 3, year = 2021, epsg = 3035)
nuts0 <- gisco_get_nuts(nuts_level = 0, year = 2021, epsg = 3035)

# ---------------------------
# Load model metrics
# ---------------------------
metrics_files <- list.files(glue("{output_path}/Metrics"), full.names = TRUE)
metrics_files <- metrics_files[!grepl("all.csv$", metrics_files)]

metrics <- metrics_files |>
  map(\(x) read.csv(x)) |>
  list_rbind()

head(metrics |> arrange(dic))

# Save combined metrics
metrics |>
  arrange(dic) |>
  write.csv(glue("{output_path}/Metrics/metrics_all.csv"), row.names = FALSE)

# ---------------------------
# Select and load the best model
# ---------------------------
best_model_name <- metrics |>
  arrange(dic) |>
  slice(1) |>
  pull(model)

best_model <- readRDS(glue("{output_path}/Models/{best_model_name}.RDS"))

# ---------------------------
# Extract random effects
# ---------------------------
re <- get_random_INLA(
  model = best_model,
  covars = names(best_model$summary.random)[
    !names(best_model$summary.random) %in% c("month","nuts_index","year")
  ],
  transform = TRUE
)

unique(re$covar)

# ---------------------------
# Plot smooth terms
# ---------------------------
re_plot <- re |>
  ggplot(aes(x = value, y = mean)) +
  geom_line(linewidth = 0.5, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, fill = "grey", size = 0.3) +
  geom_hline(yintercept = 1, lty = 2) +
  facet_wrap(~covar, scales = "free", labeller = as_labeller(c(
    "fhn_g" = "Forest human nexus",
    "spei_lag_g" = "3-month SPEI\n(1-month lag)",
    "t2m_lag_g" = "Tmean K (1-month lag)",
    "mammal_richness_g" = "Mammal richness"
  ))) +
  theme_matplotlib(base_size = 12) +
  labs(x = "", y = "Relative risk") +
  theme(
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "white", size = 0.5)
  )

print(re_plot)

ggsave(glue("{output_path}/Figures/RR_smooth_terms_free_scales.tiff"), 
       plot = last_plot(), width = 10, height = 8, dpi = 300)

# ---------------------------
# Seasonal effect
# ---------------------------
# Extract seasonality from best model
se_best <- get_random_INLA(
  model = best_model,
  covars = names(best_model$summary.random)[names(best_model$summary.random) == "month"],
  transform = TRUE
)

# Load alternative models
model_without_t2m <- readRDS(glue("{output_path}/Models/model_smooth_14.RDS"))
model_without_spei <- readRDS(glue("{output_path}/Models/model_smooth_13.RDS"))
model_only_seasonal <- readRDS(glue("{output_path}/Models/model_only_seasonal.RDS"))

# Extract seasonality effects
se_without_t2m <- get_random_INLA(model_without_t2m, covars = "month", transform = TRUE)
se_without_spei <- get_random_INLA(model_without_spei, covars = "month", transform = TRUE)
se_only_seasonal <- get_random_INLA(model_only_seasonal, covars = "month", transform = TRUE)

# Combine all seasonality effects
se_all <- se_best |> mutate(model = "Best model") |>
  rbind(se_without_t2m |> mutate(model = "Model without Tmean (1-month lag)")) |>
  rbind(se_without_spei |> mutate(model = "Model without 3-month SPEI (1-month lag)")) |>
  rbind(se_only_seasonal |> mutate(model = "Model with only seasonality"))

# Plot seasonal effect with facets
se_all |>
  ggplot(aes(x = value, y = mean)) +
  geom_line(size = 0.5, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "grey", size = 0.3) +
  scale_x_continuous(breaks = 1:12, labels = month.name) +
  facet_wrap(~model) +
  theme_matplotlib(base_size = 12) +
  labs(x = "Month", y = "Relative risk") +
  theme(
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "white", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(glue("{output_path}/Figures/seasonal_effect_facets.tiff"),
       plot = last_plot(), width = 12, height = 7.2, dpi = 300)

# Plot seasonal effect overlayed by model
se_all |>
  ggplot(aes(x = value, y = mean, color = model)) +
  geom_line(size = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.2, size = 0.3) +
  scale_x_continuous(breaks = 1:12, labels = month.name) +
  scale_color_manual(values = c(
    "Best model" = "#1b9e77",
    "Model without Tmean (1-month lag)" = "#d95f02",
    "Model without 3-month SPEI (1-month lag)" = "#7570b3",
    "Model with only seasonality" = "#e7298a"
  )) +
  scale_fill_manual(values = c(
    "Best model" = "#1b9e77",
    "Model without Tmean (1-month lag)" = "#d95f02",
    "Model without 3-month SPEI (1-month lag)" = "#7570b3",
    "Model with only seasonality" = "#e7298a"
  )) +
  theme_matplotlib(base_size = 12) +
  labs(x = "Month", y = "Relative risk", color = "", fill = "") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

ggsave(glue("{output_path}/Figures/seasonal_effect.tiff"),
       plot = last_plot(), width = 12, height = 7.2, dpi = 300)


# ---------------------------
# Inter-annual variability
# ---------------------------
model_only_inter_annual <- readRDS(glue("{output_path}/Models/model_only_inter_annual.RDS"))
model_without_fhn <- readRDS(glue("{output_path}/Models/model_smooth_12.RDS"))

# Extract random effects for inter-annual
inter_annual_best <- get_random_INLA(best_model, covars = "year", transform = TRUE)
inter_annual_without_fhn <- get_random_INLA(model_without_fhn, covars = "year", transform = TRUE)
inter_annual_without_t2m <- get_random_INLA(model_without_t2m, covars = "year", transform = TRUE)
inter_annual_without_spei <- get_random_INLA(model_without_spei, covars = "year", transform = TRUE)
only_inter_annual <- get_random_INLA(model_only_inter_annual, covars = "year", transform = TRUE)

# Combine inter-annual effects
inter_annual_all <- inter_annual_best |> mutate(model = "Best model") |>
  rbind(inter_annual_without_t2m |> mutate(model = "Model without Tmean (1-month lag)")) |>
  rbind(inter_annual_without_spei |> mutate(model = "Model without 3-month SPEI (1-month lag)")) |>
  rbind(inter_annual_without_fhn |> mutate(model = "Model without forest human nexus")) |>
  rbind(only_inter_annual |> mutate(model = "Model with only inter-annual effect"))

# Plot inter-annual effects with facets
inter_annual_all |>
  ggplot(aes(x = value, y = mean)) +
  geom_line(size = 0.5, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "grey", size = 0.3) +
  scale_x_continuous(breaks = 2010:2023) +
  facet_wrap(~model) +
  theme_matplotlib(base_size = 12) +
  labs(x = "Year", y = "Relative risk") +
  theme(
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "white", size = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(glue("{output_path}/Figures/inter_annual_effect_facets.tiff"),
       plot = last_plot(), width = 15, height = 7.2, dpi = 300)

# Plot inter-annual effects overlayed by model
inter_annual_all |>
  ggplot(aes(x = value, y = mean, color = model)) +
  geom_line(size = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.2, size = 0.3) +
  scale_x_continuous(breaks = 2010:2023) +
  scale_color_manual(values = c(
    "Best model" = "#1b9e77",
    "Model without Tmean (1-month lag)" = "#d95f02",
    "Model without 3-month SPEI (1-month lag)" = "#7570b3",
    "Model with only inter-annual effect" = "#e7298a",
    "Model without forest human nexus" = "#66a61e"
  )) +
  scale_fill_manual(values = c(
    "Best model" = "#1b9e77",
    "Model without Tmean (1-month lag)" = "#d95f02",
    "Model without 3-month SPEI (1-month lag)" = "#7570b3",
    "Model with only inter-annual effect" = "#e7298a",
    "Model without forest human nexus" = "#66a61e"
  )) +
  theme_matplotlib(base_size = 12) +
  labs(x = "Year", y = "Relative risk", color = "", fill = "") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

ggsave(glue("{output_path}/Figures/inter_annual_effect.tiff"),
       plot = last_plot(), width = 18, height = 7.2, dpi = 300)


# ---------------------------
# Leave-one-out cross-validation (LOOCV) & leave-group-out cross-validation (LGOCV)
# ---------------------------
alpha <- 0.01
u <- 1
precision_prior <- list(prec = list(prior = "pc.prec", param = c(u, alpha)))

loocv_res <- inla.group.cv(best_model)
cv <- loocv_res$cv[!is.na(loocv_res$cv)]
ULOOCV <- mean(log(cv))
exp(ULOOCV)  # Probability of predicting exactly what we observed

# Automatic leave-group-out CV
lgocv_auto_res <- inla.group.cv(result = best_model, num.level.sets = 3)
ULGOCV_auto <- mean(log(lgocv_auto_res$cv), na.rm = TRUE)
print(paste0("ULGOCV_auto is ", round(exp(ULGOCV_auto), 5), "."))

# ---------------------------
# Monthly risk mapping
# ---------------------------
pred <- fread(glue("{output_path}/Predictions/pred_{best_model_name}.csv"), header = TRUE)

# Compute mean by month across years for each NUTS3
monthly_incidence <- pred |>
  group_by(month, nuts_id) |>
  summarise(
    n = n(),
    monthly_sd = sd(mean),
    monthly_incidence = mean(mean, na.rm = TRUE),
    monthly_cv = monthly_sd / monthly_incidence,
    monthly_se = monthly_sd / sqrt(n),
    monthly_upper = monthly_incidence + (1.96 * monthly_se),
    monthly_lower = monthly_incidence - (1.96 * monthly_se)
  ) |>
  ungroup() |>
  left_join(nuts3, by = c("nuts_id" = "NUTS_ID")) |>
  st_as_sf()

# Define risk palette
risk_palette <- c("#f1eef6", "#d7b5d8", "#df65b0", "#dd1c77", "#980043")
extent <- st_bbox(st_crop(nuts0, c(xmin = 2377294, ymin = 1313597, xmax = 7453440, ymax = 8028510)))
osm_eur_carto <- basemap(extent, map_service = "carto", map_type = "light_no_labels", class = "terra")
osm_eur_carto <- osm_eur_carto |> project(crs(monthly_incidence)) |> crop(ext(monthly_incidence) * 1.02)

# Monthly incidence map
monthly_incidence_map <- tm_shape(osm_eur_carto) +
  tm_rgb(col_alpha = 0.8) +
  tm_shape(monthly_incidence) +
  tm_polygons(
    fill = "monthly_incidence",
    lwd = 0,
    fill.scale = tm_scale_intervals(n = 5, style = "fisher", interval.closure = "left", values = risk_palette),
    fill.legend = tm_legend(title = "", frame = FALSE, orientation = "landscape",
                            text.size = 1.2, ticks.col = "black",
                            position = tm_pos_out("center", "bottom", pos.h = "center"), show = TRUE,
                            width = 15)
  ) +
  tm_facets(by = "month", nrow = 3) +
  tm_shape(nuts0) +
  tm_borders("black", lwd = 0.2) +
  tm_layout(bg.color = "white", space.color = "white", asp = 0,
            outer.margins = c(0, 0, 0, 0), inner.margins = c(0, 0, 0, 0),
            panel.label.bg.color = "white",
            panel.labels = setNames(month.name, seq_along(month.name)),
            panel.label.size = 1.1,
            panel.label.height = 0.8)

tmap_save(monthly_incidence_map, filename = glue("{output_path}/Figures/estimated_monthly_incidence_map.tiff"),
          width = 4.5, height = 4.5, dpi = 300, asp = 0)

# ---------------------------
# SD as a measure of uncertainty
# ---------------------------
monthly_sd_map <- 
  tm_shape(osm_eur_carto) +
  tm_rgb(col_alpha = 0.8) +
  
  tm_shape(monthly_incidence) +
  tm_polygons(
    fill = "monthly_sd",
    lwd = 0,
    fill.scale = tm_scale_continuous(
      values = "scico.roma",
      limits = c(0, max(monthly_incidence$monthly_sd))
    ),
    fill.legend = tm_legend(
      title = "",
      frame = FALSE,
      orientation = "landscape",
      text.size = 0.7,
      ticks.col = "black",
      position = tm_pos_out("center", "bottom", pos.h = "center"),
      show = TRUE,
      width = 10
    )
  ) +
  
  tm_facets(by = "month", nrow = 3) +
  
  # tm_title(
  #   text = "Standard deviation of estimated incidence per 100,000 people", 
  #   position = tm_pos_out("center", "top", pos.h = "center"),
  #   frame = FALSE, padding = c(0, 0, 0, 0), size = 0.8
  # ) +
  
  tm_shape(nuts0) +
  tm_borders("black", lwd = 0.2) +
  
  tm_layout(
    bg.color = "white",
    space.color = "white",
    asp = 0,
    outer.margins = c(0, 0, 0, 0),
    inner.margins = c(0, 0, 0, 0),
    panel.label.bg.color = "white",
    panel.labels = setNames(month.name, seq_along(month.name)),
    panel.label.size = 1.1,
    panel.label.height = 0.8
  )

tmap_save(
  monthly_sd_map,
  filename = glue("{output_path}/Figures/estimated_monthly_sd_map.tiff"),
  width = 9 * 0.5,
  height = 9 * 0.5,
  dpi = 300,
  asp = 0
)


# ---------------------------
# Coefficient of variation (extent of prediction variability)
# ---------------------------
monthly_cv_map <- 
  tm_shape(osm_eur_carto) +
  tm_rgb(col_alpha = 0.8) +
  
  tm_shape(monthly_incidence) +
  tm_polygons(
    fill = "monthly_cv",
    lwd = 0,
    fill.scale = tm_scale_continuous(
      values = "scico.roma",
      # limits = c(0, round(max(monthly_incidence$monthly_cv), digits = 1))
      limits = c(0.10, 0.44)
    ),
    fill.legend = tm_legend(
      title = "",
      frame = FALSE,
      orientation = "landscape",
      text.size = 0.7,
      ticks.col = "black",
      position = tm_pos_out("center", "bottom", pos.h = "center"),
      show = TRUE,
      width = 10
    )
  ) +
  
  tm_facets(by = "month", nrow = 3) +
  
  # tm_title(
  #   text = "Coefficient of variation of estimated incidence per 100,000 people", 
  #   position = tm_pos_out("center", "top", pos.h = "center"),
  #   frame = FALSE, padding = c(0, 0, 0, 0), size = 0.8
  # ) +
  
  tm_shape(nuts0) +
  tm_borders("black", lwd = 0.2) +
  
  tm_layout(
    bg.color = "white",
    space.color = "white",
    asp = 0,
    outer.margins = c(0, 0, 0, 0),
    inner.margins = c(0, 0, 0, 0),
    panel.label.bg.color = "white",
    panel.labels = setNames(month.name, seq_along(month.name)),
    panel.label.size = 1.1,
    panel.label.height = 0.8
  )

tmap_save(
  monthly_cv_map,
  filename = glue("{output_path}/Figures/estimated_monthly_cv_map.tiff"),
  width = 9 * 0.5,
  height = 9 * 0.5,
  dpi = 300,
  asp = 0
)

# ---------------------------
# Convert incidence to risk probability
# ---------------------------
monthly_risk <- monthly_incidence |>
  mutate(monthly_risk = 1 - exp(-monthly_incidence)) |>
  relocate(monthly_risk, .after = monthly_incidence)

breaks_cat <- c(0, 0.15, 0.30, 0.45, 0.6, 1)
monthly_risk$risk_interval <- cut(monthly_risk$monthly_risk, breaks = breaks_cat, right = FALSE)
monthly_risk$risk_class <- cut(monthly_risk$monthly_risk, breaks = breaks_cat,
                               labels = c("Very low", "Low", "Medium", "High", "Very high"), right = FALSE)
monthly_risk$risk_interval <- case_when(
  monthly_risk$risk_interval == "[0.6,1)" ~ "[0.6,1]",
  TRUE ~ monthly_risk$risk_interval
)

# Stacked bar plot of NUTS3 by risk class
monthly_risk_nuts3_summary <- monthly_risk |>
  group_by(month, risk_class) |>
  summarise(n_nuts3 = n()) |>
  ungroup()

risk_order <- c("Very low", "Low", "Medium", "High", "Very high")
monthly_risk_nuts3_summary$risk_class <- factor(monthly_risk_nuts3_summary$risk_class, levels = risk_order)

ggplot(monthly_risk_nuts3_summary, aes(x = factor(month), y = n_nuts3, fill = risk_class)) +
  geom_col(position = position_stack(reverse = TRUE), color = "black", size = 0.1) +
  scale_fill_manual(values = risk_palette, breaks = risk_order) +
  scale_x_discrete(labels = month.name, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1500, by = 200), expand = c(0, 0)) +
  labs(x = "Month", y = "Number of NUTS 3", fill = "Risk Class") +
  theme_matplotlib(base_size = 12) +
  theme(axis.ticks.length = unit(0, "pt"),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white", color = "black", size = 0.5),
        panel.spacing = unit(30, "pt"),
        legend.position = "bottom",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))

ggsave(glue("{output_path}/Figures/n_nuts3_monthly_risk.tiff"),
       plot = last_plot(), width = 12, height = 8, dpi = 300)

# ---------------------------
# Approximate population at risk
# ---------------------------

# Calculate average population per NUTS3 region (scaled by 10^5)
pop_nuts3 <- data_model |> 
  group_by(nuts_id) |> 
  summarise(pop_mean = mean(population) * 10^5) |> 
  ungroup()

# Combine with monthly risk data and estimate total population per risk class and month
pop_risk <- monthly_risk |> 
  st_drop_geometry() |> 
  select(month, nuts_id, risk_class) |> 
  left_join(pop_nuts3, by = "nuts_id") |> 
  group_by(month, risk_class) |> 
  summarise(pop_sum = sum(pop_mean)) |> 
  ungroup()

# Add readable month names before plotting
pop_risk <- pop_risk |>
  mutate(month_name = factor(month.name[month], levels = month.name))

ggplot(pop_risk, aes(x = risk_class, y = pop_sum, fill = risk_class)) +
  geom_col() +
  facet_wrap(~ month_name, ncol = 3) +
  scale_fill_manual(values = risk_palette, breaks = risk_order) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Risk Class", y = "Population", fill = "Risk Class") +
  theme_matplotlib(base_size = 12) +
  theme(
    axis.ticks.length = unit(0, "pt"),
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "white", size = 0.5),
    panel.spacing = unit(30, "pt"),
    legend.position = "bottom",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13)
  )

ggsave(
  glue("{output_path}/Figures/pop_risk.tiff"),
  plot = last_plot(),
  width = 12,
  height = 8,
  dpi = 300
)

