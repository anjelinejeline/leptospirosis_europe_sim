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
# Define output paths and create folders if not exist
# ---------------------------
output_path <- "Output"
# dirs_to_create <- c(
#   output_path,
#   glue("{output_path}/Figures")
# )

#walk(dirs_to_create, ~ if(!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# ---------------------------
# Load NUTS shapefiles
# ---------------------------
nuts3 <- gisco_get_nuts(nuts_level = 3, year = 2021, epsg = 3035)
nuts0 <- gisco_get_nuts(nuts_level = 0, year = 2021, epsg = 3035)

# ---------------------------
# List linear model results
# ---------------------------
lm_future_scenarios <- list.files(
  glue("{output_path}/NEX-GDDP/Future_scenarios/Linear_models"),
  full.names = TRUE, pattern = "\\.RDS$"
)

# ---------------------------
# Map significant % change
# ---------------------------
walk(lm_future_scenarios, function(x) {
  
  incidence_change <- readRDS(x)
  filename <- basename(x)
  parts <- str_match(filename, "^lm_results_(ssp[0-9]+)_(.*)\\.RDS$")
  ssp <- parts[,2]
  timeframe <- parts[,3]
  
  ssp_title <- ifelse(ssp == "ssp245", "SSP2-4.5", "SSP5-8.5")
  timeframe_title <- ifelse(timeframe == "future_short_term", "short-term future", "long-term future")
  
  incidence_change_significant <- incidence_change$monthly_lm |>
    filter(pvalue < 0.05) |>
    left_join(nuts3, by = c("nuts_id" = "NUTS_ID")) |>
    st_as_sf()
  
  extent <- st_bbox(st_crop(nuts0, c(xmin = 2377294, ymin = 1313597, xmax = 7453440, ymax = 8028510)))
  
  osm_eur_carto <- basemap(extent, map_service = "carto", map_type = "light_no_labels", class = "terra") |>
    project(crs(incidence_change_significant)) |>
    crop(ext(incidence_change_significant) * 1.02)
  
  incidence_change_significant_map <- tm_shape(osm_eur_carto) +
    tm_rgb(col_alpha = 0.8) +
    tm_shape(incidence_change_significant) +
    tm_polygons(
      fill = "perc_diff", lwd = 0,
      tm_scale_continuous(values = "seismic", midpoint = 0, limits = c(-42, 140), ticks = seq(-42, 140, by = 30)),
      fill.legend = tm_legend(title = "", frame = FALSE, orientation = "landscape",
                              text.size = 0.7, ticks.col = "black",
                              position = tm_pos_out("center", "bottom", pos.h = "center"),
                              show = TRUE, width = 13)
    ) +
    tm_facets(by = "month", nrow = 3) +
    tm_title(glue("Expected % change in incidence\nunder {ssp_title} ({timeframe_title})"),
             position = tm_pos_out("center", "top", pos.h = "center"),
             frame = FALSE, padding = c(0,0,0,0), size = 0.8) +
    tm_shape(nuts0) + tm_borders("black", lwd = 0.2) +
    tm_layout(bg.color = "white", space.color = "white", asp = 0,
              outer.margins = c(0,0,0,0), inner.margins = c(0,0,0,0),
              panel.label.bg.color = "white",
              panel.labels = setNames(month.name, seq_along(month.name)),
              panel.label.size = 1.1, panel.label.height = 0.8)
  
  tmap_save(incidence_change_significant_map,
            glue("{output_path}/Figures/incidence_change_significant_{ssp}_{timeframe}.tiff"),
            width = 9*0.5, height = 9*0.5, dpi = 300, asp = 0)
  
}, .progress = TRUE)

# ---------------------------
# Map future estimates and uncertainty
# ---------------------------
estimates_future_scenarios <- list.files(
  glue("{output_path}/NEX-GDDP/Future_scenarios/Predictions"),
  full.names = TRUE, pattern = "across_years\\.RDS$"
)

walk(estimates_future_scenarios, function(x) {
  
  estimates <- readRDS(x) |> filter(period_group == 2)
  filename <- basename(x)
  parts <- str_match(filename, "^monthly_incidence_by_period_(ssp[0-9]+)_(.*)\\.RDS$")
  ssp <- parts[,2]
  timeframe <- gsub("_across_years", "", parts[,3])
  
  ssp_title <- ifelse(ssp == "ssp245", "SSP2-4.5", "SSP5-8.5")
  timeframe_title <- ifelse(timeframe == "future_short_term", "short-term future", "long-term future")
  
  incidence_future <- estimates |> left_join(nuts3, by = c("nuts_id" = "NUTS_ID")) |> st_as_sf()
  extent <- st_bbox(st_crop(nuts0, c(xmin = 2377294, ymin = 1313597, xmax = 7453440, ymax = 8028510)))
  osm_eur_carto <- basemap(extent, map_service = "carto", map_type = "light_no_labels", class = "terra") |>
    project(crs(incidence_future)) |> crop(ext(incidence_future)*1.02)
  
  risk_palette <- c("#f1eef6", "#d7b5d8", "#df65b0", "#dd1c77", "#980043")
  
  # Monthly incidence map
  incidence_future_map <- tm_shape(osm_eur_carto) +
    tm_rgb(col_alpha = 0.8) +
    tm_shape(incidence_future) +
    tm_polygons(fill = "monthly_incidence", lwd = 0,
                fill.scale = tm_scale_intervals(n = 5, style = "fisher", interval.closure = "left", values = risk_palette),
                fill.legend = tm_legend(title = "", frame = FALSE, orientation = "landscape",
                                        text.size = 1.2, ticks.col = "black",
                                        position = tm_pos_out("center","bottom",pos.h="center"), show = TRUE, width = 15)) +
    tm_facets(by = "month", nrow = 3) +
    tm_title(glue("Estimated incidence per 100,000 people\nunder {ssp_title} ({timeframe_title})"),
             position = tm_pos_out("center","top",pos.h="center"), frame = FALSE, padding = c(0,0,0,0), size = 0.8) +
    tm_shape(nuts0) + tm_borders("black", lwd = 0.2) +
    tm_layout(bg.color = "white", space.color = "white", asp = 0,
              outer.margins = c(0,0,0,0), inner.margins = c(0,0,0,0),
              panel.label.bg.color = "white",
              panel.labels = setNames(month.name, seq_along(month.name)),
              panel.label.size = 1.1, panel.label.height = 0.8)
  
  tmap_save(incidence_future_map,
            glue("{output_path}/Figures/incidence_future_{ssp}_{timeframe}.tiff"),
            width = 9*0.5, height = 9*0.5, dpi = 300, asp = 0)
  
  # SD map
  monthly_sd_map <- tm_shape(osm_eur_carto) +
    tm_rgb(col_alpha = 0.8) +
    tm_shape(incidence_future) +
    tm_polygons(fill = "monthly_sd", lwd = 0,
                fill.scale = tm_scale_continuous(values = "scico.roma", limits = c(0, max(incidence_future$monthly_sd))),
                fill.legend = tm_legend(title = "", frame = FALSE, orientation = "landscape",
                                        text.size = 0.7, ticks.col = "black",
                                        position = tm_pos_out("center","bottom",pos.h="center"), show = TRUE, width = 10)) +
    tm_facets(by = "month", nrow = 3) +
    tm_title(glue("Standard deviation of estimated incidence per 100,000 people\nunder {ssp_title} ({timeframe_title})"),
             position = tm_pos_out("center","top",pos.h="center"), frame = FALSE, padding = c(0,0,0,0), size = 0.8) +
    tm_shape(nuts0) + tm_borders("black", lwd = 0.2) +
    tm_layout(bg.color = "white", space.color = "white", asp = 0,
              outer.margins = c(0,0,0,0), inner.margins = c(0,0,0,0),
              panel.label.bg.color = "white",
              panel.labels = setNames(month.name, seq_along(month.name)),
              panel.label.size = 1.1, panel.label.height = 0.8)
  
  tmap_save(monthly_sd_map,
            glue("{output_path}/Figures/incidence_future_monthly_sd_map_{ssp}_{timeframe}.tiff"),
            width = 9*0.5, height = 9*0.5, dpi = 300, asp = 0)
  
  # CV map
  monthly_cv_map <- tm_shape(osm_eur_carto) +
    tm_rgb(col_alpha = 0.8) +
    tm_shape(incidence_future) +
    tm_polygons(fill = "monthly_cv", lwd = 0,
                fill.scale = tm_scale_continuous(values = "scico.roma", limits = c(0, round(max(incidence_future$monthly_cv), 1))),
                fill.legend = tm_legend(title = "", frame = FALSE, orientation = "landscape",
                                        text.size = 0.7, ticks.col = "black",
                                        position = tm_pos_out("center","bottom",pos.h="center"), show = TRUE, width = 10)) +
    tm_facets(by = "month", nrow = 3) +
    tm_title(glue("Coefficient of variation of estimated incidence per 100,000 people\nunder {ssp_title} ({timeframe_title})"),
             position = tm_pos_out("center","top",pos.h="center"), frame = FALSE, padding = c(0,0,0,0), size = 0.8) +
    tm_shape(nuts0) + tm_borders("black", lwd = 0.2) +
    tm_layout(bg.color = "white", space.color = "white", asp = 0,
              outer.margins = c(0,0,0,0), inner.margins = c(0,0,0,0),
              panel.label.bg.color = "white",
              panel.labels = setNames(month.name, seq_along(month.name)),
              panel.label.size = 1.1, panel.label.height = 0.8)
  
  tmap_save(monthly_cv_map,
            glue("{output_path}/Figures/incidence_future_monthly_cv_map_{ssp}_{timeframe}.tiff"),
            width = 9*0.5, height = 9*0.5, dpi = 300, asp = 0)
  
}, .progress = TRUE)

# ---------------------------
# Sample NUTS 3 plots
# ---------------------------
sample_nuts3 <- c("ES300", "ITI43", "FR101", "DE300")
scenario_labels <- c("ssp245" = "under SSP2-4.5", "ssp585" = "under SSP5-8.5")

walk(c("ssp245", "ssp585"), function(scenario) {
  
  # Load short-term and long-term predictions
  short_term <- readRDS(glue("{output_path}/NEX-GDDP/Future_scenarios/Predictions/monthly_incidence_by_period_{scenario}_future_short_term_across_years.RDS")) |>
    mutate(period_group = case_when(period_group == 1 ~ "Historical-reference",
                                    period_group == 2 ~ "Short-term future"))
  
  long_term <- readRDS(glue("{output_path}/NEX-GDDP/Future_scenarios/Predictions/monthly_incidence_by_period_{scenario}_future_long_term_across_years.RDS")) |>
    mutate(period_group = case_when(period_group == 1 ~ "Historical-reference",
                                    period_group == 2 ~ "Long-term future"))
  
  pred_all <- bind_rows(short_term, long_term) |>
    mutate(period_group = factor(period_group, levels = c("Historical-reference", "Short-term future", "Long-term future")))
  
  incidence_sample_plot <- pred_all |>
    filter(nuts_id %in% sample_nuts3) |>
    ggplot() +
    geom_line(aes(month, monthly_incidence, color = period_group), size = 0.5) +
    geom_ribbon(aes(x = month, ymin = monthly_lower, ymax = monthly_upper, fill = period_group), alpha = 0.2) +
    scale_x_continuous(breaks = seq(1, 12, by = 1)) +
    scale_color_manual(values = c("Historical-reference" = "#1b9e77", "Short-term future" = "#d95f02", "Long-term future" = "#7570b3")) +
    scale_fill_manual(values = c("Historical-reference" = "#1b9e77", "Short-term future" = "#d95f02", "Long-term future" = "#7570b3")) +
    labs(x = "Month", y = "Estimated incidence per 100,000 people",
         title = paste("Estimated incidence", scenario_labels[[scenario]])) +
    facet_wrap(~ nuts_id, nrow = 2, scales = "free", labeller = labeller(
      nuts_id = c("ES300" = "Madrid", "ITI43" = "Rome", "FR101" = "Paris", "DE300" = "Berlin")
    )) +
    theme_matplotlib(base_size = 12) +
    theme(strip.text = element_text(size = 14), strip.background = element_rect(fill = "white", color = "white", size = 0.5),
          legend.position = "bottom", legend.text = element_text(size = 12), axis.text = element_text(size = 12),
          axis.title = element_text(size = 12))
  
  ggsave(glue("{output_path}/Figures/incidence_sample_plot_{scenario}.tiff"), plot = incidence_sample_plot,
         width = 10, height = 8, dpi = 300)
  
})

# ---------------------------
# Save sample NUTS 3 estimates
# ---------------------------
sample_nuts3_estimates <- map(lm_future_scenarios, function(x) {
  filename <- basename(x)
  parts <- str_match(filename, "^lm_results_(ssp[0-9]+)_(.*)\\.RDS$")
  ssp <- parts[,2]
  incidence_change <- readRDS(x)
  incidence_change$monthly_lm |> 
    filter(nuts_id %in% sample_nuts3) |> 
    mutate(nuts_name = case_when(nuts_id=="ES300" ~"Madrid",
                                 nuts_id=="ITI43" ~ "Rome",
                                 nuts_id=="FR101" ~"Paris",
                                 nuts_id=="DE300" ~ "Berlin"),
           scenario = ssp)
}) |> list_rbind()

write.csv(sample_nuts3_estimates,
          glue("{output_path}/NEX-GDDP/Future_scenarios/Linear_models/sample_nuts3_estimates.csv"),
          row.names = FALSE)

# ---------------------------
# Plot number of significant NUTS
# ---------------------------
nuts_change_significant_summary <- read.csv(glue("{output_path}/NEX-GDDP/Future_scenarios/Linear_models/nuts_change_significant_summary_grouped_by_month.csv"))

# Transform for positive/negative plotting
nuts_change_significant_summary_long <- nuts_change_significant_summary |>
  mutate(sig_neg = -sig_neg, period = factor(period, levels = c("future_short_term", "future_long_term"))) |>
  pivot_longer(cols = c(sig_pos, sig_neg), names_to = "sign", values_to = "count")

ggplot(nuts_change_significant_summary_long, aes(x = factor(month), y = count, fill = sign)) +
  geom_col(color = "black", size = 0.1) +
  scale_fill_manual(values = c(sig_pos = "darkred", sig_neg = "darkblue"),
                    labels = c(sig_pos = "Significant increase", sig_neg = "Significant decrease"), name = "") +
  scale_x_discrete(labels = month.name, expand = c(0,0)) +
  scale_y_continuous(labels = abs, expand = c(0,0)) +
  labs(x = "Month", y = "Number of NUTS 3") +
  coord_flip() +
  facet_wrap(~ ssp + period, labeller = labeller(
    ssp = c("ssp245" = "SSP2-4.5", "ssp585" = "SSP5-8.5"),
    period = c("future_short_term" = "(short-term future)", "future_long_term" = "(long-term future)")
  )) +
  theme_matplotlib(base_size = 12) +
  theme(axis.ticks.length = unit(0, "pt"), strip.text = element_text(size = 13),
        strip.background = element_rect(fill = "white", color = "white", size = 0.5),
        panel.spacing = unit(30, "pt"), legend.position = "bottom",
        legend.text = element_text(size = 13), legend.title = element_text(size = 13))

ggsave(glue("{output_path}/Figures/summary_n_nuts_future_scenarios.tiff"), width = 10, height = 8, dpi = 300)

# ---------------------------
# Plot uncertainty of all future estimates
# ---------------------------
files_all_years <- list.files("Output/NEX-GDDP/Future_scenarios/Predictions", pattern = "all_years\\.RDS$", full.names = TRUE)

format_scenario <- function(scenario) {
  case_when(
    scenario == "ssp245" ~ "SSP2-4.5",
    scenario == "ssp585" ~ "SSP5-8.5",
    TRUE ~ scenario
  )
}

uncertainty_table <- map(files_all_years, function(f) {
  scenario <- str_extract(f, "ssp[0-9]+")
  horizon <- str_extract(f, "short_term|long_term")
  
  df <- readRDS(f) |> filter(year > 2023)
  
  data.frame(
    scenario = scenario,
    horizon = horizon,
    avg_month_incidence = mean(df$monthly_incidence),
    avg_monthly_sd = mean(df$monthly_sd),
    avg_monthly_cv = mean(df$monthly_cv)
  )
}) |> list_rbind()

write.csv(uncertainty_table,
          "Output/NEX-GDDP/Future_scenarios/Predictions/average_uncertainty.csv",
          row.names = FALSE)
