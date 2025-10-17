#####################
### Author: Angela Fanelli
### Role: Epidemiologist
### Institution: JRC
#####################

# ---------------------------
# Source required scripts
# ---------------------------
source("Scripts/00_Load_packages.R")

# ---------------------------
# Create adjacency matrix if it does not exist
# ---------------------------
nuts3_graph <- "Input/Data_model/nuts3_graph"

if (!file.exists(nuts3_graph)) {
  
  print(glue("map_graph does not exist... creating map_graph at {Sys.time()}"))
  
  # Load NUTS3 shapefile for Europe
  nuts3 <- gisco_get_nuts(nuts_level = 3, year = 2021, epsg = 3035)
  
  # Load modeling data
  data <- fread("Input/Data_model/data_2010_2023.csv", header = TRUE)
  
  # Filter nuts3 to match data and convert to sf
  nuts3 <- nuts3 |>
    filter(NUTS_ID %in% data$nuts_id) |>
    st_as_sf()
  
  # Reorder rows of nuts3 based on the order of NUTS_ID in data
  nuts3 <- nuts3 |>
    arrange(match(NUTS_ID, data$nuts_id))
  
  # Quick check
  # head(nuts3)
  # head(data)
  
  # Create adjacency matrix
  nuts3_nb <- poly2nb(as_Spatial(nuts3$geometry))
  nb2INLA(nuts3_graph, nuts3_nb)
  
  # Load the graph into INLA
  g <- inla.read.graph(filename = "Input/Data_model/nuts3_graph")  
  
} else {
  
  print(glue("map_graph exists .. loading it"))
  
  # Load existing graph
  g <- inla.read.graph(filename = "Input/Data_model/nuts3_graph")
}

# ---------------------------
# Create data for modeling if it does not exist
# ---------------------------
data_model <- "Input/Data_model/data_model.csv"

if (!file.exists(data_model)) {
  
  print(glue("data_model does not exist... creating data_model at {Sys.time()}"))
  
  # Load raw data
  data <- fread("Input/Data_model/data_2010_2023.csv", header = TRUE)
  
  # Prepare data for modeling
  data_model <- data |>
    # Reclassify variables for smooth terms in the model
    mutate(
      t2m_lag_g = inla.group(data$t2m_lag, method = "cut", n = 20),
      spei_lag_g = inla.group(data$spei_lag, method = "cut", n = 20),
      fhn_g = inla.group(data$fhn, method = "cut", n = 20),
      mammal_richness_g = inla.group(data$mammal_richness, method = "cut", n = 20)
    ) |>
    # Rescale population for offset
    mutate(population = population / 10^5) |>
    # Rename outcome variable
    rename(Y = cases) |>
    # Select relevant columns
    select(
      country_code, nuts_index, nuts_id, nuts_name, 
      disease, month, year, Y, 
      t2m_lag_g, spei_lag_g, fhn_g, mammal_richness_g,
      population
    )
  
  # Save processed data for future use
  fwrite(data_model, file = "Input/Data_model/data_model.csv", row.names = FALSE)
  
} else {
  
  print(glue("data_model exists.. loading data at {Sys.time()}"))
  
  # Load existing data_model
  data_model <- fread(data_model, header = TRUE)
}

