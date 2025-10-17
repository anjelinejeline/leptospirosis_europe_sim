#####################
### Author: Angela Fanelli
### Role: Epidemiologist
### Institution: JRC
#####################

# ---------------------------
# Load the required packages
# ---------------------------

libraries <- c(
  # Data manipulation and reading
  "tidyverse", "dplyr", "lubridate", "purrr", "readxl", "ISOweek", "data.table",
  
  # Visualization
  "corrplot", "ggsci", "ggplot2", "gridExtra", "classInt",
  "ggridges", "colorspace", "scales", "flextable",
  
  # GIS and spatial
  "giscoR", "eurostat", "sf", "tmap", "basemaps", "terra", "sp", "spdep",
  
  # Modeling
  "INLA", "INLAOutputs",
  
  # Utility
  "glue"
)

# ---------------------------
# Install and load packages
# ---------------------------
for (lib in libraries) {
  if (!(lib %in% installed.packages())) {
    
    # Special installation for INLA
    if (lib == "INLA") {
      install.packages(
        "INLA", 
        repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), 
        dep = TRUE
      )
    } else {
      install.packages(lib)
    }
  }
  
  library(lib, character.only = TRUE)
}


