# Code and data to accompany the article "Spatiotemporal dynamics of leptospirosis in Europe and projected climate change impacts"

## Description

This repository contains the data and code used to perform the analyses, described in the article *Spatiotemporal dynamics of leptospirosis in Europe and projected climate change impacts*.

The repository is divided as follows:

-   Input
-   Scripts

It assumes that the results will be saved in a folder called Output and relative subfolders. For any issues with the code please contact [Angela Fanelli](Angela.FANELLI@ec.europa.eu).

## Data Availability

-   **Human cases**: Available upon request from [The European Surveillance System (TESSy)](https://www.ecdc.europa.eu/en/publications-data/european-surveillance-system-tessy).
-   **Monthly temperature**: Retrieved from the [Copernicus Climate Data Store – ERA5](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means).
-   **Forest–human nexus**: Accessible via the code provided by [Massaro et al. (2025)](https://www.nature.com/articles/s43247-025-02514-8), available on [GitHub](https://github.com/emanuelemassaro/Forest_Human_Nexus).
-   **Mammal species richness**: Available from the [IUCN Red List Spatial Dataset](https://www.iucnredlist.org/resources/spatial-data-download), processed following the approach developed by [the Knowledge Centre for Biodiversity – Global Biodiversity Data (KCBD-GBD; ex-DOPA) unit of the JRC](https://www.mdpi.com/2073-445X/13/9/1506).
-   **Population data**: Sourced from the [Global Human Settlement Layer (GHSL)](https://human-settlement.emergency.copernicus.eu/).
-   **SPEI data**: Obtained from the [SPEIbase v.2.10](https://digital.csic.es/handle/10261/364137) developed by the Climatology and Climate Services Laboratory.
-   **Future climate projections**:
    -   **SPEI**: From the [NEX-GDDP-CMIP6 dataset](https://www.ciesin.columbia.edu/data/globaldrought/) via NASA’s Socioeconomic Data and Applications Center (SEDAC).
    -   **Temperature**: From the [NASA Center for Climate Simulation (NCCS)](https://www.nccs.nasa.gov/services/data-collections/land-based-products/nex-gddp).
-   **Healthcare capacity**: Number of hospital beds per country obtained from [EUROSTAT](https://ec.europa.eu/eurostat/databrowser/view/HLTH_RS_BDS1/default/table?lang=en).

## Files in this repository

The *Scripts* folder contains all R scripts, while the *Input* folder holds the data files sourced by these scripts. Please note that the scripts are applied to a simulated dataset, and the results will slightly differ from those presented in the manuscript.

### Input

This folder contains all data used as input for the project.

#### Input/Data_model/data_2010_2023.csv

This CSV contains simulated leptospirosis cases and associated covariates for the reference period (2010–2023).

| Column | Description |
|---------------------|---------------------------------------------------|
| `country_code` | Country code of the European country |
| `nuts_index` | Unique index of the NUTS region, used for the spatial adjacency matrix |
| `nuts_id` | NUTS 3 code (version 2021) |
| `nuts_name` | Full name of the NUTS 3 region |
| `disease` | Disease name |
| `month` | Month of observation |
| `year` | Year of observation |
| `cases` | Number of human cases (simulated) |
| `t2m_lag` | 2-metre air temperature lagged by 1 month (Tmean in the article) |
| `spei_lag` | 3-month Standardised Precipitation–Evapotranspiration Index (SPEI), lagged by 1 month |
| `population` | Population of the NUTS 3 region |
| `fhn` | Forest human nexus |
| `mammal_richness` | Mammal species richness |

#### Input/NEX-GDPP

Contains the 3-month SPEI and Tmean projections for all selected climate models.

#### Input/Historical_climate

Contains the 3-month SPEI and Tmean data for the period 2004–2009.

### Session info

```         
─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.5.0 (2025-04-11)
 os       Ubuntu 22.04.5 LTS
 system   x86_64, linux-gnu
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Rome
 date     2025-10-17
 rstudio  2025.05.0+496 Mariposa Orchid (desktop)
 pandoc   3.4 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64/ (via rmarkdown)

─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package     * version  date (UTC) lib source
 base        * 4.5.0    2025-05-04 [4] local
 basemaps    * 0.0.8    2024-11-01 [1] CRAN (R 4.5.0)
 classInt    * 0.4-11   2025-01-08 [1] CRAN (R 4.5.0)
 colorspace  * 2.1-1    2024-07-26 [1] CRAN (R 4.5.0)
 corrplot    * 0.92     2021-11-18 [3] CRAN (R 4.2.0)
 data.table  * 1.17.2   2025-05-12 [1] CRAN (R 4.5.0)
 datasets    * 4.5.0    2025-05-04 [4] local
 dplyr       * 1.1.4    2023-11-17 [3] CRAN (R 4.3.2)
 eurostat    * 4.0.0    2023-12-19 [1] CRAN (R 4.5.0)
 flextable   * 0.9.7    2024-10-27 [1] CRAN (R 4.5.0)
 forcats     * 1.0.0    2023-01-29 [3] CRAN (R 4.2.2)
 ggplot2     * 3.5.2    2025-04-09 [2] CRAN (R 4.5.0)
 ggridges    * 0.5.6    2024-01-23 [3] CRAN (R 4.3.2)
 ggsci       * 3.0.0    2023-03-08 [3] CRAN (R 4.2.2)
 giscoR      * 0.6.1    2025-01-27 [1] CRAN (R 4.5.0)
 glue        * 1.8.0    2024-09-30 [2] CRAN (R 4.5.0)
 graphics    * 4.5.0    2025-05-04 [4] local
 grDevices   * 4.5.0    2025-05-04 [4] local
 gridExtra   * 2.3      2017-09-09 [3] CRAN (R 4.0.1)
 INLA        * 25.09.04 2025-09-04 [1] local
 INLAOutputs * 1.4.11   2025-08-11 [1] Github (oswaldosantos/INLAOutputs@276db00)
 ISOweek     * 0.6-2    2011-09-07 [1] CRAN (R 4.5.0)
 lubridate   * 1.9.3    2023-09-27 [3] CRAN (R 4.3.1)
 Matrix      * 1.7-3    2025-03-11 [4] CRAN (R 4.4.3)
 methods     * 4.5.0    2025-05-04 [4] local
 purrr       * 1.0.4    2025-02-05 [2] CRAN (R 4.5.0)
 readr       * 2.1.5    2024-01-10 [3] CRAN (R 4.3.2)
 readxl      * 1.4.3    2023-07-06 [3] CRAN (R 4.3.1)
 scales      * 1.4.0    2025-04-24 [1] CRAN (R 4.5.0)
 sf          * 1.0-21   2025-05-15 [1] CRAN (R 4.5.0)
 sp          * 2.2-0    2025-02-01 [1] CRAN (R 4.5.0)
 spData      * 2.3.0    2023-07-06 [3] CRAN (R 4.3.1)
 spdep       * 1.3-1    2023-11-23 [3] CRAN (R 4.3.2)
 stats       * 4.5.0    2025-05-04 [4] local
 stringr     * 1.5.1    2023-11-14 [3] CRAN (R 4.3.2)
 terra       * 1.8-50   2025-05-09 [1] CRAN (R 4.5.0)
 tibble      * 3.2.1    2023-03-20 [3] CRAN (R 4.3.1)
 tidyr       * 1.3.1    2024-01-24 [3] CRAN (R 4.3.2)
 tidyverse   * 2.0.0    2023-02-22 [3] CRAN (R 4.2.2)
 tmap        * 4.1      2025-05-26 [1] Github (r-tmap/tmap@e24262a)
 utils       * 4.5.0    2025-05-04 [4] local

 [1] /home/panelan/R/x86_64-pc-linux-gnu-library/4.5
 [2] /usr/local/lib/R/site-library
 [3] /usr/lib/R/site-library
 [4] /usr/lib/R/library

───────────────────────────────────────────────────────────────────────────────────
```
