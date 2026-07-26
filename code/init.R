
library(sf)
library(raster)
#library(rgdal)

library(viridisLite)
library(RColorBrewer)
library(scales)
library(plotrix)

library(dplyr)
library(tidyr)
library(readr)

setwd("~/Desktop/jev/pathogen_paper")

# load Australia mask
# throw this out? I don't end up using it? It doesn't match other stuff?
# But it's used in pathogen_model.R? Gah!
aus_mask <- raster('data/aus_only_mask.tif')

# load aus state boundaries
aus_simple <- st_read('data/aus_simplified.shp')

AGG_FACTOR <- 10

reso <- 300
