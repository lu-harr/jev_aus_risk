library(raster)
library(dismo)

thresholdify <- function(suitability_raster, point_coords, inclusion_quantile = 0.9) {
  
  # get the threshold value
  vals <- raster::extract(suitability_raster, point_coords)
  threshold <- quantile(vals, 1 - inclusion_quantile, na.rm=TRUE)
  
  # binarise the raster
  binary_raster <- suitability_raster >= threshold
  binary_raster
  
}

blobify <- function(binary_raster, blob_distance = 0.1, agg_factor=1) {
  
  # blobbyify the raster
  binary_raster_agg <- aggregate(binary_raster, agg_factor, fun=max, na.rm=TRUE)
  buffer_window <- focalWeight(binary_raster_agg, blob_distance, "circle")
  buffer_window <- buffer_window / max(buffer_window)
  blobby_binary_raster <- focal(binary_raster_agg, buffer_window, fun = max, na.rm = TRUE)
  blobby_binary_raster <- mask(blobby_binary_raster, binary_raster_agg)
  
  blobby_binary_raster
}

# threshold and blobbify a raster 
threshy_blobs <- function(suitability_raster, point_coords, inclusion_quantile = 0.9, blob_distance = 0.1) {
  
  binary_raster <- thresholdify(suitability_raster, point_coords, inclusion_quantile)
  blobby_binary_raster <- blobify(binary_raster, blob_distance)
  blobby_binary_raster
  
}
