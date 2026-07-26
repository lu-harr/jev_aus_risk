
make_legend <- function(title,
                        pal=cols){
  
  if (title == 'Habitat suitability index') {
    
    upper_label <- '1'
    lower_label <- '0'
    
  } else {
    
    upper_label <- 'High'
    lower_label <- 'Low'
    
  }
    
  # add legend
  long_origin <- 116
  lat_origin <- -37.5 
  
  # title text
  text(long_origin, 
       lat_origin, 
       title, 
       cex=1.5, 
       col="grey20", 
       adj=0)
  
  # vertical legend
  gradient.rect(long_origin, 
                lat_origin-2, 
                long_origin+1.5, 
                lat_origin-5-0.2, 
                col=rev(pal), 
                border='grey20',
                gradient='y')
  
  text(
    #long_origin+1.5+0.5,
    long_origin+1.5+2,
    lat_origin-2, 
    upper_label, 
    cex=1.2, 
    col="grey20", 
    adj=0.5)
  
  text(
    #long_origin+1.5+0.5,
    long_origin+1.5+2,
    lat_origin-5, 
    lower_label, 
    cex=1.2, 
    col="grey20", 
    adj=0.5)
}

sort_outputs <- function(outputs, mask_layer){
  # take a list of raster outputs
  # and raster mask layer
  
  # crop
  outputs_cropped <- lapply(outputs,
                            crop,
                            mask_layer)
  
  # make sure the extents are the same
  outputs_cropped <- lapply(outputs_cropped,
                            setExtent, 
                            ext = mask_layer)
  
  # mask
  outputs_masked <- lapply(outputs_cropped, 
                           mask,
                           mask_layer)
  
  
  return(outputs_masked)
  
}


plot_continuous_output <- function(output,
                                   cols,
                                   title,
                                   region,
                                   legend = TRUE,
                                   base.legend = FALSE,
                                   lgas = NULL,
                                   zlim = NULL,
                                   axis.args = NULL,
                                   labels = NULL,
                                   legend.mar = 5.1) { 
  # adding in legend.mar and leaving it as default
  # removing borders (this function only ends up used for potential niche fig)
  
  
  if (region == 'national'){
    
    borders <- shapefile('data/aus_simplified.shp')
    
  } else {
    
    borders <- shapefile(paste0('data/state_shp/', region, '_simple.shp'))
    
  }  
  
  plot(output, 
       col = cols,
       legend = base.legend,
       axes = FALSE,
       box = FALSE,
       main = title,
       cex.main = 1.5,
       zlim = zlim,
       axis.args = axis.args,
       labels = labels,
       legend.mar = legend.mar)
  
  # if (is.null(lgas)) {
  #   
  #   # add borders
  #   plot(borders, 
  #        lwd = 0.5, 
  #        col = "transparent",
  #        border = 'grey20',
  #        add = TRUE)
  #   
  # } else {
  #   
  #   # add borders
  #   plot(lgas, 
  #        lwd = 0.6, 
  #        col = "transparent",
  #        border = 'grey30',
  #        add = TRUE)
  #   
  # }  
  
  #add a legend 
  if (legend == TRUE) {
    
    make_legend(region,
                "Suitability",
                cols)
  }
  
}

plot_blobby_output <- function(output,
                               cols,
                               i,
                               region){
  
  if (region == 'national'){
    
    borders <- shapefile('data/aus_simplified.shp')
    
  } else {
    
    borders <- shapefile(paste0(path, 'data/state_shp/', region, '_simple.shp'))
    
  }  
  
  plot(output, 
       col = cols,
       legend = FALSE,
       axes = FALSE,
       box = FALSE,
       main = paste0("Thresholded suitability of\n", i, " hosts\nfor JEV persistence"),
       cex.main = 1.5)
  
  # add borders
  plot(borders, 
       lwd = 1.2, 
       col = "transparent",
       border = 'grey20',
       add = TRUE)
  
  legend("bottomleft", legend = c("50%","80%","90%"),
         fill = rev(cols[2:4]), cex=1.3)
  
} 


subfigure_label = function(plot_region, x_displacement, y_displacement, label,
                           cex.label=1){
  # gives me an (a)
  text(plot_region[1] + (plot_region[2] - plot_region[1])*x_displacement, 
       plot_region[3] + (plot_region[4] - plot_region[3])*y_displacement,
       label, cex=cex.label)
}


empty_plot_for_legend = function(){
  plot(0, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
}


make_legend_for_empty_plot = function(title, pal=cols){
  # fills in fourth panel in transmission suitability figure with legend
  if (title == 'Habitat suitability index') {
    upper_label <- '1'
    lower_label <- '0'
  } else {
    upper_label <- 'High'
    lower_label <- 'Low'
  }
  
  # add legend
  long_origin <- 0.7
  lat_origin <- 0.8
  
  # title text
  text(long_origin, lat_origin, 
       title, cex=1.6, 
       col="grey20", adj=0)
  
  # vertical legend
  gradient.rect(long_origin, lat_origin-0.15, 
                long_origin+0.04, lat_origin-0.8, 
                col=rev(pal), 
                border='grey20',
                gradient='y')
  
  # low
  text(long_origin+0.08, lat_origin-0.15, 
       upper_label, cex=1.4, col="grey20", adj=0.5)
  
  # high
  text(long_origin+0.08, lat_origin-0.8, 
       lower_label, cex=1.4, col="grey20", adj=0.5)
}

