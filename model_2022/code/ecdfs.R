# Only difference here is pathogen_presence_cut?
# Using birds/vectors from 2023 model?
setwd("..")
source("code/init.R")
source("code/plotting_functions.R")
set.seed(834903)

# read all of the separate species maps in:
birds <- brick('data/egrons_brts_new.tif')
birds.feb.22 <- birds[[which(names(birds)=="feb_2022")]]

vectors.feb.22 <- raster('data/cx_oi_wa_reup_mxnt__feb_22.tif')

abund_feral_pig <- raster('data/Feral_Pig_density_5k.tif') %>%
  projectRaster(vectors.feb.22)

# fill in feral pig land NAs with zeroes - use same mask as birds/vectors
fpigs_inv <- abund_feral_pig
values(fpigs_inv) <- ifelse(
  is.na(values(abund_feral_pig)), 1, NA)
values(abund_feral_pig)[is.na(values(abund_feral_pig))] <- ifelse(
  is.na(raster::extract(vectors.feb.22, rasterToPoints(fpigs_inv)[,1:2])), NA, 0)

dom_pig_density <- raster('data/national_piggery_density.tif') %>% 
  disaggregate(fac=10, method='bilinear')/100

# stack reservoir surfaces
suitability <- stack(
  vectors.feb.22,
  birds.feb.22,
  abund_feral_pig
) %>% 
  setExtent(dom_pig_density)

suitability <- stack(
  suitability,
  dom_pig_density)

# name them
names(suitability)  <- c(
  "vectors", 
  "birds", 
  "feral_pigs", 
  "domestic_pigs"
)

#################
# trim and mask them

suitability_aus <- mask(
  setExtent(
    crop(
      suitability, 
      aus_mask), 
    aus_mask), 
  aus_mask
)

#######################
# load and organise JEV infection data
aus_dat <- read_csv('model_2022/data/pathogen_presence_national_cut.csv') %>% 
  dplyr::select(
    'longitude', 
    'latitude'
  )


######################
# compute the thresholds

suitability <- suitability_aus
pos_points <- aus_dat

# extract suitabilities at presence locations to build ECDF objects
# extractions for the presence should technically occur according to the date of specimen collection
# but there's no temporal treatment here

# extract suitabilities at locations of detection
bird_suitability_at_pos <- raster::extract(suitability[[which(names(suitability)=='birds')]], pos_points)
vector_suitability_at_pos <- raster::extract(suitability[[which(names(suitability)=='vectors')]], pos_points)
feral_pig_suitability_at_pos <- raster::extract(suitability[[which(names(suitability)=='feral_pigs')]], pos_points)
domestic_pig_suitability_at_pos <- raster::extract(suitability[[which(names(suitability)=='domestic_pigs')]], pos_points)

# build empirical cumulative density function
bird_ecdf <- ecdf(bird_suitability_at_pos)
vector_ecdf <- ecdf(vector_suitability_at_pos)
feral_pig_ecdf <- ecdf(feral_pig_suitability_at_pos)
domestic_pig_ecdf <- ecdf(domestic_pig_suitability_at_pos)

# re-scale suitability layer subject to ECDF
bird_suitability_ecdf <- vector_suitability_ecdf <- feral_suitability_ecdf <- domestic_suitability_ecdf <- suitability[[which(names(suitability)=='birds')]]
bird_suitability_ecdf[] <- bird_ecdf(suitability[[which(names(suitability)=='birds')]][])
vector_suitability_ecdf[] <- vector_ecdf(suitability[[which(names(suitability)=='vectors')]][])
feral_suitability_ecdf[] <- feral_pig_ecdf(suitability[[which(names(suitability)=='feral_pigs')]][])
domestic_suitability_ecdf[] <- domestic_pig_ecdf(suitability[[which(names(suitability)=='domestic_pigs')]][])

# now consider combinations of groups
host_composition <- c('vectors and avian',
                      'vectors, avian, and feral pigs',
                      'vectors, avian, and domestic pigs')

cont_outputs <- list()

for (i_hosts in host_composition) {
  
  if (i_hosts == 'vectors and avian') {
    # vectors and birds only
    combined_stack <- stack(bird_suitability_ecdf, vector_suitability_ecdf)
    
  } else if (i_hosts == 'vectors, avian, and feral pigs') {
    # vectors and birds and feral pigs
    combined_stack <- stack(bird_suitability_ecdf, vector_suitability_ecdf, feral_suitability_ecdf)
    
    
  } else {
    # vectors and birds and all pigs
    combined_stack <- stack(bird_suitability_ecdf, vector_suitability_ecdf, feral_suitability_ecdf, domestic_suitability_ecdf)
    
  }
  
  combined_suitability_ecdf <- calc(combined_stack, min)
  
  writeRaster(combined_suitability_ecdf,
              paste0('output/continuous suit ', i_hosts, '.tif'),
              overwrite=TRUE)
  
  cont_outputs[[i_hosts]] <- combined_suitability_ecdf
  
}  


# MAKE THE FIGURE

###############################################################################
# TRANSMISSION SUITABILITY
# edit mfrow, mar, remove figure titles

# specify colours
pn_cols <- colorRampPalette(colors = c("#f7f7f7", "#c23375"))(1000)
# plot png  
reso <- 300

cont_outputs = lapply(cont_outputs, trim) # let's get more map for our map

{png('model_2022/figures/potential_niche.png',
     width = 1000*reso/72, 
     height = 900*reso/72, 
     res = reso,
     bg = 'white')
  
  par(mfrow = c(2,2),
      #mar = c(0.5, 0.5, 0.5, 0.5), 
      mar=c(0.5,0.5,0.5,0),
      bty="n",
      #oma=c(6.5,2.5,0.5,0)
      oma = c(0,0,0,0))
  
  # continuous prediction
  plot_continuous_output(cont_outputs[[1]], pn_cols, '', region='national', legend=FALSE)
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1.6)
  plot_continuous_output(cont_outputs[[2]], pn_cols, '', region='national', legend=FALSE)
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1.6)
  plot_continuous_output(cont_outputs[[3]], pn_cols, '', region='national', legend=FALSE)
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", 1.6)
  empty_plot_for_legend()
  make_legend_for_empty_plot("Suitability", pn_cols)
  
  dev.off()}


# with labelled map of aus ....
source("code/aus_labelled_figure.R")

{png('model_2022/figures/potential_niche_with_aus_labelled_koppen.png',
     width = 1000*reso/72, 
     height = 1000*reso/72, 
     res = reso,
     bg = 'white')
  
  par(mfrow = c(2,2),
      #mar = c(0.5, 0.5, 0.5, 0.5), 
      mar=c(0.5,0.5,0.5,0),
      bty="n", xpd=TRUE,
      #oma=c(6.5,2.5,0.5,0)
      oma = c(8,0,0,0))
  
  # plot(aus_elevation, breaks=c(seq(0, maxValue(aus_elevation), 200), maxValue(aus_elevation)), col=purps,
  #      legend=FALSE, 
  #      xaxt="n", yaxt="n")
  # gradient.rect(120, -40, 135, -39,
  #               col=purps[1:7],
  #               border='grey20',
  #               gradient='x')
  # axis(1, at=c(120,135), labels=c("Low", "High"), line=-5, cex.axis=1.2)
  # text(x=(135+120)/2, y=-41, labels="Elevation", cex=1.4)
  
  plot(clim_zones, col=alpha(clim_pal, 0.4), legend=FALSE, xaxt="n", yaxt="n")
  legend(117,-37, rev(c("Temperate", "Grassland", "Desert", "Subtropical", "Tropical", "Equatorial")),
         fill=rev(clim_pal), title="Koppen classification", bty="n")
  lines(c(116.8, 137, 137, 116.8, 116.8), c(-37, -37, -44.7, -44.7, -37))
  axis(1, at=c(128, 130.5, 133), labels=c(0, 250, 500), line=-2.3) # distHaversine(c(128,-25), c(133, -25))/1000
  text(130.5, -41.5, "km")
  # plot(aus_elevation, breaks=c(seq(0, maxValue(aus_elevation), 200), maxValue(aus_elevation)), col=purps,
  #      legend=FALSE,
  #      xaxt="n", yaxt="n")
  # plot(aus_elevation, col=purps[1:7], legend.only=TRUE, border="grey80", lwd=2,
  #      legend.args=list(text="Elevation", side=2, line=0.5), legend.width=1.2,
  #      axis.args=list(at=c(minValue(aus_elevation) + 100, maxValue(aus_elevation) - 100), c("Low", "High")),
  #      bty="o")
  plot_mdb("#D780A9", 0.2, lwd=1.2)
  plot(aus_simple, add=TRUE, border="grey20", lwd=1.2)
  plot_rivers(lwd_fac=0.6)
  plot_state_labels(cex=1.5, alph=0.6)
  plot_great_dividing_range(cex=1.2)
  plot_murray()
  par(xpd=TRUE) # Don't let city labels fall off the plot!
  plot_cities(cex=1.2)
  text(143.5, -16, "Cape\nYork", adj=c(0.5, 0.5))
  text(142.5, -10, "Torres Strait", adj=c(0.5,0.5))
  
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1.6)
  
  # continuous prediction
  plot_continuous_output(cont_outputs[[1]], pn_cols, '', region='national', legend=FALSE)
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1.6)
  plot_continuous_output(cont_outputs[[2]], pn_cols, '', region='national', legend=FALSE)
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", 1.6)
  plot_continuous_output(cont_outputs[[3]], pn_cols, '', region='national', legend=FALSE)
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", 1.6)
  #empty_plot_for_legend()
  make_legend_for_empty_plot("Suitability", pn_cols)
  
  par(new=TRUE, mfrow=c(1,1), mar=c(4.5,0,0,0), oma=c(0.5,0,0,0))
  plot(0, type="n", new=TRUE, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", xlab="Suitability", 
       cex.lab=1.6)
  gradient.rect(0.4, -0.04, 0.6, -0.02,
                col=pn_cols,
                border='grey20',
                gradient='x')
  axis(1, at=c(0.4,0.6), line=0, labels=c("Low", "High"), cex.axis=1.4)
  
  dev.off()}
