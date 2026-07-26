# TODO
# - CHECK I have the right birds/vectors from DHD (I didn't notice the difference ...)
# - pathogen_presence_cut.csv

# ade <- c(137, 140, -37.5, -34.5)

source("code/init.R")
source("code/plotting_functions.R")
set.seed(834903)

# read all of the separate species maps in:
birds.feb.22 <- raster('data/enms/brt_egrons_halfwt_tc4__brt_feb_22.tif')
#birds.feb.22 <- birds[[which(names(birds)=="feb_2022")]]

#vectors.feb.221 <- raster('data/enms/cx_oi_bg1x_noAccess__mxnt_feb_22.tif')
# it's now August 2025 and DHD has submitted MS:
vectors.feb.22 <- raster("data/enms/Cx_oi_bg_tailor_rm3_mxnt_feb_22.tif")

# ooft I mighta broken it .. for some reason piping doesn't work on last step?
abund_feral_pig <- raster('data/Feral_Pig_density_5k.tif') %>% # abund 2-83 per pixel, at 5kmsq
  disaggregate(fac=5, method='')/25 #%>% # don't want this to be interpolated .. abund now 0.08 - 3.32 per pixel
abund_feral_pig <- projectRaster(abund_feral_pig, vectors.feb.22) # also 2-83

# fill in feral pig land NAs with zeroes - use same mask as birds/vectors
fpigs_inv <- abund_feral_pig
values(fpigs_inv) <- ifelse(
  is.na(values(abund_feral_pig)), 1, NA)
values(abund_feral_pig)[is.na(values(abund_feral_pig))] <- ifelse(
  is.na(raster::extract(vectors.feb.22, rasterToPoints(fpigs_inv)[,1:2])), NA, 0)

dom_pig_density <- raster('data/national_piggery_density.tif') %>% # this is at 11km
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
  dom_pig_density) %>%
  mask(vectors.feb.22)

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
# LH WATCH OUT HERE
# aus_dat <- read_csv('data/pathogen_1/pathogen_presence_national_cut.csv') %>% 
#   dplyr::select(
#     'longitude', 
#     'latitude'
#   )

aus_dat <- read_csv('data/secure/phase_2/clean/pathogen/pathogen_presence_national_cut.csv') %>% 
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

writeRaster(bird_suitability_ecdf, "output/scaled_birds.tif", overwrite = TRUE)
writeRaster(vector_suitability_ecdf, "output/scaled_vecs.tif", overwrite = TRUE)

# now consider combinations of groups
host_composition <- c('vectors and avian',
                      'vectors, avian, and feral pigs',
                      'vectors, avian, and domestic pigs')

cont_outputs <- list()

# AND-OR loop: min(vectors, max(hosts))
for (i_hosts in host_composition) {
  
  if (i_hosts == 'vectors and avian') {
    # vectors and birds only
    combined_stack <- stack(bird_suitability_ecdf, vector_suitability_ecdf)
    combined_suitability_ecdf <- calc(combined_stack, min)
    
  } else if (i_hosts == 'vectors, avian, and feral pigs') {
    # vectors and birds and feral pigs
    host_ras <- calc(stack(bird_suitability_ecdf, feral_suitability_ecdf), max)
    combined_stack <- stack(host_ras, vector_suitability_ecdf)
    #combined_stack <- stack(bird_suitability_ecdf, vector_suitability_ecdf, feral_suitability_ecdf)
    combined_suitability_ecdf <- calc(combined_stack, min)
    
    
  } else {
    # vectors and birds and all pigs
    host_ras <- calc(stack(bird_suitability_ecdf, feral_suitability_ecdf, domestic_suitability_ecdf), max)
    combined_stack <- stack(host_ras, vector_suitability_ecdf)
    #combined_stack <- stack(bird_suitability_ecdf, vector_suitability_ecdf, feral_suitability_ecdf, domestic_suitability_ecdf)
    combined_suitability_ecdf <- calc(combined_stack, min)
  }
  
  # combined_suitability_ecdf <- calc(combined_stack, min)
  
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

{png('figures/potential_niche.png',
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

{png('figures/potential_niche_with_aus_labelled.png',
     width = 1000*reso/72, 
     height = 1100*reso/72, 
     res = reso,
     bg = "transparent")
  
  par(mfrow = c(2,2),
      #mar = c(0.5, 0.5, 0.5, 0.5), 
      mar=c(0.5,0.5,2.5,0),
      bty="n", xpd=NA,
      #oma=c(6.5,2.5,0.5,0)
      oma = c(6,0,0,0))
  
  # plot(clim_zones, col=alpha(clim_pal, 0.4), legend=FALSE, xaxt="n", yaxt="n")
  # legend(117,-37, rev(c("Temperate", "Grassland", "Desert", "Subtropical", "Tropical", "Equatorial")),
  #        fill=alpha(rev(clim_pal), 0.8), title="Koppen classification", bty="n")
  # lines(c(116.8, 137, 137, 116.8, 116.8), c(-37, -37, -44.7, -44.7, -37))
  # axis(1, at=c(128, 130.5, 133), labels=c(0, 250, 500), line=-4.5) # distHaversine(c(128,-25), c(133, -25))/1000
  # text(130.5, -41.5, "km")
  # plot_mdb("#D780A9", 0.2, lwd=1.2)
  # plot(st_geometry(aus_simple), add=TRUE, border="grey20", lwd=1.2)
  # plot_rivers(lwd_fac=0.6)
  # plot_state_labels(cex=1.5, alph=0.6)
  # plot_great_dividing_range(cex=1.2)
  # plot_murray()
  # par(xpd=TRUE) # Don't let city labels fall off the plot!
  # plot_cities(cex=1.2)
  # text(143.5, -16, "Cape\nYork", adj=c(0.5, 0.5))
  # text(142.5, -10, "Torres Strait", adj=c(0.5,0.5))
  
  plot(raster(clim_zones), col=alpha(clim_pal, 0.4), legend=FALSE, xaxt="n", yaxt="n")
  legend(117,-37, rev(c("Summer Dominant", "Arid", "Summer", "Uniform", "Winter", "Winter Dominant")),
         fill=rev(clim_pal), title="Long-term Rainfall Pattern", bty="n")
  lines(c(116.8, 137, 137, 116.8, 116.8), c(-37, -37, -44.7, -44.7, -37))
  axis(1, at=c(129, 131.5, 134), labels=c(0, 250, 500), line=-4.5) # distHaversine(c(128,-25), c(133, -25))/1000
  text(131.5, -41, "km")
  plot_mdb("#D780A9", points = TRUE, alph=0.7, lwd=1)
  plot(st_geometry(aus_simple), add=TRUE, border="grey50", lwd=3)
  plot_rivers(lwd_fac=1.5)
  plot_state_labels()
  plot_great_dividing_range()
  plot_murray()
  par(xpd = TRUE)
  plot_cities()
  text(143.5, -16, "Cape\nYork", adj=c(0.5, 0.5))
  text(126.5, -17, "Kimberley", adj=c(0.5, 0.5))
  
  par(xpd=TRUE)
  subfigure_label(par("usr"), 0.1, 1, "(a)", 1.6)
  
  # continuous prediction
  
  plot_continuous_output(cont_outputs[[1]], pn_cols, 
                         "", region='national', legend=FALSE)
  par(xpd=NA)
  mtext("Predicted transmission suitability:\n mosquitoes and wading birds", 3, 
        font=2, cex=1.5, line=-1)
  subfigure_label(par("usr"), 0.1, 1, "(b)", 1.6)
  plot_continuous_output(cont_outputs[[2]], pn_cols, 
                         "", region='national', legend=FALSE)
  par(xpd=NA)
  mtext("Predicted transmission suitability:\n mosquitoes, wading birds and feral pigs", 3, 
        font=2, cex=1.5, line=-1)
  subfigure_label(par("usr"), 0.1, 1, "(c)", 1.6)
  plot_continuous_output(cont_outputs[[3]], pn_cols, 
                         "", region='national', legend=FALSE)
  par(xpd=NA)
  mtext("Predicted transmission suitability:\n mosquitoes, wading birds,\n feral and domestic pigs", 3, 
        font=2, cex=1.5, line=-1)
  subfigure_label(par("usr"), 0.1, 1, "(d)", 1.6)
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


###############################################################################
# SUPP FIGURE: ECDFS ON HISTOGRAMS OF UNSCALED SURFACES

pal=viridis(5)
cexlab = 1.4
ecdfcol = "#c23375"

ecdf_subplot_chunky <- function(all, 
                                jev, 
                                logx=FALSE, 
                                epsilon=1, 
                                breaks=50,
                                title="Vectors",
                                xlab="Unscaled suitability"){
  if (logx){
    #all <- log10(all + epsilon)
    #jev <- log10(jev + epsilon)

    all <- log(all + epsilon)
    jev <- log(jev + epsilon)
  }
  
  all_dens <- hist(all, breaks=breaks, plot=FALSE)
  jev_dens <- hist(jev, breaks=all_dens$breaks, plot=FALSE)
  
  max_density <- max(all_dens$density, jev_dens$density)
  xlim <- c(min(0, min(all, na.rm=TRUE)),
            ifelse(logx, max(all, na.rm=TRUE), 1))
  
  plot(0,
       type="n",
       xlim = xlim,
       ylim = c(0, max_density),
       axes = FALSE,
       xlab = xlab,
       ylab = "",
       main = title,
       cex.main=2,
       cex.lab=1.5)
  
  rect(xleft = all_dens$breaks[1:length(all_dens$breaks) - 1],
        xright = all_dens$breaks[2:length(all_dens$breaks)],
        ybottom = 0,
        ytop = all_dens$density,
       border = grey(0.8),
       col = grey(0.9))
  
  rect(xleft = jev_dens$breaks[1:length(jev_dens$breaks) - 1],
       xright = jev_dens$breaks[2:length(jev_dens$breaks)],
       ybottom = 0,
       ytop = jev_dens$density,
       border = ecdfcol,
       col = alpha(ecdfcol, 0.5))
  
  
  ecdf_x <- seq(xlim[1], xlim[2], length.out = 1000)
  ecdf_y <- ecdf(jev)(ecdf_x) * max_density
  
  lines(ecdf_x,
        ecdf_y,
        col = ecdfcol,
        lwd = 2)
  if (!logx){
    axis(1)
  }
  #axis(1)
  
}


{png('figures/ecdfs_chunky.png',
     width = 1100*reso/72,
     height = 700*reso/72,
     res = reso,
     bg = 'white')

par(mfrow=c(2,2))

par(mar = c(5, 1, 4, 1))

ecdf_subplot_chunky(values(suitability$vectors), 
                    vector_suitability_at_pos)

ecdf_subplot_chunky(values(suitability$birds), 
                    bird_suitability_at_pos,
                    title="Birds")

ecdf_subplot_chunky(values(suitability$feral_pigs), 
                    feral_pig_suitability_at_pos, 
                    TRUE,
                    title="Feral pigs",
                    xlab="Unscaled abundance")
axis(1, at=log1p(0:3), labels=0:3)
# log1p becomes increasingly linear the closer we get to 0 you duffer ... do some more thinking
ecdf_subplot_chunky(values(suitability$domestic_pigs), 
                    domestic_pig_suitability_at_pos, 
                    TRUE, 1e-3,
                    title="Domestic pigs",
                    xlab="Unscaled abundance of premises")
axis(1, at=log(seq(0, 0.004, length.out=5) + 0.001), labels=seq(0, 0.004, length.out=5))
dev.off()}


# vectors_hist <- hist(vectors.feb.22, breaks=50, plot=FALSE)
# vectors_hist$counts <- cumsum(vectors_hist$counts)
# maxvectors <- max(vectors_hist$counts)
# vectors_hist$counts <- vectors_hist$counts / max(vectors_hist$counts) #* 1.1
# 
# birds_hist <- hist(birds.feb.22, breaks=50, plot=FALSE) # more detail in xax?
# birds_hist$counts <- cumsum(birds_hist$counts)
# maxbirds <- max(birds_hist$counts)
# birds_hist$counts <- birds_hist$counts / max(birds_hist$counts) #* 1.1
# 
# dpig_hist <- hist(dom_pig_density, breaks=50, plot=FALSE)
# dpig_hist$counts <- cumsum(dpig_hist$counts)
# maxdpig <- max(dpig_hist$counts)
# dpig_hist$counts <- dpig_hist$counts / max(dpig_hist$counts) #* 1.1
# 
# fpig_hist <- hist(abund_feral_pig, breaks=50, plot=FALSE)
# fpig_hist$counts <- cumsum(fpig_hist$counts)
# maxfpig <- max(fpig_hist$counts)
# fpig_hist$counts <- fpig_hist$counts / max(fpig_hist$counts) #* 1.1
# 
# {png('figures/ecdf_hists_cumul.png',
#      width = 1100*reso/72, 
#      height = 700*reso/72, 
#      res = reso,
#      bg = 'white')
# par(mfrow = c(2,2),
#     mar = c(4.1,6.1,4.1,6.1),
#     bty="n")
# plot(vectors_hist, main="", xlab="Unscaled vector surface", #yaxt = "n",
#      cex.lab=cexlab, ylab="Empirical cumulative distribution:\n vector surface")
# axis(4, at=c(0,0.25,0.5,0.75,1), col=ecdfcol, col.axis=ecdfcol)
# #axis(2,  at=seq(0, 4000, length.out = 5) / maxvectors, labels=seq(0, 4000, length.out = 5))
# mtext("Empirical cumulative distribution: JEV detections", 4, col=ecdfcol, cex=cexlab*0.9, line=3)
# plot(vector_ecdf, add=TRUE, pch=1, col=ecdfcol)
# 
# plot(birds_hist, main="", xlab="Unscaled bird surface", #yaxt = "n",
#      cex.lab=cexlab, ylab="Empirical cumulative distribution:\n bird surface")
# axis(4, at=c(0,0.25,0.5,0.75,1), col=ecdfcol, col.axis=ecdfcol)
# #axis(2,  at=seq(0, 15000, length.out = 4) / maxbirds, labels=seq(0, 15000, length.out = 4))
# mtext("Empirical cumulative distribution: JEV detections", 4, col=ecdfcol, cex=cexlab*0.9, line=3)
# plot(bird_ecdf, add=TRUE, pch=1, col=ecdfcol)
# 
# plot(dpig_hist, main="", xlab="Unscaled domestic pig surface", #yaxt = "n",
#      cex.lab=cexlab, ylab="Empirical cumulative distribution:\n domestic pig surface")
# axis(4, at=c(0,0.25,0.5,0.75,1), col=ecdfcol, col.axis=ecdfcol)
# #axis(2,  at=seq(0, 8e6, length.out = 5) / maxdpig, labels=seq(0, 8e6, length.out = 5))
# mtext("Empirical cumulative density", 4, col=ecdfcol, cex=cexlab*0.9, line=3)
# plot(domestic_pig_ecdf, add=TRUE, pch=1, col=ecdfcol)
# 
# plot(fpig_hist, main="", xlab="Unscaled feral pig surface", #yaxt = "n",
#      cex.lab=cexlab, ylab="Empirical cumulative distribution:\n feral pig surface")
# axis(4, at=c(0,0.25,0.5,0.75,1), col=ecdfcol, col.axis=ecdfcol)
# #axis(2,  at=seq(0, 6e6, length.out = 5) / maxfpig, labels=seq(0, 6e6, length.out = 5))
# mtext("Empirical cumulative density", 4, col=ecdfcol, cex=cexlab*0.9, line=3)
# plot(feral_pig_ecdf, add=TRUE, pch=1, col=ecdfcol)
# 
# 
# par(mfrow=c(1,1), new=TRUE, mar=c(0.5,0.5,0.5,0.5))
# plot(0, type="n", axes=FALSE, bty="n")
# subfigure_label(par("usr"), 0.55, 0.45, "(d)", 1.4)
# subfigure_label(par("usr"), 0.04, 0.45, "(c)", 1.4)
# subfigure_label(par("usr"), 0.55, 0.95, "(b)", 1.4)
# subfigure_label(par("usr"), 0.04, 0.95, "(a)", 1.4)
# dev.off()}


############################################################################
# PIG DENSITIES FOR SUPP
purps = colorRampPalette(brewer.pal(9,"Purples"))(100)
# harsher palette of purples:
extreme_purps = colorRampPalette(purps[c(seq(1,30,5),30:100)])(100)
state.lwd <- 0.4

blues <- colorRampPalette(brewer.pal(9,"Blues"))(100)
greens <- colorRampPalette(brewer.pal(9,"Greens"))(100)

{png("figures/pig_density.png",
     width = 2000, 
     height = 3300, 
     res = 300,
     type = 'cairo')
  
  par(mfrow = c(2,1), mar=c(1.1,1.1,1.1,2.1))
  
  plot(suitability_aus$domestic_pigs*100, axes=FALSE, box=FALSE, col=extreme_purps, 
       legend.args = list(text='Pig premises\nper square km', cex=1, line=0.9))
  plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1.2)
  
  
  plot(suitability_aus$feral_pigs, axes=FALSE, box=FALSE, col=extreme_purps, 
       legend.args = list(text='Pigs per\nsquare km', cex=1, line=0.9))
  plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1.2)
  
  dev.off()}

{png("figures/birds_vectors.png",
     width = 2000, 
     height = 3300, 
     res = 300,
     type = 'cairo')
  
  par(mfrow = c(2,1), mar=c(1.1,1.1,1.1,2.1))
  
  plot(suitability_aus$vectors, axes=FALSE, box=FALSE, col=blues, 
       legend.args = list(text='Habitat suitability\n index', cex=1, line=0.9))
  plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1.2)
  
  
  plot(suitability_aus$birds, axes=FALSE, box=FALSE, col=greens, 
       legend.args = list(text='Habitat suitability\n index', cex=1, line=0.9))
  plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1.2)
  
  dev.off()}

{png("figures/birds_vectors.png",
     width = 3300, 
     height = 3000, 
     res = 300,
     type = 'cairo')
  
  par(mfrow = c(2,2), mar=c(1.1,1.1,1.1,2.1), oma = rep(2,4))
  
  plot(vectors.feb.22, axes=FALSE, box=FALSE, col=blues, 
       legend.args = list(text='Habitat suitability\n index', cex=1, line=0.9))
  plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1.2)
  
  plot(birds.feb.22, axes=FALSE, box=FALSE, col=greens, 
       legend.args = list(text='Habitat suitability\n index', cex=1, line=0.9))
  plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1.2)
  
  plot(vector_suitability_ecdf, axes=FALSE, box=FALSE, col=blues, 
       legend.args = list(text='Habitat suitability\n index', cex=1, line=0.9))
  plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", 1.2)
  
  plot(bird_suitability_ecdf, axes=FALSE, box=FALSE, col=greens, 
       legend.args = list(text='Habitat suitability\n index', cex=1, line=0.9))
  plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", 1.2)
  
  mtext(outer = TRUE, adj = 0.22, "Vectors", font = 2)
  mtext(outer = TRUE, adj = 0.75, "Birds", font = 2)
  mtext(side = 2, outer = TRUE, adj = 0.25, "Scaled", font = 2)
  mtext(side = 2, outer = TRUE, adj = 0.8, "Unscaled", font = 2)
  
  dev.off()}


###############################################################################
# all <- rbeta(1e4, 5, 20)
# jev <- rbeta(100, 30, 40)

ecdf_subplot_continuous <- function(all, jev, logx=FALSE, epsilon=0){
  if (logx){
    all <- log10(all + epsilon)
    jev <- log10(jev + epsilon)
  }
  
  all_dens <- density(all, na.rm=TRUE)
  jev_dens <- density(jev, na.rm=TRUE)
  
  max_density <- max(all_dens$y, jev_dens$y)
  xlim <- c(min(0, min(all, na.rm=TRUE)),
            max(1, max(all, na.rm=TRUE)))
  
  plot(all_dens,
       type = "n",
       xlim = xlim,
       ylim = c(0, max_density),
       axes = FALSE,
       xlab = "Unscaled suitability",
       ylab = "",
       main = "Vectors")
  
  polygon(x = all_dens$x,
          y = all_dens$y,
          border = grey(0.8),
          col = grey(0.9))
  
  polygon(x = jev_dens$x,
          y = jev_dens$y,
          border = ecdfcol,
          col = scales::alpha(ecdfcol, 0.5))
  
  ecdf_x <- seq(xlim[1], xlim[2], length.out = 1000)
  ecdf_y <- ecdf(jev)(ecdf_x) * max_density
  
  lines(ecdf_x,
        ecdf_y,
        col = ecdfcol,
        lwd = 2)
  
  axis(1)
}

# host control map?
critters <- list(mozz=1, bird=2, fpig=3, dpig=4)
bvfd_df <- data.frame(bird=getValues(bird_suitability_ecdf),
                      mozz=getValues(vector_suitability_ecdf),
                      dpig=getValues(domestic_suitability_ecdf),
                      fpig=getValues(feral_suitability_ecdf)) %>%
  drop_na() %>%
  mutate(which_mozz_bird = ifelse(bird > mozz, critters[["mozz"]], critters[["bird"]]),
         which_bird_fpig = ifelse(bird < fpig, critters[["fpig"]], critters[["bird"]]),
         bird_fpig = pmax(bird, fpig),
         which_bird_fpig_dpig = ifelse(which_bird_fpig == critters[["bird"]], # we want the maximum of the three
                                       ifelse(bird < dpig, critters[["dpig"]], critters[["bird"]]), # test for bird
                                       ifelse(fpig < dpig, critters[["dpig"]], critters[["fpig"]])), # test for fpig
         bird_fpig_dpig = pmax(bird_fpig, dpig),
         which_mozz_bird_fpig = ifelse(which_bird_fpig == critters[["bird"]], # we want the minimum of the maximum hosts
                                       ifelse(bird > mozz, critters[["mozz"]], critters[["bird"]]), # test for bird
                                       ifelse(fpig > mozz, critters[["mozz"]], critters[["fpig"]])), # test for fpig
         mozz_bird_fpig = pmin(bird_fpig, mozz),
         which_mozz_bird_fpig_dpig = ifelse(which_bird_fpig_dpig == critters[["bird"]],
                                            ifelse(bird > mozz, critters[["mozz"]], critters[["bird"]]), # test for bird
                                            ifelse(which_bird_fpig_dpig == critters[["fpig"]],
                                                   ifelse(fpig > mozz, critters[["mozz"]], critters[["fpig"]]),
                                                   ifelse(dpig > mozz, critters[["mozz"]], critters[["dpig"]]))),
         mozz_bird_fpig_dpig = pmin(bird_fpig_dpig, mozz)
         ) 


tmp <- bird_suitability_ecdf
values(tmp)[!is.na(values(tmp))] <- bvfd_df$which_mozz_bird
plot(tmp)


tmp <- trim(bird_suitability_ecdf)
vectors_birds <- tmp
values(vectors_birds)[!is.na(values(vectors_birds))] <- bvfd_df$which_mozz_bird

vectors_birds_fpigs <- tmp
values(vectors_birds_fpigs)[!is.na(values(vectors_birds_fpigs))] <- bvfd_df$which_mozz_bird_fpig

vectors_birds_fpigs_dpigs <- tmp
values(vectors_birds_fpigs_dpigs)[!is.na(values(vectors_birds_fpigs_dpigs))] <- bvfd_df$which_mozz_bird_fpig_dpig

birds_fpigs_dpigs <- tmp
values(birds_fpigs_dpigs)[!is.na(values(birds_fpigs_dpigs))] <- bvfd_df$which_bird_fpig_dpig


brewer.pal(9,"RdPu")
# "#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5" "#08519C" "#08306B"
# "#F7FCF5" "#E5F5E0" "#C7E9C0" "#A1D99B" "#74C476" "#41AB5D" "#238B45" "#006D2C" "#00441B"
# "#FCFBFD" "#EFEDF5" "#DADAEB" "#BCBDDC" "#9E9AC8" "#807DBA" "#6A51A3" "#54278F" "#3F007D"
# "#FFF7F3" "#FDE0DD" "#FCC5C0" "#FA9FB5" "#F768A1" "#DD3497" "#AE017E" "#7A0177" "#49006A"

# re-label, re-calculate - should increasingly resemble mosquito surface

{png('figures/potential_niche_unpack.png',
     width = 1000*reso/72, 
     height = 1000*reso/72, 
     res = reso,
     bg = 'white')
  
my_blue <- "#e0effd" # this can be lighter :)
#my_green <- "#74C476"
my_green <- "#3dc440"
#my_purp <- "#9E9AC8"
my_purp <- "#917ef6"
my_pink <- "#F768A1"

main.cex <- 1.8
legend.cex <- 1.6
  
par(mfrow=c(2,2), mar=c(0.1,0.1,4,0.1), bty="n", oma=c(0,0,0,0), xpd=NA)
plot(vectors_birds, col=c(my_blue, my_green), 
     legend=FALSE, bty="n", axes=FALSE,
     main="Mosquitoes AND birds", cex.main=main.cex)
legend("bottomleft", fill=c(my_blue, my_green), legend=c("Mosquito", "Bird"), cex=legend.cex)
par(xpd=NA)
subfigure_label(par()$usr, 0.1,1.05,"(a)",main.cex)

plot(vectors_birds_fpigs, col=c(my_blue, my_green, my_purp), 
     legend=FALSE, bty="n", axes=FALSE,
     main="Mosquitoes AND (birds OR feral pigs)", cex.main=main.cex)
legend("bottomleft", fill=c(my_blue, my_green, my_purp), 
       legend=c("Mosquito", "Bird", "Feral pig"), cex=legend.cex)
par(xpd=NA)
subfigure_label(par()$usr, 0.1,1.05,"(b)",main.cex)

plot(vectors_birds_fpigs_dpigs, col=c(my_blue, my_green, my_purp, my_pink),
     legend=FALSE, bty="n", axes=FALSE,
     main="Mosquitoes AND\n (birds OR feral pigs OR domestic pigs)", cex.main=main.cex)
legend("bottomleft", fill=c(my_blue, my_green, my_purp, my_pink), 
       legend=c("Mosquito", "Bird", "Feral pig", "Domestic pig"), cex=legend.cex)
par(xpd=NA)
subfigure_label(par()$usr, 0.1,1.05,"(c)",main.cex)

plot(birds_fpigs_dpigs, col=c(my_green, my_purp, my_pink),
     legend=FALSE, bty="n", axes=FALSE,
     main="Birds OR feral pigs OR domestic pigs", cex.main=main.cex)
legend("bottomleft", fill=c(my_green, my_purp, my_pink), 
       legend=c("Bird", "Feral pig", "Domestic pig"), cex=legend.cex)
par(xpd=NA)
subfigure_label(par()$usr, 0.1,1.05,"(d)",main.cex)
#mtext("Bonus plot: where does the value in the potential niche fig come from?", side=3, outer=TRUE, 
#      font=3, cex=1.5)
dev.off()}

# {png('figures/ecdfs_cont.png',
#      width = 1100*reso/72,
#      height = 700*reso/72,
#      res = reso,
#      bg = 'white')
#   
# par(mfrow=c(2,2))
# 
# par(mar = c(5, 1, 4, 1))
# 
# ecdf_subplot_continuous(values(suitability$vectors), vector_suitability_at_pos)
# ecdf_subplot_continuous(values(suitability$birds), bird_suitability_at_pos)
# ecdf_subplot_continuous(values(suitability$feral_pigs), feral_pig_suitability_at_pos, TRUE, 1e-1)
# ecdf_subplot_continuous(values(suitability$domestic_pigs), domestic_pig_suitability_at_pos, TRUE, 1e-15)
# 
# dev.off()}

