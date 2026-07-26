# A copy of `ecdfs.R` with February 2022 bird and vector surfaces substituted
# for median predictions from December-March 2022

source("code/init.R")
source("code/plotting_functions.R")
set.seed(834903)

# read all of the separate species maps in:
birds.med <- raster('data/enms/median_surfaces/egrons_q0.5_DtM.tif')
# updating post egron/cx MS submission
# vectors.med <- raster('data/enms/median_surfaces/cx_q0.5_DtM.tif')
vectors.med <- raster("data/enms/median_surfaces/cx_oi_q0.5_Dec21tMar22.tif")

birds.feb.22 <- raster('data/enms/brt_egrons_halfwt_tc4__brt_feb_22.tif')
# vectors.feb.22 <- raster('data/enms/cx_oi_bg1x_noAccess__mxnt_feb_22.tif')
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
  vectors.med,
  birds.med,
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
  
  cont_outputs[[i_hosts]] <- combined_suitability_ecdf
  
}  



# MAKE THE FIGURE
existing_pink <- raster("output/continuous suit vectors and avian.tif")
pal = colorRampPalette(c("blue", "lightyellow", "red"))
blues <- colorRampPalette(brewer.pal(9,"Blues"))(100)
greens <- colorRampPalette(brewer.pal(9,"Greens"))(100)

{png("figures/dec_mar_median.png",
    height=2500, width=2000, pointsize=35)
par(mfrow=c(3,2), mar=c(1.1,2.1,1.1,6.1), oma=c(0,2,0,2))

plot(suitability_aus$vectors, axes=FALSE, box=FALSE, col=blues, 
     legend.args = list(text='Habitat suitability\n index', cex=1, line=0.9))
plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
mtext("Vectors", font=2)
subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1.2)
mtext("Median prediction:\nDecember 2021-March 2022", side=2, font=2, line=1)

plot(suitability_aus$birds, axes=FALSE, box=FALSE, col=greens, 
     legend.args = list(text='Habitat suitability\n index', cex=1, line=0.9))
plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1.2)
mtext("Ardeid birds", font=2)

leg_end = -0.25
leg_pos = 160
plot(vectors.feb.22 - vectors.med, breaks=seq(leg_end, leg_end*-1, length.out=100), 
     col=pal(100), box=FALSE, axes = FALSE,
     axis.args = list(at = seq(leg_end, leg_end*-1, length.out=5)))
plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
mtext("Feb 22\nhigher", at=leg_pos, line=-6, cex=0.9)
mtext("Dec 21-Mar 22\nmedian\nhigher", at=leg_pos, line=-2.5, side = 1, cex=0.9)
subfigure_label(par("usr"), 0.1, 0.9, "(c)", 1.2)
mtext("Difference between February 2022\nand December-March median", side=2, font=2, line=1)

leg_end = -0.215
plot(birds.feb.22 - birds.med, breaks=seq(leg_end, leg_end*-1, length.out=100), 
     col=pal(100), box = FALSE, axes = FALSE,
     axis.args = list(at = seq(-0.2, 0.2, length.out=5)))
plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
mtext("Feb 22\nhigher", at=leg_pos, line=-6, cex=0.9)
mtext("Dec 21-Mar 22\nmedian\nhigher", at=leg_pos, line=-2.5, side = 1, cex=0.9)
subfigure_label(par("usr"), 0.1, 0.9, "(d)", 1.2)

par(xpd=NA)
plot(cont_outputs[[1]], col=pn_cols, box=FALSE, axes=FALSE, 
     legend.args=list(text="JEV\ntransmission\nsuitability", cex=0.85, line=1.5),
     cex.main=1.4,
     main="JEV transmission suitability\n(December 2021-March 2022 median birds and vectors)")
plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
subfigure_label(par("usr"), 0.1, 0.9, "(e)", 1.2)

leg_end = -0.46
par(xpd=NA)
plot(existing_pink - cont_outputs[[1]], breaks=seq(leg_end, leg_end*-1, length.out=100), 
     col=pal(100), box = FALSE, axes = FALSE,
     axis.args = list(at = seq(-0.4, 0.4, length.out=5)),
     cex.main=1.4,
     main = "Comparison between transmission suitability outputs\n(birds and vectors only)")
plot(st_geometry(aus_simple), lwd=state.lwd, border='grey30', add=TRUE)
mtext("Feb 22\nhigher", at=leg_pos, line=-6, cex=0.9)
mtext("Dec 21-Mar 22\nmedian\nhigher", at=leg_pos, line=-2.5, side = 1, cex=0.9)
subfigure_label(par("usr"), 0.1, 0.9, "(f)", 1.2)

par(xpd=NA)
abline(h=-3, lty=2, lwd=2)
dev.off()}



###############################################################################

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

