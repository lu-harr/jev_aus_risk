# 2022 model!

source("code/init.R")
source("model_2022/code/pathogen_functions.R")
source("code/plotting_functions.R")
# fit transmission intensity model, predict, and do power calculations

# jumping to 11km res
AGG_FACTOR <- 10
aus_agg_mask <- aggregate(aus_mask,
                          fact = AGG_FACTOR)
set.seed(834903)

# load blurred human pop layer
hpop_blur <- raster('data/hpop_blur_aus_0.2.tif')

# load human, mozzie and feral pig, and domestic pig infection data
all_pres <- read_csv('model_2022/data/pathogen_presence_national_cut.csv')
ftable(all_pres$host_species)

all_humans <- all_pres %>%
  dplyr::filter(host_species == 'human') %>%
  mutate(virus_detection = "JEV positive")

# load feral pig testing data
all_ferals <- read_csv('model_2022/data/all_feral_pigs_clean.csv')

# load domestic pig testing data
all_dompigs <- read_csv('model_2022/data/all_domestic_pigs_clean.csv') %>%
  mutate(virus_detection = case_when(
    virus_detection == "Negative" ~ "JEV negative",
    virus_detection == "Positive" ~ "JEV positive"
  ))

# load mozzie testing data
all_mozzies <- read.csv('model_2022/data/mosquito_detections_all_w_qld_sa.csv') %>%
  dplyr::filter(virus_detection %in% c('JEV positive',
                                       'JEV negative'))




# restrict the model-fitting region to only those pixels within some distance
# of presence points - clip all rasters used in model fitting
mask_clipped <- clip_to_distance(aus_agg_mask, as.matrix(all_pres[,c('longitude', 'latitude')]),
                                 max_distance = 1300000)
hpop_blur_clipped <- mask(hpop_blur, mask_clipped)

# find mozzie points within clipped mask
vals <- raster::extract(mask_clipped,
                        all_mozzies[, c('longitude',
                                        'latitude')])

all_mozzies_clipped <- all_mozzies[!(is.na(vals)),]

# organise all of the data: need surveillance effort surfaces, etc.
dat <- prepare_pathogen_inputs_national(all_humans,
                                        all_mozzies_clipped,
                                        all_dompigs,
                                        all_ferals,
                                        mask_clipped,
                                        hpop_blur_clipped,
                                        mozz_bandwidth_scalar = 0.8,
                                        dpig_bandwidth_scalar = 0.8,
                                        feral_bandwidth_scalar = 0.8)

# fit the model!!!
m <- fit_pathogen_spatial_model(dat$df, dat$coordinates)

saveRDS(m, "model_2022/output/pathogen_model_object.rds")
m <- readRDS("model_2022/output/pathogen_model_object.rds")

################################

# pull out intercepts for the host species ids
int <- m$model$coefficients["(Intercept)"]
contrasts_idx <- grepl("host_species_id", names(m$model$coefficients))
contrasts <- c(0, m$model$coefficients[contrasts_idx])
intercepts <- int + contrasts

# for each species, pull out the appropriate intercept (will need to check
# against earlier code), specify the required intensity level, and use them to
# construct the expected rate of detection

# we need to specify this - can read it off the plot, or better to get it as a threshold
# value for the overall pathogen threshold (e.g. the 50% inclusion level)
# calculate threshold on pred_inc returned with model object

m_pred_inc <- m$pred_inc

r <- m_pred_inc
r[is.na(r)] <- cellStats(r, min)
r <- mask(r, hpop_blur)

writeRaster(r,
            'model_2022/output/continuous_pathogen_national_cut.tif',
            overwrite=TRUE)

# find threshold: median of intensities at locations of detection
pos_locs = dat$df[dat$df$cases > 0, c("x", "y")]
vals <- raster::extract(m_pred_inc, pos_locs)
threshold <- quantile(vals, 1 - 0.5)

species_id_lookup <- c(human = 1, mossie = 2, dpig = 3, fpig = 4)

# make mozzie density plot: this step is also in prepare_pathogen_inputs, but want unclipped
mean_sample <- round(mean(all_mozzies$n_individuals, na.rm=TRUE), 0)

all_mozzies_expanded <- all_mozzies %>%
  mutate(n_individuals=case_when(
    is.na(n_individuals) ~ mean_sample,
    TRUE ~ round(n_individuals, 0)
  )) %>% 
  uncount(n_individuals)

tmp <- rasterize(all_mozzies_expanded[,c('longitude', 'latitude')], 
                 aus_agg_mask, background=0, fun='count')
gf <- focalWeight(aus_agg_mask, 0.8, "Gauss")
nat_mozz_density <- raster::focal(tmp, gf, pad=TRUE, na.rm=TRUE)
nat_mozz_density = mask(nat_mozz_density, aus_agg_mask)

#human_expected_detections <- exp(log(rel_intensity_cutoff) + intercepts[species_id_lookup["human"]] + log(dat$densities$human))
human_expected_detections <- exp(log(threshold) + intercepts[species_id_lookup["human"]] + log(hpop_blur))
#mossie_expected_detections <- exp(log(threshold) + intercepts[species_id_lookup["mossie"]] + log(dat$densities$mosquito_sampling))
mossie_expected_detections <- exp(log(threshold) + intercepts[species_id_lookup["mossie"]] + log(nat_mozz_density))

dpig_expected_detections <- exp(log(threshold) + intercepts[species_id_lookup["dpig"]] + log(dat$densities$domestic_pig_sampling))
dpig_expected_detections <- extend_raster(dpig_expected_detections, hpop_blur)

fpig_expected_detections <- exp(log(threshold) + intercepts[species_id_lookup["fpig"]] + log(dat$densities$feral_pig_sampling))
fpig_expected_detections <- extend_raster(fpig_expected_detections, hpop_blur)


# probabilities of observing one or more detections in each host species at each
# location
detection_prob <- function(expected_detections) {
  1 - dpois(0, expected_detections)
}

human_detection_probability <- calc(human_expected_detections, detection_prob)
mossie_detection_probability <- calc(mossie_expected_detections, detection_prob)
dpig_detection_probability <- calc(dpig_expected_detections, detection_prob)
fpig_detection_probability <- calc(fpig_expected_detections, detection_prob)


writeRaster(human_detection_probability,
            'model_2022/output/human_detection_probability.tif',
            overwrite = TRUE)

writeRaster(mossie_detection_probability,
            'model_2022/output/vector_detection_probability.tif',
            overwrite = TRUE)

writeRaster(dpig_detection_probability,
            'model_2022/output/domestic_pig_detection_probability.tif',
            overwrite = TRUE)

writeRaster(fpig_detection_probability,
            'model_2022/output/feral_pig_detection_probability.tif',
            overwrite = TRUE)


# combine these to get the overall probability of detection per pixel
overall_detection_probability <- 1 - (1 - human_detection_probability) * (1 - mossie_detection_probability) * (1 - dpig_detection_probability) * (1 - fpig_detection_probability)

# par(mfrow = c(1, 1))
plot(overall_detection_probability)

# save for later
writeRaster(overall_detection_probability,
            'output/national_detection_probability_cut.tif',
            overwrite = TRUE)


#######################################################################
# Play around with thresholds a little more ...


all_detections <- all_pres %>% 
  rename(x='longitude', y='latitude')

# set up thresholds for plot colouring
threshes <- quantile(raster::extract(m$pred_rel_inc, all_detections[,c("x","y")]), na.rm=TRUE, 
                     probs=c(0, 0.01, 0.05,0.25,0.5))
threshes_full <- c(minValue(m$pred_rel_inc), threshes, maxValue(m$pred_rel_inc)) # with min/max for plot

# smooth and re-mask intensity surface
pred_rel_inc_2dp_ish <- m$pred_rel_inc %>%
  disaggregate(fact=10) %>%
  focal(w=matrix(1, 11, 11), mean, na.rm=TRUE) %>% #%>% # ignore NAs or we lose boundary pixels
  mask(aus_mask)

###############################################################################
# PATHOGEN MODEL INPUTS FIGURE (detections + surveillance effort surfaces)

# extend clipped matrices ..
domestic_pig_sampling <- dat$densities$domestic_pig_sampling
feral_pig_sampling <- dat$densities$feral_pig_sampling
domestic_pig_sampling[is.na(domestic_pig_sampling) & 
                        !is.na(nat_mozz_density)] = minValue(domestic_pig_sampling)
feral_pig_sampling[is.na(feral_pig_sampling) & 
                        !is.na(nat_mozz_density)] = minValue(feral_pig_sampling)


purps = colorRampPalette(brewer.pal(9,"Purples"))(100)

# harsher palette of purples for surveillance effort:
extreme_purps = colorRampPalette(purps[c(seq(1,30,5),30:100)])(100)
reso <- 300

surveillance_effort_points = alpha("blue", 1)
surv_cex = 0.6
surv_pch = 16

{png("figures/pathogen_inputs.png",
     width = 1000*reso/72, 
     height = 1100*reso/72, 
     res = reso)
  #type = 'cairo')
  
  par(mfrow = c(2, 2),
      mar = c(3.5, 0.5, 1, 0.5), bty="n",
      oma=c(2,0,0,0))
  
  # humans
  plot(sqrt(hpop_blur/100), col=extreme_purps,
       box = FALSE, axes=FALSE, legend=FALSE)
  plot(sqrt(hpop_blur/100), col=extreme_purps, legend.width=1.4, horizontal=TRUE,
       legend.args = list(text='Human population\nper square km', cex=1.2, line=0.9),
       axis.args = list(at=sqrt(c(10,100,400,900)), label=c(10,100,400,900)),
       legend.only=TRUE)
  plot(aus_simple, add=TRUE, lwd=1.2, border='grey20')
  points(dat$human_points, cex=surv_cex, col="red", pch=surv_pch)
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1.6)
  
  # mosquitoes
  plot(nat_mozz_density/100, box = FALSE, axes=FALSE, col=extreme_purps, 
       legend.width=1.4,horizontal=TRUE,
       legend.args = list(text='Mosquito catches\nscreened\nper square km', 
                          cex=1.2, line=0.9))
  plot(aus_simple, lwd=1.2, border='grey20', add=TRUE)
  points(
    all_mozzies[,c("longitude","latitude")],
    cex = surv_cex, pch=surv_pch, lwd=0.4, col=surveillance_effort_points)
  points(dat$mozz_points[which(dat$mozz_points$cases==1),c("x","y")], 
         col="red", pch=surv_pch, cex=surv_cex)
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1.6)
  
  # domestic pigs
  #plot(aus_simple, lwd=0.01, border='grey20')
  plot(domestic_pig_sampling/100, box = FALSE, axes=FALSE, 
       col=extreme_purps, horizontal=TRUE,
       legend.args = list(text='Sampled domestic pigs\nper square km', cex=1.2, line=0.9), 
       #add=TRUE, 
       legend.width=1.4)
  plot(aus_simple, lwd=1.2, border='grey20', add=TRUE)
  points(
    dat$dpig_points[,c("x","y")],
    lwd=0.4, col=surveillance_effort_points,pch=surv_pch, cex=surv_cex)
  points(dat$dpig_points[which(dat$dpig_points$cases==1),c("x","y")], 
         col="red", pch=surv_pch, cex=surv_cex)
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", 1.6)
  
  # feral pigs
  #plot(aus_simple, lwd=0.01, border='grey20')
  plot(feral_pig_sampling/100, box = FALSE, axes=FALSE, 
       col=extreme_purps, horizontal=TRUE,
       legend.args = list(text='Sampled feral pigs\nper square km', cex=1.2, line=0.9), 
       #add=TRUE, 
       legend.width=1.4)
  plot(aus_simple, lwd=1.2, border='grey20', add=TRUE)
  points(
    dat$feral_points[,c("x","y")],
    pch=surv_pch, cex=surv_cex, lwd=0.4, col=surveillance_effort_points)
  points(dat$feral_points[which(dat$feral_points$cases==1),c("x","y")], 
         col="red", pch=surv_pch, cex=surv_cex)
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", 1.6)
  
  dev.off()}


###############################################################################
# TRANSMISSION INTENSITY IN RESERVOIR FIGURE

pal <- brewer.pal(9,"Reds")[2:9]
pal_cont <- colorRampPalette(pal)(100)
reso <- 300

darkred = rgb(colorRamp(c(pal[8], "black"))(0.7)/100)
purps = colorRampPalette(brewer.pal(9,"Purples"))(100)
extreme_purps = colorRampPalette(purps[c(seq(1,30,5),30:100)])(100)

{png("figures/pathogen_intensity_national.png",
     width = 1000*reso/72, 
     height = 500*reso/72, 
     res = reso)
  
  par(mfrow = c(1,2),
      mar = c(1, 0.1, 1, 1), bty="n",
      oma=c(0,0,0,0))
  
  plot(aus_simple, lwd=0.0001, border='grey20', col="#f7f7f7")
  plot(pred_rel_inc_2dp_ish, breaks=threshes_full, col=c("#f7f7f7", pal[(c(1,2,3,4,5))]), legend=FALSE,
       xaxt="n", yaxt="n", legend.mar=0, add=TRUE)
  plot(rasterToContour(m$pred_rel_inc, levels=threshes), add=TRUE, col=darkred)
  plot(aus_simple, lwd=1, border='grey20', add=TRUE)
  text(c(135, 141, 144, 148.5, 144.7), 
       c(-19, -21, -23.5,-27.9, -27.9), 
       c("100%", "99%", "95%", "75%", "50%"), cex=1)
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1.6)
  
  # plot overall detection probability
  detect_breaks <- c(0, 0.01, 0.1, 0.25, 0.5)
  
  plot(aus_simple, lwd=0.01, border='grey20')
  plot(sqrt(overall_detection_probability), xaxt="n", yaxt="n", 
       col=extreme_purps, axis.args = list(at = sqrt((detect_breaks)), labels = detect_breaks), 
       add=TRUE, legend.width=1, horizontal=TRUE)
  plot(aus_simple, lwd=1.2, border='grey20', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1.6)
  
  dev.off()}

# Here's a check of the statement about power in the discussion
tmp = overall_detection_probability
tmp[tmp < 1/4] = NA
plot(tmp, col="red")



