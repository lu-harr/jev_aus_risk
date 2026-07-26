# 2023 model! (i.e. Phase 2: chickens and diagnostic weights)
setwd("~/Desktop/jev/pathogen_paper")

source("code/init.R")
source("code/pathogen_functions.R")
source("code/plotting_functions.R")
source("code/diagnostic_weights_table.R")
set.seed(834903)

states <- st_read("~/Desktop/jev/data/admin/STE_2021_AUST_SHP_GDA2020/STE_2021_AUST_GDA2020.shp")

# let's jump to 10km res
AGG_FACTOR <- 10
aus_agg_mask <- aggregate(aus_mask, 
                          fact = AGG_FACTOR)

# load blurred human pop layer
hpop_blur <- raster('data/hpop_blur_aus_0.2.tif')

# all presences
all_pres <- read_csv('data/secure/phase_2/clean/pathogen/pathogen_presence_national_cut.csv')
# somehow there are non-distinct combinations of lon, lat, date in mozzies?
# also same for domestic pigs ...
# I have once again made more work for myself :(

# tmp = all_pres %>%
# #  mutate(date = ifelse(date < "2020-01-01 UTC", as.POSIXct("2022-04-26 UTC"), as.POSIXct(date))) %>%
#   filter(date < "2022-01-01 UTC" & date > "2021-01-01 UTC")
# # Victorian domestic pig
# 
# tmp = dplyr::select(filter(all_pres, host_type == "human" | host_type == "domestic_pig"), date) %>%
#   as.vector() %>%
#   unlist() %>%
#   as.POSIXct()
# hist(tmp, breaks="month")


# for first paragraph of paper
ftable(all_pres$host_type)


# all humans (presences)
all_humans <- all_pres %>% 
  dplyr::filter(host_type == 'human') %>%
  mutate(virus_detection = "JEV positive")

# tmp <- st_as_sf(all_humans,
#                 coords=c("longitude", "latitude"),
#                 crs=crs(states)) %>%
#   st_join(states)
# 
# all_humans %>%
#   mutate(state = tmp$STE_NAME21) %>%
#   group_by(state, virus_detection) %>%
#   summarise(event = n_distinct(longitude, latitude, date))

# load feral pig testing data
all_ferals <- read_csv('data/secure/phase_2/clean/pathogen/aus_feral_pigs_combined.csv') %>%
  drop_na(c(longitude, latitude)) %>%
  filter(virus_detection %in% c("JEV positive", "JEV negative"))
# LH removing inconclusives Dec 25
# could add a Pending row here?
# diagnostic_weights_table[diagnostic_weights_table$host_type == "feral_pig" &
#                            diagnostic_weights_table$diagnostic_type == "Serology" &
#                            diagnostic_weights_table$diagnostic_endpoint == "ELISA",]

# tmp <- st_as_sf(all_ferals,
#                 coords=c("longitude", "latitude"),
#                 crs=crs(states)) %>%
#   st_join(states)
# 
# all_ferals %>%
#   mutate(state = tmp$STE_NAME21) %>%
#   group_by(state, virus_detection) %>%
#   summarise(event = n_distinct(longitude, latitude, date))

# load domestic pig testing data 
all_dompigs <- read_csv('data/secure/phase_2/clean/pathogen/all_dompigs_combined.csv') %>% 
  dplyr::select(
    -c(lga,
       age_class,
       #sample_age,
       sample_unit)
  ) %>%
  drop_na(longitude:latitude) %>%
  mutate(exposure_weight = if_else(
           is.na(exposure_weight), 1, exposure_weight), # if no info on exposure weight, assume 1
         diagnostic_detail = case_when(
           is.na(diagnostic_detail) ~ "PCR", 
           diagnostic_detail == "IgG ELISA JKM" ~ "ELISA", # check this with DHD
           TRUE ~ diagnostic_detail), # should only be one thing here (300)... be sure to check !
         
   ) %>%
  filter(
    virus_detection %in% c("JEV positive", "JEV negative")
  )

# adding line to scoop up some dompigs that end up with diagnostic weight NA ...
# these have now disappeared ?
# diagnostic_weights_table[diagnostic_weights_table$host_type == "domestic_pig" &
#                            diagnostic_weights_table$diagnostic_type == "PCR",]
# all_dompigs[all_dompigs$diagnostic_type == "PCR" &
#               all_dompigs$diagnostic_detail == "Unassigned" &
#               all_dompigs$diagnostic_type == "PCR", "diagnostic_detail"] = "PCR"
all_dompigs %>%
    group_by(state, virus_detection) %>%
    #summarise(n=n())
    summarise(event = n_distinct(longitude, latitude, date))
# these numbers do not match the table at all ...

# load mozzie testing data
# all_mozzies <- read.csv(paste0(path, "data/national_mozzie_data/mosquito_detections_all.csv")) %>% 
# all_mozzies <- read.csv(paste0(path, 'data/national_mozzie_data/mosquito_detections_all_w_qld_sa.csv')) %>% 
#   dplyr::filter(!is.na(longitude),
#                 virus_detection %in% c('JEV positive',
#                                        'JEV negative')) 
all_mozzies <- read.csv('data/secure/phase_2/clean/pathogen/all_moz_surveillance.csv') %>%
  rename(host_type=host_species) %>% 
  dplyr::select(
    -c(site,
       n_individual_spp_tests,
    )
  ) %>% 
  drop_na(longitude:virus_detection)

all_mozzies %>%
  group_by(state, virus_detection) %>%
  summarise(events = n_distinct(longitude, latitude, date))
# o dea there are NA states in here ... surprise surprise they're all queenslanders

all_chickens <- read_csv(
  'data/secure/phase_2/clean/pathogen/all_chicken_combined.csv'
) %>% 
  dplyr::select(
    -c(yr_mn,
       sample_unit,
       n_tested)
  ) %>% 
  mutate(
    diagnostic_detail = case_when(
      diagnostic_detail=="IgM ELISA Ab" ~ "ELISA",
      #diagnostic_detail=="panflavi ELISA" ~ "panflavi_screen",
      # not sure what happened here ... panflavi ELISA/panflavi screen only applied to chickies
      # panflavi ELISA is what's in the diagnostic weights table ...
      TRUE ~ diagnostic_detail
    )
  ) %>% 
   filter(
     virus_detection %in% c("JEV positive", "JEV negative")
   )

all_chickens %>%
  group_by(state, virus_detection) %>%
  summarise(event = n_distinct(longitude, latitude, date))

plot(all_chickens[,c("longitude", "latitude")])
points(all_chickens %>%
         filter(state == "SA" & virus_detection=="JEV positive") %>%
         dplyr::select("longitude", "latitude"), col="red")

# SKIP MASK CLIPPING IN 2023 - NOW HAVE DETECTIONS IN WA

dat <- prepare_pathogen_inputs_national(all_humans,
                                        all_mozzies,
                                        all_dompigs,
                                        all_ferals,
                                        all_chickens,
                                        aus_agg_mask, 
                                        hpop_blur,
                                        mozz_bandwidth_scalar = 0.8,
                                        dpig_bandwidth_scalar = 0.8,
                                        feral_bandwidth_scalar = 0.8,
                                        chickens_bandwidth_scalar = 0.8)

dat$df <- left_join(
  dat$df,
  diagnostic_weights_table,
  by=c("diagnostic_type", "diagnostic_detail"="diagnostic_endpoint", 
       "host_type", "virus_detection")) %>% 
  mutate(
    likelihood_weight = exposure_weight * diagnostic_weight
  ) %>%
  filter(!is.na(likelihood_weight)) %>%
  # setting order of levels explicitly
  mutate(host_type = factor(host_type, 
                            levels = c("human", "mosquito", "domestic_pig", 
                                       "feral_pig", "chicken")),
         diagnostic_type = factor(diagnostic_type,
                                  levels = c("PCR", "Serology", "Unassigned")))

saveRDS(dat, "output/dat.rds")

m <- fit_pathogen_spatial_model(dat$df, dat$coordinates)
saveRDS(m, "output/national_pathogen_model_object_2026_04_04.rds")
#saveRDS(m, "output/national_pathogen_model_object_sero_demo.rds")
m <- readRDS("output/national_pathogen_model_object_2026_04_04.rds")

# this should give us back the number of (weight-adjusted) detections we give it
common_sense_check(dat$df, m$model)


##########################################################################

# pull out intercepts for the host species ids
int <- m$model$coefficients["(Intercept)"]
species_contrasts <- c(0, m$model$coefficients[grepl("host_type_id", names(m$model$coefficients))])
detection_contrasts <- c(0, m$model$coefficients[grepl("diagnostic_type", names(m$model$coefficients))])
intercepts <- int + species_contrasts 

# from 2022 model:
#                 host_species_id2 host_species_id3 host_species_id4 
# -15.200581       -14.788745        -4.665704        -5.201265

# with bug?
# s(host_type_id).1 s(host_type_id).2 s(host_type_id).3 s(host_type_id).4 
# -6.563939        -12.376992        -12.231256         -1.858678         -2.511064 
# s(host_type_id).5 
# -3.841707 

# with all of the pigs?
# s(host_type_id).1 s(host_type_id).2 s(host_type_id).3 s(host_type_id).4 s(host_type_id).5 
# -12.439706        -12.285886         -1.943317         -2.543849         -3.877953 

# without all of the pigs?
# s(host_type_id).1 s(host_type_id).2 s(host_type_id).3 s(host_type_id).4 s(host_type_id).5 
# -12.376992        -12.231256         -1.858678         -2.511064         -3.841707 
# so that's not where the problem is ...

message(paste0(length(species_contrasts), " species contrasts; ",
               length(detection_contrasts), " detection contrasts"))

# for each species, pull out the appropriate intercept, 
# specify the required intensity level, and use them to
# construct the expected rate of detection

# we need to specify this - can read it off the plot, or better to get it as a threshold
# value for the overall pathogen threshold (e.g. the 50% inclusion level)
# calculate threshold on pred_inc returned with model object
#rel_intensity_cutoff <- 0.1

m_pred_inc <- m$pred_inc

r <- m_pred_inc
r[is.na(r)] <- cellStats(r, min)
r <- mask(r, hpop_blur)

# save for later
writeRaster(r,
            'output/continuous_pathogen_national_cut.tif',
            overwrite=TRUE)

# get the threshold value - we're basing this on positives only
pos_locs = dat$df[dat$df$cases > 0, c("x", "y")]
vals <- raster::extract(m_pred_inc, pos_locs)
which(is.na(vals))

rel_intensity_cutoff <- quantile(vals, 1 - 0.5)

plot(m_pred_inc)
plot(m_pred_inc > rel_intensity_cutoff, main = "Areas with est. intensity above assumed 'high intensity' value")

# tried to dplyr this and we're just not getting anywhere :/
# tmp <- m_pred_inc %>%
#   values() %>%
#   cut(breaks = seq(0, 0.1, 0.02), labels = seq(0, 0.08, 0.02))
# tmp <- as.numeric(levels(tmp)[tmp]) %>%
#   setValues(m_pred_inc, .)

# we need to know the order of the species
species_id_lookup <- c(human = 1, mossie = 2, dpig = 3, fpig = 4, chicken = 5)
detection_id_lookup <- c(PCR = 1, serology = 2, unassigned = 3)

# length of species_contrasts and detection_contrasts should match the number of unique things in here:
unique(dat$df[, c("host_type", "host_type_id")])
unique(dat$df[, c("diagnostic_type", "diagnostic_type_id")])

# make mozzie density plot
mean_sample <- round(mean(all_mozzies$n_individuals, na.rm=TRUE), 0)

all_mozzies_expanded <- all_mozzies %>%
  mutate(n_individuals=case_when(
    is.na(n_individuals) ~ mean_sample,
    TRUE ~ round(n_individuals, 0)
  )) %>% 
  uncount(n_individuals)

tmp <- rasterize(all_mozzies_expanded[,c('longitude', 'latitude')], aus_agg_mask, 
                 background=0, fun='count')
gf <- focalWeight(aus_agg_mask, 0.8, "Gauss")
nat_mozz_density <- raster::focal(tmp, gf, pad=TRUE, na.rm=TRUE)
nat_mozz_density = mask(nat_mozz_density, aus_agg_mask)

# from 2022 model:
# threshold
# 50% 
# 115.7252 

# This is E[y_i,g,d | lambda_50] = exp(log(lambda_50) + log(beta_g) + log(gamma_d) + log(phi_i,g,d))
human_expected_detections <- exp(log(rel_intensity_cutoff) + 
                                   intercepts[species_id_lookup["human"]] +
                                   detection_contrasts[detection_id_lookup["PCR"]] +
                                    log(hpop_blur))
                              
mossie_expected_detections <- exp(log(rel_intensity_cutoff) + 
                                    intercepts[species_id_lookup["mossie"]] + 
                                    detection_contrasts[detection_id_lookup["PCR"]] +
                                    log(dat$densities$mosquito_sampling))

dpig_expected_detections <- exp(log(rel_intensity_cutoff) + 
                                  intercepts[species_id_lookup["dpig"]] + 
                                  detection_contrasts[detection_id_lookup["PCR"]] +
                                  log(dat$densities$domestic_pig_sampling)) +
                              exp(log(rel_intensity_cutoff) +
                                    intercepts[species_id_lookup["dpig"]] +
                                    detection_contrasts[detection_id_lookup["serology"]] +
                                    log(dat$densities$domestic_pig_sampling))+
                              exp(log(rel_intensity_cutoff) +
                                    intercepts[species_id_lookup["dpig"]] +
                                    detection_contrasts[detection_id_lookup["unassigned"]] +
                                    log(dat$densities$domestic_pig_sampling))

fpig_expected_detections <- exp(log(rel_intensity_cutoff) + 
                                  intercepts[species_id_lookup["fpig"]] + 
                                  detection_contrasts[detection_id_lookup["PCR"]] +
                                  log(dat$densities$feral_pig_sampling)) +
                            exp(log(rel_intensity_cutoff) +
                                  intercepts[species_id_lookup["fpig"]] +
                                  detection_contrasts[detection_id_lookup["serology"]] +
                                  log(dat$densities$feral_pig_sampling))+
                            exp(log(rel_intensity_cutoff) +
                                  intercepts[species_id_lookup["fpig"]] +
                                  detection_contrasts[detection_id_lookup["unassigned"]] +
                                  log(dat$densities$feral_pig_sampling))

chicken_expected_detections <- exp(log(rel_intensity_cutoff) +
                                      intercepts[species_id_lookup["chicken"]] +
                                      detection_contrasts[detection_id_lookup["serology"]] +
                                      log(dat$densities$sentinel_chicken_sampling))

sum(values(human_expected_detections), na.rm=TRUE)
sum(values(mossie_expected_detections), na.rm=TRUE)
sum(values(dpig_expected_detections), na.rm=TRUE)
sum(values(fpig_expected_detections), na.rm=TRUE)
sum(values(chicken_expected_detections), na.rm=TRUE)

plot(human_expected_detections)
plot(mossie_expected_detections)
plot(dpig_expected_detections)
plot(fpig_expected_detections)
plot(chicken_expected_detections)

# to include in text of MS:
human_expected_detections_per_test <- exp(log(rel_intensity_cutoff) + 
                                           intercepts[species_id_lookup["human"]] +
                                           detection_contrasts[detection_id_lookup["PCR"]])


mossie_expected_detections_per_test <- exp(log(rel_intensity_cutoff) + 
                                            intercepts[species_id_lookup["mossie"]] + 
                                            detection_contrasts[detection_id_lookup["PCR"]])


dpig_expected_detections_per_test <- exp(log(rel_intensity_cutoff) + 
                                          intercepts[species_id_lookup["dpig"]] + 
                                          detection_contrasts[detection_id_lookup["PCR"]]) +
                                    exp(log(rel_intensity_cutoff) +
                                          intercepts[species_id_lookup["dpig"]] +
                                          detection_contrasts[detection_id_lookup["serology"]])+
                                    exp(log(rel_intensity_cutoff) +
                                          intercepts[species_id_lookup["dpig"]] +
                                          detection_contrasts[detection_id_lookup["unassigned"]])


fpig_expected_detections_per_test <- exp(log(rel_intensity_cutoff) + 
                                          intercepts[species_id_lookup["fpig"]] + 
                                          detection_contrasts[detection_id_lookup["PCR"]]) +
                                    exp(log(rel_intensity_cutoff) +
                                          intercepts[species_id_lookup["fpig"]] +
                                          detection_contrasts[detection_id_lookup["serology"]]) +
                                    exp(log(rel_intensity_cutoff) +
                                          intercepts[species_id_lookup["fpig"]] +
                                          detection_contrasts[detection_id_lookup["unassigned"]])


chicken_expected_detections_per_test <- exp(log(rel_intensity_cutoff) +
                                     intercepts[species_id_lookup["chicken"]] +
                                     detection_contrasts[detection_id_lookup["serology"]])

human_expected_detections_per_test
# or rather, positives per human resident
mossie_expected_detections_per_test
# noting that scale here was per mosquito
dpig_expected_detections_per_test
fpig_expected_detections_per_test
chicken_expected_detections_per_test




# probabilities of observing one or more detections in each host species at each
# location
detection_prob <- function(expected_detections) {
  1 - dpois(0, expected_detections)
}

human_detection_probability <- calc(human_expected_detections, detection_prob)
mossie_detection_probability <- calc(mossie_expected_detections, detection_prob)
dpig_detection_probability <- calc(dpig_expected_detections, detection_prob)
fpig_detection_probability <- calc(fpig_expected_detections, detection_prob)
chicken_detection_probability <- calc(chicken_expected_detections, detection_prob)

# save for later
writeRaster(human_detection_probability,
            'output/human_detection_probability.tif',
            overwrite = TRUE)

writeRaster(mossie_detection_probability,
            'output/vector_detection_probability.tif',
            overwrite = TRUE)

writeRaster(dpig_detection_probability,
            'output/domestic_pig_detection_probability.tif',
            overwrite = TRUE)

writeRaster(fpig_detection_probability,
            'output/feral_pig_detection_probability.tif',
            overwrite = TRUE)

writeRaster(chicken_detection_probability,
            'output/chicken_detection_probability.tif',
            overwrite = TRUE)

maxValue(human_detection_probability)
maxValue(dpig_detection_probability)
maxValue(fpig_detection_probability)
maxValue(mossie_detection_probability)
maxValue(chicken_detection_probability)

# combine these to get the overall probability of detection per pixel
overall_detection_probability <- 1 - (1 - human_detection_probability) * 
  (1 - mossie_detection_probability) * (1 - dpig_detection_probability) * 
  (1 - fpig_detection_probability) * (1 - chicken_detection_probability)

plot(overall_detection_probability, main="power")
maxValue(overall_detection_probability)

# save for later
writeRaster(overall_detection_probability,
            'output/national_detection_probability_cut.tif',
            overwrite = TRUE)

###############################################################################
# SET UP FOR THRESHOLDED PLOT OF INTENSITY SURFACE

# set up for thresholded plot
# might need to edit colnames 
all_detections <- all_pres %>% 
  rename(x='longitude', y='latitude')

# set up thresholds for plot colouring
threshes <- quantile(raster::extract(m$pred_rel_inc, pos_locs),
                     probs=c(0, 0.01, 0.05,0.25,0.5))
threshes_full <- c(minValue(m$pred_rel_inc), threshes, maxValue(m$pred_rel_inc)) # with min/max for plot

# save for state-level plotting
write_csv(as.data.frame(threshes_full),
          'output/intensity_threshold_values.csv')

# threshes2 <- quantile(raster::extract(m$pred_rel_inc, all_detections[,c("x","y")]), na.rm=TRUE, 
#                       probs=c(0.05, 0.25,0.5)) 
# 
# threshes_cut <- c(minValue(m$pred_rel_inc), threshes2, maxValue(m$pred_rel_inc)) # with min/max for plot

# smooth and re-mask intensity surface
pred_rel_inc_2dp_ish <- m$pred_rel_inc %>%
  disaggregate(fact=10) %>%
  focal(w=matrix(1, 11, 11), mean, na.rm=TRUE) %>% #%>% # ignore NAs or we lose boundary pixels
  mask(aus_mask)

# save for later
writeRaster(pred_rel_inc_2dp_ish,
            'output/smooth_continuous_pathogen_national_cut.tif',
            overwrite=TRUE)


###############################################################################
# TRANSMISSION INTENSITY IN RESERVOIR

pal <- brewer.pal(9,"Reds")[2:9]
pal_cont <- colorRampPalette(pal)(100)

darkred = rgb(colorRamp(c(pal[8], "black"))(0.7)/100)

purps = colorRampPalette(brewer.pal(9,"Purples"))(100)

# harsher palette of purples for surveillance effort:
extreme_purps = colorRampPalette(purps[c(seq(1,30,5),30:100)])(100)

{png("figures/pathogen_intensity_national.png",
     width = 1000*reso/72, 
     height = 500*reso/72, 
     res = reso)
  #type = 'cairo')
  
  par(mfrow = c(1,2),
      mar = c(1, 0.1, 1, 1), bty="n",
      oma=c(0,0,0,0))
  
  plot(st_geometry(aus_simple), lwd=0.0001, border='grey20', col="#f7f7f7",
       main="Thresholded intensity of JEV transmission", cex.main=1.5)
  plot(pred_rel_inc_2dp_ish, breaks=threshes_full, 
       col=c("#f7f7f7", pal[(c(1,2,3,4,5))]), legend=FALSE,
       xaxt="n", yaxt="n", legend.mar=0, add=TRUE)
  plot(rasterToContour(m$pred_rel_inc, levels=threshes_full[2:(length(threshes_full - 1))]), 
       add=TRUE, col=darkred)
  plot(st_geometry(aus_simple), lwd=1, border='grey20', add=TRUE)
  #points(all_detections[,c('x','y')], pch=22, col=alpha('black', 0.2), lwd=0.8, cex=1)
  # text(c(125, 131, 135, 133, 143), 
  #      c(-30, -26.6, -24, -19.8, -27.9), 
  #      c("100%", "99%", "95%", "75%", "50%"), cex=1)
  text(c(124, 123, 131, 144, 143), 
       c(-26, -22, -23.5, -32, -27.9), 
       c("100%", "99%", "95%", "75%", "50%"), cex=1)
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1.6)
  
  # plot overall detection probability
  detect_breaks <- c(0, 0.01, 0.1, 0.25, 0.5, 0.75)
  #detect_breaks <- c(0, 0.0001, 0.001, 0.0025, 0.005)
  
  plot(st_geometry(aus_simple), lwd=0.01, border='grey20',main="Power to detect JEV", cex.main=1.5)
  plot(sqrt(overall_detection_probability), xaxt="n", yaxt="n", 
       col=extreme_purps, axis.args = list(at = sqrt((detect_breaks)), labels = detect_breaks), 
       add=TRUE, legend.width=1, horizontal=TRUE)
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1.6)
  
  dev.off()}

tmp = overall_detection_probability
tmp[tmp <= 0.002] = NA
values(tmp)[!is.na(values(tmp))] = 1
plot(st_geometry(aus_simple), main="Power > 0.002")
plot(tmp, legend=FALSE, add=TRUE, col="red")



# {png("figures/pathogen_intensity_national_sero_demo.png",
#      width = 1000*reso/72, 
#      height = 500*reso/72, 
#      res = reso)
#   #type = 'cairo')
#   
#   par(mfrow = c(1,2),
#       mar = c(1, 0.1, 1, 1), bty="n",
#       oma=c(0,0,0,0))
#   
#   plot(st_geometry(aus_simple), lwd=0.0001, border='grey20', col="#f7f7f7",
#        main="Thresholded intensity of JEV transmission", cex.main=1.5)
#   plot(pred_rel_inc_2dp_ish, breaks=threshes_full, 
#        col=c("#f7f7f7", pal[(c(1,2,3,4,5))]), legend=FALSE,
#        xaxt="n", yaxt="n", legend.mar=0, add=TRUE)
#   plot(rasterToContour(m$pred_rel_inc, levels=threshes_full[2:(length(threshes_full - 1))]), 
#        add=TRUE, col=darkred)
#   plot(st_geometry(aus_simple), lwd=1, border='grey20', add=TRUE)
#   #points(all_detections[,c('x','y')], pch=22, col=alpha('black', 0.2), lwd=0.8, cex=1)
#   # text(c(125, 131, 135, 133, 143), 
#   #      c(-30, -26.6, -24, -19.8, -27.9), 
#   #      c("100%", "99%", "95%", "75%", "50%"), cex=1)
#   text(c(125, 131, 121, 148, 143), 
#        c(-30, -27.5, -21, -27.9, -27.9), 
#        c("100%", "99%", "95%", "75%", "50%"), cex=1)
#   subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1.6)
#   
#   # plot overall detection probability
#   #detect_breaks <- c(0, 0.01, 0.1, 0.25, 0.5)
#   detect_breaks <- c(0, 0.0001, 0.001, 0.005, 0.01)
#   
#   plot(st_geometry(aus_simple), lwd=0.01, border='grey20',main="Power to detect JEV", cex.main=1.5)
#   plot(sqrt(overall_detection_probability), xaxt="n", yaxt="n", 
#        col=extreme_purps, 
#        axis.args = list(at = sqrt((detect_breaks)), labels = detect_breaks), 
#        add=TRUE, legend.width=1, horizontal=TRUE)
#   plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE)
#   subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1.6)
#   
#   dev.off()}

# {png("figures/detection_probabilities.png",
#      width = 1000*reso/72, 
#      height = 1100*reso/72, 
#      res = reso)
# par(mfrow = c(3,2))
# plot(human_detection_probability, col=viridis(100), main = "People")
# plot(st_geometry(aus_simple), add=TRUE, lwd=0.8, border='grey80')
# plot(dpig_detection_probability, col=viridis(100), main = "D pig")
# plot(st_geometry(aus_simple), add=TRUE, lwd=0.8, border='grey80')
# plot(fpig_detection_probability, col=viridis(100), main = "F pig")
# plot(st_geometry(aus_simple), add=TRUE, lwd=0.8, border='grey80')
# plot(mossie_detection_probability, col=viridis(100), main = "Mosquito")
# plot(st_geometry(aus_simple), add=TRUE, lwd=0.8, border='grey80')
# plot(chicken_detection_probability, col=viridis(100), main = "Chix")
# plot(st_geometry(aus_simple), add=TRUE, lwd=0.8, border='grey80')
# dev.off()}


sub_cex_main <- 2
leg_cex_main <- 1.3
leg_cex_axis <- 1.4
{png("figures/detection_probabilities.png",
     width = 1000*reso/72, 
     height = 1100*reso/72, 
     res = reso)
  
  par(mfrow = c(3, 2),
      mar = c(0.1, 0.1, 3.1, 1), bty="n",
      oma=c(0,0,0,3))
  
  # humans
  axis_breaks = c(0.01, 0.1, 0.25, 0.5, 0.75)
  plot(sqrt(human_detection_probability), 
       col=extreme_purps, legend.width=1.4,
       legend.args = list(text='Power to detect', 
                          line=0.9, side=2, cex=leg_cex_main), 
       axis.args = list(at = sqrt(axis_breaks),
                        labels = axis_breaks, 
                        cex.axis = leg_cex_axis),
       box = FALSE, axes = FALSE, cex.main=sub_cex_main,
       main="Humans")
  plot(st_geometry(aus_simple), add=TRUE, lwd=1.2, border='grey20')
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", sub_cex_main)
  
  # mosquitoes
  plot(mossie_detection_probability,
    box = FALSE, axes=FALSE, col=extreme_purps, 
    legend.width=1.4,
    legend.args = list(text='Power to detect', 
                       line=0.9, cex=leg_cex_main, side=2), 
    axis.args = list(#at=c(0, 5e-6, 1e-5, 1.5e-5, 2e-5),
                     #labels=c("0", "5e-6", "1e-5", "1.5e-5", "2e-5"),
                     cex.axis=leg_cex_axis),
    cex.main=sub_cex_main,
    main="Mosquitoes")
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", sub_cex_main)
  
  # domestic pigs
  plot(dpig_detection_probability, box = FALSE, axes=FALSE, 
       col=extreme_purps,
       legend.args = list(text='Power to detect', 
                          cex=leg_cex_main, line=0.9, side=2), 
       axis.args = list(#at=seq(0, 6e-04, length.out=4),
                        #labels = c("0", "2e-4", "4e-4", "6e-4"),
                        cex.axis=leg_cex_axis),
       legend.width=1.4, cex.main=sub_cex_main,
       main="Domestic pigs")
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", sub_cex_main)
  
  # feral pigs
  plot(fpig_detection_probability, box = FALSE, axes=FALSE, 
       col=extreme_purps,
       legend.args = list(text='Power to detect', cex=leg_cex_main, 
                          line=0.9, side=2),
       axis.args = list(#at=seq(0, 0.002, length.out = 5),
                        #labels=c("0","2e-3","4e-3","6e-3","8e-3"),
                        cex.axis=leg_cex_axis),
       legend.width=1.4, cex.main=sub_cex_main,
       main="Feral pigs")
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE) 
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", sub_cex_main)
  
  plot(chicken_detection_probability, box=FALSE, axes=FALSE,
       col=extreme_purps,
       legend.args = list(text='Power to detect', cex=leg_cex_main, 
                          line=0.9, side=2),
       axis.args = list(# at=c(0,5e-4,1e-3,1.5e-3,2e-3), 
                       # at = seq(0, 0.00002, length.out=5),
                       # labels=c("0","5e-6","1e-5","1.5e-5","2e-5"),
         cex.axis=leg_cex_axis),
       legend.width=1.4, cex.main=sub_cex_main,
       main="Sentinel chickens")
  plot(st_geometry(aus_simple), lwd=1.2, border='grey20', add=TRUE)
  subfigure_label(par("usr"), 0.1, 0.9, "(e)", sub_cex_main)
  
  dev.off()}





################################################################################
# mozzie surveillance against detection: coarser window and sum detections?

# nrow(dat$df %>% filter(host_type == "mosquito"))
# nrow(dat$mozz_points)
# 2501+4196+589+1372+1026+747 # numbers in paper

# plot(as.factor(all_mozzies$state), all_mozzies$n_individuals)
# points(as.factor(all_mozzies$state)[all_mozzies$virus_detection == "JEV positive"], 
#        all_mozzies$n_individuals[all_mozzies$virus_detection == "JEV positive"], col="red")

# all_mozzies %>%
#   group_by(state, virus_detection) %>%
#   summarise(events = n_distinct(longitude, latitude, date)) %>% 
#   ungroup() %>% 
#   dplyr::select(events) %>% 
#   sum()

# things are likely to be weird in the NT as we imputed dates...
# all_mozzies %>%
#   group_by(state, virus_detection, longitude, latitude, date) %>%
#   summarise(n = n()) %>%
#   filter(n > 1) %>%
#   group_by(state, virus_detection) %>%
#   summarise(n = n())

# tests <- dat$mozz_points %>%
#   dplyr::select(x, y) %>%
#   rasterize(hpop_blur, background = 0, fun = "count") %>%
#   aggregate(fact = 10, sum)
# 
# tests <- dat$mozz_points %>%
#   dplyr::select(x, y) %>%
#   rasterize(hpop_blur %>% aggregate(fact = 10), 
#             background = 0, field = dat$mozz_points$exposure_weight, fun = sum)

# calculate total number of mozzies tested (irrespective of number of tests)
# exposure_weight doesn't do much as it's only < 1 for those NT +ves
# tests <- dat$mozz_points %>%
#   dplyr::select(x, y) %>%
#   rasterize(hpop_blur %>% aggregate(fact = 10), 
#             background = 0, field = dat$mozz_points$density*dat$mozz_points$exposure_weight, fun = sum)
# 
# tests[tests == 0] <- NA
# 
# # dets <- dat$mozz_points %>%
# #   filter(virus_detection == "JEV positive") %>%
# #   dplyr::select(x, y) %>%
# #   rasterize(hpop_blur, background = 0, fun = "count") %>%
# #   aggregate(fact = 10, sum)
# 
# # calculate total number of +ve tests
# # (exposure_weight is important in the NT)
# dets <- dat$mozz_points %>%
#   filter(virus_detection == "JEV positive")
# dets <- dets %>% 
#   dplyr::select(x, y) %>%
#   rasterize(hpop_blur %>% aggregate(fact = 10), 
#             background = 0, field = dets$exposure_weight, fun = sum)
# 
# dets[dets == 0] <- NA

# {png("~/Desktop/mosquito_power.png", height = 1800, width = 1000, pointsize = 30)
# par(mfrow = c(3,1), oma = c(0,0,0,2), bty="n", mar = rep(1,4))
# plot(tests, col = viridis(100), main = "Total mosquitoes tested", axes=FALSE)
# subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1)
# plot(st_geometry(aus_simple), add=TRUE)
# plot(dets, col = viridis(3), breaks = seq(0.5, 3.5), axis.args = list(at = seq(1,3)),
#      main = "Total detections", axes=FALSE)
# subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1)
# plot(st_geometry(aus_simple), add=TRUE)
# leg_breaks = c(1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3)
# plot(log10(dets / tests), col = viridis(100), main="Detections / mosquitoes tested", 
#      #breaks = log10(leg_breaks), 
#      legend=FALSE, axes=FALSE)
# subfigure_label(par("usr"), 0.1, 0.9, "(c)", 1)
# plot(log10(dets/tests), col=viridis(100), legend.only = TRUE, 
#      #breaks = log10(leg_breaks),
#      axis.args = list(at = log10(leg_breaks), labels = leg_breaks))
# plot(rasterToContour(m$pred_rel_inc, levels=threshes_full["50%"]), 
#      add=TRUE, col=darkred) #, fill = alpha(brewer.pal(9,"Reds")[3], 0.5))
# plot(st_geometry(aus_simple), add=TRUE)
# dev.off()}

# just vic:
# par(mfrow = c(1,3), oma = c(0,0,0,2))
# plot(tests, col = viridis(100), xlim = c(140, 150), ylim = c(-40, -30), main = "Total mozzies tested")
# plot(st_geometry(aus_simple), add=TRUE)
# plot(dets, col = viridis(100), xlim = c(140, 150), ylim = c(-40, -30), main = "Total mozzies in +ve tests")
# plot(st_geometry(aus_simple), add=TRUE)
# plot(dets / tests, col = viridis(100), xlim = c(140, 150), ylim = c(-40, -30), main="+ve / tested")
# plot(st_geometry(aus_simple), add=TRUE)

# plot(mossie_detection_probability, col=viridis(100),
#      main="Mozzie detection probability")
# plot(st_geometry(aus_simple), add=TRUE, border="white")
# plot(rasterToContour(m$pred_rel_inc, levels=threshes_full["50%"]), 
#      add=TRUE, col="lightyellow")

# have an equivalent look at chickens:
# different to mozzies as `density` is from KDE ...

################################################################################
# now having a look at chickens:
# chix <- dat$df %>%
#   filter(host_type == "chicken")
# tests <- chix %>%
#   dplyr::select(x, y) %>%
#   rasterize(hpop_blur %>% aggregate(fact = 10), 
#             background = 0, field = chix$diagnostic_weight, fun = sum)
# tests[tests == 0] <- NA
# 
# # calculate total number of +ve tests
# dets <- chix %>%
#   filter(virus_detection == "JEV positive")
# dets <- dets %>% 
#   dplyr::select(x, y) %>%
#   rasterize(hpop_blur %>% aggregate(fact = 10), 
#             background = 0, field = dets$likelihood_weight, fun = sum)
# dets[dets == 0] <- NA
# 
# par(mfrow = c(2,2), oma = c(0,0,0,2), mar=rep(3,4))
# plot(tests, col = viridis(100), main = "Total chicken tests \n(adjusted for diagnostic weight)")
# plot(st_geometry(aus_simple), add=TRUE)
# plot(dets, col = viridis(10),
#      main = "Total chicken detections \n(adjusted for diagnostic weight)")
# plot(st_geometry(aus_simple), add=TRUE)
# plot(dets/tests, col=viridis(100), main="detections / tests")
# plot(rasterToContour(m$pred_rel_inc, levels=threshes_full["50%"]), 
#      add=TRUE, col=darkred) #, fill = alpha(brewer.pal(9,"Reds")[3], 0.5))
# plot(st_geometry(aus_simple), add=TRUE)
# 
# 
# {png("~/Desktop/chicken_power.png", height = 1800, width = 1000, pointsize = 30)
#   par(mfrow = c(3,1), oma = c(0,0,0,2), bty="n", mar = rep(1,4))
#   plot(tests, col = viridis(100), main = "Total chicken tests*", axes=FALSE)
#   subfigure_label(par("usr"), 0.1, 0.9, "(a)", 1)
#   plot(st_geometry(aus_simple), add=TRUE)
#   plot(dets, col = viridis(10),
#        main = "Total chicken detections*", axes=FALSE)
#   subfigure_label(par("usr"), 0.1, 0.9, "(b)", 1)
#   plot(st_geometry(aus_simple), add=TRUE)
#   leg_breaks = c(1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3)
#   plot(dets/tests, col=viridis(100), main="Detections / tests", axes=FALSE)
#   subfigure_label(par("usr"), 0.1, 0.9, "(c)", 1)
#   plot(rasterToContour(m$pred_rel_inc, levels=threshes_full["50%"]), 
#        add=TRUE, col=darkred) #, fill = alpha(brewer.pal(9,"Reds")[3], 0.5))
#   plot(st_geometry(aus_simple), add=TRUE)
#   dev.off()}

###############################################################################
# LH TODO: raw ratio within thresholded region (lambda_50) of total number of 
# tests, and total number of positives


# how many mozzies tested and how many detections within thresholded region?
# high_intens_region <- m_pred_inc
# high_intens_region[high_intens_region < rel_intensity_cutoff] <- NA
# 
# rat <- function(df, ras){
#   df <- df %>% filter(virus_detection %in% c("JEV positive", "JEV negative"))
#   
#   tests <- df %>%
#     dplyr::select(x, y) %>%
#     rasterize(ras, 
#               background = 0, 
#               field = df$density *
#                 df$exposure_weight, fun = sum) %>%
#     mask(ras)
#   
#   
#   dets <- df %>% filter(virus_detection == "JEV positive")
#   
#   dets <- dets %>%
#     dplyr::select(x, y) %>%
#     rasterize(ras, 
#               background = 0, 
#               field = df$exposure_weight, fun = sum) %>%
#     mask(ras)
#   
#   plot(tests)
#   plot(dets)
#     
#   sum(values(dets), na.rm=TRUE) /
#     sum(values(tests), na.rm=TRUE)
# }
# 
# rat(dat$mozz_points, high_intens_region)
# rat(dat$dpig_points, high_intens_region)
# rat(dat$feral_points, high_intens_region)
# rat(dat$chicken_points, high_intens_region)
# this doesn't work for humans - we didn't assign outcomes/exposure weights as they were all pos



# dat$df %>% 
#   filter(host_type == "mosquito") %>%
#   st_as_sf(coords = c("x", "y"), crs = st_crs(oz_outline)) %>%
#   st_join(oz_outline %>% dplyr::select(name)) %>%
#   st_drop_geometry() %>%
#   filter(name == "Queensland")

