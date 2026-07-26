# power-to-detect: can we make sense out of this?
library(terra)

lgas <- st_read(paste0("/Users/harrisonl2/Desktop/jev/data/admin/", 
                       '/LGA_2021_AUST_GDA2020_SHP/LGA_2021_AUST_GDA2020.shp')) %>%
  rename(lga_name = LGA_NAME21) %>%
  filter(STE_NAME21 != "Other Territories" & STE_NAME21 != "Outside Australia") %>%
  # dplyr::select(lga_name) %>%
  # mutate(alt_lga_name = gsub(" Regional", "", lga_name)) %>%
  # mutate(alt_lga_name = gsub(" Shire", "", alt_lga_name)) %>%
  # filter(alt_lga_name %in% hpop_df$lga_name)
  st_simplify(dTolerance = 1000) %>%
  filter(!st_is_empty(st_geometry(.)))

###############################################################################
# bring in bits from pathogen_model.R

m <- readRDS("output/national_pathogen_model_object_2025_12_11.rds")
dat <- readRDS("output/dat.rds")


# rel_intensity_cutoff
# intercepts
# detection_contrasts

###############################################################################
# realised surveillance effort - discarding KDEs for power estimates
hpop_blur <- raster('data/hpop_blur_aus_0.2.tif')

# rasterize tests and extract for LGAs
rasterize_tests_and_dets <- function(df, lgas, ras){
  df <- df %>% 
    filter(virus_detection %in% c("JEV positive", "JEV negative"))
  
  tests <- df %>%
    dplyr::select(x, y) %>%
    terra::rasterize(ras, 
              background = 0, 
              field = df$density * df$exposure_weight, 
              fun = sum)
  
  dets <- df %>% filter(virus_detection == "JEV positive")
  
  dets <- dets %>%
    dplyr::select(x, y) %>%
    terra::rasterize(ras, 
              background = 0, 
              field = df$exposure_weight, fun = sum)
  
  c(rast(tests), rast(dets))
}

tests_and_dets <- c(rasterize_tests_and_dets(dat$mozz_points, lgas, hpop_blur),
                   rasterize_tests_and_dets(dat$chicken_points, lgas, hpop_blur),
                   rasterize_tests_and_dets(dat$feral_points, lgas, hpop_blur),
                   rasterize_tests_and_dets(dat$dpig_points, lgas, hpop_blur))
names(tests_and_dets) <- apply(expand.grid(c("tests", "dets"),
                                           c("mozz", "chix", "feral", "dom")),
                               1, function(x) paste(x[2], x[1], sep = "_"))

# humans are a little annoying
tests_and_dets$human_tests <- rast(hpop_blur)

# (it seems that exposure weights have gone missing here but they were extant
# in earlier model iterations?)
tests_and_dets$human_dets <- dat$human_points %>%
  as.data.frame() %>%
  rename(x = longitude, y = latitude) %>%
  rasterize(hpop_blur, 
            background = 0, field = 1, fun = sum) %>%
  rast()

# now to extract tests and dats
test_det_lga = terra::extract(tests_and_dets, vect(lgas))

test_det_lga <- group_by(test_det_lga, ID) %>%
  summarise(across(everything(), sum), npix = n())

lgas_tests = bind_cols(lgas, test_det_lga)

plot(lgas_tests[, "human_dets"], pal = viridis, nbreaks = 99, 
     border = NA)
plot(lgas_tests[, "mozz_tests"], pal = viridis, nbreaks = 99, 
     border = NA)

###############################################################################
# calculate expected detections per LGA given tests per LGA
# (maths as in pathogen_model.R, substituting LGA tests)

human_expected_detections <- exp(log(rel_intensity_cutoff) + 
                                   intercepts[species_id_lookup["human"]] +
                                   detection_contrasts[detection_id_lookup["PCR"]] +
                                   log(lgas_tests$human_tests))

mossie_expected_detections <- exp(log(rel_intensity_cutoff) + 
                                    intercepts[species_id_lookup["mossie"]] + 
                                    detection_contrasts[detection_id_lookup["PCR"]] +
                                    log(lgas_tests$mozz_tests))

dpig_expected_detections <- exp(log(rel_intensity_cutoff) + 
                                  intercepts[species_id_lookup["dpig"]] + 
                                  detection_contrasts[detection_id_lookup["PCR"]] +
                                  log(lgas_tests$dom_tests)) +
  exp(log(rel_intensity_cutoff) +
        intercepts[species_id_lookup["dpig"]] +
        detection_contrasts[detection_id_lookup["serology"]] +
        log(lgas_tests$dom_tests))+
  exp(log(rel_intensity_cutoff) +
        intercepts[species_id_lookup["dpig"]] +
        detection_contrasts[detection_id_lookup["unassigned"]] +
        log(lgas_tests$dom_tests))

fpig_expected_detections <- exp(log(rel_intensity_cutoff) + 
                                  intercepts[species_id_lookup["fpig"]] + 
                                  detection_contrasts[detection_id_lookup["PCR"]] +
                                  log(lgas_tests$feral_tests)) +
  exp(log(rel_intensity_cutoff) +
        intercepts[species_id_lookup["fpig"]] +
        detection_contrasts[detection_id_lookup["serology"]] +
        log(lgas_tests$feral_tests))+
  exp(log(rel_intensity_cutoff) +
        intercepts[species_id_lookup["fpig"]] +
        detection_contrasts[detection_id_lookup["unassigned"]] +
        log(lgas_tests$feral_tests))

chicken_expected_detections <- exp(log(rel_intensity_cutoff) +
                                     intercepts[species_id_lookup["chicken"]] +
                                     detection_contrasts[detection_id_lookup["serology"]] +
                                     log(lgas_tests$chix_tests))

lgas_expected_dets <- lgas_tests %>%
  mutate(chix_expected_dets = chicken_expected_detections,
         feral_expected_dets = fpig_expected_detections,
         dom_expected_dets = dpig_expected_detections,
         mozz_expected_dets = mossie_expected_detections,
         human_expected_dets = human_expected_detections)

###############################################################################
# populating down detection probability
detection_prob <- function(expected_detections) {
  1 - dpois(0, expected_detections)
}

lgas_det_probs <- lgas_expected_dets %>%
  mutate(across(ends_with("expected_dets"), 
                detection_prob,
                .names = "{sub('expected_dets', 'det_prob', .col)}"))

plot(lgas_det_probs[, "mozz_det_prob"], pal = viridis, nbreaks = 99, border = NA)
plot(lgas_tests[, "mozz_tests"], pal = viridis, nbreaks = 99, border = NA)




purps = colorRampPalette(brewer.pal(9,"Purples"))(100)

# harsher palette of purples for surveillance effort:
extreme_purps = colorRampPalette(purps[c(seq(1,30,5),30:100)])

pal = c("blue", colorRampPalette(c("white", "red"))(99))
pal = viridis(100)

sub_cex_main <- 2
leg_cex_main <- 1.3
leg_cex_axis <- 1.4
ncol = 100

det_prob_range <- lgas_det_probs %>%
  st_drop_geometry() %>%
  dplyr::select(ends_with("det_prob")) %>%
  as.vector() %>%
  range(na.rm = TRUE)
#breaks <- seq(sqrt(min(det_prob_range)), sqrt(max(det_prob_range)), length.out = 100)**2
breaks = seq(0, 1, length.out = 101)

# this should potentially be expected detections if there's going to be lots of ones ..?
{png("figures/detection_probabilities_lga.png",
     width = 1000*reso/72, 
     height = 1100*reso/72, 
     res = reso)
  
  # can't stop this from being reset after each panel, but now have to match the reset
  # settings
  par(mfrow = c(3, 2),
      mar = c(0.1, 0.1, 1.2, 2.2), bty="n",
      oma=c(0,0,0,2))
  # 
  # layout(matrix(1:12, ncol = 4, byrow = TRUE), widths = c(1, 0.08, 1, 0.08))
  
  # humans
  # breaks = seq(min(lgas_det_probs$human_det_prob, na.rm=TRUE),
  #              max(lgas_det_probs$human_det_prob, na.rm=TRUE),
  #              length.out = 100)
  plot(st_geometry(lgas_det_probs[,"human_det_prob"]), 
       col = pal[cut(lgas_det_probs$human_det_prob, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Humans", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", sub_cex_main)
  
  # .image_scale(breaks, col = extreme_purps(ncol), 
  #              key.length = lcm(20), 
  #              key.width = lcm(0.5), # erm I checked the src and this doesn't appear
  #              key.pos = 4, cex.axis = leg_cex_axis,
  #              at = seq(1e-2, 4e-2, length.out = 5))
  # mtext("Probability of detection", 4, line = -4)
  
  
  
  # mosquitoes
  # breaks = seq(min(lgas_det_probs$mozz_det_prob, na.rm=TRUE),
  #              max(lgas_det_probs$mozz_det_prob, na.rm=TRUE),
  #              length.out = 100)
  plot(st_geometry(lgas_det_probs[,"mozz_det_prob"]), 
       col = pal[cut(lgas_det_probs$mozz_det_prob, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Mosquitoes", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", sub_cex_main)
  
  # .image_scale(breaks, col = extreme_purps(ncol), 
  #              key.length = lcm(20), 
  #              key.width = lcm(0.5), # erm I checked the src and this doesn't appear
  #              key.pos = 4, cex.axis = leg_cex_axis,
  #              at = seq(0, 5e-03, length.out = 6))
  # mtext("Probability of detection", 4, line = -4)
  
  # domestic pigs
  # breaks = seq(min(lgas_det_probs$dom_det_prob, na.rm=TRUE),
  #              max(lgas_det_probs$dom_det_prob, na.rm=TRUE),
  #              length.out = 100)
  plot(st_geometry(lgas_det_probs[,"dom_det_prob"]), 
       col = pal[cut(lgas_det_probs$dom_det_prob, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Domestic pigs", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", sub_cex_main)
  
  # .image_scale(breaks, col = extreme_purps(ncol), 
  #              key.length = lcm(20), 
  #              key.width = lcm(0.5), # erm I checked the src and this doesn't appear
  #              key.pos = 4, cex.axis = leg_cex_axis,
  #              at = seq(0, 0.125, length.out = 6))
  # mtext("Probability of detection", 4, line = -4)
  
  # feral pigs
  # breaks = seq(min(lgas_det_probs$feral_det_prob, na.rm=TRUE),
  #              max(lgas_det_probs$feral_det_prob, na.rm=TRUE),
  #              length.out = 100)
  plot(st_geometry(lgas_det_probs[,"feral_det_prob"]), 
       col = pal[cut(lgas_det_probs$feral_det_prob, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Feral pigs", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", sub_cex_main)
  
  # chickens
  plot(st_geometry(lgas_det_probs[,"chix_det_prob"]), 
       col = pal[cut(lgas_det_probs$chix_det_prob, breaks=breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Sentinel chickens", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(e)", sub_cex_main)
  
  empty_plot_for_legend()
  plot(raster(matrix(breaks)), 
       col = pal, 
       legend.only = TRUE, 
       horizontal = TRUE,
       smallplot=c(0.25, 0.75, 0.48,0.52),
       axis.args = list(cex.axis = leg_cex_axis),
       legend.args = list("Probability of detection", 
                          line = 1, cex = leg_cex_main))
  legend(0.82, 0.5, fill = "white", "No surveillance recorded", 
         bty = "n", cex = 2)
  dev.off()}

exp_det_range <- lgas_det_probs %>%
  st_drop_geometry() %>%
  dplyr::select(ends_with("expected_dets")) %>%
  as.vector() %>%
  range(na.rm = TRUE)
breaks <- seq(sqrt(min(exp_det_range)), sqrt(max(exp_det_range)), length.out = 100)**2

# expected detections - sqrt this a go for a more colourful palette ..
{png("figures/expected_detections_lga.png",
     width = 1000*reso/72, 
     height = 1100*reso/72, 
     res = reso)
  # can't stop this from being reset after each panel, but now have to match the reset
  # settings
  par(mfrow = c(3, 2),
      mar = c(0.1, 0.1, 1.2, 2.2), bty="n",
      oma=c(0,0,0,2))
  
  # humans
  plot(st_geometry(lgas_det_probs[,"human_expected_dets"]), 
       col = pal[cut(lgas_det_probs$human_expected_dets, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Humans", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", sub_cex_main)
  
  # .image_scale(breaks, col = extreme_purps(ncol), 
  #              key.length = lcm(20), 
  #              key.width = lcm(0.5), # erm I checked the src and this doesn't appear
  #              key.pos = 4, 
  #              cex.axis = leg_cex_axis#,
  #              #at = seq(1e-2, 5e-2, length.out = 5)
  #              )
  # mtext("Expected detections", 4, line = -4)
  
  
  # mosquitoes
  # breaks = seq(min(lgas_det_probs$mozz_expected_dets, na.rm=TRUE),
  #              max(lgas_det_probs$mozz_expected_dets, na.rm=TRUE),
  #              length.out = 100)
  plot(st_geometry(lgas_det_probs[,"mozz_expected_dets"]), 
       col = pal[cut(lgas_det_probs$mozz_expected_dets, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Mosquitoes", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", sub_cex_main)
  
  # .image_scale(breaks, col = extreme_purps(ncol), 
  #              key.length = lcm(20), 
  #              key.width = lcm(0.5), # erm I checked the src and this doesn't appear
  #              key.pos = 4, cex.axis = leg_cex_axis#,
  #              #at = seq(0, 5e-03, length.out = 6)
  #              )
  # mtext("Expected detections", 4, line = -4)
  
  # domestic pigs
  # breaks = seq(min(lgas_det_probs$dom_expected_dets, na.rm=TRUE),
  #              max(lgas_det_probs$dom_expected_dets, na.rm=TRUE),
  #              length.out = 100)
  plot(st_geometry(lgas_det_probs[,"dom_expected_dets"]), 
       col = pal[cut(lgas_det_probs$dom_expected_dets, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Domestic pigs", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", sub_cex_main)
  
  # .image_scale(breaks, 
  #              col = extreme_purps(ncol), 
  #              key.length = lcm(20), 
  #              key.width = lcm(0.5), # erm I checked the src and this doesn't appear
  #              key.pos = 4, 
  #              cex.axis = leg_cex_axis#,
  #              #at = seq(0, 0.125, length.out = 6)
  #              )
  # mtext("Expected detections", 4, line = -4)
  
  
  # feral pigs
  # breaks = seq(min(lgas_det_probs$feral_expected_dets, na.rm=TRUE),
  #              max(lgas_det_probs$feral_expected_dets, na.rm=TRUE),
  #              length.out = 100)
  plot(st_geometry(lgas_det_probs[,"feral_expected_dets"]), 
       col = pal[cut(lgas_det_probs$feral_expected_dets, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Feral pigs", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", sub_cex_main)
  
  # .image_scale(breaks, 
  #              col = extreme_purps(ncol), 
  #              key.length = lcm(20), 
  #              key.width = lcm(0.5), # erm I checked the src and this doesn't appear
  #              key.pos = 4, cex.axis = leg_cex_axis#,
  #              #at = seq(0, 0.8, length.out = 5)
  #              )
  # mtext("Expected detections", 4, line = -4)
  
  # chickens
  # breaks = seq(min(lgas_det_probs$chix_expected_dets, na.rm=TRUE),
  #              max(lgas_det_probs$chix_expected_dets, na.rm=TRUE),
  #              length.out = 100)
  plot(st_geometry(lgas_det_probs[,"chix_expected_dets"]), 
       col = pal[cut(lgas_det_probs$chix_expected_dets, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Sentinel chickens", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(e)", sub_cex_main)
  
  # .image_scale(breaks, 
  #              col = extreme_purps(ncol), 
  #              key.length = lcm(20), 
  #              key.width = lcm(0.5), # erm I checked the src and this doesn't appear
  #              key.pos = 4, cex.axis = leg_cex_axis#,
  #              #at = seq(0, 5e-3, length.out = 5)
  #              )
  # mtext("Expected detections", 4, line = -4)
  
  #empty_plot_for_legend()
  # .image_scale(breaks,
  #              col = extreme_purps(ncol),
  #              key.length = lcm(20),
  #              key.width = lcm(0.5), # erm I checked the src and this doesn't appear
  #              key.pos = 4, cex.axis = leg_cex_axis#,
  #              #at = seq(0, 5e-3, length.out = 5)
  #              )
  
  empty_plot_for_legend()
  plot(raster(matrix(breaks)), 
       col = extreme_purps(ncol), 
       legend.only = TRUE, 
       horizontal = TRUE,
       smallplot=c(0.25, 0.75, 0.48,0.52),
       axis.args = list(cex.axis = leg_cex_axis),
       legend.args = list("Expected Detections", line = 1, cex = leg_cex_main))
  dev.off()}



# total tests
{png("figures/tests_lga.png",
     width = 1000*reso/72, 
     height = 1100*reso/72, 
     res = reso)
  
  # can't stop this from being reset after each panel, but now have to match the reset
  # settings
  par(mfrow = c(3, 2),
      mar = c(0.1, 0.1, 1.2, 2.2), bty="n",
      oma=c(0,0,0,2))
  
  layout(matrix(1:12, ncol = 4, byrow = TRUE), widths = c(1, 0.08, 1, 0.08))
  
  # humans
  breaks = seq(min(lgas_det_probs$human_tests, na.rm=TRUE),
               max(lgas_det_probs$human_tests, na.rm=TRUE),
               length.out = 100)
  plot(st_geometry(lgas_det_probs[,"human_expected_dets"]), 
       col = extreme_purps(ncol)[cut(lgas_det_probs$human_tests, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Humans", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", sub_cex_main)
  
  .image_scale(breaks, col = extreme_purps(ncol), 
               key.length = lcm(20), 
               key.width = lcm(0.5), # erm I checked the src and this doesn't appear
               key.pos = 4, cex.axis = leg_cex_axis,
               at = seq(1e5, 5e5, length.out = 5))
  mtext("Tests", 4, line = -4)
  
  
  # mosquitoes
  breaks = seq(min(lgas_det_probs$mozz_tests, na.rm=TRUE),
               max(lgas_det_probs$mozz_tests, na.rm=TRUE),
               length.out = 100)
  plot(st_geometry(lgas_det_probs[,"mozz_expected_dets"]), 
       col = extreme_purps(ncol)[cut(lgas_det_probs$mozz_tests, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Mosquitoes", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", sub_cex_main)
  
  .image_scale(breaks, col = extreme_purps(ncol), 
               key.length = lcm(20), 
               key.width = lcm(0.5), # erm I checked the src and this doesn't appear
               key.pos = 4, cex.axis = leg_cex_axis,
               at = seq(0, 8e4, length.out = 5))
  mtext("Tests", 4, line = -4)
  
  # domestic pigs
  breaks = seq(min(lgas_det_probs$dom_tests, na.rm=TRUE),
               max(lgas_det_probs$dom_tests, na.rm=TRUE),
               length.out = 100)
  plot(st_geometry(lgas_det_probs[,"dom_expected_dets"]), 
       col = extreme_purps(ncol)[cut(lgas_det_probs$dom_tests, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Domestic pigs", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", sub_cex_main)
  
  .image_scale(breaks, col = extreme_purps(ncol), 
               key.length = lcm(20), 
               key.width = lcm(0.5), # erm I checked the src and this doesn't appear
               key.pos = 4, cex.axis = leg_cex_axis,
               at = seq(0, 40, length.out = 5))
  mtext("Tests", 4, line = -4)
  
  
  # feral pigs
  breaks = seq(min(lgas_det_probs$feral_tests, na.rm=TRUE),
               max(lgas_det_probs$feral_tests, na.rm=TRUE),
               length.out = 100)
  plot(st_geometry(lgas_det_probs[,"feral_expected_dets"]), 
       col = extreme_purps(ncol)[cut(lgas_det_probs$feral_tests, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Feral pigs", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", sub_cex_main)
  
  .image_scale(breaks, col = extreme_purps(ncol), 
               key.length = lcm(20), 
               key.width = lcm(0.5), # erm I checked the src and this doesn't appear
               key.pos = 4, cex.axis = leg_cex_axis,
               at = seq(0, 300, length.out = 7))
  mtext("Tests", 4, line = -4)
  
  # chickens
  breaks = seq(min(lgas_det_probs$chix_tests, na.rm=TRUE),
               max(lgas_det_probs$chix_tests, na.rm=TRUE),
               length.out = 100)
  plot(st_geometry(lgas_det_probs[,"chix_expected_dets"]), 
       col = extreme_purps(ncol)[cut(lgas_det_probs$chix_tests, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Sentinel chickens", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(e)", sub_cex_main)
  
  .image_scale(breaks, col = extreme_purps(ncol), 
               key.length = lcm(20), 
               key.width = lcm(0.5), # erm I checked the src and this doesn't appear
               key.pos = 4, cex.axis = leg_cex_axis,
               at = seq(0, 60, length.out = 7))
  mtext("Tests", 4, line = -4)
  
  
  dev.off()}



# total detections
{png("figures/detections_lga.png",
     width = 1000*reso/72, 
     height = 1100*reso/72, 
     res = reso)
  
  # can't stop this from being reset after each panel, but now have to match the reset
  # settings
  par(mfrow = c(3, 2),
      mar = c(0.1, 0.1, 1.2, 2.2), bty="n",
      oma=c(0,0,0,2))
  
  layout(matrix(1:12, ncol = 4, byrow = TRUE), widths = c(1, 0.08, 1, 0.08))
  
  # humans
  breaks = 0:12
  plot(st_geometry(lgas_det_probs[,"human_expected_dets"]), 
       col = extreme_purps(length(breaks))[cut(lgas_det_probs$human_dets, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Humans", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(a)", sub_cex_main)
  
  .image_scale(breaks, col = extreme_purps(length(breaks)), 
               key.length = lcm(20), 
               key.width = lcm(0.5), # erm I checked the src and this doesn't appear
               key.pos = 4, cex.axis = leg_cex_axis,
               at = seq(0, 12, length.out = 7))
  mtext("Detections", 4, line = -4)
  
  
  # mosquitoes
  breaks = 0:6
  plot(st_geometry(lgas_det_probs[,"mozz_expected_dets"]), 
       col = extreme_purps(length(breaks))[cut(lgas_det_probs$mozz_dets, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Mosquitoes", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(b)", sub_cex_main)
  
  .image_scale(breaks, col = extreme_purps(length(breaks)), 
               key.length = lcm(20), 
               key.width = lcm(0.5), # erm I checked the src and this doesn't appear
               key.pos = 4, cex.axis = leg_cex_axis,
               at = 0:6)
  mtext("Detections", 4, line = -4)
  
  # domestic pigs
  breaks = 0:14
  plot(st_geometry(lgas_det_probs[,"dom_expected_dets"]), 
       col = extreme_purps(length(breaks))[cut(lgas_det_probs$dom_dets, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Domestic pigs", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(c)", sub_cex_main)
  
  .image_scale(breaks, col = extreme_purps(length(breaks)), 
               key.length = lcm(20), 
               key.width = lcm(0.5), # erm I checked the src and this doesn't appear
               key.pos = 4, cex.axis = leg_cex_axis,
               at = seq(0, 14, length.out = 8))
  mtext("Detections", 4, line = -4)
  
  
  # feral pigs
  breaks = 0:62
  plot(st_geometry(lgas_det_probs[,"feral_expected_dets"]), 
       col = extreme_purps(length(breaks))[cut(lgas_det_probs$feral_dets, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Feral pigs", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(d)", sub_cex_main)
  
  .image_scale(breaks, col = extreme_purps(length(breaks)), 
               key.length = lcm(20), 
               key.width = lcm(0.5), # erm I checked the src and this doesn't appear
               key.pos = 4, cex.axis = leg_cex_axis,
               at = seq(0, 60, length.out = 4))
  mtext("Detections", 4, line = -4)
  
  # chickens
  breaks = 0:4
  plot(st_geometry(lgas_det_probs[,"chix_expected_dets"]), 
       col = extreme_purps(length(breaks))[cut(lgas_det_probs$chix_dets, breaks= breaks)],
       key.pos = NULL,
       reset = FALSE,
       cex.main=sub_cex_main,
       main="Sentinel chickens", lwd = 0.1, xlim = c(113, 154))
  plot(st_geometry(aus_simple), add=TRUE)
  
  subfigure_label(par("usr"), 0.1, 0.9, "(e)", sub_cex_main)
  
  .image_scale(breaks, col = extreme_purps(length(breaks)), 
               key.length = lcm(20), 
               key.width = lcm(0.5), # erm I checked the src and this doesn't appear
               key.pos = 4, cex.axis = leg_cex_axis,
               at = 0:4)
  mtext("Detections", 4, line = -4)
  
  
  dev.off()}




