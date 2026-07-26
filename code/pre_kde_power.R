# run the top of pathogen_model.R
# rel_intensity_cutoff <- 0.029*50

# probabilities of observing one or more detections in each host species at each
# location
detection_prob <- function(expected_detections) {
  1 - dpois(0, expected_detections)
}

mozz_tests <- dat$mozz_points %>%
  dplyr::select(x, y) %>%
  rasterize(hpop_blur, 
            background = 0, 
            field = dat$mozz_points$density*dat$mozz_points$exposure_weight, 
            fun = sum)

chix <- dat$df %>%
  filter(host_type == "chicken")
chix_tests <- chix %>%
  dplyr::select(x, y) %>%
  rasterize(hpop_blur,
            background = 0, field = chix$diagnostic_weight, fun = sum)

dpigs <- dat$df %>%
  filter(host_type == "domestic_pig")
dpig_tests <- dpigs %>%
  dplyr::select(x, y) %>%
  rasterize(hpop_blur,
            background = 0, field = dpigs$diagnostic_weight*dpigs$exposure_weight, 
            fun = sum)

fpigs <- dat$df %>%
  filter(host_type == "feral_pig")
fpig_tests <- fpigs %>%
  dplyr::select(x, y) %>%
  rasterize(hpop_blur,
            background = 0, field = fpigs$diagnostic_weight*fpigs$exposure_weight, 
            fun = sum)

#human_tests == hpop

tests <- stack(chix_tests, mozz_tests, fpig_tests, dpig_tests, hpop_blur)
names(tests) <- c("chix", "mozz", "fpig", "dpig", "hpop")

tests[tests == 0] <- NA



human_expected_detections <- exp(log(rel_intensity_cutoff) + 
                                   intercepts[species_id_lookup["human"]] +
                                   detection_contrasts[detection_id_lookup["PCR"]] +
                                   log(hpop_blur))

sum(values(human_expected_detections), na.rm=TRUE)

mossie_expected_detections <- exp(log(rel_intensity_cutoff) + 
                                    intercepts[species_id_lookup["mossie"]] + 
                                    detection_contrasts[detection_id_lookup["PCR"]] +
                                    log(tests$mozz))

sum(values(mossie_expected_detections), na.rm=TRUE)

dpig_expected_detections <- exp(log(rel_intensity_cutoff) + 
                                  intercepts[species_id_lookup["dpig"]] + 
                                  detection_contrasts[detection_id_lookup["PCR"]] +
                                  log(tests$dpig)) +
  exp(log(rel_intensity_cutoff) +
        intercepts[species_id_lookup["dpig"]] +
        detection_contrasts[detection_id_lookup["serology"]] +
        log(tests$dpig))+
  exp(log(rel_intensity_cutoff) +
        intercepts[species_id_lookup["dpig"]] +
        detection_contrasts[detection_id_lookup["unassigned"]] +
        log(tests$dpig))

sum(values(dpig_expected_detections), na.rm=TRUE)

fpig_expected_detections <- exp(log(rel_intensity_cutoff) + 
                                  intercepts[species_id_lookup["fpig"]] + 
                                  detection_contrasts[detection_id_lookup["PCR"]] +
                                  log(tests$fpig)) +
  exp(log(rel_intensity_cutoff) +
        intercepts[species_id_lookup["fpig"]] +
        detection_contrasts[detection_id_lookup["serology"]] +
        log(tests$fpig))+
  exp(log(rel_intensity_cutoff) +
        intercepts[species_id_lookup["fpig"]] +
        detection_contrasts[detection_id_lookup["unassigned"]] +
        log(tests$fpig))

sum(values(fpig_expected_detections), na.rm=TRUE)

chicken_expected_detections <- exp(log(rel_intensity_cutoff) +
                                     intercepts[species_id_lookup["chicken"]] +
                                     detection_contrasts[detection_id_lookup["serology"]] +
                                     log(tests$chix))

sum(values(chicken_expected_detections), na.rm=TRUE)

human_detection_probability <- calc(human_expected_detections, detection_prob)
mossie_detection_probability <- calc(mossie_expected_detections, detection_prob)
dpig_detection_probability <- calc(dpig_expected_detections, detection_prob)
fpig_detection_probability <- calc(fpig_expected_detections, detection_prob)
chicken_detection_probability <- calc(chicken_expected_detections, detection_prob)

maxValue(human_detection_probability)
maxValue(dpig_detection_probability)
maxValue(fpig_detection_probability)
maxValue(mossie_detection_probability)
maxValue(chicken_detection_probability)








