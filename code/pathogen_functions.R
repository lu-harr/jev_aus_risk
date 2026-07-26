# PAPER VERSION

# FUNCTION BANK FOR PATHOGEN MODEL
# INCLUDES DATA HANDLING AND FITTING
library(dismo)
library(mgcv)
library(MASS)


ranking_to_prob <- function(ranking, max_rank = NULL, increasing = TRUE) {
  if (is.null(max_rank)) {
    max_rank = max(ranking)
  }
  
  if (!increasing) {
    ranking <- max_rank - ranking + 1
  }
  
  ranking / max_rank
}


# clip a raster to only the region within max_distance of point_coords
clip_to_distance <- function(raster, point_coords, max_distance) {
  dist <- raster::distanceFromPoints(raster, point_coords)
  dist[dist > max_distance] <- NA
  dist <- dist * 0
  dist <- mask(dist, raster)
  mask(raster, dist)
}


prepare_pathogen_inputs_national <- function(human_df,
                                             mosquito_df,
                                             dompigs_df,
                                             feral_df,
                                             chickens_df,
                                             ras,
                                             hpop,
                                             hums_lon="longitude", 
                                             hums_lat="latitude",
                                             mozz_bandwidth_scalar=1,
                                             dpig_bandwidth_scalar=1,
                                             feral_bandwidth_scalar=1,
                                             chickens_bandwidth_scalar=1
){
  
  # HUMANS - Lucy's version passed truncated coords and sampled points first
  # probably somewhat ruins any prospect of bootstrapping
  # I have already done this step in sort_pathogen_data
  human_case_points <- human_df %>% 
    dplyr::select('longitude',
                  'latitude') %>% 
    as.matrix(.)
  
  # need to assign the weights according to this...
  idx_human_cases <- cellFromXY(ras, human_case_points) 
  
  human_cases_df <- rbind(
    data.frame(
      idx = idx_human_cases,
      count = 1,
      weight = human_df$exposure_weight),
    data.frame(
      idx = 1:ncell(ras),
      count = 0,
      weight = 0)) %>% 
    group_by(idx) %>% 
    summarise(count=sum(count),
              weight=sum(weight)) %>%
    mutate(weight=case_when(
      count==0 ~ 1,
      TRUE ~ weight/count
    )) %>% 
    arrange(idx) 
  
  # MOSQUITOES - given coords -------------------------------------------
  # duplicates coordinates according to sample size - i.e. mosquito per row
  # first calculate average because sample size is not available for many
  mean_sample <- round(mean(mosquito_df$n_individuals, na.rm=TRUE), 0)
  
  mosquito_df_expanded <- all_mozzies %>%
    mutate(n_individuals=case_when(
      is.na(n_individuals) ~ mean_sample,
      TRUE ~ round(n_individuals, 0)
    )) %>% 
    uncount(n_individuals)
  
  # mosquito_df_expanded <- data.frame()
  # 
  # # there's a much faster way to do this with uncount
  # # see nat pathogen spatial model script!!
  # for (i_row in 1:nrow(mosquito_df)) {
  #   
  #   # get no. of rows
  #   n_rows <- mosquito_df[i_row,]$n_individuals
  #   
  #   if(is.na(n_rows)){n_rows <- mean_sample}
  #   
  #   info_mat <- mosquito_df[rep(i_row, n_rows),]
  #   
  #   mosquito_df_expanded <- rbind(mosquito_df_expanded, info_mat)
  #   
  # }
  
  tmp <- rasterize(mosquito_df_expanded[,c('longitude', 'latitude')], ras, background=0, fun='count')
  gf <- focalWeight(ras, mozz_bandwidth_scalar, "Gauss")
  mozz_density <- raster::focal(tmp, gf, pad=TRUE, na.rm=TRUE)
  mozz_density = mask(mozz_density, ras)
  
  mosquito_pos_vec = ifelse(!is.na(mosquito_df$virus_detection) & 
                              mosquito_df$virus_detection == "JEV positive", 1, 0)
  
  # DOMESTIC PIGS -------------------------------------------------------
  # do the same for the domestic pigs as mosquito
  tmp <- rasterize(dompigs_df[,c('longitude', 'latitude')], ras, background=0, fun='count')
  gf <- focalWeight(ras, dpig_bandwidth_scalar, "Gauss")
  dpig_density <- raster::focal(tmp, gf, pad=TRUE, na.rm=TRUE)
  dpig_density = mask(dpig_density, ras)
  
  dpig_pos_vec = ifelse(dompigs_df$virus_detection == "JEV positive", 1, 0)
  
  # FERALS --------------------------------------------------------------
  # do the same for the ferals as mosquito
  tmp <- rasterize(feral_df[,c('longitude', 'latitude')], ras, background=0, fun='count')
  gf <- focalWeight(ras, feral_bandwidth_scalar, "Gauss")
  feral_density <- raster::focal(tmp, gf, pad=TRUE, na.rm=TRUE)
  feral_density = mask(feral_density, ras)
  
  feral_pos_vec = ifelse(feral_df$virus_detection == "JEV positive", 1, 0)
  
  # Chickens --------------------------------------------------------------
  # do the same for the chickens as ferals
  tmp <- rasterize(chickens_df[,c('longitude', 'latitude')], ras, background=0, fun='count')
  gf <- focalWeight(ras, chickens_bandwidth_scalar, "Gauss")
  chicken_density <- raster::focal(tmp, gf, pad=TRUE, na.rm=TRUE)
  chicken_density = mask(chicken_density, ras)
  
  chicken_pos_vec = ifelse(chickens_df$virus_detection == "JEV positive", 1, 0)
  
  # combine densities
  
  coordinates <- stack(ras, ras)
  names(coordinates) <- c("x", "y")
  coordinates[] <- coordinates(ras)
  coordinates <- mask(coordinates, ras)
  
  densities <- stack(hpop,
                     mozz_density,
                     dpig_density,
                     feral_density,
                     chicken_density) 
  
  names(densities) <- c("human", "mosquito_sampling", "domestic_pig_sampling", "feral_pig_sampling", "sentinel_chicken_sampling")
  
  # PUTTING IT ALL TOGETHER!
  non_na_idx <- which(!is.na(getValues(ras)))  # index to non-NA cells
  coordinates_df <- as.data.frame(coordinates[non_na_idx])
  
  human_final_df <- cbind(
    cases = human_cases_df$count[non_na_idx], 
    density = densities$human[non_na_idx],
    coordinates_df,
    host_type = "human",
    exposure_weight = human_cases_df$weight[non_na_idx],
    diagnostic_type = "PCR",
    diagnostic_detail = "PCR",
    virus_detection = "JEV positive"
  )
  
  # mosquitoes: we have sampling points: only supply coords at sample/trap locs!
  # To do: Add weight-based imputation for WA mozzie samples size
  coordinates_df_mosquitoes = mosquito_df[,c("longitude","latitude")]
  names(coordinates_df_mosquitoes) = c("x","y")
  mosquito_sample_sizes = mosquito_df$n_individuals
  mosquito_sample_sizes[mosquito_sample_sizes == 0] = NA
  mosquito_sample_sizes[is.na(mosquito_sample_sizes)] = mean(mosquito_sample_sizes, na.rm=TRUE) # ish ...
  
  mosquito_final_df <- cbind(
    cases = mosquito_pos_vec,
    density = mosquito_sample_sizes, 
    coordinates_df_mosquitoes,
    host_type = "mosquito",
    exposure_weight = mosquito_df$exposure_weight,
    diagnostic_type = mosquito_df$diagnostic_type,
    diagnostic_detail = mosquito_df$diagnostic_detail,
    virus_detection = mosquito_df$virus_detection
  )
  
  ####
  
  # domestic pigs
  coordinates_df_dpig <- dompigs_df[,c("longitude","latitude")]
  names(coordinates_df_dpig) <- c("x","y")
  
  dpig_sample_sizes = rep(1, nrow(coordinates_df_dpig))
  
  # sample size change to actual number of pigs at sampled locations but not necessary
  
  dpig_final_df <- cbind(
    cases = dpig_pos_vec,
    density = dpig_sample_sizes,
    coordinates_df_dpig,
    host_type = "domestic_pig",
    exposure_weight = dompigs_df$exposure_weight,
    diagnostic_type = dompigs_df$diagnostic_type,
    diagnostic_detail = dompigs_df$diagnostic_detail,
    virus_detection = dompigs_df$virus_detection
  )
  
  # feral pigs  
  coordinates_df_feral <- feral_df[,c("longitude","latitude")]
  names(coordinates_df_feral) <- c("x","y")
  
  feral_sample_sizes = rep(1, nrow(coordinates_df_feral))
  
  # feral pigs
  feral_final_df <- cbind(
    cases = feral_pos_vec,
    density = feral_sample_sizes,
    coordinates_df_feral,
    host_type = "feral_pig",
    exposure_weight = feral_df$exposure_weight,
    diagnostic_type = feral_df$diagnostic_type,
    diagnostic_detail = feral_df$diagnostic_detail,
    virus_detection = feral_df$virus_detection
  )
  
  # sentinel chickens
  
  coordinates_df_chickens <- chickens_df[,c("longitude","latitude")]
  names(coordinates_df_chickens) <- c("x","y")
  
  chicken_sample_sizes = rep(1, nrow(coordinates_df_chickens))
  
  chicken_final_df <- cbind(
    cases=chicken_pos_vec,
    density=chicken_sample_sizes,
    coordinates_df_chickens,
    host_type = "chicken",
    exposure_weight = chickens_df$exposure_weight,
    diagnostic_type = chickens_df$diagnostic_type,
    diagnostic_detail = chickens_df$diagnostic_detail,
    virus_detection = chickens_df$virus_detection
  )
  
  # combine into 1 dataframe, with the host species as a covariate
  df <- rbind(human_final_df, mosquito_final_df, dpig_final_df, 
              feral_final_df, chicken_final_df)
  host_type_id_order <- unique(df$host_type)
  df$host_type_id <- match(df$host_type, host_type_id_order)
  df$host_type_id <- factor(df$host_type_id)
  
  # 2024: diagnostic_type_id is missing?
  diagnostic_type_id_order <- unique(df$diagnostic_type)
  df$diagnostic_type_id <- match(df$diagnostic_type, diagnostic_type_id_order)
  df$diagnostic_type_id <- factor(df$diagnostic_type_id)
  
  return(list(df=df,
              human_points = human_case_points,
              mozz_points = mosquito_final_df,
              dpig_points = dpig_final_df,
              feral_points = feral_final_df,
              chicken_points = chicken_final_df,
              coordinates = coordinates,
              densities = densities,
              host_type_order = host_type_id_order))
}



fit_pathogen_spatial_model = function(df, pred_coordinates){
  # (should be able to rbind runs of prepare_spatial_pathogen_inputs to fit
  # model to multiple states at once)
  
  # fit a spatial GAM to all three host species data simultaneously
  # with difference intercepts for PCR and serology by species where relevant
  m <- bam(cases ~ s(x, y, k=100) + 
             host_type_id  +
             #s(host_type_id, bs="re") + 
             diagnostic_type_id +
             offset(log(density)), 
           data = df,
           family = poisson,
           weights = likelihood_weight)
  
  # predict back ()
  pred_covs <- stack(
    pred_coordinates, 
    pred_coordinates[[1]] * 0 + 1, #adding in human incidence
    pred_coordinates[[1]] * 0 + 1, # add PCR as the diagnostic type
    pred_coordinates[[1]] * 0 + 1 # density?
  )
  names(pred_covs) <- c(names(pred_coordinates), "host_type_id", "diagnostic_type_id", "density")
  pred_inc_human <- predict(pred_covs, m, type = "response")
  pred_inc <- exp(log(pred_inc_human)-m$coefficients[1])
  pred_rel_inc <- pred_inc / cellStats(pred_inc, max)
  
  return(list(pred_inc=pred_inc,
              pred_rel_inc=pred_rel_inc,
              model=m))
}


extend_raster <- function(r, to){
  
  r[is.na(r)] <- cellStats(r, min)
  r <- mask(r, to)
  return(r)
  
}

common_sense_check <- function(df, m){
  # predict the model back on the data we gave it and check it gives us back the
  # right number of detections for species groups with only one diagnostic type
  # (i.e. a free parameter under fixed effect model)
  df %>%
    mutate(pred = as.vector(predict(m, newdata = ., type = "response"))) %>%
    group_by(host_type, diagnostic_type) %>%
    summarise(ncases = sum(cases),
              ncases_wt = sum(cases*likelihood_weight),
              pred_cases = sum(pred),
              pred_cases_wt = sum(pred*likelihood_weight),
              total_weights = sum(likelihood_weight),
              total_density = sum(density),
              nsamples = n(),
              rat = pred_cases_wt / ncases_wt) %>%
    mutate(across(ncases_wt:pred_cases_wt, round))
}

##############################################################
# POPULATION AT RISK CODE

# given continuous rasters of the potential niche and the pathogen distribution,
# rasters of the populations of people and pigs in each cell, a *raster* of the
# states to which each cell belongs, a two-column matrix of the coordinates of
# the pathogen presence points, the vector of inclusion quantile thresholds to
# consider and the blob distance, return a list containing a dataframe of the
# populations (people and pigs) at risk in each state, in the potential and
# realised niches at each threshold, and the binary rasters of the potential and
# realised niches, and the pathogen distribution, at each thrseshold.
all_populations_at_risk <- function(potential,
                                    pathogen,
                                    people,
                                    pigs,
                                    states_raster,
                                    pres_points,
                                    inclusion_quantiles = c(0.5, 0.9, 1),
                                    blob_distance = 0.5) {
  
  res <- list()
  
  for (quantile in inclusion_quantiles) {
    
    threshold_name <- paste0(100 * quantile, "%")
    
    par_people <- population_at_risk(potential = potential,
                                     pathogen = pathogen,
                                     population = people,
                                     states_raster = states_raster,
                                     pres_points = pres_points,
                                     inclusion_quantile = quantile,
                                     blob_distance = blob_distance)
    
    par_pigs <- population_at_risk(potential = potential,
                                   pathogen = pathogen,
                                   population = pigs,
                                   states_raster = states_raster,
                                   pres_points = pres_points,
                                   inclusion_quantile = quantile,
                                   blob_distance = blob_distance)
    
    par_table_people <- par_people$populations_at_risk %>%
      mutate(species = "people")
    par_table_pigs <- par_pigs$populations_at_risk %>%
      mutate(species = "pigs")
    
    par_table <- bind_rows(
      par_table_people,
      par_table_pigs
    ) %>%
      relocate(
        species,
        .after = state
      ) %>%
      mutate(
        threshold = threshold_name,
        .after = niche
      )
    
    res[[threshold_name]] <- list(
      population_at_risk = par_table,
      potential_niche = par_people$potential_niche,
      realised_niche = par_people$realised_niche,
      pathogen_intensity = par_people$pathogen_intensity
    )
    
  }
  
  # reorder results
  tables <- lapply(res, `[[`, "population_at_risk")
  table <- do.call(rbind, tables)
  
  # combine rasters
  potential_niches <- lapply(res, `[[`, "potential_niche")
  realised_niches <- lapply(res, `[[`, "realised_niche")
  pathogen_intensities <- lapply(res, `[[`, "pathogen_intensity")
  
  list(
    population_at_risk = table,
    potential_niches = potential_niches,
    realised_niches = realised_niches,
    pathogen_intensities = pathogen_intensities
  )
  
}



# given continuous rasters of the potential niche and the pathogen distribution,
# a raster of the population in each cell, a *raster* of the states to which
# each cell belongs, a two-column matrix of the coordinates of the pathogen
# presence points, the inclusion quantile threshold and blob distance, return a
# list containing a dataframe of the populations at risk in each state in the
# potential and realised niches at this threshold, and the rasters of the
# potential and realised niches, and the pathogen distribution.
population_at_risk <- function(potential,
                               pathogen,
                               population,
                               states_raster,
                               pres_points,
                               inclusion_quantile = 0.5,
                               blob_distance = 0.5) {
  
  # blobbify the potential and pathogen layers at different thresholds, and intersect them
  potential_thresh <- threshy_blobs(potential,
                                    pres_points,
                                    inclusion_quantile = inclusion_quantile,
                                    blob_distance = blob_distance)
  pathogen_thresh <- threshy_blobs(pathogen,
                                   pres_points,
                                   inclusion_quantile = inclusion_quantile,
                                   blob_distance = blob_distance)
  realised_thresh <- potential_thresh * pathogen_thresh
  
  potential_par <- zonal(potential_thresh * population, states, fun = "sum")
  realised_par <- zonal(realised_thresh * population, states, fun = "sum") 
  total_par <- zonal(population, states, fun='sum')
  
  par_table <- data.frame(
    state = potential_par[, "zone"],
    potential = potential_par[, "sum"],
    realised = realised_par[, "sum"],
    total = total_par[, 'sum']
  ) %>%
    pivot_longer(
      cols = c("potential", "realised", "total"),
      names_to = "niche",
      values_to = "population_at_risk"
    )
  
  list(
    populations_at_risk = par_table,
    potential_niche = potential_thresh,
    realised_niche = realised_thresh,
    pathogen_intensity = pathogen_thresh
  )
  
}

# add factor levels to a raster layer (by default, those for potential and realised niche)
add_levels <- function(raster, names = c("outside niche", "potential niche", "realised niche")) {
  raster <- as.factor(raster)
  lvl <- levels(raster)[[1]]
  lvl$code <- names
  levels(raster) <- lvl
  raster
}

