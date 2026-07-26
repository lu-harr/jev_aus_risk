# PAPER VERSION

# FUNCTION BANK FOR PATHOGEN MODEL
# INCLUDES DATA HANDLING AND FITTING
library(dismo)
library(mgcv)
library(MASS)

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
                                             ras,
                                             hpop,
                                             mozz_bandwidth_scalar=1,
                                             dpig_bandwidth_scalar=1,
                                             feral_bandwidth_scalar=1){
  
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
      weight = human_df$weight),
    data.frame(
      idx = 1:ncell(ras),
      count = 0,
      weight = 0)) %>% 
    group_by(idx) %>% 
    summarise(count=sum(count),
              weight=sum(weight)) %>%
    mutate(weights=case_when(
      count==0 ~ 1,
      TRUE ~ weight/count
    )) %>% 
    arrange(idx) 
  
  # MOSQUITOES - given coords
  # duplicates coordinates according to sample size - i.e. mosquito per row
  # first calculate average because sample size is not available for many
  mean_sample <- round(mean(mosquito_df$n_individuals, na.rm=TRUE), 0)
  
  # swapped in uncount method ...
  mosquito_df_expanded <- mosquito_df %>%
    mutate(n_individuals=case_when(
      is.na(n_individuals) ~ mean_sample,
      TRUE ~ round(n_individuals, 0)
    )) %>% 
    uncount(n_individuals)
  
  tmp <- rasterize(mosquito_df_expanded[,c('longitude', 'latitude')], ras, background=0, fun='count')
  gf <- focalWeight(ras, mozz_bandwidth_scalar, "Gauss")
  mozz_density <- raster::focal(tmp, gf, pad=TRUE, na.rm=TRUE)
  mozz_density = mask(mozz_density, ras)
  
  mosquito_pos_vec = ifelse(!is.na(mosquito_df$virus_detection) & 
                              mosquito_df$virus_detection == "JEV positive", 1, 0)
  
  # DOMESTIC PIGS
  # do the same for the domestic pigs as mosquito
  tmp <- rasterize(dompigs_df[,c('longitude', 'latitude')], ras, background=0, fun='count')
  gf <- focalWeight(ras, dpig_bandwidth_scalar, "Gauss")
  dpig_density <- raster::focal(tmp, gf, pad=TRUE, na.rm=TRUE)
  dpig_density = mask(dpig_density, ras)
  
  dpig_pos_vec = ifelse(dompigs_df$virus_detection == "JEV positive", 1, 0)
  
  # FERALS
  # do the same for the ferals as mosquito
  tmp <- rasterize(feral_df[,c('longitude', 'latitude')], ras, background=0, fun='count')
  gf <- focalWeight(ras, feral_bandwidth_scalar, "Gauss")
  feral_density <- raster::focal(tmp, gf, pad=TRUE, na.rm=TRUE)
  feral_density = mask(feral_density, ras)
  
  feral_pos_vec = ifelse(feral_df$virus_detection == "JEV positive", 1, 0)
  
  coordinates <- stack(ras, ras)
  names(coordinates) <- c("x", "y")
  coordinates[] <- coordinates(ras)
  coordinates <- mask(coordinates, ras)
  
  densities <- stack(hpop,
                     mozz_density,
                     dpig_density,
                     feral_density)
  
  names(densities) <- c("human", "mosquito_sampling", "domestic_pig_sampling", "feral_pig_sampling")
  
  # PUTTING IT ALL TOGETHER!
  non_na_idx <- which(!is.na(getValues(ras)))  # index to non-NA cells
  coordinates_df <- as.data.frame(coordinates[non_na_idx])
  
  human_final_df <- cbind(
    cases = human_cases_df$count[non_na_idx], # check this line....I'm not sure what it's doing
    density = densities$human[non_na_idx],
    coordinates_df,
    host_species = "human",
    weight = human_cases_df$weights[non_na_idx]
  )
  
  # mosquitoes: we have sampling points: only supply coords at sample/trap locs!
  coordinates_df_mosquitoes = mosquito_df[,c("longitude","latitude")]
  names(coordinates_df_mosquitoes) = c("x","y")
  mosquito_sample_sizes = mosquito_df$n_individuals
  mosquito_sample_sizes[mosquito_sample_sizes == 0] = NA
  mosquito_sample_sizes[is.na(mosquito_sample_sizes)] = mean(mosquito_sample_sizes, na.rm=TRUE) # ish ...
  
  mosquito_final_df <- cbind(
    cases = mosquito_pos_vec,
    density = mosquito_sample_sizes, 
    coordinates_df_mosquitoes,
    host_species = "mosquito",
    weight = 1
  )
  
  ####
  
  # domestic pigs
  coordinates_df_dpig <- dompigs_df[,c("longitude","latitude")]
  names(coordinates_df_dpig) <- c("x","y")
  
  dpig_sample_sizes = rep(1, nrow(coordinates_df_dpig))
  
  dpig_final_df <- cbind(
    cases = dpig_pos_vec,
    density = dpig_sample_sizes,
    coordinates_df_dpig,
    host_species = "domestic_pig",
    weight = 1
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
    host_species = "feral_pig",
    weight = 1
  )
  
  # combine into 1 dataframe, with the host species as a covariate
  df <- rbind(human_final_df, mosquito_final_df, dpig_final_df, feral_final_df)
  host_species_id_order <- unique(df$host_species)
  df$host_species_id <- match(df$host_species, host_species_id_order)
  df$host_species_id <- factor(df$host_species_id)
  
  return(list(df=df,
              human_points = human_case_points,
              mozz_points = mosquito_final_df,
              dpig_points = dpig_final_df,
              feral_points = feral_final_df,
              coordinates = coordinates,
              densities = densities,
              host_species_order = host_species_id_order))
}



fit_pathogen_spatial_model = function(df, pred_coordinates){
  
  # fit a spatial GAM to all ~three~ four host species simultaneously
  m <- bam(cases ~ s(x, y, k=100) + host_species_id,
           data = df,
           offset = log(density),
           family = poisson,
           weights = weight)
  
  # predict back (adding in human incidence as a covariate for prediction, since it needs one)
  pred_covs <- stack(pred_coordinates, pred_coordinates[[1]] * 0 + 1)
  names(pred_covs) <- c(names(pred_coordinates), "host_species_id")
  pred_inc_human <- predict(pred_covs, m, type = "response")
  pred_inc <- exp(log(pred_inc_human)-m$coefficients[1]) 
  # why is the intercept subtracted off here? Because that's beta_human
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



