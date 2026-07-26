library(xtable)

source("init.R")
source("pathogen_functions.R")

states <- raster('data/state_shp/state_ras_agg.tif')
lgas <- st_read("~/Desktop/jev/data/admin/")

potential_continuous <- raster('output/continuous suit vectors and avian.tif') %>%
  aggregate(fact=AGG_FACTOR)
realised_continuous <- raster('output/continuous_pathogen_national_cut.tif')

people <- raster('data/hpop_agg_10km.tif')
pigs <- raster('data/national_piggery_density.tif')

people_corrected <- mask(people, states)

pres_points <- read_csv('data/pathogen_presence_national_cut.csv') %>% 
  dplyr::select(
    'longitude', 
    'latitude'
  )

# calculate populations at risk
blob_distance <- 0.1

# this is the ABS estimated national population for 31/03/22 ..
people_corrected <- people_corrected*25890773/cellStats(people_corrected, 'sum')

uncorrected_totals <- zonal(people_corrected, states, fun='sum') 

# compute niches and populations at risk at different thresholds
par <- all_populations_at_risk(potential = potential_continuous,
                               pathogen = realised_continuous,
                               people = people_corrected,
                               pigs = pigs,
                               states_raster = states,
                               pres_points = pres_points,
                               inclusion_quantiles = c(0.5, 0.75, 0.95),
                               blob_distance = blob_distance)

tab <- par$population_at_risk %>%
  filter(species=='people') %>% 
  relocate(
    species,
    .before = everything()
  ) %>%
  # mutate(
  #   population_at_risk = format(
  #     round(population_at_risk/1000),
  #     big.mark = ",",
  #     scientific = FALSE
  #   )
  # ) %>%
  arrange(state,
          species) %>%
  pivot_wider(
    names_from = c(niche, threshold),
    values_from = population_at_risk
  ) %>% 
  mutate(state=case_when(
    state==1 ~ 'NSW',
    state==2 ~ 'VIC',
    state==3 ~ 'QLD',
    state==4 ~ 'SA',
    state==5 ~ 'WA',
    state==6 ~ 'TAS',
    state==7 ~ 'NT',
    state==8 ~ 'ACT',
  )) %>% 
  # mutate(species=case_when(
  #   species=='pigs' ~ 'domestic pigs',
  #   TRUE ~ species
  # )) %>% 
  mutate('outside_50%'=`total_95%` - `potential_50%`,
         'outside_75%'=`total_95%` - `potential_75%`,
         'outside_95%'=`total_95%` - `potential_95%`) %>%
  dplyr::select(-'total_50%', -'total_75%', -'species') %>% 
  # relocate('total_95%',
  #          .after = state) %>%
  # relocate('outside_50%',
  #          .after = `total_95%`) %>%
  # relocate('outside_75%',
  #          .after = `realised_50%`) %>%
  # relocate('outside_95%',
  #          .after = `realised_75%`) %>%
  dplyr::select("state", "total_95%", "realised_50%", "potential_50%", "outside_50%",
                "realised_75%", "potential_75%", "outside_75%",
                "realised_95%", "potential_95%", "outside_95%") %>%
  mutate(across(!state, function(x){
    format(
      round(x/1000),
      big.mark = ",",
      scientific = FALSE
    )
  })) # this is a bit cheeky :)

print(xtable::xtable(tab),
      include.rownames=FALSE)

# get the intersections of these niches, codes as 0 for outside niche, 1 for
# potential (but unrealised), and 2 for realised
intersections <- mapply(`+`, par$potential_niches, par$realised_niches)

# code these up as factor rasters
intersections <- lapply(
  intersections,
  add_levels,
  names = c("outside niche", "potential niche", "realised niche")
)

###############################################################################
# NICHE EXTENTS FIGURE
# edit mfrow, only one legend, mar

{png("figures/niches_thresholds_avian_vectors.png",
     width = 900*reso/72, 
     height = 900*reso/72, 
     res = reso
)
  
  par(mfrow = c(2, 2),
      mar = c(0.5, 0.5, 0.5, 0.5), bty="n",
      oma=c(0,0,0,0))
  
  
  for(i in 1:length(names(intersections))) {
    plot(intersections[[names(intersections)[i]]], col=cols, legend=FALSE,
         xaxt="n", yaxt="n", legend.mar=0)
    plot(aus_simple, lwd=1, border='grey20', add=TRUE)
    subfigure_label(par("usr"), 0.1, 0.9, c("(a)","(b)","(c)")[i], 1.6)
    
  }
  
  empty_plot_for_legend()
  legend("center", legend = rev(c("realised niche","potential niche","outside niche")),
         fill = cols, cex=1.4)
  
  dev.off()}
