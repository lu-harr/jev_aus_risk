# run things top to bottom:
# (intended more as a road map than for actually running)

# leaving this script out of github repo as there are raw data in there:
# (e.g. where it was sent as raw text in an email)
source("secure_code/pathogen_p2_sort_surveillance_data.R")
# ensure PROMOTE_SEROS is set to TRUE
source("code/surveillance_table.R")
# makes up a nice .tex table

source("code/pathogen_model.R")
# check whether intensity model fitting is turned off:
# the fitted model used in the manuscript is provided as an .rds in output/

source("code/pathogen_model_inputs_figures.R")
# this is making five pictures, which is more than you need
# requires top half of `pathogen_model` to have been run: 
# uses `dat` object of prepared data for model

source("code/ecdfs.R")
# does the pink maps and supporting supp figures

source("code/niche_extents_population_at_risk.R")
# does the neapolitan map and supporting table

# there are lots of other scripts in the repo for supplementary bits and bobs

