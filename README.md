Geospatial assessment of Japanese encephalitis virus infection risk in
Australia
================

A repository of all code relating to the manuscript *Geospatial
assessment of Japanese encephalitis virus infection risk in Australia*.
This code was written by Lucinda Harrison, in collaboration with Freya
Shearer, David Duncan, Nick Golding, Jennifer Flegg, and David Price.
The data associated with this manuscript are not for circulation and are
not attached in this repository.

Start at `top_to_bottom.R`:

``` r
source("code/pathogen_p2_sort_surveillance_data.R")
# ensure PROMOTE_SEROS is set to TRUE

source("code/pathogen_model.R")
# check whether intensity model fitting is turned off:
# the fitted model used in the manuscript is provided as an .rds in output/

source("code/pathogen_model_inputs_figure.R")
# this is making five pictures, which is more than you need
# requires top half of `pathogen_model` to have been run: 
# uses `dat` object of prepared data for model

source("code/ecdfs.R")
# does the pink maps and supporting supp figures

source("code/niche_extents_population_at_risk.R")
# does the neapolitan map and supporting table
```

(Probably don’t run this script but do use it as a guide to navigate the
repo)

An older iteration of the model (prior to the incorporation of sentinel
chicken data and diagnostic uncertainty) is included in `model_2022`.

#### Acknowledgements

Populate when they’re in the manuscript.
