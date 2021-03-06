![ZEN logo](http://zenscience.org/wp-content/uploads/2011/09/Zen_header_logo_50_pct.png)
# ZEN_geography
###### this repository contains data, code, and various outputs associated with the _Zostera_ Experimental Network's evolutionary legacy project

## Please contact us if you'd like to use the data

This project considers the ecosystem-wide consequences of evolutionary legacies associated with the establishment and history of _Zostera marina_ (eelgrass) in the Atlantic Ocean relative to its Pacific origins. Atlantic eelgrass is smaller and shorter than Pacific eelgrass, forming "meadows" rather than "forests", and these morphological features are a legacy of population genetic bottlenecks and genetic differentiation. In turn, differences in the stature and biomass of eelgrass influence the biomass of associated algae and invertebrate animals
# 

##### citation:
Duffy, J.E., et al. 2022. A Pleistocene legacy structures variation in modern seagrass ecosystems. Proceedings of the National Academy of Sciences of the United States of America. https://https://doi.org/10.1073/pnas.2121425119

contact: DuffyE@si.edu

### Guide to this repository
`data/` contains 
- `input/` contains data from the project and previous projects used as inputs in scripts located in `code/`
- `output/` contains data and results from the R scripts located in `code/`
- `popgen/` contains eelgrass microsatellite data
- `evol_envir_path_calculations/` contains an excel file (all sheets saved as .csv files) used to calculate and organize direct and indirect paths for comparing effects of evolutionary history and environment on eelgrass ecosystem components

`code/` contains R code associated with the project listed in order of running:
- `environmental_vars_precip.R` - gathers WorldClim precipitation data.
- `environmental_biooracle_merge.R` - gathers BioOracle data variables (n=24) and merges with precipipation data
- `data_assembly.R` - prepares data for investigation (summaries, figures, statistical models)
- `model_comparison.R` - constructs and compares linear models of site-level data
- `figures.R` - generates figures of data and model outputs
- `impute_missing.R` - imputes missing data using random forest models. This was primarily used to impute _Zostera_ morphological variables prior to a principal components analysis

`ZEN_geography.Rproj` is an R project file through which all R code should be run, noting that file paths called in `code` begin in the parent directory

### NOTES
Genetic analyses are not included in this repository. For analysis of genetic diversity spectra and isolation by distance, see https://zenodo.org/record/3660013#.YWF8C0bMLpJ
