# ZEN_geography
###### this repository contains data, code, and various outputs associated with the _Zostera_ Experimental Network's evolutionary legacy project
This project considers the ecosystem-wide consequences of evolutionary legacies associated with the establishment and history of _Zostera marina_ (eelgrass) in the Atlantic Ocean relative to its Pacific origins. Atlantic eelgrass is smaller and shorter than Pacific eelgrass, forming "meadows" rather than "forests", and these morphological features are a legacy of population genetic bottlenecks and genetic differentiation. In turn, differences in the stature and biomass of eelgrass influence the biomass of associated algae and invertebrate animals
# 
![ZEN logo](http://zenscience.org/wp-content/uploads/2011/09/Zen_header_logo_50_pct.png)
##### citation:
Duffy, J.E., et al. 2022. A Pleistocene legacy structures variation in modern seagrass ecosystems. Proceedings of the National Academy of Sciences of the United States of America. **DOI**

### Guide to this repository
`data/` contains 
- `input/` contains data from the project and previous projects
- `output/` contains data the result from the R scripts located in `code/`
`code/` contains R code associated with the project listed in order of running:
- `environmental_vars_precip` - gathers WorldClim precipitation data.
- `environmental_vars_precip` - gathers BioOracle data variables (n=24)
- `data_assembly` - prepares data for investigation (summaries, figures, statistical models)
- `model_comparison` - constructs and compares linear models of site-level data
- `figures` - generates figures of data and model outputs

### NOTES
Genetic analyses are not included in this repository. For analysis of genetic diversity spectra and isolation by distance, see https://zenodo.org/record/3660013#.YWF8C0bMLpJ
