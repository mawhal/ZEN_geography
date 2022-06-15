# To-do list for ZEN_geography repository

- [X] get someone to provide content for the README
- [X] check outputs of model comparison script
- [X] check outputs of data assembly script
- [X] check outputs of figures script
- [X] consider removing "2014" from script file names - year does not matter relative to purpose
- [X] add zenodo link to genetic analysis - see email from Marlene Jahnke
- [X] add pop gen data - see email from Marlene -  we made it sound like dryad would be the location, but maybe github is okay?
- [ ] check minor differences in output from lm's in `ZEN_2014_model_comparison_site_means_range_standardized.R`


### Zipped data file sent by Emmett
- [X] One file with Zostera Relative Growth rates does not seem to be used, or maybe I am missing something?  `ZEN_2011_ZRG_AllSites_Edit141102 copy.csv` -- NEVERMIND, FOUND IT
- [ ] *Should additional scripts be added? for example, PCA scripts, FCA scripts, other genetic analyses?*



### Model Comparison script
###### question about comment in code
> NOTE: INTEGRATE THIS FILE WITH "DEFINITIVE" SCRIPT. RENAME THAT ONE TOO.
 PUT ALL SCRIPTS INTO A LOGICAL NAMING CONVENTION MAKING CLEAR THEIR TEMPORAL SEQUENCE

- [X] What are definitive scripts? 
- [X] Does anything need to be renamed at this stage?
###### naming conventions
How are we feeling about naming conventions now that we are using git and these things can be tracked?
- [X] remove date suffixes from file names
- [X] modify plotting routines so they are reproducible - saving to .svg or pdf?
- [X] fix code that cannot run on its own, see below
> Missing data? No.
nrow(ZEN_2014_plot_49) # 980
sum(is.na(ZEN_2014_plot_49$ocean.code)) # 0
sum(is.na(ZEN_2014_plot_49$zPC1.env.global)) # 0
sum(is.na(ZEN_2014_plot_49$zPC2.env.global)) # 0
sum(is.na(ZEN_2014_plot_49$zPC3.env.global)) # 0
sum(is.na(ZEN_2014_plot_49$zFC1)) # 0
sum(is.na(ZEN_2014_plot_49$zFC2)) # 0
sum(is.na(ZEN_2014_plot_49$zPC1.zos)) # 0
sum(is.na(ZEN_2014_plot_49$zPC2.zos)) # 0
sum(is.na(ZEN_2014_plot_49$zperiphyton.perg)) # 0
sum(is.na(ZEN_2014_plot_49$zmeso.mass.perg)) # 0

- [X] Note that model comparison script used `ZEN_2014_site_means_Atlantic` instead of `ZEN_2014_site_means_49_Atlantic` which was not part of the zipped directory Whalen received. Code did not run, so Whalen updated to `ZEN_2014_site_means_49_Atlantic`
- [X] Whalen does not have column of data frame for `leaf.CN.ratio.site`

### Data Assembly script
###### PACKAGES
We might want shift away from ddply in plyr to just using tidyverse
- [X] replace `read.csv()` with `read_csv()` - DECIDED NOT TO DO THIS
- [X] randomForest --- was data actually imputed? make this clear early in the code

###### DATA SOURCES
- [ ] **Clarify input file origins.** Links to other repos (github/doi), papers (doi), or people? this could all be listed in README.md
- [X] * for Whalen - add a repo for the code that made that environmental data

