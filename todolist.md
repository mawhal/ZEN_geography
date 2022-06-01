# To-do list for ZEN_geography repository

- [ ] get someone to provide content for the README
- [X] check outputs of model comparison script
- [ ] check outputs of data assembly script
- [ ] check outputs of figures script
- [ ] consider removing "2014" from script file names - year does not matter relative to purpose
- [ ] add data/results tables for genetic analysis - see email from Marlene Jahnke



### Zipped data file sent by Emmett
- [X] One file with Zostera Relative Growth rates does not seem to be used, or maybe I am missing something?  `ZEN_2011_ZRG_AllSites_Edit141102 copy.csv` -- NEVERMIND, FOUND IT
- [ ] Should additional scripts be added? for example, PCA scripts, FCA scripts, other genetic analyses?



### Model Comparison script
###### question about comment in code
> NOTE: INTEGRATE THIS FILE WITH "DEFINITIVE" SCRIPT. RENAME THAT ONE TOO.
 PUT ALL SCRIPTS INTO A LOGICAL NAMING CONVENTION MAKING CLEAR THEIR TEMPORAL SEQUENCE

- [ ] What are definitive scripts? 
- [ ] Does anything need to be renamed at this stage?
###### naming conventions
How are we feeling about naming conventions now that we are using git and these things can be tracked?
- [ ] remove date suffixes from file names
- [ ] modify plotting routines so they are reproducible - saving to .svg or pdf?
- [ ] fix code that cannot run on its own, see below
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
- [ ] replace `read.csv()` with `read_csv()`
- [X] randomForest --- was data actually imputed? make this clear early in the code

###### DATA SOURCES
- [ ] Clarify input file origins. Links to other repos (github/doi), papers (doi), or people?
- [ ] * for Whalen - add a repo for the code that made that environmental data

