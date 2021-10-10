# To-do list for ZEN_geography repository

- [ ] check outputs of model comparison script
- [ ] check outputs of figures script
- [ ] check outputs of data assembly script
- [ ] consider removing "2014" from script file names - year does not matter relative to purpose



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


### Data Assembly script
###### PACKAGES
We might want shift away from ddply in plyr to just using tidyverse
- [ ] replace `read.csv()` with `read_csv()`
- [ ] randomForest --- was data actually imputed? make this clear early in the code

###### DATA SOURCES
- [ ] Clarify input file origins. Links to other repos (github/doi), papers (doi), or people?
- [ ] * for Whalen - add a repo for the code that made that environmental data

