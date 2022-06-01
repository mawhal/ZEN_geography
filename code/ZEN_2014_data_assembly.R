###################################################################################
#                                                                                ##
# ZEN 2014 Global eelgrass ecosystem structure: Data assembly                    ##
# RAW data are current as of 2017-04-24                                          ##
# Emmett Duffy (duffye@si.edu)                                                   ##  
# Last updated 2022-05-29                                                        ##
#                                                                                ##
###################################################################################

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# METADATA                                                                        #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# CREATE DERIVED VARIABLES                                                        #
# EXPLORE DISTRIBUTIONS OF VARIABLES (PLOT LEVEL)                                 #
# LOG TRANSFORMS                                                                  #
# OBTAIN SITE MEANS                                                               #
# PCA - ENVIRONMENTAL VARIABLES (GLOBAL)                                          #
# PCA - ENVIRONMENTAL VARIABLES (ATLANTIC)                                        #
# PCA - ENVIRONMENTAL VARIABLES (PACIFIC)                                         #
# EXPLORE DATA COMPLETENESS                                                       #
# IMPUTE MISSING DATA                                                             #
# PCA - EELGRASS VARIABLES (GLOBAL)                                               #
# CREATE SCALED VARIABLES                                                         #
# INTEGRATE DATA FRAMES                                                           #
# CREATE REDUCED DATA FRAMES WITHOUT MISSING SITES OR VALUES                      #
# SUBSET DATA SETS BY GEOGRAPHY                                                   #
# OUTPUT CURATED DATA SETS                                                        #
#                                                                                 #
###################################################################################

###################################################################################
# METADATA                                                                        #
###################################################################################

# This script assembles raw data from the ZEN 2014 global eelgrass ecosystem sampling 
# project, and outputs data files for use in modeling and other applications. See also:
#   ZEN_2014_model_comparison_20210227.R: for data exploration and model building
#   ZEN_2014_figures_20210227.R series: for building figures for the MS

# Source data: For most of the history of this script I was using 
# ZEN_2014_Site&PlotData_2016_05_17_Released.xlsx.


###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages:
library(tidyverse) # for reformatting epibiota data
library(plyr) # to use ddply below in fixing richness values
library(randomForest) # needed for data imputation
library(car) # needed or vif analysis
library(psych) # to visualize relationshiops in pairs panels


###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# To skip ahead TO AFTER imputation, but before other steps, use the following file:
# zen2014_imputed <- read.csv("zen2014_imputed_20210227.csv", header = TRUE)
# zen2014_49_imputed <- read.csv("zen2014_imputed_49_20210227.csv", header = TRUE)

# Summary of data files required:
# 2017-03-16 ZEN Clean MASTER V3 JSL copy.csv
# ...

# MAIN ZEN 2014 DATA SET
# Read in summary data set for ZEN 2014:
# IS THIS THE MOST RECENT, COMPLETE DATA FILE? PROBABLY SHOULD RE-EXTACT FROM DEFINITIVE MASTER EXCEL
zen2014 <- read.csv("data/input/2017-03-16 ZEN Clean MASTER V3 JSL copy.csv", header = TRUE)
names(zen2014)

# BIO-ORACLE CLIMATE AND ENVIRONMENTAL DATA
# Read in Bio-ORACLE and WorldClim environmental data for ZEN sites from Matt Whalen's script:
zen2014.env <- read.csv("data/input/ZEN_2014_data_environmental_20210124_copy.csv", header = TRUE)
names(zen2014.env)

# PERCENT COVER
# Read in data on percent cover:
zen2014.cover <- read.csv("data/input/ZEN_2014_R_PercentCover_by_plot_2016_01_06_modified_copy.csv", header = TRUE)
names(zen2014.cover)
# NOTE: WA.A was not able to collect percent cover data - NO DATA for this site. 

# PREDATION INTENSITY
# Read in data on predation intensity from amphipod tethering assay (Reynolds et al. 2017):
zen2014.pred <- read.csv("data/input/ZEN_2014_R_Predation_by_site_2016_01_21_copy.csv", header = TRUE)
names(zen2014.pred)

# EPIFAUNA AND ASSOCIATED MACROPHYTE SAMPLES
# Load epifauna data
epifauna = read.csv("data/input/ZEN_2014_abund_traits_merge_2016-08-19_copy.csv")

################################################################################
# ----------  WHALEN UNABLE TO FIND THIS FILE
# ----------  BUT, OBJECT `macro` NOT USED IN THE REST OF THE SCRIPT
# # Load macrophyte data
# macro = read.csv("ZEN_2014_DATA_MASTER_2016_08_24_EpifaunaMacrophytes_copy.csv")
################################################################################

# EELGRASS GENETICS
zen2014_gen_fca <- read.csv("data/input/ZEN 2014 FCA scores 20201016.csv", header = TRUE)
zen2014_gen_fca_atlantic <- read.csv("data/input/ZEN_2014_fca_scores_atlantic_20210125_copy.csv", header = TRUE)
zen2014_gen_fca_pacific <- read.csv("data/input/ZEN_2014_fca_scores_pacific_20210125_copy.csv", header = TRUE)

# EPIBIOTA 
# Epibiota (periphyton) data were not calculated correctly in the summary 'Master" data file, 
# specifically filtered material ("Epibiota filter") was not divided by the mass of Zostera scraped. 
# So we have to regenerate these numbers from scratch. 

################################################################################
# ----------  WHALEN UNABLE TO FIND THIS FILE
# # Read in the epibiota data
# zen2014.epibiota <- read.csv("data/input/ZEN_2014_epibiota_mass_2016_01_06.csv", header = TRUE)
# ----------  REPLACING WITH FILE THAT WAS INCLUDED IN REPO
zen2014.epibiota <- read.csv("data/input/ZEN_2014_epibiota_mass_2016_01_06 copy.csv", header = TRUE)
################################################################################

# First, Japan B separated "Chl large epiphytes filter" from "Epibiota filter", the latter of which are all zero. 
# We need to  recode those labeled "Chl large epiphytes filter" as "Epibiota filter", so we can add them together. 
zen2014.epibiota$Species <- as.factor(zen2014.epibiota$Species)
levels(zen2014.epibiota$Species)[levels(zen2014.epibiota$Species) == "Chl Large Epiphytes Filter"] <- "Epibiota filter"
levels(zen2014.epibiota$Species) #This appears to have worked ...

# Second, Oregon A misnamed the epibiota filters ...
levels(zen2014.epibiota$Species)[levels(zen2014.epibiota$Species) == "Epibiota Packet"] <- "Epibiota filter"

# Subset to only the variables of interest here
zen2014.epibiota <- subset(zen2014.epibiota, select = c(Site, Site.Code, Subsite, 
  Sampling.Time, Plot.ID, Unique.ID, Species, Taxa, Group, Type, Dry.Mass.g.))

# # Read in the epibiota data
# zen2014.epibiota <- read.csv("ZEN_2014_epibiota_mass_2016_01_06_JSBfixed.csv", header = TRUE)
# names(zen2014.epibiota)

# Recode second sampling time as 1 for LI
zen2014.epibiota$Site.Code <- as.factor(zen2014.epibiota$Site.Code)
zen2014.epibiota[zen2014.epibiota$Site.Code == "LI" & zen2014.epibiota$Sampling.Time == "1", "Sampling.Time"] = 3
zen2014.epibiota[zen2014.epibiota$Site.Code == "LI" & zen2014.epibiota$Sampling.Time == "2", "Sampling.Time"] = 1

# Regenerate Unique IDs
zen2014.epibiota$Unique.ID = as.character(zen2014.epibiota$Unique.ID)
zen2014.epibiota[zen2014.epibiota$Site.Code == "LI", "Unique.ID"] = 
  paste(zen2014.epibiota[zen2014.epibiota$Site.Code == "LI", "Site.Code"], 
        zen2014.epibiota[zen2014.epibiota$Site.Code == "LI", "Subsite"], 
        zen2014.epibiota[zen2014.epibiota$Site.Code == "LI", "Sampling.Time"], 
        zen2014.epibiota[zen2014.epibiota$Site.Code == "LI", "Plot.ID"],
        sep = ".")

zen2014.epibiota$Unique.ID = as.factor(zen2014.epibiota$Unique.ID)

# Include only first sampling time
zen2014.epibiota = subset(zen2014.epibiota, Sampling.Time == "1")

names(zen2014.epibiota)

# create new wide-form data frame containing only the dry mass data
epibiota.temp <- zen2014.epibiota %>%
  # only take along columns that are unique, otherwise output is staggered in chunks
  select(Unique.ID, Species, Dry.Mass.g.) %>%
  # group the data by species and sample
  group_by(Unique.ID, Species) %>%
  # sum the dry mass for each species in each sample (i.e., sum measurements from the same unique.ID)
  summarize(Dry.mass.g = sum(Dry.Mass.g.)) %>%
  # cast species as columns
  spread(Species, Dry.mass.g, fill = 0)
# No idea how or why 'Dry.Mass.g.' changed to 'Dry.mass.g' but the code somehow did this and won't work without it ...

names(epibiota.temp)
# rename some variables 
names(epibiota.temp)[names(epibiota.temp)=="Epibiota filter"] <- "epibiota.filter"
names(epibiota.temp)[names(epibiota.temp)=="Zostera marina"] <- "epibiota.zostera.marina"

# Subset to only the variables of interest
epibiota.temp <- subset(epibiota.temp, select = c(Unique.ID, epibiota.filter, epibiota.zostera.marina))

# Now the prize: Calculate periphyton mass per g Zostera marina:
epibiota.temp$periphyton.mass.per.g.zostera <- epibiota.temp$epibiota.filter / epibiota.temp$epibiota.zostera.marina

# Add raw and normalized periphyton data back into main data frame
zen2014$epibiota.filter <- epibiota.temp$epibiota.filter[match(zen2014$Unique.ID, epibiota.temp$Unique.ID)]
zen2014$epibiota.zostera.marina <- epibiota.temp$epibiota.zostera.marina[match(zen2014$Unique.ID, epibiota.temp$Unique.ID)]
zen2014$periphyton.mass.per.g.zostera <- epibiota.temp$periphyton.mass.per.g.zostera[match(zen2014$Unique.ID, epibiota.temp$Unique.ID)]

# Remove miscalculated periphyton variable from summary data set
zen2014 <- subset(zen2014, select = -c(Epibiota.Periphyton))

# CLEAN UP AND CONSOLIDATE

# Convert categorical variables to factors
zen2014$Site.Code <- as.factor(zen2014$Site.Code)
zen2014$Ocean <- as.factor(zen2014$Ocean)

# Rename Long Island sites (c'mon guys, read the handbook ...)
# re-leveling apparently only works for factor variables, so convert
zen2014$Site <- as.factor(zen2014$Site)
levels(zen2014$Site)[levels(zen2014$Site)=="LI.1"] <- "LI.A"
levels(zen2014$Site)[levels(zen2014$Site)=="LI.2"] <- "LI.B"


# Rename misspelled or confusing variables
names(zen2014)[names(zen2014)=="Mean.Sheath.Width.cm."] <- "Zostera.sheath.width"
names(zen2014)[names(zen2014)=="Mean.Shealth.Length.cm."] <- "Zostera.sheath.length"
names(zen2014)[names(zen2014)=="Mean.Longest.Leaft.Length.cm."] <- "Zostera.longest.leaf.length"
names(zen2014)[names(zen2014)=="Mean.Above.Zmarina.g"] <- "Zostera.aboveground.mean.mass"
names(zen2014)[names(zen2014)=="Mean.Below.Zmarina.g"] <- "Zostera.belowground.mean.mass"
names(zen2014)[names(zen2014)=="Shoots.Zmarina.per.m2"] <- "Zostera.shoots.per.m2.core"
names(zen2014)[names(zen2014)=="Mean.Fetch"] <- "mean.fetch"
names(zen2014)[names(zen2014)=="PopDens2"] <- "pop.density.2015"
names(zen2014)[names(zen2014)=="mesograzer.total.site.richness"] <- "grazer.richness.site"

names(zen2014.cover)[names(zen2014.cover)=="PercBare"] <- "pct.cover.bare"
names(zen2014.cover)[names(zen2014.cover)=="PercMacroalgae"] <- "pct.cover.macroalgae"
names(zen2014.cover)[names(zen2014.cover)=="PercSeagrass"] <- "pct.cover.seagrass"


# MESOGRAZER SITE RICHNESS: FIX MISSING VALUES 
# Create vector of plots with missing values to see what is missing:
missing.richness <- zen2014[is.na(zen2014$grazer.richness.site), c(3,7)] # columns 3 and 7 are Site, Unique.ID
# We will replace all site richness values with "mean" for that site. First, create vector of means:
temp <- ddply(zen2014, c("Site"), summarize, grazer.richness.site = mean(grazer.richness.site, na.rm = T))
# But CR.A has NO mesograzers at all so returns NaN.  Assume species pool is same as for CR.B (S = 3) and replace:
temp[temp == "NaN"] = NA 
temp$grazer.richness.site[is.na(temp$grazer.richness.site)] <- 3 # CR.A grazer richness now  = 3
zen2014$grazer.richness.site <- temp$grazer.richness.site[match(zen2014$Site, temp$Site)]
sum(is.na(zen2014$grazer.richness.site)) # 0. Voila.  


# Add cover and predation data to main ZEN dataframe:
zen2014$pct.cover.bare <- zen2014.cover$pct.cover.bare[match(zen2014$Unique.ID, zen2014.cover$Unique.ID)]
zen2014$pct.cover.macroalgae <- zen2014.cover$pct.cover.macroalgae[match(zen2014$Unique.ID, zen2014.cover$Unique.ID)]
zen2014$pct.cover.seagrass <- zen2014.cover$pct.cover.seagrass[match(zen2014$Unique.ID, zen2014.cover$Unique.ID)]

zen2014$predation.amphipods <- zen2014.pred$Mean.Pred.Amphipod[match(zen2014$Site, zen2014.pred$Site)]
zen2014$predation.gastropods <- zen2014.pred$Mean.Pred.Gastropod[match(zen2014$Site, zen2014.pred$Site)]
zen2014$predation.squidpops <- zen2014.pred$Mean.Pred.Squid[match(zen2014$Site, zen2014.pred$Site)]


# Add BioOracle environmental data to main ZEN dataframe:
zen2014$sst.min <- zen2014.env$sstmin[match(zen2014$Site, zen2014.env$Site)]
zen2014$sst.mean <- zen2014.env$sstmean[match(zen2014$Site, zen2014.env$Site)]
zen2014$sst.max <- zen2014.env$sstmax[match(zen2014$Site, zen2014.env$Site)]
zen2014$sst.range <- zen2014.env$sstrange[match(zen2014$Site, zen2014.env$Site)]
zen2014$chlomean <- zen2014.env$chlomean[match(zen2014$Site, zen2014.env$Site)]
zen2014$nitrate <- zen2014.env$nitrate[match(zen2014$Site, zen2014.env$Site)]
zen2014$parmean <- zen2014.env$parmean[match(zen2014$Site, zen2014.env$Site)]
zen2014$cloudmean <- zen2014.env$cloudmean[match(zen2014$Site, zen2014.env$Site)]
zen2014$day.length <- zen2014.env$Day.length.hours[match(zen2014$Site, zen2014.env$Site)]
zen2014$ph <- zen2014.env$ph[match(zen2014$Site, zen2014.env$Site)]
zen2014$phosphate <- zen2014.env$phosphate[match(zen2014$Site, zen2014.env$Site)]
zen2014$salinity <- zen2014.env$salinity[match(zen2014$Site, zen2014.env$Site)]
zen2014$precipitation <- zen2014.env$precip[match(zen2014$Site, zen2014.env$Site)]

# Reorder variables 'Coast': WP to EA 
zen2014$Coast <- as.factor(zen2014$Coast)
zen2014$Coast <- factor(zen2014$Coast, levels = c("West Pacific", "East Pacific", "West Atlantic", "East Atlantic"))

names(zen2014)


###################################################################################
# CREATE DERIVED VARIABLES                                                        #
###################################################################################

# # Create variables for water column stiochiometry
# zen2014$NP.ratio <-  zen2014$nitrate / zen2014$phosphate
# hist(zen2014$NP.ratio) # not too far from normal
# zen2014$PARP.ratio <-  zen2014$parmean / zen2014$phosphate

# Percentage of crustaceans and gastropods among the mesograzers
zen2014$crust.pct.mass <-  zen2014$Malacostraca.mesograzer.plot.biomass.std.mg.g / zen2014$mesograzer.total.plot.biomass.std.mg.g
zen2014$gast.pct.mass <-  zen2014$Gastropoda.mesograzer.plot.biomass.std.mg.g / zen2014$mesograzer.total.plot.biomass.std.mg.g

# grazer and periphyton nunmbers per unit bottom area (i.e., core)
zen2014$mesograzer.abund.per.area <-  zen2014$mesograzer.total.plot.abund.std.g * zen2014$Zostera.aboveground.mean.mass
zen2014$crustacean.mass.per.area <-  zen2014$Malacostraca.mesograzer.plot.biomass.std.mg.g * zen2014$Zostera.aboveground.mean.mass
zen2014$gastropod.mass.per.area <-  zen2014$Gastropoda.mesograzer.plot.biomass.std.mg.g * zen2014$Zostera.aboveground.mean.mass
zen2014$mesograzer.mass.per.area <-  zen2014$mesograzer.total.plot.biomass.std.mg.g * zen2014$Zostera.aboveground.mean.mass
zen2014$periphyton.mass.per.area <-  zen2014$periphyton.mass.per.g.zostera * zen2014$Zostera.aboveground.mean.mass

# Leaf C:N ratio
zen2014$leaf.CN.ratio <-  zen2014$Leaf.PercC / zen2014$Leaf.PercN


###################################################################################
# EXPLORE DISTRIBUTIONS OF VARIABLES (PLOT LEVEL)                                 #
###################################################################################

# Examine frequency distribution of sites by environmental factor
# par(mfrow = c(1,1))
# par(mfrow = c(2,4))
hist(zen2014$Latitude, col = "cyan", main = "Surveys by latitude")    
hist(zen2014$Longitude, col = "cyan", main = "Surveys by longitude")    
hist(zen2014$Temperature.C, col = "cyan", main = "Surveys by temperature")    
hist(zen2014$Salinity.ppt, col = "cyan", main = "Surveys by salinity")    
hist(zen2014$pop.density.2015, col = "cyan", main = "Surveys by population density")    
hist(zen2014$day.length, col = "cyan", main = "Surveys by day length")    
hist(zen2014$mean.fetch, col = "cyan", main = "Surveys by mean fetch")    

hist(zen2014$Zostera.aboveground.mean.mass, col = "cyan", main = "Surveys by Zostera AG biomass")    
hist(zen2014$periphyton.mass.per.g.zostera, col = "cyan", main = "Surveys by periphyton biomass")    
hist(zen2014$Malacostraca.mesograzer.plot.abund.std.g, col = "cyan", main = "Surveys by crustacean biomass")    
hist(zen2014$Gastropoda.mesograzer.plot.biomass.std.mg.g, col = "cyan", main = "Surveys by gastropod biomass")    
hist(zen2014$grazer.richness.site, col = "cyan", main = "Surveys by mesograzer richness")    

hist(zen2014$mesograzer.total.plot.biomass.std.mg.g, col = "cyan", main = "Surveys by mesograzer biomass")    
hist(zen2014$epifauna.total.plot.biomass.std.mg.g, col = "cyan", main = "Surveys by mobile epifauna biomass")    

hist(zen2014$pct.cover.seagrass, col = "cyan", main = "Surveys by seagrass % cover")    
hist(zen2014$pct.cover.macroalgae, col = "cyan", main = "Surveys by macroalgal % cover")    

hist(zen2014$predation.amphipods, col = "cyan", main = "Surveys by amphipod loss to predation")    
hist(zen2014$predation.gastropods, col = "cyan", main = "Surveys by gastropod loss to predation")    
hist(zen2014$predation.squidpops, col = "cyan", main = "Surveys by squidpop loss to predation")    


###################################################################################
# LOG TRANSFORMS                                                                  #
###################################################################################

# NOTE: For many variables I add a constant roughly equal to the smallest value recorded 

zen2014$log10.Zostera.AG.mass <- log10(zen2014$Zostera.aboveground.mean.mass + 1) 
zen2014$log10.Zostera.BG.mass <- log10(zen2014$Zostera.belowground.mean.mass + 1) 
zen2014$log10.Zostera.shoots.core <- log10(zen2014$Zostera.shoots.per.m2.core) 
zen2014$log10.Zostera.sheath.width <- log10(zen2014$Zostera.sheath.width) 
zen2014$log10.Zostera.sheath.length <- log10(zen2014$Zostera.sheath.length) 
zen2014$log10.Zostera.longest.leaf.length <- log10(zen2014$Zostera.longest.leaf.length) 

zen2014$log10.epibiota.filter <- log10(zen2014$epibiota.filter) 
zen2014$log10.epibiota.zostera.marina <- log10(zen2014$epibiota.zostera.marina) 
zen2014$log10.periphyton.mass.per.g.zostera <- log10(zen2014$periphyton.mass.per.g.zostera + 0.001) 
zen2014$log10.periphyton.mass.per.area <- log10(zen2014$periphyton.mass.per.area + 0.1) 

zen2014$log10.mesograzer.abund.per.g.plant <- log10(zen2014$mesograzer.total.plot.abund.std.g + 0.01) 
zen2014$log10.crustacean.abund.per.g.plant <- log10(zen2014$Malacostraca.mesograzer.plot.abund.std.g + 0.01) 
zen2014$log10.gastropod.abund.per.g.plant <- log10(zen2014$Gastropoda.mesograzer.plot.abund.std.g + 0.01) 

zen2014$log10.mesograzer.mass.per.g.plant <- log10(zen2014$mesograzer.total.plot.biomass.std.mg.g + 0.01) 
zen2014$log10.crustacean.mass.per.g.plant <- log10(zen2014$Malacostraca.mesograzer.plot.biomass.std.mg.g + 0.01) 
zen2014$log10.gastropod.mass.per.g.plant <- log10(zen2014$Gastropoda.mesograzer.plot.biomass.std.mg.g + 0.01) 

zen2014$log10.mesograzer.abund.per.area <-  log10(zen2014$mesograzer.abund.per.area + 1) 

zen2014$log10.crustacean.mass.per.area <-  log10(zen2014$crustacean.mass.per.area + 1) 
zen2014$log10.gastropod.mass.per.area <-  log10(zen2014$gastropod.mass.per.area + 1) 
zen2014$log10.mesograzer.mass.per.area <-  log10(zen2014$mesograzer.mass.per.area + 1) 

zen2014$log10.grazer.richness.site <- log10(zen2014$grazer.richness.site + 1) 

zen2014$log10.day.length <- log10(zen2014$day.length) 
zen2014$log10.Leaf.PercN <- log10(zen2014$Leaf.PercN) 
zen2014$sqrt.nitrate <- sqrt(zen2014$nitrate) 
zen2014$log10.phosphate <- log10(zen2014$phosphate) 
zen2014$log10.chlomean <- log10(zen2014$chlomean) 
zen2014$log10.mean.fetch <- log10(zen2014$mean.fetch) 

zen2014$log10.pct.cover.bare <- log10(zen2014$pct.cover.bare)
zen2014$log10.pct.cover.macroalgae <- log10(zen2014$pct.cover.macroalgae + 0.1)
zen2014$log10.pct.cover.seagrass <- log10(zen2014$pct.cover.seagrass)
zen2014$sqrt.pct.cover.seagrass <- sqrt(zen2014$pct.cover.seagrass)

zen2014$log10.predation.amphipods <- log10(zen2014$predation.amphipods)
zen2014$log10.predation.gastropods <- log10(zen2014$predation.gastropods)
zen2014$log10.predation.squidpops <- log10(zen2014$predation.squidpops)


hist(zen2014$nitrate)
hist(zen2014$sqrt.nitrate)

hist(zen2014$log10.Zostera.AG.mass)

hist(zen2014$log10.pct.cover.bare)
hist(zen2014$log10.pct.cover.macroalgae)
hist(zen2014$log10.pct.cover.seagrass) # Crap. Raw data are better

hist(zen2014$log10.predation.amphipods)    
hist(zen2014$log10.predation.gastropods) # Good - better than raw  
hist(zen2014$log10.predation.squidpops) # worse than raw 

# Change values of NaN to NA:
zen2014[zen2014 == "NaN"] = NA 
zen2014[zen2014 == "Inf"] = NA 
names(zen2014)


###################################################################################
# OBTAIN SITE MEANS                                                               #
###################################################################################

# CAN THIS GO AFTER IMPUTATION SECTION? SHOULD IT? 

# Obtain mean values per site
ZEN_2014_site_means <- ddply(zen2014, c("Site"), summarize, 
                             Seagrass.pct.cover.site = mean(pct.cover.seagrass, na.rm = T),
                             Macroalgae.pct.cover.site = mean(pct.cover.macroalgae, na.rm = T),
                             Zostera.AG.mass.site = mean(Zostera.aboveground.mean.mass, na.rm = T), 
                             Zostera.BG.mass.site = mean(Zostera.belowground.mean.mass, na.rm = T),                      
                             Zostera.shoots.core.site = mean(Zostera.shoots.per.m2.core, na.rm = T),                      
                             Zostera.sheath.width.site = mean(Zostera.sheath.width, na.rm = T),                      
                             Zostera.sheath.length.site = mean(Zostera.sheath.length, na.rm = T),                      
                             Zostera.longest.leaf.length.site = mean(Zostera.longest.leaf.length, na.rm = T), 
                             epibiota.filter.site = mean(epibiota.filter, na.rm = T), 
                             epibiota.zostera.marina.site = mean(epibiota.zostera.marina, na.rm = T), 
                             periphyton.mass.per.g.zostera.site = mean(periphyton.mass.per.g.zostera, na.rm = T),                      
                             
                             mesograzer.abund.per.g.plant.site = mean(mesograzer.total.plot.abund.std.g, na.rm = T),                      
                             crustacean.abund.per.g.plant.site = mean(Malacostraca.mesograzer.plot.abund.std.g, na.rm = T), 
                             gastropod.abund.per.g.plant.site = mean(Gastropoda.mesograzer.plot.abund.std.g, na.rm = T), 
                             mesograzer.mass.per.g.plant.site = mean(mesograzer.total.plot.biomass.std.mg.g, na.rm = T),                      
                             crustacean.mass.per.g.plant.site = mean(Malacostraca.mesograzer.plot.biomass.std.mg.g, na.rm = T), 
                             gastropod.mass.per.g.plant.site = mean(Gastropoda.mesograzer.plot.biomass.std.mg.g, na.rm = T), 
                             mesograzer.mass.per.area.site = mean(mesograzer.mass.per.area, na.rm = T),                      
                             crustacean.mass.per.area.site = mean(crustacean.mass.per.area, na.rm = T),                      
                             gastropod.mass.per.area.site = mean(gastropod.mass.per.area, na.rm = T),  
                             periphyton.mass.per.area.site = mean(periphyton.mass.per.area, na.rm = T),  
                             log10.grazer.richness.site = mean(log10.grazer.richness.site, na.rm = T), 
                             
                             crust.pct.mass.site = mean(crust.pct.mass, na.rm = T), 
                             gast.pct.mass.site = mean(gast.pct.mass, na.rm = T), 
                             
                             Leaf.PercN.site = mean(Leaf.PercN, na.rm = T), 
                             leaf.CN.ratio.site = mean(leaf.CN.ratio, na.rm = T), 
                             
                             log10.Macroalgae.pct.cover.site = mean(log10.pct.cover.macroalgae, na.rm = T),
                             log10.Zostera.AG.mass.site = mean(log10.Zostera.AG.mass, na.rm = T), 
                             log10.Zostera.BG.mass.site = mean(log10.Zostera.BG.mass, na.rm = T), 
                             log10.Zostera.shoots.core.site = mean(log10.Zostera.shoots.core, na.rm = T), 
                             log10.Zostera.sheath.width.site = mean(log10.Zostera.sheath.width, na.rm = T), 
                             log10.Zostera.sheath.length.site = mean(log10.Zostera.sheath.length, na.rm = T), 
                             log10.Zostera.longest.leaf.length.cm.site = mean(log10.Zostera.longest.leaf.length, na.rm = T), 
                             log10.periphyton.mass.per.g.zostera.site = mean(log10.periphyton.mass.per.g.zostera, na.rm = T), 
                             
                             log10.mesograzer.abund.per.g.plant.site = mean(log10.mesograzer.abund.per.g.plant, na.rm = T), 
                             log10.crustacean.abund.per.g.plant.site = mean(log10.crustacean.abund.per.g.plant, na.rm = T), 
                             log10.gastropod.abund.per.g.plant.site = mean(log10.gastropod.abund.per.g.plant, na.rm = T), 
                             
                             log10.mesograzer.mass.per.g.plant.site = mean(log10.mesograzer.mass.per.g.plant, na.rm = T), 
                             log10.crustacean.mass.per.g.plant.site = mean(log10.crustacean.mass.per.g.plant, na.rm = T), 
                             log10.gastropod.mass.per.g.plant.site = mean(log10.gastropod.mass.per.g.plant, na.rm = T), 
                             
                             log10.mesograzer.abund.per.area.site = mean(log10.mesograzer.abund.per.area, na.rm = T), # check this - did it work?                   
                             
                             log10.mesograzer.mass.per.area.site = mean(log10.mesograzer.mass.per.area, na.rm = T),                      
                             log10.crustacean.mass.per.area.site = mean(log10.crustacean.mass.per.area, na.rm = T),                      
                             log10.gastropod.mass.per.area.site = mean(log10.gastropod.mass.per.area, na.rm = T),                      
                             log10.periphyton.mass.per.area.site = mean(log10.periphyton.mass.per.area, na.rm = T),                      
                             
                             log10.Leaf.PercN.site = mean(log10.Leaf.PercN, na.rm = T), 
                             
                             log10.predation.gastropods.site = mean(log10.predation.gastropods, na.rm = T)
)


# Add in predation data (already site level)
ZEN_2014_site_means$predation.amphipods <- zen2014$predation.amphipods[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$predation.squidpops <- zen2014$predation.squidpops[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$predation.gastropods <- zen2014$predation.gastropods[match(ZEN_2014_site_means$Site, zen2014$Site)]

ZEN_2014_site_means$grazer.richness.site <- zen2014$grazer.richness.site[match(ZEN_2014_site_means$Site, zen2014$Site)]

# Change values of NaN to NA:
ZEN_2014_site_means[ZEN_2014_site_means == "NaN"] = NA 

# Add site-level environmental (and other) variables back in
ZEN_2014_site_means$Ocean <- zen2014$Ocean[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$Coast <- zen2014$Coast[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$Latitude <- zen2014$Latitude[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$Longitude <- zen2014$Longitude[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$Temperature.C <- zen2014$Temperature.C[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$Salinity.ppt <- zen2014$Salinity.ppt[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$log10.mean.fetch <- zen2014$log10.mean.fetch[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$day.length <- zen2014$day.length[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$log10.day.length <- zen2014$log10.day.length[match(ZEN_2014_site_means$Site, zen2014$Site)]

ZEN_2014_site_means$sst.min <- zen2014$sst.min[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$sst.mean <- zen2014$sst.mean[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$sst.max <- zen2014$sst.max[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$sst.range <- zen2014$sst.range[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$salinity <- zen2014$salinity[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$parmean <- zen2014$parmean[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$cloudmean <- zen2014$cloudmean[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$precipitation <- zen2014$precipitation[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$nitrate <- zen2014$nitrate[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$sqrt.nitrate <- zen2014$sqrt.nitrate[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$ph <- zen2014$ph[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$phosphate <- zen2014$phosphate[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$log10.phosphate <- zen2014$log10.phosphate[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$NP.ratio <- zen2014$NP.ratio[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$chlomean <- zen2014$chlomean[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$log10.chlomean <- zen2014$log10.chlomean[match(ZEN_2014_site_means$Site, zen2014$Site)]
ZEN_2014_site_means$pop.density.2015 <- zen2014$pop.density.2015[match(ZEN_2014_site_means$Site, zen2014$Site)]

names(ZEN_2014_site_means)

# Add genetic data to site means data frame
ZEN_2014_site_means$FC1 <- zen2014_gen_fca$FC1[match(ZEN_2014_site_means$Site, zen2014_gen_fca$Site)]
ZEN_2014_site_means$FC2 <- zen2014_gen_fca$FC2[match(ZEN_2014_site_means$Site, zen2014_gen_fca$Site)]

# For boxplots, reorder variable 'Coast': WP to EA
ZEN_2014_site_means$Coast <- factor(ZEN_2014_site_means$Coast, levels = c("West Pacific", "East Pacific", "West Atlantic", "East Atlantic"))

# Create separate data sets by Ocean - SITE level
ZEN_2014_site_means_Atlantic <- droplevels(subset(ZEN_2014_site_means, Ocean == "Atlantic"))
ZEN_2014_site_means_Pacific <- droplevels(subset(ZEN_2014_site_means, Ocean == "Pacific"))
ZEN_2014_site_means_49_Atlantic <- droplevels(subset(ZEN_2014_site_means_Atlantic, Site != "SW.A"))


# Incorporate genetic FCA axes done separately by ocean
ZEN_2014_site_means_Atlantic$FC1.atl <- zen2014_gen_fca_atlantic$fca1.atl[match(ZEN_2014_site_means_Atlantic$Site, zen2014_gen_fca_atlantic$Site)]
ZEN_2014_site_means_Atlantic$FC2.atl <- zen2014_gen_fca_atlantic$fca2.atl[match(ZEN_2014_site_means_Atlantic$Site, zen2014_gen_fca_atlantic$Site)]
ZEN_2014_site_means_Pacific$FC1.pac <- zen2014_gen_fca_pacific$fca1.pac[match(ZEN_2014_site_means_Pacific$Site, zen2014_gen_fca_pacific$Site)]
ZEN_2014_site_means_Pacific$FC2.pac <- zen2014_gen_fca_pacific$fca2.pac[match(ZEN_2014_site_means_Pacific$Site, zen2014_gen_fca_pacific$Site)]


###################################################################################
# PCA - ENVIRONMENTAL VARIABLES (GLOBAL)                                          #
###################################################################################

# Explore correlations among environmental drivers
pairs.panels(ZEN_2014_site_means[,c("Latitude", "sst.mean", "sst.range", "sst.min", "sst.max", "Salinity.ppt", 
  "parmean", "log10.day.length", "cloudmean", "precipitation", "sqrt.nitrate", "log10.phosphate", "log10.chlomean", 
  "Leaf.PercN.site", "log10.mean.fetch")], 
  smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 8)

# Create data frame containing the ZEN 2014 environmental variables for PCA
# Note: Some exploration shows that nitrate is closely correlated with several other
# variables, and taking it out results in first 3 PC axes explaining ~75% of variation. This 
# is parsimonious and simplifies the analysis. 
ZEN.env <- ZEN_2014_site_means[c("sst.mean", "sst.range", "Salinity.ppt", "parmean", 
   "cloudmean", "log10.phosphate", "log10.chlomean", "Leaf.PercN.site"
  # , "precipitation", "log10.day.length",
  )]
ZEN.sites <- ZEN_2014_site_means[c("Site")]


# Compute PCAs
ZEN.env.pca <- prcomp(ZEN.env, center = TRUE, scale. = TRUE) 

print(ZEN.env.pca)
#                        PC1         PC2         PC3        PC4         PC5         PC6         PC7         PC8
# sst.mean         0.5344090 -0.04221968  0.12650153 -0.2221002  0.11595693 -0.56707288  0.49861640  0.25230424
# sst.range       -0.1607624 -0.40262794  0.45918615 -0.4862507  0.41315371  0.25966358 -0.15719348  0.31925476
# Salinity.ppt     0.3702257  0.16135868 -0.48106388 -0.4651378  0.05646463 -0.08442206 -0.61172656  0.06779392
# parmean          0.4076216  0.22572201  0.39507514  0.3928616 -0.25219684  0.21903419 -0.29892746  0.52108800
# cloudmean       -0.4937825 -0.21507910 -0.27382435  0.1300389 -0.18748290 -0.44075941 -0.12798127  0.61010822
# log10.phosphate -0.2101797  0.54450089 -0.13760560 -0.4243534 -0.22277173  0.36170941  0.41340358  0.33010411
# log10.chlomean  -0.2566312  0.34762747  0.53996106 -0.2846051 -0.31346195 -0.45082306 -0.26740350 -0.26025590
# Leaf.PercN.site -0.1774368  0.54363232  0.01286878  0.2560322  0.75235033 -0.16600039 -0.06571552  0.09672818

# Interpretation:
# PCe1: latitude/climate: high = warmer, brighter, less cloudy (lower latitude)
# PCe2: nutrient status: high = high PO4, leaf N  
# PCe3: estuarine: low salinity, variable temp, high chl

# plot cumulative proportion of variance explained by PC axes
plot(ZEN.env.pca, type = "l")

# Calculate proportion of variance explained by each PC
summary(ZEN.env.pca)
#                           PC1    PC2    PC3    PC4     PC5     PC6     PC7    PC8
# Standard deviation     1.6849 1.4240 1.1552 0.9516 0.65646 0.48125 0.36494 0.3124
# Proportion of Variance 0.3549 0.2535 0.1668 0.1132 0.05387 0.02895 0.01665 0.0122
# Cumulative Proportion  0.3549 0.6083 0.7751 0.8883 0.94220 0.97115 0.98780 1.0000


# Combine PCA scores with SITE-level data frame
site.env.pca.scores <- ZEN.env.pca$x
site.env.pca.scores <- cbind(ZEN.sites, site.env.pca.scores) 
ZEN_2014_site_means <- cbind(ZEN_2014_site_means, site.env.pca.scores) 

# Rename PCA variables 1-3 and cull PC4-7
names(ZEN_2014_site_means)[names(ZEN_2014_site_means)=="PC1"] <- "PC1.env.global"
names(ZEN_2014_site_means)[names(ZEN_2014_site_means)=="PC2"] <- "PC2.env.global"
names(ZEN_2014_site_means)[names(ZEN_2014_site_means)=="PC3"] <- "PC3.env.global"
ZEN_2014_site_means <- subset(ZEN_2014_site_means, select = -c(PC4,PC5,PC6, PC7, PC8))

names(ZEN_2014_site_means)


###################################################################################
# PCA - ENVIRONMENTAL VARIABLES (ATLANTIC)                                        #
###################################################################################

# Explore correlations among environmental drivers
pairs.panels(ZEN_2014_site_means_Atlantic[,c("sst.mean", "sst.range", "Salinity.ppt", "parmean", 
  "cloudmean", "log10.phosphate", "log10.chlomean", "Leaf.PercN.site"
  # , "precipitation", "log10.day.length"
  )], 
  smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 8)

# Create data frame containing the ZEN 2014 environmental variables for PCA
# Note: Some exploration shows that nitrate is closely corrtelated with several other
# variables, and taking it out results in first 3 PC axes explaining ~75% of variation. This 
# is parsimonious and simplifies the analysis. 
ZEN.env.atl <- ZEN_2014_site_means_Atlantic[c("sst.mean", "sst.range", "Salinity.ppt", "parmean", 
  "cloudmean", "log10.phosphate", "log10.chlomean", "Leaf.PercN.site"
  # , "precipitation", "log10.day.length"
  )]
ZEN.sites.atl <- ZEN_2014_site_means_Atlantic[c("Site")]

# Compute PCAs
ZEN.env.pca.atl <- prcomp(ZEN.env.atl, center = TRUE, scale. = TRUE) 

print(ZEN.env.pca.atl)
#                          PC1         PC2         PC3         PC4         PC5         PC6         PC7         PC8
# sst.mean        -0.550319750  0.07256028 -0.14266055  0.01964309 -0.26247919  0.50693440  0.34783455  0.47358063
# sst.range        0.008028728  0.53243059  0.13905815  0.56739502 -0.43108238 -0.27143385 -0.27518892  0.19985403
# Salinity.ppt    -0.312338254 -0.33887929 -0.52503373  0.18367192 -0.01847826  0.09370635 -0.67915383 -0.08853019
# parmean         -0.307553079  0.44084782  0.04824442 -0.51027750  0.40436656 -0.23475705 -0.34111047  0.33671071
# cloudmean        0.486920633 -0.28474069  0.27891671  0.01450237  0.05975327  0.32618638 -0.33098360  0.61992583
# log10.phosphate  0.294237976  0.02478199 -0.66842063  0.18661880  0.25670169 -0.33170669  0.31530953  0.39478092
# log10.chlomean   0.265024764  0.54345377 -0.20872625  0.13627880  0.32327853  0.62268122 -0.08243877 -0.27063675
# Leaf.PercN.site  0.333217372  0.15821912 -0.33789872 -0.57441251 -0.63831315  0.03592546 -0.09954655 -0.03411444

# Interpretation:
# PCe1: latitude/climate: high = cooler, cloudier
# PCe2: estuarine/eutrophic: high = high phytoplankton, variable temperature, bright, lowish salinity
# PCe3: arid watershed? oligotrophic Baltic?: high = low salinity, low PO4

# plot cumulative proportion of variance explained by PC axes
plot(ZEN.env.pca.atl, type = "l")

# Calculate proportion of variance explained by each PC
summary(ZEN.env.pca.atl)
#                           PC1    PC2    PC3    PC4     PC5     PC6     PC7     PC8
# Standard deviation     1.6778 1.4182 1.2097 0.9687 0.62444 0.46015 0.35168 0.21673
# Proportion of Variance 0.3519 0.2514 0.1829 0.1173 0.04874 0.02647 0.01546 0.00587
# Cumulative Proportion  0.3519 0.6033 0.7862 0.9035 0.95220 0.97867 0.99413 1.00000

# Output PCA scores for each site and combine with site means data frame
site.env.pca.scores.atl <- ZEN.env.pca.atl$x
site.env.pca.scores.atl <- cbind(ZEN.sites.atl, site.env.pca.scores.atl) 
ZEN_2014_site_means_Atlantic <- cbind(ZEN_2014_site_means_Atlantic, site.env.pca.scores.atl) 

# Rename PCA variables 1-3 and cull PC4-7
names(ZEN_2014_site_means_Atlantic)[names(ZEN_2014_site_means_Atlantic)=="PC1"] <- "PC1.env.atl"
names(ZEN_2014_site_means_Atlantic)[names(ZEN_2014_site_means_Atlantic)=="PC2"] <- "PC2.env.atl"
names(ZEN_2014_site_means_Atlantic)[names(ZEN_2014_site_means_Atlantic)=="PC3"] <- "PC3.env.atl"
ZEN_2014_site_means_Atlantic <- subset(ZEN_2014_site_means_Atlantic, select = -c(PC4,PC5,PC6, PC7, PC8))
names(ZEN_2014_site_means_Atlantic)


###################################################################################
# PCA - ENVIRONMENTAL VARIABLES (PACIFIC)                                         #
###################################################################################

# Explore correlations among environmental drivers
pairs.panels(ZEN_2014_site_means_Pacific[,c("Latitude", "sst.mean", "sst.range", "sst.min", "sst.max", "Salinity.ppt", 
  "parmean", "log10.day.length", "cloudmean", "precipitation", "sqrt.nitrate", "log10.phosphate", "log10.chlomean", 
  "Leaf.PercN.site", "log10.mean.fetch")], 
  smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 8)

# Crap, nearly everything is strongly correlated with latitude ...

# Create data frame containing the ZEN 2014 environmental variables for PCA
# Note: Some exploration shows that nitrate is closely correlated with several other
# variables, and taking it out results in first 3 PC axes explaining ~75% of variation. This 
# is parsimonious and simplifies the analysis. 
ZEN.env.pac <- ZEN_2014_site_means_Pacific[c("sst.mean", "sst.range", "Salinity.ppt", "parmean", 
  "cloudmean", "log10.phosphate", "log10.chlomean", "Leaf.PercN.site"
  # , "precipitation", "log10.day.length"
  )]
ZEN.sites.pac <- ZEN_2014_site_means_Pacific[c("Site")]

# Compute PCAs
ZEN.env.pca.pac <- prcomp(ZEN.env.pac, center = TRUE, scale. = TRUE) 

print(ZEN.env.pca.pac)
#                        PC1         PC2         PC3         PC4         PC5         PC6         PC7         PC8
# sst.mean         0.4416493 -0.14998580  0.38471592 -0.09795308  0.11434105 -0.46072831  0.20973408  0.59625174
# sst.range       -0.1192591 -0.58280840  0.13287760 -0.61360069 -0.05905439  0.33920457  0.35256555 -0.09539264
# Salinity.ppt     0.4002213  0.04551641 -0.50374668  0.10047825  0.56645013  0.33137386  0.37126335  0.07337350
# parmean          0.4058142  0.32386570  0.34788599  0.04747739 -0.21892434 -0.03638956  0.46415902 -0.58519351
# cloudmean       -0.3739858 -0.36483629 -0.19281131  0.31505264  0.17158648 -0.57553674  0.40395943 -0.25831575
# log10.phosphate -0.4215990  0.32143191  0.04272324  0.20357878 -0.29958265  0.25871215  0.55282804  0.46191518
# log10.chlomean  -0.3422080  0.18681817  0.58063017  0.01610370  0.69496221  0.13777490 -0.03697576 -0.08532484
# Leaf.PercN.site -0.1764750  0.50946333 -0.28882340 -0.67866604  0.11171810 -0.37997422  0.09088341 -0.01326676

# Interpretation:
# PCe1: latitude/climate: high = warmer, brighter, higher salinity, lower PO4
# PCe2: nutrient status: high = high nutrients (especially leaf N), more stable temperature
# PCe3: estuarine/eutrophic: high = low salinity, high chl

# plot cumulative proportion of variance explained by PC axes
plot(ZEN.env.pca.pac, type = "l")

# Calculate proportion of variance explained by each PC
summary(ZEN.env.pca.pac)
#                           PC1    PC2    PC3     PC4     PC5     PC6     PC7     PC8
# Standard deviation     1.9641 1.4390 0.9141 0.71060 0.62592 0.49046 0.24605 0.19570
# Proportion of Variance 0.4822 0.2588 0.1045 0.06312 0.04897 0.03007 0.00757 0.00479
# Cumulative Proportion  0.4822 0.7410 0.8455 0.90860 0.95758 0.98765 0.99521 1.00000

# Output PCA scores for each site and combine with site means data frame
site.env.pca.scores.pac <- ZEN.env.pca.pac$x
site.env.pca.scores.pac <- cbind(ZEN.sites.pac, site.env.pca.scores.pac) 
ZEN_2014_site_means_Pacific <- cbind(ZEN_2014_site_means_Pacific, site.env.pca.scores.pac) 

# Rename PCA variables 1-3 and cull PC4-7
names(ZEN_2014_site_means_Pacific)[names(ZEN_2014_site_means_Pacific)=="PC1"] <- "PC1.env.pac"
names(ZEN_2014_site_means_Pacific)[names(ZEN_2014_site_means_Pacific)=="PC2"] <- "PC2.env.pac"
names(ZEN_2014_site_means_Pacific)[names(ZEN_2014_site_means_Pacific)=="PC3"] <- "PC3.env.pac"
ZEN_2014_site_means_Pacific <- subset(ZEN_2014_site_means_Pacific, select = -c(PC4,PC5,PC6, PC7, PC8))

names(ZEN_2014_site_means_Pacific)


###################################################################################
# EXPLORE DATA COMPLETENESS                                                       #
###################################################################################

# NOTE: AIC comparisons among models are invalid unless exactly the same number of plots
# are used in each comparison, because the DF influences calculation of the AIC score. 
# This means that we need data on all plots and need to impute missing data for 
# valid AIC model comparisons. 

# How many observations are missing for each variable?
sum(is.na(zen2014$pct.cover.seagrass)) # 60
sum(is.na(zen2014$pct.cover.macroalgae)) # 60
sum(is.na(zen2014$log10.Zostera.AG.mass)) # 24
sum(is.na(zen2014$log10.Zostera.shoots.core)) # 15
sum(is.na(zen2014$Zostera.longest.leaf.length)) # 0
sum(is.na(zen2014$Leaf.PercN)) # 14
sum(is.na(zen2014$Temperature.C)) # 0
sum(is.na(zen2014$Salinity.ppt)) # 0
sum(is.na(zen2014$pop.density.2015)) # 20 huh?
sum(is.na(zen2014$GenotypicRichness)) # 0
sum(is.na(zen2014$AllelicRichness)) # 0 
sum(is.na(zen2014$grazer.richness.site)) # 0
sum(is.na(zen2014$log10.periphyton.mass.per.g.zostera)) # 4
sum(is.na(zen2014$log10.mesograzer.abund.per.g.plant)) # 9
sum(is.na(zen2014$log10.crustacean.abund.per.g.plant)) # 9 
sum(is.na(zen2014$log10.gastropod.abund.per.g.plant)) # 9

# Look at percentage of values missing for each variable
# First create function to calculate % of missing values infor each variable in a data frame… 
pMiss <- function(x){sum(is.na(x))/length(x)*100}
# Now apply it to the data frame: 
apply(zen2014,2,pMiss)

# Results: Most variables have fewer < 4% missing. Exceptions are cover bare (42%), 
# predation (16%), the mean body mass variables, which are missing ~25% of values because 
# body mass can't be calculated for samples with no animals. 


###################################################################################
# IMPUTE MISSING DATA                                                             #
###################################################################################

# NOTE: To sidestep the long imputation process, read in the most recent versions of the 
# imputed data sets (see intro section above)

# First the rationale: Model comparisons requires that alternative models use exactly the 
# same source data set. Missing cells result in slightly different data sets for models
# that include, versus do not include, that variable, with different DF such that 
# resulting AIC scores used to compare them will be invalid. Solving this requires 
# either (1) throwing out all rows that have a missing cell, or (2) imputing (modeling) 
# the missing values. The former is far the worse alternative since it will end up 
# discarding a substantial part of the entire data set. 

# A few other points: we cannot include the following predictors because biased by 
# missing from entire site: periphyton (missing from SW.A). Also BC.A is missing 12 
# of 20 samples for Zostera biomass and shoot density. I am going ahead to impute the 
# missing Zostera data for BC.A but we may want to exclude this site from any analysis 
# that uses Zostera mass or shoot density as either a predictor or  response. 

# Note that a data set used for RF imputation model cannot include any rows that lack 
# data for one of the predictor variables.

# We first impute missing values for Zostera variables (AG biomass, shoot density, %N), then
# proceed to the grazer variables that depend on these grass variables

# Remind me which variables are missing data:
apply(zen2014,2,pMiss)

# create a temporary dataframe used to model (impute) the missing values
y <- zen2014

# LEAF % NITROGEN
# Use random forest to impute missing values. First build predictive model:
sum(is.na(zen2014$log10.Leaf.PercN)) # 14

leafN.rf = randomForest(log10.Leaf.PercN ~ Ocean + Coast + Latitude + Longitude 
                        + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                        + log10.chlomean
                        + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                        + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                        + log10.periphyton.mass.per.g.zostera
                        + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                        ## + log10.Leaf.PercN
                        # + pop.density.2015 + AllelicRichness
                        ,
                        na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.Leaf.PercN), 
  "log10.Leaf.PercN"] = predict(leafN.rf, y[is.na(y$log10.Leaf.PercN), ])
sum(is.na(y$log10.Leaf.PercN)) # 0


# ZOSTERA SHOOT DENSITY
# Use random forest to impute missing values. First build predictive model:
sum(is.na(zen2014$log10.Zostera.shoots.core)) # 15
shootdensity.rf = randomForest(log10.Zostera.shoots.core ~ Ocean + Coast + Latitude + Longitude 
                               + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                               + log10.chlomean
                               # + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                               + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                               # + log10.periphyton.mass.per.g.zostera  
                               # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                               + log10.Leaf.PercN 
                               # + pop.density.2015 
                               + AllelicRichness
                               ,
                               na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.Zostera.shoots.core), 
  "log10.Zostera.shoots.core"] = predict(shootdensity.rf, y[is.na(y$log10.Zostera.shoots.core), ])
sum(is.na(y$log10.Zostera.shoots.core)) # 0 


# ZOSTERA ABOVE-GROUND BIOMASS
# Use random forest to impute missing values. First build predictive model:
sum(is.na(zen2014$log10.Zostera.AG.mass)) # 24
ZAG.rf = randomForest(log10.Zostera.AG.mass ~ Ocean + Coast + Latitude + Longitude 
                      # + sst.mean + Salinity.ppt + parmean + log10.nitrate + log10.phosphate
                      + log10.chlomean 
                      # + log10.Zostera.AG.mass
                      + log10.Zostera.shoots.core
                      + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                      # + log10.periphyton.mass.per.g.zostera  
                      # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                      + log10.Leaf.PercN 
                      # + pop.density.2015 
                      + AllelicRichness
                      ,
                      na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.Zostera.AG.mass), 
  "log10.Zostera.AG.mass"] = predict(ZAG.rf, y[is.na(y$log10.Zostera.AG.mass), ])
sum(is.na(y$log10.Zostera.AG.mass)) # 0


# ZOSTERA BELOW-GROUND BIOMASS
# Use random forest to impute missing values. First build predictive model:
sum(is.na(zen2014$log10.Zostera.BG.mass)) # 29
ZBG.rf = randomForest(log10.Zostera.BG.mass ~ Ocean + Coast + Latitude + Longitude 
                      # + sst.mean + Salinity.ppt + parmean + log10.nitrate + log10.phosphate
                      + log10.chlomean 
                      # + log10.Zostera.BG.mass
                      + log10.Zostera.shoots.core
                      + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                      # + log10.periphyton.mass.per.g.zostera  
                      # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                      + log10.Leaf.PercN 
                      # + pop.density.2015 
                      + AllelicRichness
                      ,
                      na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.Zostera.BG.mass), 
  "log10.Zostera.BG.mass"] = predict(ZBG.rf, y[is.na(y$log10.Zostera.BG.mass), ])
sum(is.na(y$log10.Zostera.BG.mass)) # 0


# CRUSTACEAN mesograzer biomass
# Use random forest to impute missing values. First build predictive model:
sum(is.na(zen2014$log10.crustacean.mass.per.g.plant)) # 9
crust.rf = randomForest(log10.crustacean.mass.per.g.plant ~ Ocean + Coast + Latitude + Longitude 
                        + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                        + log10.chlomean
                        + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                        + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                        # + log10.periphyton.mass.per.g.zostera
                        # # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                        + log10.Leaf.PercN
                        # # + pop.density.2015 
                        + AllelicRichness
                        ,
                        na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.crustacean.mass.per.g.plant), 
  "log10.crustacean.mass.per.g.plant"] = predict(crust.rf, y[is.na(y$log10.crustacean.mass.per.g.plant), ])
sum(is.na(y$log10.crustacean.mass.per.g.plant)) # 0 


# GASTROPOD mesograzer biomass
# Use random forest to impute missing values. First build predictive model:
sum(is.na(zen2014$log10.gastropod.mass.per.g.plant)) # 9
gast.rf = randomForest(log10.gastropod.mass.per.g.plant ~ Ocean + Coast + Latitude + Longitude 
                       + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                       + log10.chlomean
                       + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                       + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                       # + log10.periphyton.mass.per.g.zostera
                       # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                       + log10.Leaf.PercN 
                       # + pop.density.2015 
                       + AllelicRichness
                       ,
                       na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.gastropod.mass.per.g.plant), 
  "log10.gastropod.mass.per.g.plant"] = predict(gast.rf, y[is.na(y$log10.gastropod.mass.per.g.plant), ])
sum(is.na(y$log10.gastropod.mass.per.g.plant)) # 0


# total MESOGRAZER biomass
# Use random forest to impute missing values. First build predictive model:
sum(is.na(zen2014$log10.mesograzer.mass.per.g.plant)) # 9
meso.rf = randomForest(log10.mesograzer.mass.per.g.plant ~ Ocean + Coast + Latitude + Longitude 
                       + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                       + log10.chlomean 
                       + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                       + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                       # + log10.periphyton.mass.per.g.zostera
                       # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                       + log10.Leaf.PercN
                       # + pop.density.2015
                       + AllelicRichness
                       ,
                       na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.mesograzer.mass.per.g.plant), 
  "log10.mesograzer.mass.per.g.plant"] = predict(meso.rf, y[is.na(y$log10.mesograzer.mass.per.g.plant), ])
sum(is.na(y$log10.mesograzer.mass.per.g.plant)) # 0


# total MESOGRAZER abundance
# Use random forest to impute missing values. First build predictive model:
sum(is.na(zen2014$log10.mesograzer.abund.per.g.plant)) # 9
meso.abund.rf = randomForest(log10.mesograzer.abund.per.g.plant ~ Ocean + Coast + Latitude + Longitude 
                             + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                             + log10.chlomean 
                             + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                             + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                             # + log10.periphyton.mass.per.g.zostera
                             # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                             + log10.Leaf.PercN
                             # + pop.density.2015
                             + AllelicRichness
                             ,
                             na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.mesograzer.abund.per.g.plant), 
  "log10.mesograzer.abund.per.g.plant"] = predict(meso.abund.rf, y[is.na(y$log10.mesograzer.abund.per.g.plant), ])
sum(is.na(y$log10.mesograzer.abund.per.g.plant)) # 0


# CRUSTACEAN mesograzer biomass PER AREA
# Use random forest to impute missing values. First build predictive model:
sum(is.na(zen2014$log10.crustacean.mass.per.area)) # 33
crust.area.rf = randomForest(log10.crustacean.mass.per.area ~ Ocean + Coast + Latitude + Longitude 
                             + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                             + log10.chlomean 
                             + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                             + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                             # + log10.periphyton.mass.per.g.zostera
                             # # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                             + log10.Leaf.PercN
                             # # + pop.density.2015 
                             + AllelicRichness
                             ,
                             na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.crustacean.mass.per.area), 
  "log10.crustacean.mass.per.area"] = predict(crust.area.rf, y[is.na(y$log10.crustacean.mass.per.area), ])
sum(is.na(y$log10.crustacean.mass.per.area)) # 0 


# GASTROPOD mesograzer biomass PER AREA
# Use random forest to impute missing values. First build predictive model:
sum(is.na(zen2014$log10.gastropod.mass.per.area)) # 33
gast.area.rf = randomForest(log10.gastropod.mass.per.area ~ Ocean + Coast + Latitude + Longitude 
                            + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                            + log10.chlomean
                            + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                            + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                            # + log10.periphyton.mass.per.g.zostera
                            # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                            + log10.Leaf.PercN 
                            # + pop.density.2015 
                            + AllelicRichness
                            ,
                            na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.gastropod.mass.per.area), 
  "log10.gastropod.mass.per.area"] = predict(gast.area.rf, y[is.na(y$log10.gastropod.mass.per.area), ])
sum(is.na(y$log10.gastropod.mass.per.area)) # 0


# total MESOGRAZER biomass PER AREA
# Use random forest to impute missing values. First build predictive model:
sum(is.na(zen2014$log10.mesograzer.mass.per.area)) # 33
meso.area.rf = randomForest(log10.mesograzer.mass.per.area ~ Ocean + Coast + Latitude + Longitude 
                            + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                            + log10.chlomean 
                            + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                            + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                            # + log10.periphyton.mass.per.g.zostera
                            # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                            + log10.Leaf.PercN
                            # + pop.density.2015
                            + AllelicRichness
                            ,
                            na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.mesograzer.mass.per.area), 
  "log10.mesograzer.mass.per.area"] = predict(meso.area.rf, y[is.na(y$log10.mesograzer.mass.per.area), ])
sum(is.na(y$log10.mesograzer.mass.per.area)) # 0


# total MESOGRAZER abundance PER AREA
# Use random forest to impute missing values. First build predictive model:
sum(is.na(zen2014$log10.mesograzer.abund.per.area)) # 33
meso.abund.area.rf = randomForest(log10.mesograzer.abund.per.area ~ Ocean + Coast + Latitude + Longitude 
                                  + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                                  + log10.chlomean 
                                  + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                                  + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                                  # + log10.periphyton.mass.per.g.zostera
                                  # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                                  + log10.Leaf.PercN
                                  # + pop.density.2015
                                  + AllelicRichness
                                  ,
                                  na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = y)

# Impute missing values
y[is.na(y$log10.mesograzer.abund.per.area), 
  "log10.mesograzer.abund.per.area"] = predict(meso.abund.area.rf, y[is.na(y$log10.mesograzer.abund.per.area), ])
sum(is.na(y$log10.mesograzer.abund.per.area)) # 0


# Create data frame containing the imputed values and add them to the master data frame
names(y)
imputed.values.y <- y[c("Unique.ID",  "log10.Zostera.shoots.core", "log10.Zostera.AG.mass", 
                        "log10.Zostera.BG.mass", "log10.Leaf.PercN",  
                        "log10.crustacean.mass.per.g.plant", "log10.crustacean.mass.per.area", 
                        "log10.gastropod.mass.per.g.plant", "log10.gastropod.mass.per.area",
                        "log10.mesograzer.mass.per.g.plant", "log10.mesograzer.mass.per.area", 
                        "log10.mesograzer.abund.per.g.plant", "log10.mesograzer.abund.per.area" 
)]
names(imputed.values.y)

# Rename imputed values
colnames(imputed.values.y)[2:13] <- c("log10.Zostera.shoots.core.imputed", 
                                      "log10.Zostera.AG.mass.imputed", 
                                      "log10.Zostera.BG.mass.imputed", 
                                      "log10.Leaf.PercN.imputed", 
                                      "log10.crustacean.mass.per.g.plant.imputed", 
                                      "log10.crustacean.mass.per.area.imputed", 
                                      "log10.gastropod.mass.per.g.plant.imputed", 
                                      "log10.gastropod.mass.per.area.imputed",
                                      "log10.mesograzer.mass.per.g.plant.imputed", 
                                      "log10.mesograzer.mass.per.area.imputed",
                                      "log10.mesograzer.abund.per.g.plant.imputed",
                                      "log10.mesograzer.abund.per.area.imputed"
) 


# Integrate the imputed values back into master data set 
zen2014_imputed <- zen2014

zen2014_imputed$log10.Zostera.shoots.core.imputed <- 
  imputed.values.y$log10.Zostera.shoots.core.imputed[match(zen2014_imputed$Unique.ID, imputed.values.y$Unique.ID)]

zen2014_imputed$log10.Zostera.AG.mass.imputed <- 
  imputed.values.y$log10.Zostera.AG.mass.imputed[match(zen2014_imputed$Unique.ID, imputed.values.y$Unique.ID)]

zen2014_imputed$log10.Zostera.BG.mass.imputed <- 
  imputed.values.y$log10.Zostera.BG.mass.imputed[match(zen2014_imputed$Unique.ID, imputed.values.y$Unique.ID)]

zen2014_imputed$log10.Leaf.PercN.imputed <- 
  imputed.values.y$log10.Leaf.PercN.imputed[match(zen2014_imputed$Unique.ID, imputed.values.y$Unique.ID)]

zen2014_imputed$log10.crustacean.mass.per.g.plant.imputed <- 
  imputed.values.y$log10.crustacean.mass.per.g.plant.imputed[match(zen2014_imputed$Unique.ID, imputed.values.y$Unique.ID)]

zen2014_imputed$log10.crustacean.mass.per.area.imputed <- 
  imputed.values.y$log10.crustacean.mass.per.area.imputed[match(zen2014_imputed$Unique.ID, imputed.values.y$Unique.ID)]

zen2014_imputed$log10.gastropod.mass.per.g.plant.imputed <- 
  imputed.values.y$log10.gastropod.mass.per.g.plant.imputed[match(zen2014_imputed$Unique.ID, imputed.values.y$Unique.ID)]

zen2014_imputed$log10.gastropod.mass.per.area.imputed <- 
  imputed.values.y$log10.gastropod.mass.per.area.imputed[match(zen2014_imputed$Unique.ID, imputed.values.y$Unique.ID)]

zen2014_imputed$log10.mesograzer.mass.per.g.plant.imputed <- 
  imputed.values.y$log10.mesograzer.mass.per.g.plant.imputed[match(zen2014_imputed$Unique.ID, imputed.values.y$Unique.ID)]

zen2014_imputed$log10.mesograzer.mass.per.area.imputed <- 
  imputed.values.y$log10.mesograzer.mass.per.area.imputed[match(zen2014_imputed$Unique.ID, imputed.values.y$Unique.ID)]

zen2014_imputed$log10.mesograzer.abund.per.g.plant.imputed <- 
  imputed.values.y$log10.mesograzer.abund.per.g.plant.imputed[match(zen2014_imputed$Unique.ID, imputed.values.y$Unique.ID)]

zen2014_imputed$log10.mesograzer.abund.per.area.imputed <- 
  imputed.values.y$log10.mesograzer.abund.per.area.imputed[match(zen2014_imputed$Unique.ID, imputed.values.y$Unique.ID)]

names(zen2014_imputed)
nrow(zen2014_imputed) # 1000 - good


# PERIPHYTON

# NOTE: Here we need a reduced data set of 49 sites (i.e., excluding SW.A) for predictive model and imputation,
# because ALL plots from SW.A. had no periphyton values so it s not valid to impute values for that site.
# This requires two steps: 

# First, we subset the dataframe of imputed values created above. After we derive imputed values for 
# periphyton, we will paste them back inot this dataframe: 
zen2014_49_imputed <- droplevels(subset(zen2014_imputed, Site != "SW.A"))

# Second: To rigorously estimate imputed values for a variable (in this case periphyton), we should use only 
# empirical data, so we next subset the original (pre-imputation) dataframe for this purpose:
zen2014_49 <- droplevels(subset(zen2014, Site != "SW.A"))
x <- zen2014_49

# PERIPHYTON mass per g Zostera
# Use random forest to impute missing values for periphyton in the reduced dataframe (49 sites). 
# First build predictive model:
sum(is.na(zen2014_49$log10.periphyton.mass.per.g.zostera)) # 4
peri.rf = randomForest(log10.periphyton.mass.per.g.zostera ~ Ocean + Coast + Latitude + Longitude 
                       + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                       + log10.chlomean
                       + log10.Zostera.AG.mass + log10.Zostera.shoots.core
                       + log10.Zostera.sheath.length + log10.Zostera.longest.leaf.length
                       # # + log10.periphyton.mass.per.g.zostera
                       # + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                       + log10.Leaf.PercN
                       # + pop.density.2015
                       + AllelicRichness
                       ,
                       na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = x)

# Impute missing values
x[is.na(x$log10.periphyton.mass.per.g.zostera),
  "log10.periphyton.mass.per.g.zostera"] = predict(peri.rf, x[is.na(x$log10.periphyton.mass.per.g.zostera), ])
sum(is.na(x$log10.periphyton.mass.per.g.zostera)) # 0


# PERIPHYTON mass per AREA
# Use random forest to impute missing values for periphyton in the reduced dataframe (49 sites). 
# First build predictive model:
sum(is.na(zen2014_49$log10.periphyton.mass.per.area)) # 28
peri.area.rf = randomForest(log10.periphyton.mass.per.area ~ Ocean + Coast + Latitude 
                            + Longitude
                            + sst.mean + Salinity.ppt + parmean + sqrt.nitrate + log10.phosphate
                            + log10.chlomean
                            + log10.Zostera.shoots.core
                            + log10.Zostera.sheath.length
                            + log10.Zostera.longest.leaf.length
                            # # + log10.periphyton.mass.per.g.zostera
                            + log10.mesograzer.mass.per.g.plant + log10.crustacean.mass.per.g.plant + log10.gastropod.mass.per.g.plant
                            + log10.Leaf.PercN
                            # + pop.density.2015
                            + AllelicRichness
                            ,
                            na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = x)

# Impute missing values
x[is.na(x$log10.periphyton.mass.per.area),
  "log10.periphyton.mass.per.area"] = predict(peri.area.rf, x[is.na(x$log10.periphyton.mass.per.area), ])
sum(is.na(x$log10.periphyton.mass.per.area)) # 0 (if it returns 3, then take out most predictors, and add back gradually)


# Create data frame containing the imputed value for periphyton and add to the 49-site data frame
names(x)
imputed.values.x <- x[c("Unique.ID", "log10.periphyton.mass.per.g.zostera", "log10.periphyton.mass.per.area")]
names(imputed.values.x)

# Rename imputed values
colnames(imputed.values.x)[2-3] <- c("log10.periphyton.mass.per.g.zostera.imputed", "log10.periphyton.mass.per.area.imputed" ) 
names(imputed.values.x)

# Now paste the imputed values for periphyton back into the 49-site dataframe with the other 
# imputed variables created above:
zen2014_49_imputed$log10.periphyton.mass.per.g.zostera.imputed <- 
  imputed.values.x$log10.periphyton.mass.per.g.zostera.imputed[match(zen2014_49_imputed$Unique.ID, imputed.values.x$Unique.ID)]
sum(is.na(zen2014_49_imputed$log10.periphyton.mass.per.g.zostera.imputed)) # 0

zen2014_49_imputed$log10.periphyton.mass.per.area.imputed <- 
  imputed.values.x$log10.periphyton.mass.per.area.imputed[match(zen2014_49_imputed$Unique.ID, imputed.values.x$Unique.ID)]
sum(is.na(zen2014_49_imputed$log10.periphyton.mass.per.area.imputed)) # 0

names(zen2014_imputed)
nrow(zen2014_imputed) # 1000 - good

names(zen2014_49_imputed)
nrow(zen2014_49_imputed) # 980 - good

# Summary: We can now build and compare models that will have same number of observations
# BUT can't include periphyton as predictor because biased by missing from entire site SW.A.

# Export imputed data frames
write.csv(zen2014_imputed, "data/output/zen2014_imputed_20210227.csv", row.names = F)
write.csv(zen2014_49_imputed, "data/output/zen2014_imputed_49_20210227.csv", row.names = F)


###################################################################################
# PCA - EELGRASS VARIABLES (GLOBAL)                                               #
###################################################################################

# NOTE: This includes all available ZEN eelgrass morphological variables. We use the 
# first two axes, which together explain 83% of the variation in input variables, under 
# the (arbitrary) criterion of using those PC axes necessary to capture 75% of the variation. 

# PCA - EELGRASS VARIABLES (PLOT LEVEL)                                           

# Create data frame containing the ZEN 2014 eelgrass morphological variables
zos.morph.plot.2 <- zen2014_imputed[c("log10.Zostera.AG.mass.imputed", "log10.Zostera.BG.mass.imputed", 
  "log10.Zostera.shoots.core.imputed", "log10.Zostera.sheath.length", "log10.Zostera.sheath.width", "log10.Zostera.longest.leaf.length")]

# Compute PCAs
zos.morph.plot.2.pca <- prcomp(zos.morph.plot.2, center = TRUE, scale. = TRUE) 

print(zos.morph.plot.2.pca)
#                                           PC1         PC2         PC3        PC4         PC5         PC6
# log10.Zostera.AG.mass.imputed     -0.29772190 -0.58976969  0.16131419 -0.7076165  0.12385514 -0.14645813
# log10.Zostera.BG.mass.imputed      0.08114321 -0.67078182 -0.63774621  0.3664483 -0.03986877  0.02955342
# log10.Zostera.shoots.core.imputed  0.34930322 -0.42578505  0.70199747  0.3770211  0.20963800  0.13341998
# log10.Zostera.sheath.length       -0.51441226 -0.05711932  0.21262143  0.4040899 -0.27044926 -0.67117666
# log10.Zostera.sheath.width        -0.50068037  0.09723378 -0.08264182  0.2209389  0.81254579  0.15488847
# log10.Zostera.longest.leaf.length -0.51716912 -0.09062856  0.14973149  0.1036680 -0.45359545  0.69671169

# Interpretation:
# PCz1: growth form: high = short canopy, denser shoots
# PCz2: biomass: high values = low AG and especially BG biomass


# plot cumulative proportion of variance explained by PC axes
plot(zos.morph.plot.2.pca, type = "l")

# Calculate proportion of variance explained by each PC
summary(zos.morph.plot.2.pca)
#                           PC1    PC2     PC3     PC4     PC5     PC6
# Standard deviation     1.8230 1.2796 0.71769 0.48452 0.45114 0.29318
# Proportion of Variance 0.5539 0.2729 0.08585 0.03913 0.03392 0.01433
# Cumulative Proportion  0.5539 0.8268 0.91263 0.95175 0.98567 1.00000

# RESULT: First two PC axes explain 83% of variation in eelgrass morphology with ALL input variables. 

# Output PCA scores and combine with plot data frame
zos.morph.plot.2.pca.scores <- zos.morph.plot.2.pca$x
ZEN_2014_plot <- cbind(zen2014_imputed, zos.morph.plot.2.pca.scores) 

# Rename PCA variables 1-2 and cull PC3-4
names(ZEN_2014_plot)[names(ZEN_2014_plot)=="PC1"] <- "PC1.zos"
names(ZEN_2014_plot)[names(ZEN_2014_plot)=="PC2"] <- "PC2.zos"
ZEN_2014_plot <- subset(ZEN_2014_plot, select = -c(PC3,PC4,PC5,PC6))
names(ZEN_2014_plot)


# NOTE: PROBABLY NEED NEW SECTION HEADING HERE ...

# NOTE: IS THIS WHERE THIS SHOULD BE?
# Add environmental and genetic PC scores to plot-level data frame
ZEN_2014_plot$PC1.env.global <- ZEN_2014_site_means$PC1.env.global[match(ZEN_2014_plot$Site, ZEN_2014_site_means$Site)]
ZEN_2014_plot$PC2.env.global <- ZEN_2014_site_means$PC2.env.global[match(ZEN_2014_plot$Site, ZEN_2014_site_means$Site)]
ZEN_2014_plot$PC3.env.global <- ZEN_2014_site_means$PC3.env.global[match(ZEN_2014_plot$Site, ZEN_2014_site_means$Site)]
ZEN_2014_plot$FC1 <- ZEN_2014_site_means$FC1[match(ZEN_2014_plot$Site, ZEN_2014_site_means$Site)]
ZEN_2014_plot$FC2 <- ZEN_2014_site_means$FC2[match(ZEN_2014_plot$Site, ZEN_2014_site_means$Site)]


# NOTE: IS THIS WHERE THIS SHOULD BE?
# Obtain mean values per site: Eelgrass growth form PCz1 and PCz2 
add_means <- ddply(ZEN_2014_plot, c("Site"), summarize, 
                   PC1.zos.site = mean(PC1.zos, na.rm = T),
                   PC2.zos.site = mean(PC2.zos, na.rm = T)
)

# Add to site means data frame
ZEN_2014_site_means <- merge(ZEN_2014_site_means, add_means)
names(ZEN_2014_site_means)

# Add to ocean data frames
ZEN_2014_site_means_Atlantic$PC1.zos.site <- ZEN_2014_site_means$PC1.zos.site[match(ZEN_2014_site_means_Atlantic$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_Atlantic$PC2.zos.site <- ZEN_2014_site_means$PC2.zos.site[match(ZEN_2014_site_means_Atlantic$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_Pacific$PC1.zos.site <- ZEN_2014_site_means$PC1.zos.site[match(ZEN_2014_site_means_Pacific$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_Pacific$PC2.zos.site <- ZEN_2014_site_means$PC2.zos.site[match(ZEN_2014_site_means_Pacific$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_49_Atlantic$PC1.zos.site <- ZEN_2014_site_means$PC1.zos.site[match(ZEN_2014_site_means_49_Atlantic$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_49_Atlantic$PC2.zos.site <- ZEN_2014_site_means$PC2.zos.site[match(ZEN_2014_site_means_49_Atlantic$Site, ZEN_2014_site_means$Site)]


###################################################################################
# CREATE SCALED VARIABLES                                                         #
###################################################################################

# Create function to standardize and center a variable by its range of observed values. 
# The '...' allows it to work with NAs.
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}


# Combine PCA scores with PLOT-level data frame
ZEN_2014_site_means_49_Atlantic$PC1.env.global <- ZEN_2014_site_means$PC1.env.global[match(ZEN_2014_site_means_49_Atlantic$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_49_Atlantic$PC2.env.global <- ZEN_2014_site_means$PC2.env.global[match(ZEN_2014_site_means_49_Atlantic$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_49_Atlantic$PC3.env.global <- ZEN_2014_site_means$PC3.env.global[match(ZEN_2014_site_means_49_Atlantic$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_49_Atlantic$FC1 <- ZEN_2014_site_means$FC1[match(ZEN_2014_site_means_49_Atlantic$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_49_Atlantic$FC2 <- ZEN_2014_site_means$FC2[match(ZEN_2014_site_means_49_Atlantic$Site, ZEN_2014_site_means$Site)]

ZEN_2014_site_means_Pacific$PC1.env.global <- ZEN_2014_site_means$PC1.env.global[match(ZEN_2014_site_means_Pacific$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_Pacific$PC2.env.global <- ZEN_2014_site_means$PC2.env.global[match(ZEN_2014_site_means_Pacific$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_Pacific$PC3.env.global <- ZEN_2014_site_means$PC3.env.global[match(ZEN_2014_site_means_Pacific$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_Pacific$FC1 <- ZEN_2014_site_means$FC1[match(ZEN_2014_site_means_Pacific$Site, ZEN_2014_site_means$Site)]
ZEN_2014_site_means_Pacific$FC2 <- ZEN_2014_site_means$FC2[match(ZEN_2014_site_means_Pacific$Site, ZEN_2014_site_means$Site)]

# Create z-scaled variables: SITE level (GLOBAL)
ZEN_2014_site_means$zLatitude <- scale(ZEN_2014_site_means$Latitude)
ZEN_2014_site_means$zPC1.zos.site <- scale(ZEN_2014_site_means$PC1.zos.site)
ZEN_2014_site_means$zPC2.zos.site <- scale(ZEN_2014_site_means$PC2.zos.site)
ZEN_2014_site_means$zPC1.env.global <- scale(ZEN_2014_site_means$PC1.env.global)
ZEN_2014_site_means$zPC2.env.global <- scale(ZEN_2014_site_means$PC2.env.global)
ZEN_2014_site_means$zPC3.env.global <- scale(ZEN_2014_site_means$PC3.env.global)
ZEN_2014_site_means$zFC1 <- scale(ZEN_2014_site_means$FC1)
ZEN_2014_site_means$zFC2 <- scale(ZEN_2014_site_means$FC2)
ZEN_2014_site_means$zcanopy <- scale(ZEN_2014_site_means$log10.Zostera.longest.leaf.length.cm.site)
ZEN_2014_site_means$zshoots <- scale(ZEN_2014_site_means$log10.Zostera.shoots.core.site)
ZEN_2014_site_means$zagbiomass <- scale(ZEN_2014_site_means$log10.Zostera.AG.mass.site)
ZEN_2014_site_means$zbgbiomass <- scale(ZEN_2014_site_means$log10.Zostera.BG.mass.site)
ZEN_2014_site_means$zperiphyton <- scale(ZEN_2014_site_means$log10.periphyton.mass.per.area.site)
ZEN_2014_site_means$zperiphyton.perg <- scale(ZEN_2014_site_means$log10.periphyton.mass.per.g.zostera.site)
ZEN_2014_site_means$zmesograzer.mass <- scale(ZEN_2014_site_means$log10.mesograzer.mass.per.area.site)
ZEN_2014_site_means$zmesograzer.mass.perg <- scale(ZEN_2014_site_means$log10.mesograzer.mass.per.g.plant.site)
ZEN_2014_site_means$zmesograzer.abund <- scale(ZEN_2014_site_means$log10.mesograzer.abund.per.area.site)
ZEN_2014_site_means$zmesograzer.abund.perg <- scale(ZEN_2014_site_means$log10.mesograzer.abund.per.g.plant.site)

# Create RANGE-scaled variables: SITE level (GLOBAL)
ZEN_2014_site_means$rLatitude <- range01(ZEN_2014_site_means$Latitude)
ZEN_2014_site_means$rPC1.zos.site <- range01(ZEN_2014_site_means$PC1.zos.site)
ZEN_2014_site_means$rPC2.zos.site <- range01(ZEN_2014_site_means$PC2.zos.site)
ZEN_2014_site_means$rPC1.env.global <- range01(ZEN_2014_site_means$PC1.env.global)
ZEN_2014_site_means$rPC2.env.global <- range01(ZEN_2014_site_means$PC2.env.global)
ZEN_2014_site_means$rPC3.env.global <- range01(ZEN_2014_site_means$PC3.env.global)
ZEN_2014_site_means$rFC1 <- range01(ZEN_2014_site_means$FC1)
ZEN_2014_site_means$rFC2 <- range01(ZEN_2014_site_means$FC2)
ZEN_2014_site_means$rcanopy <- range01(ZEN_2014_site_means$log10.Zostera.longest.leaf.length.cm.site)
ZEN_2014_site_means$rshoots <- range01(ZEN_2014_site_means$log10.Zostera.shoots.core.site)
ZEN_2014_site_means$ragbiomass <- range01(ZEN_2014_site_means$log10.Zostera.AG.mass.site)
ZEN_2014_site_means$rbgbiomass <- range01(ZEN_2014_site_means$log10.Zostera.BG.mass.site)
ZEN_2014_site_means$rperiphyton <- range01(ZEN_2014_site_means$log10.periphyton.mass.per.area.site)
ZEN_2014_site_means$rperiphyton.perg <- range01(ZEN_2014_site_means$log10.periphyton.mass.per.g.zostera.site)
ZEN_2014_site_means$rmesograzer.mass <- range01(ZEN_2014_site_means$log10.mesograzer.mass.per.area.site)
ZEN_2014_site_means$rmesograzer.mass.perg <- range01(ZEN_2014_site_means$log10.mesograzer.mass.per.g.plant.site)
ZEN_2014_site_means$rmesograzer.abund <- range01(ZEN_2014_site_means$log10.mesograzer.abund.per.area.site)
ZEN_2014_site_means$rmesograzer.abund.perg <- range01(ZEN_2014_site_means$log10.mesograzer.abund.per.g.plant.site)


# Create z-scaled variables: SITE level (ATLANTIC 49)
# This data set scales the variables using only Atlantic values. Omit SW.A as the plot-level data set does. 
ZEN_2014_site_means_49_Atlantic$zLatitude.atl <- scale(ZEN_2014_site_means_49_Atlantic$Latitude, scale = TRUE, center = TRUE)
ZEN_2014_site_means_49_Atlantic$zPC1.zos.atl <- scale(ZEN_2014_site_means_49_Atlantic$PC1.zos.site)
ZEN_2014_site_means_49_Atlantic$zPC2.zos.atl <- scale(ZEN_2014_site_means_49_Atlantic$PC2.zos.site)
ZEN_2014_site_means_49_Atlantic$zPC1.env.global.atl <- scale(ZEN_2014_site_means_49_Atlantic$PC1.env.global)
ZEN_2014_site_means_49_Atlantic$zPC2.env.global.atl <- scale(ZEN_2014_site_means_49_Atlantic$PC2.env.global)
ZEN_2014_site_means_49_Atlantic$zPC3.env.global.atl <- scale(ZEN_2014_site_means_49_Atlantic$PC3.env.global)
ZEN_2014_site_means_49_Atlantic$zFC1.global.atl <- scale(ZEN_2014_site_means_49_Atlantic$FC1)
ZEN_2014_site_means_49_Atlantic$zFC2.global.atl <- scale(ZEN_2014_site_means_49_Atlantic$FC2)
################################################################################
# ------------------  WHALEN FOUND ERRORS IN SEVERAL OF THE FOLLOWING LINES. 
# ------------------  NEEDS TO BE DEBUGGED
ZEN_2014_site_means_49_Atlantic$zPC1.env.atl <- scale(ZEN_2014_site_means_49_Atlantic$PC1.env.atl)
ZEN_2014_site_means_49_Atlantic$zPC2.env.atl <- scale(ZEN_2014_site_means_49_Atlantic$PC2.env.atl)
ZEN_2014_site_means_49_Atlantic$zPC3.env.atl <- scale(ZEN_2014_site_means_49_Atlantic$PC3.env.atl)
ZEN_2014_site_means_49_Atlantic$zFC1.atl <- scale(ZEN_2014_site_means_49_Atlantic$FC1.atl)
ZEN_2014_site_means_49_Atlantic$zFC2.atl <- scale(ZEN_2014_site_means_49_Atlantic$FC2.atl)
ZEN_2014_site_means_49_Atlantic$zperiphyton.area.atl <- scale(ZEN_2014_site_means_49_Atlantic$log10.periphyton.mass.per.area.site)
ZEN_2014_site_means_49_Atlantic$zperiphyton.perg.atl <- scale(ZEN_2014_site_means_49_Atlantic$log10.periphyton.mass.per.g.zostera.site)
ZEN_2014_site_means_49_Atlantic$zmesograzer.mass.area.atl <- scale(ZEN_2014_site_means_49_Atlantic$log10.mesograzer.mass.per.area.site)
ZEN_2014_site_means_49_Atlantic$zmesograzer.mass.perg.atl <- scale(ZEN_2014_site_means_49_Atlantic$log10.mesograzer.mass.per.g.plant.site)
################################################################################



# Create RANGE-scaled variables: SITE level (ATLANTIC 49)
# This data set scales the variables using only Atlantic values. Omit SW.A as the plot-level data set does. 
ZEN_2014_site_means_49_Atlantic$rLatitude.atl <- range01(ZEN_2014_site_means_49_Atlantic$Latitude)
ZEN_2014_site_means_49_Atlantic$rPC1.zos.atl <- range01(ZEN_2014_site_means_49_Atlantic$PC1.zos.site)
ZEN_2014_site_means_49_Atlantic$rPC2.zos.atl <- range01(ZEN_2014_site_means_49_Atlantic$PC2.zos.site)
ZEN_2014_site_means_49_Atlantic$rPC1.env.global.atl <- range01(ZEN_2014_site_means_49_Atlantic$PC1.env.global)
ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl <- range01(ZEN_2014_site_means_49_Atlantic$PC2.env.global)
ZEN_2014_site_means_49_Atlantic$rPC3.env.global.atl <- range01(ZEN_2014_site_means_49_Atlantic$PC3.env.global)
ZEN_2014_site_means_49_Atlantic$rFC1.global.atl <- range01(ZEN_2014_site_means_49_Atlantic$FC1)
ZEN_2014_site_means_49_Atlantic$rFC2.global.atl <- range01(ZEN_2014_site_means_49_Atlantic$FC2)
ZEN_2014_site_means_49_Atlantic$rPC1.env.atl <- range01(ZEN_2014_site_means_49_Atlantic$PC1.env.atl)
ZEN_2014_site_means_49_Atlantic$rPC2.env.atl <- range01(ZEN_2014_site_means_49_Atlantic$PC2.env.atl)
ZEN_2014_site_means_49_Atlantic$rPC3.env.atl <- range01(ZEN_2014_site_means_49_Atlantic$PC3.env.atl)
ZEN_2014_site_means_49_Atlantic$rFC1.atl <- range01(ZEN_2014_site_means_49_Atlantic$FC1.atl)
ZEN_2014_site_means_49_Atlantic$rFC2.atl <- range01(ZEN_2014_site_means_49_Atlantic$FC2.atl)
ZEN_2014_site_means_49_Atlantic$rperiphyton.area.atl <- range01(ZEN_2014_site_means_49_Atlantic$log10.periphyton.mass.per.area.site)
ZEN_2014_site_means_49_Atlantic$rperiphyton.perg.atl <- range01(ZEN_2014_site_means_49_Atlantic$log10.periphyton.mass.per.g.zostera.site)
ZEN_2014_site_means_49_Atlantic$rmesograzer.mass.area.atl <- range01(ZEN_2014_site_means_49_Atlantic$log10.mesograzer.mass.per.area.site)
ZEN_2014_site_means_49_Atlantic$rmesograzer.mass.perg.atl <- range01(ZEN_2014_site_means_49_Atlantic$log10.mesograzer.mass.per.g.plant.site)

# Create z-scaled variables: SITE level (PACIFIC)
# This data set scales the variables using only Pacific values.  
ZEN_2014_site_means_Pacific$zLatitude.pac <- scale(ZEN_2014_site_means_Pacific$Latitude, scale = TRUE, center = TRUE)
ZEN_2014_site_means_Pacific$zPC1.zos.pac <- scale(ZEN_2014_site_means_Pacific$PC1.zos.site)
ZEN_2014_site_means_Pacific$zPC2.zos.pac <- scale(ZEN_2014_site_means_Pacific$PC2.zos.site)
ZEN_2014_site_means_Pacific$zPC1.env.global.pac <- scale(ZEN_2014_site_means_Pacific$PC1.env.global)
ZEN_2014_site_means_Pacific$zPC2.env.global.pac <- scale(ZEN_2014_site_means_Pacific$PC2.env.global)
ZEN_2014_site_means_Pacific$zPC3.env.global.pac <- scale(ZEN_2014_site_means_Pacific$PC3.env.global)
ZEN_2014_site_means_Pacific$zFC1.global.pac <- scale(ZEN_2014_site_means_Pacific$FC1)
ZEN_2014_site_means_Pacific$zFC2.global.pac <- scale(ZEN_2014_site_means_Pacific$FC2)
ZEN_2014_site_means_Pacific$zPC1.env.pac <- scale(ZEN_2014_site_means_Pacific$PC1.env.pac)
ZEN_2014_site_means_Pacific$zPC2.env.pac <- scale(ZEN_2014_site_means_Pacific$PC2.env.pac)
ZEN_2014_site_means_Pacific$zPC3.env.pac <- scale(ZEN_2014_site_means_Pacific$PC3.env.pac)
ZEN_2014_site_means_Pacific$zFC1.pac <- scale(ZEN_2014_site_means_Pacific$FC1.pac)
ZEN_2014_site_means_Pacific$zFC2.pac <- scale(ZEN_2014_site_means_Pacific$FC2.pac)
ZEN_2014_site_means_Pacific$zperiphyton.area.pac <- scale(ZEN_2014_site_means_Pacific$log10.periphyton.mass.per.area.site)
ZEN_2014_site_means_Pacific$zperiphyton.perg.pac <- scale(ZEN_2014_site_means_Pacific$log10.periphyton.mass.per.g.zostera.site)
ZEN_2014_site_means_Pacific$zmesograzer.mass.area.pac <- scale(ZEN_2014_site_means_Pacific$log10.mesograzer.mass.per.area.site)
ZEN_2014_site_means_Pacific$zmesograzer.mass.perg.pac <- scale(ZEN_2014_site_means_Pacific$log10.mesograzer.mass.per.g.plant.site)

# Create RANGE-scaled variables: SITE level (PACIFIC)
# This data set scales the variables using only Pacific values.  
ZEN_2014_site_means_Pacific$rLatitude.pac <- range01(ZEN_2014_site_means_Pacific$Latitude)
ZEN_2014_site_means_Pacific$rPC1.zos.pac <- range01(ZEN_2014_site_means_Pacific$PC1.zos.site)
ZEN_2014_site_means_Pacific$rPC2.zos.pac <- range01(ZEN_2014_site_means_Pacific$PC2.zos.site)
ZEN_2014_site_means_Pacific$rPC1.env.global.pac <- range01(ZEN_2014_site_means_Pacific$PC1.env.global)
ZEN_2014_site_means_Pacific$rPC2.env.global.pac <- range01(ZEN_2014_site_means_Pacific$PC2.env.global)
ZEN_2014_site_means_Pacific$rPC3.env.global.pac <- range01(ZEN_2014_site_means_Pacific$PC3.env.global)
ZEN_2014_site_means_Pacific$rFC1.global.pac <- range01(ZEN_2014_site_means_Pacific$FC1)
ZEN_2014_site_means_Pacific$rFC2.global.pac <- range01(ZEN_2014_site_means_Pacific$FC2)
ZEN_2014_site_means_Pacific$rPC1.env.pac <- range01(ZEN_2014_site_means_Pacific$PC1.env.pac)
ZEN_2014_site_means_Pacific$rPC2.env.pac <- range01(ZEN_2014_site_means_Pacific$PC2.env.pac)
ZEN_2014_site_means_Pacific$rPC3.env.pac <- range01(ZEN_2014_site_means_Pacific$PC3.env.pac)
ZEN_2014_site_means_Pacific$rFC1.pac <- range01(ZEN_2014_site_means_Pacific$FC1.pac)
ZEN_2014_site_means_Pacific$rFC2.pac <- range01(ZEN_2014_site_means_Pacific$FC2.pac)
ZEN_2014_site_means_Pacific$rperiphyton.area.pac <- range01(ZEN_2014_site_means_Pacific$log10.periphyton.mass.per.area.site)
ZEN_2014_site_means_Pacific$rperiphyton.perg.pac <- range01(ZEN_2014_site_means_Pacific$log10.periphyton.mass.per.g.zostera.site)
ZEN_2014_site_means_Pacific$rmesograzer.mass.area.pac <- range01(ZEN_2014_site_means_Pacific$log10.mesograzer.mass.per.area.site)
ZEN_2014_site_means_Pacific$rmesograzer.mass.perg.pac <- range01(ZEN_2014_site_means_Pacific$log10.mesograzer.mass.per.g.plant.site)


# Subset PLOT data to 49 sites for modeling (SW.A has no periphyton data)
ZEN_2014_plot_49 <- droplevels(subset(ZEN_2014_plot, Site != "SW.A"))
ZEN_2014_plot_49$log10.periphyton.mass.per.g.zostera.imputed <- 
  zen2014_49_imputed$log10.periphyton.mass.per.g.zostera.imputed[match(ZEN_2014_plot_49$Unique.ID, zen2014_49_imputed$Unique.ID)]
ZEN_2014_plot_49$log10.periphyton.mass.per.area.imputed <- 
  zen2014_49_imputed$log10.periphyton.mass.per.area.imputed[match(ZEN_2014_plot_49$Unique.ID, zen2014_49_imputed$Unique.ID)]

# Create scaled variables: PLOT level (49 sites: minus SW.A with no periphyton data)
ZEN_2014_plot_49$zLatitude <- scale(ZEN_2014_plot_49$Latitude)
ZEN_2014_plot_49$zPC1.zos <- scale(ZEN_2014_plot_49$PC1.zos)
ZEN_2014_plot_49$zPC2.zos <- scale(ZEN_2014_plot_49$PC2.zos)
ZEN_2014_plot_49$zPC1.env.global <- scale(ZEN_2014_plot_49$PC1.env.global)
ZEN_2014_plot_49$zPC2.env.global <- scale(ZEN_2014_plot_49$PC2.env.global)
ZEN_2014_plot_49$zPC3.env.global <- scale(ZEN_2014_plot_49$PC3.env.global)
ZEN_2014_plot_49$zFC1 <- scale(ZEN_2014_plot_49$FC1)
ZEN_2014_plot_49$zFC2 <- scale(ZEN_2014_plot_49$FC2)
ZEN_2014_plot_49$zmeso.mass.perg <- scale(ZEN_2014_plot_49$log10.mesograzer.mass.per.g.plant.imputed)
ZEN_2014_plot_49$zmeso.mass.area <- scale(ZEN_2014_plot_49$log10.mesograzer.mass.per.area.imputed)
ZEN_2014_plot_49$zmeso.abund.perg <- scale(ZEN_2014_plot_49$log10.mesograzer.abund.per.g.plant.imputed)
ZEN_2014_plot_49$zmeso.abund.area <- scale(ZEN_2014_plot_49$log10.mesograzer.abund.per.area.imputed)
ZEN_2014_plot_49$zperiphyton.area <- scale(ZEN_2014_plot_49$log10.periphyton.mass.per.area.imputed)
ZEN_2014_plot_49$zperiphyton.perg <- scale(ZEN_2014_plot_49$log10.periphyton.mass.per.g.zostera.imputed)

# Create RANGE-scaled variables: PLOT level (49 sites: minus SW.A with no periphyton data)
ZEN_2014_plot_49$rLatitude <- range01(ZEN_2014_plot_49$Latitude)
ZEN_2014_plot_49$rPC1.zos <- range01(ZEN_2014_plot_49$PC1.zos)
ZEN_2014_plot_49$rPC2.zos <- range01(ZEN_2014_plot_49$PC2.zos)
ZEN_2014_plot_49$rPC1.env.global <- range01(ZEN_2014_plot_49$PC1.env.global)
ZEN_2014_plot_49$rPC2.env.global <- range01(ZEN_2014_plot_49$PC2.env.global)
ZEN_2014_plot_49$rPC3.env.global <- range01(ZEN_2014_plot_49$PC3.env.global)
ZEN_2014_plot_49$rFC1 <- range01(ZEN_2014_plot_49$FC1)
ZEN_2014_plot_49$rFC2 <- range01(ZEN_2014_plot_49$FC2)
ZEN_2014_plot_49$rmeso.mass.perg <- range01(ZEN_2014_plot_49$log10.mesograzer.mass.per.g.plant.imputed)
ZEN_2014_plot_49$rmeso.mass.area <- range01(ZEN_2014_plot_49$log10.mesograzer.mass.per.area.imputed)
ZEN_2014_plot_49$rmeso.abund.perg <- range01(ZEN_2014_plot_49$log10.mesograzer.abund.per.g.plant.imputed)
ZEN_2014_plot_49$rmeso.abund.area <- range01(ZEN_2014_plot_49$log10.mesograzer.abund.per.area.imputed)
ZEN_2014_plot_49$rperiphyton.area <- range01(ZEN_2014_plot_49$log10.periphyton.mass.per.area.imputed)
ZEN_2014_plot_49$rperiphyton.perg <- range01(ZEN_2014_plot_49$log10.periphyton.mass.per.g.zostera.imputed)

names(ZEN_2014_plot_49)

# NOTE: for following, add FC1 and FC2 to the plot data set above. Then add (and scalke) the ocean-specific genetc scores below.
# Create separate PLOT-level data sets by Ocean (49) and add ocean PC scores 
ZEN_2014_plot_49_Atlantic <- droplevels(subset(ZEN_2014_plot_49, Ocean == "Atlantic"))
ZEN_2014_plot_49_Atlantic$FC1.atl <- zen2014_gen_fca_atlantic$fca1.atl[match(ZEN_2014_plot_49_Atlantic$Site, zen2014_gen_fca_atlantic$Site)]
ZEN_2014_plot_49_Atlantic$FC2.atl <- zen2014_gen_fca_atlantic$fca2.atl[match(ZEN_2014_plot_49_Atlantic$Site, zen2014_gen_fca_atlantic$Site)]
ZEN_2014_plot_49_Atlantic$PC1.env.atl <- ZEN_2014_site_means_Atlantic$PC1.env.atl[match(ZEN_2014_plot_49_Atlantic$Site, ZEN_2014_site_means_Atlantic$Site)]
ZEN_2014_plot_49_Atlantic$PC2.env.atl <- ZEN_2014_site_means_Atlantic$PC2.env.atl[match(ZEN_2014_plot_49_Atlantic$Site, ZEN_2014_site_means_Atlantic$Site)]
ZEN_2014_plot_49_Atlantic$PC3.env.atl <- ZEN_2014_site_means_Atlantic$PC3.env.atl[match(ZEN_2014_plot_49_Atlantic$Site, ZEN_2014_site_means_Atlantic$Site)]

ZEN_2014_plot_49_Pacific <- droplevels(subset(ZEN_2014_plot_49, Ocean == "Pacific"))
ZEN_2014_plot_49_Pacific$FC1.pac <- zen2014_gen_fca_pacific$fca1.pac[match(ZEN_2014_plot_49_Pacific$Site, zen2014_gen_fca_pacific$Site)]
ZEN_2014_plot_49_Pacific$FC2.pac <- zen2014_gen_fca_pacific$fca2.pac[match(ZEN_2014_plot_49_Pacific$Site, zen2014_gen_fca_pacific$Site)]
ZEN_2014_plot_49_Pacific$PC1.env.pac <- ZEN_2014_site_means_Pacific$PC1.env.pac[match(ZEN_2014_plot_49_Pacific$Site, ZEN_2014_site_means_Pacific$Site)]
ZEN_2014_plot_49_Pacific$PC2.env.pac <- ZEN_2014_site_means_Pacific$PC2.env.pac[match(ZEN_2014_plot_49_Pacific$Site, ZEN_2014_site_means_Pacific$Site)]
ZEN_2014_plot_49_Pacific$PC3.env.pac <- ZEN_2014_site_means_Pacific$PC3.env.pac[match(ZEN_2014_plot_49_Pacific$Site, ZEN_2014_site_means_Pacific$Site)]


###################################################################################
# INTEGRATE DATA FRAMES                                                           #
###################################################################################

# NOTE: ADD IN RANGE-SCALED VARIABLES AS ABOVE

# Create z-scaled variables - PLOT level: ATLANTIC
ZEN_2014_plot_49_Atlantic$zFC1.atl <- scale(ZEN_2014_plot_49_Atlantic$FC1.atl,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zFC2.atl <- scale(ZEN_2014_plot_49_Atlantic$FC2.atl,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zFC1.global.atl <- scale(ZEN_2014_plot_49_Atlantic$FC1,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zFC2.global.atl <- scale(ZEN_2014_plot_49_Atlantic$FC2,scale=TRUE,center=TRUE)

ZEN_2014_plot_49_Atlantic$zPC1.env.global.atl <- scale(ZEN_2014_plot_49_Atlantic$PC1.env.global,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zPC2.env.global.atl <- scale(ZEN_2014_plot_49_Atlantic$PC2.env.global,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zPC3.env.global.atl <- scale(ZEN_2014_plot_49_Atlantic$PC3.env.global,scale=TRUE,center=TRUE)

ZEN_2014_plot_49_Atlantic$zPC1.env.atl <- scale(ZEN_2014_plot_49_Atlantic$PC1.env.atl,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zPC2.env.atl <- scale(ZEN_2014_plot_49_Atlantic$PC2.env.atl,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zPC3.env.atl <- scale(ZEN_2014_plot_49_Atlantic$PC3.env.atl,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zPC1.zos.global.atl <- scale(ZEN_2014_plot_49_Atlantic$PC1.zos,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zPC2.zos.global.atl <- scale(ZEN_2014_plot_49_Atlantic$PC2.zos,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zperi.area.atl <- scale(ZEN_2014_plot_49_Atlantic$log10.periphyton.mass.per.area.imputed,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zperi.perg.atl <- scale(ZEN_2014_plot_49_Atlantic$log10.periphyton.mass.per.g.zostera.imputed,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zmeso.area.atl <- scale(ZEN_2014_plot_49_Atlantic$log10.mesograzer.mass.per.area.imputed,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Atlantic$zmeso.perg.atl <- scale(ZEN_2014_plot_49_Atlantic$log10.mesograzer.mass.per.g.plant.imputed,scale=TRUE,center=TRUE)

# Create z-scaled variables - PLOT level: PACIFIC
ZEN_2014_plot_49_Pacific$zFC1.pac <- scale(ZEN_2014_plot_49_Pacific$FC1.pac,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zFC2.pac <- scale(ZEN_2014_plot_49_Pacific$FC2.pac,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zFC1.global.pac <- scale(ZEN_2014_plot_49_Pacific$FC1,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zFC2.global.pac <- scale(ZEN_2014_plot_49_Pacific$FC2,scale=TRUE,center=TRUE)

ZEN_2014_plot_49_Pacific$zPC1.env.global.pac <- scale(ZEN_2014_plot_49_Pacific$PC1.env.global,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zPC2.env.global.pac <- scale(ZEN_2014_plot_49_Pacific$PC2.env.global,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zPC3.env.global.pac <- scale(ZEN_2014_plot_49_Pacific$PC3.env.global,scale=TRUE,center=TRUE)

ZEN_2014_plot_49_Pacific$zPC1.env.pac <- scale(ZEN_2014_plot_49_Pacific$PC1.env.pac,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zPC2.env.pac <- scale(ZEN_2014_plot_49_Pacific$PC2.env.pac,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zPC3.env.pac <- scale(ZEN_2014_plot_49_Pacific$PC3.env.pac,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zPC1.zos.global.pac <- scale(ZEN_2014_plot_49_Pacific$PC1.zos,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zPC2.zos.global.pac <- scale(ZEN_2014_plot_49_Pacific$PC2.zos,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zperi.area.pac <- scale(ZEN_2014_plot_49_Pacific$log10.periphyton.mass.per.area.imputed,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zperi.perg.pac <- scale(ZEN_2014_plot_49_Pacific$log10.periphyton.mass.per.g.zostera.imputed,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zmeso.area.pac <- scale(ZEN_2014_plot_49_Pacific$log10.mesograzer.mass.per.area.imputed,scale=TRUE,center=TRUE)
ZEN_2014_plot_49_Pacific$zmeso.perg.pac <- scale(ZEN_2014_plot_49_Pacific$log10.mesograzer.mass.per.g.plant.imputed,scale=TRUE,center=TRUE)


# # Add periphyton data back to the PCA data frame (must be a a better way to do this ...)
# ZEN_2014_plot_49_PCA$log10.periphyton.mass.per.area.imputed <- zen2014_49_imputed$log10.periphyton.mass.per.area.imputed[match(ZEN_2014_plot_49_PCA$Unique.ID, zen2014_49_imputed$Unique.ID)]
# ZEN_2014_plot_49_PCA$log10.periphyton.mass.per.g.zostera.imputed <- zen2014_49_imputed$log10.periphyton.mass.per.g.zostera.imputed[match(ZEN_2014_plot_49_PCA$Unique.ID, zen2014_49_imputed$Unique.ID)]


###################################################################################
# CREATE REDUCED DATA FRAMES WITHOUT MISSING SITES OR VALUES                      #
###################################################################################

# NOTE: DO THIS AFTER EVERTYTHING ELSE, IMMEDIATELY BEFORE SUBSEETING GEOGRAPHICALLY AND WRITING OUT FINAL DATA FRAMES

# NOTE: piecewiseSEM apparently kills model fitting when there are NAs in dataframe, 
# even if they are not included in the model. Therefore, I create new dataframes that 
# include ONLY the variables that will be used in modeling:

# 49-SITE DATA FRAME (deletes SW.A, which has no periphyton data)
ZEN_2014_plot_49_noNA <- subset(ZEN_2014_plot_49, select = c(Site, Latitude, Longitude, 
  Coast, Ocean, sst.mean, sst.range, Salinity.ppt, ph, parmean, cloudmean, log10.day.length, 
  sqrt.nitrate, log10.phosphate, log10.chlomean, log10.Leaf.PercN.imputed, log10.mean.fetch, 
  AllelicRichness, log10.Zostera.AG.mass.imputed, log10.Zostera.BG.mass.imputed,
  log10.Zostera.shoots.core.imputed, log10.Zostera.longest.leaf.length, log10.Zostera.sheath.length, log10.Zostera.sheath.width, 
  log10.mesograzer.abund.per.g.plant.imputed, log10.mesograzer.abund.per.area.imputed, 
  log10.crustacean.mass.per.g.plant.imputed,  log10.gastropod.mass.per.g.plant.imputed,
  log10.mesograzer.mass.per.g.plant.imputed, log10.periphyton.mass.per.g.zostera.imputed, 
  log10.crustacean.mass.per.area.imputed, log10.gastropod.mass.per.area.imputed, 
  log10.mesograzer.mass.per.area.imputed, log10.periphyton.mass.per.area.imputed,
  log10.grazer.richness.site, predation.squidpops,
  PC1.env.global, PC2.env.global, PC3.env.global, FC1, FC2, PC1.zos, PC2.zos,
  zPC1.zos, zPC2.zos, zPC1.env.global, zPC2.env.global,zPC3.env.global, zFC1, zFC2,
  zperiphyton.perg, zmeso.abund.perg, zmeso.mass.perg, 
  zperiphyton.area, zmeso.abund.area, zmeso.mass.area
  ))

# Any NAs in this dataframe?
apply(ZEN_2014_plot_49_noNA, 2, function(x) any(is.na(x))) # predation.squidpops = "TRUE"


###################################################################################
# SUBSET DATA SETS BY GEOGRAPHY                                                   #
###################################################################################

# # Include only first sampling time
# zen2014.epibiota = subset(zen2014.epibiota, Sampling.Time == "1")

# # Create separate dataframes for each coast
# zen2014_imputed_WA <- subset(zen2014_49_imputed, Coast == "West Atlantic")
# zen2014_imputed_EA <- subset(zen2014_49_imputed, Coast == "East Atlantic")
# zen2014_imputed_WP <- subset(zen2014_49_imputed, Coast == "West Pacific")
# zen2014_imputed_EP <- subset(zen2014_49_imputed, Coast == "East Pacific")


# Create reduced data sets

# # Create separate data sets by Ocean
# ZEN_2014_site_means_PCA_Atlantic <- droplevels(subset(ZEN_2014_site_means_PCA, Ocean == "Atlantic"))
# ZEN_2014_site_means_PCA_Pacific <- droplevels(subset(ZEN_2014_site_means_PCA, Ocean == "Pacific"))

# # Delete SW.A because no periphyton data
# ZEN_2014_site_means_49_PCA <- droplevels(subset(ZEN_2014_site_means_PCA, Site != "SW.A"))

# # Create separate data sets by Ocean (49)
# ZEN_2014_site_means_49_PCA_Atlantic <- droplevels(subset(ZEN_2014_site_means_49_PCA, Ocean == "Atlantic"))
# ZEN_2014_site_means_49_PCA_Pacific <- droplevels(subset(ZEN_2014_site_means_49_PCA, Ocean == "Pacific"))

# # Create separate data set excluding SW.A (no periphyton data)
ZEN_2014_site_means_49 <- droplevels(subset(ZEN_2014_site_means, Site != "SW.A"))


###################################################################################
# OUTPUT CURATED DATA SETS                                                        #
###################################################################################

# Export zen2014 PLOT-level data set for modeling (MINUS SW.A): missing data imputed
write.csv(ZEN_2014_plot, "ZEN_2014_plot_20210430.csv", row.names = F)
write.csv(ZEN_2014_plot_49_noNA, "ZEN_2014_plot_49_noNA_20210227.csv", row.names = F)

# Export zen2014 PLOT-level data set, separately by ocean, for modeling (MINUS SW.A): missing data imputed
write.csv(ZEN_2014_plot_49_Atlantic, "ZEN_2014_plot_49_Atlantic_20210227.csv", row.names = F)
write.csv(ZEN_2014_plot_49_Pacific, "ZEN_2014_plot_49_Pacific_20210227.csv", row.names = F)

# Export zen2014 PLOT-level data set for all 50 sites (includes NAs) 
write.csv(zen2014_imputed, "zen2014_imputed_20210227.csv", row.names = F)

# Export zen2014 SITE-level data set
write.csv(ZEN_2014_site_means, "ZEN_2014_site_means_20220529.csv", row.names = F)

write.csv(ZEN_2014_site_means_Atlantic, "ZEN_2014_site_means_Atlantic_20220529.csv", row.names = F)
write.csv(ZEN_2014_site_means_49_Atlantic, "ZEN_2014_site_means_49_Atlantic_20210314.csv", row.names = F)
write.csv(ZEN_2014_site_means_Pacific, "ZEN_2014_site_means_Pacific_20210314.csv", row.names = F)






