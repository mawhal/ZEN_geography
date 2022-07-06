###################################################################################
#                                                                                ##
# ZEN 2014 Global eelgrass ecosystem structure: Data assembly                    ##
# RAW data are current as of 2017-04-24                                          ##
# Emmett Duffy (duffye@si.edu)                                                   ##  
# updated 2022-06-28 by Matt Whalen                                                        ##
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
# PCA - EELGRASS VARIABLES (GLOBAL)                                               #
# CREATE SCALED VARIABLES                                                         #
# SUBSET DATA SETS BY GEOGRAPHY                                                   #
# OUTPUT CURATED DATA SETS                                                        #
#                                                                                 #
###################################################################################

###################################################################################
# METADATA                                                                        #
###################################################################################

# This script assembles raw data from the ZEN 2014 global eelgrass ecosystem sampling 
# project, and outputs data files for use in modeling and other applications. See also:
#   ZEN_2014_model_comparison.R: for data exploration and model building
#   ZEN_2014_figures.R series: for building figures for the MS

# Source data: For most of the history of this script I was using 
# ZEN_2014_Site&PlotData_2016_05_17_Released.xlsx.


###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages:
library(tidyverse) # for reformatting epibiota data
library(randomForest) # needed for data imputation
library(car) # needed or vif analysis
library(psych) # to visualize relationshiops in pairs panels
library(plyr) # to use ddply below in fixing richness values



###################################################################################
# READ AND PREPARE DATA                                                           #
###################################################################################

# MAIN ZEN 2014 DATA SET
# Read in summary data set for ZEN 2014:
d <- read.csv("data/input/ZEN_2014_main_data.csv", header = TRUE)

# General site data
sites <- read.csv("data/input/ZEN_2014_site_metadata.csv", header = TRUE)

# BIO-ORACLE CLIMATE AND ENVIRONMENTAL DATA
# Read in Bio-ORACLE and WorldClim environmental data for ZEN sites from Matt Whalen's script:
env <- read.csv("data/output/ZEN_2014_environmental.csv", header = TRUE)
# add in situ data
env.insitu <- read.csv("data/input/ZEN_2014_environmental_in_situ.csv") %>% 
  mutate(site=Site)
env <- left_join(env, env.insitu)



# EELGRASS GENETICS
d.gen_fca <- read.csv("data/input/ZEN_2014_FCA_scores.csv", header = TRUE)
# d.gen_fca_atlantic <- read.csv("data/input/ZEN_2014_fca_scores_atlantic_20210125_copy.csv", header = TRUE)
# d.gen_fca_pacific <- read.csv("data/input/ZEN_2014_fca_scores_pacific_20210125_copy.csv", header = TRUE)



#### CLEAN UP AND CONSOLIDATE
# Convert categorical variables to factors
d$Site.Code <- as.factor(d$Site.Code)
d$Ocean <- as.factor(d$Ocean)

# Rename Long Island sites 
d$Site <- as.factor(d$Site)
levels(d$Site)[levels(d$Site)=="LI.1"] <- "LI.A"
levels(d$Site)[levels(d$Site)=="LI.2"] <- "LI.B"


# Rename misspelled or confusing variables
names(d)[names(d)=="Mean.Sheath.Width.cm."] <- "Zostera.sheath.width"
names(d)[names(d)=="Mean.Shealth.Length.cm."] <- "Zostera.sheath.length"
names(d)[names(d)=="Mean.Longest.Leaft.Length.cm."] <- "Zostera.longest.leaf.length"
names(d)[names(d)=="Mean.Above.Zmarina.g"] <- "Zostera.aboveground.mean.mass"
names(d)[names(d)=="Mean.Below.Zmarina.g"] <- "Zostera.belowground.mean.mass"
names(d)[names(d)=="Shoots.Zmarina.per.m2"] <- "Zostera.shoots.per.m2.core"
names(d)[names(d)=="Mean.Fetch"] <- "mean.fetch"
names(d)[names(d)=="PopDens2"] <- "pop.density.2015"
names(d)[names(d)=="mesograzer.total.site.richness"] <- "grazer.richness.site"





# MESOGRAZER SITE RICHNESS: FIX MISSING VALUES 
# Create vector of plots with missing values to see what is missing:
missing.richness <- d[is.na(d$grazer.richness.site), c(3,7)] # columns 3 and 7 are Site, Unique.ID
# replace all site richness values with "mean" for that site. First, create vector of means:
temp <- d %>% 
  group_by( Site) %>% 
  summarize( grazer.richness.site = mean(grazer.richness.site, na.rm = T))
# But CR.A has NO mesograzers at all so returns NaN.  Assume species pool is same as for CR.B (S = 3) and replace:
# temp$grazer.richness.site[is.na(temp$grazer.richness.site)] <- 3 # CR.A grazer richness now  = 3
temp$grazer.richness.site[temp$Site == "CR.A" ] <- 3 # CR.A grazer richness now  = 3

d$grazer.richness.site <- temp$grazer.richness.site[match(d$Site, temp$Site)]



# Add BioOracle environmental data to main ZEN dataframe:
d$sst.min <- env$sstmin[match(d$Site, env$Site)]
d$sst.mean <- env$sstmean[match(d$Site, env$Site)]
d$sst.max <- env$sstmax[match(d$Site, env$Site)]
d$sst.range <- env$sstrange[match(d$Site, env$Site)]
d$chlomean <- env$chlomean[match(d$Site, env$Site)]
d$nitrate <- env$nitrate[match(d$Site, env$Site)]
d$parmean <- env$parmean[match(d$Site, env$Site)]
d$cloudmean <- env$cloudmean[match(d$Site, env$Site)]
d$day.length <- env$Day.length.hours[match(d$Site, env$Site)]
d$ph <- env$ph[match(d$Site, env$Site)]
d$phosphate <- env$phosphate[match(d$Site, env$Site)]
d$salinity <- env$salinity[match(d$Site, env$Site)]
d$precipitation <- env$precip[match(d$Site, env$Site)]

# Reorder variables 'Coast': WP to EA 
d$Coast <- as.factor(d$Coast)
d$Coast <- factor(d$Coast, levels = c("West Pacific", "East Pacific", "West Atlantic", "East Atlantic"))



###################################################################################
# CREATE DERIVED VARIABLES                                                        #
###################################################################################


# Percentage of crustaceans and gastropods among the mesograzers
d$crust.pct.mass <-  d$Malacostraca.mesograzer.plot.biomass.std.mg.g / d$mesograzer.total.plot.biomass.std.mg.g
d$gast.pct.mass <-  d$Gastropoda.mesograzer.plot.biomass.std.mg.g / d$mesograzer.total.plot.biomass.std.mg.g

# grazer and periphyton nunmbers per unit bottom area (i.e., core)
d$mesograzer.abund.per.area <-  d$mesograzer.total.plot.abund.std.g * d$Zostera.aboveground.mean.mass
d$crustacean.mass.per.area <-  d$Malacostraca.mesograzer.plot.biomass.std.mg.g * d$Zostera.aboveground.mean.mass
d$gastropod.mass.per.area <-  d$Gastropoda.mesograzer.plot.biomass.std.mg.g * d$Zostera.aboveground.mean.mass
d$mesograzer.mass.per.area <-  d$mesograzer.total.plot.biomass.std.mg.g * d$Zostera.aboveground.mean.mass
d$periphyton.mass.per.area <-  d$periphyton.mass.per.g.zostera * d$Zostera.aboveground.mean.mass

# Leaf C:N ratio
d$leaf.CN.ratio <-  d$Leaf.PercC / d$Leaf.PercN


###################################################################################
# EXPLORE DISTRIBUTIONS OF VARIABLES (PLOT LEVEL)                                 #
###################################################################################

# Examine frequency distribution of sites by environmental factor
# par(mfrow = c(1,1))
# par(mfrow = c(2,4))
# hist(d$Latitude, col = "cyan", main = "Surveys by latitude")    
# hist(d$Longitude, col = "cyan", main = "Surveys by longitude")    
# hist(d$Temperature.C, col = "cyan", main = "Surveys by temperature")    
# hist(d$Salinity.ppt, col = "cyan", main = "Surveys by salinity")    
# hist(d$pop.density.2015, col = "cyan", main = "Surveys by population density")    
# hist(d$day.length, col = "cyan", main = "Surveys by day length")    
# hist(d$mean.fetch, col = "cyan", main = "Surveys by mean fetch")    
# 
# hist(d$Zostera.aboveground.mean.mass, col = "cyan", main = "Surveys by Zostera AG biomass")    
# hist(d$periphyton.mass.per.g.zostera, col = "cyan", main = "Surveys by periphyton biomass")    
# hist(d$Malacostraca.mesograzer.plot.abund.std.g, col = "cyan", main = "Surveys by crustacean biomass")    
# hist(d$Gastropoda.mesograzer.plot.biomass.std.mg.g, col = "cyan", main = "Surveys by gastropod biomass")    
# hist(d$grazer.richness.site, col = "cyan", main = "Surveys by mesograzer richness")    
# 
# hist(d$mesograzer.total.plot.biomass.std.mg.g, col = "cyan", main = "Surveys by mesograzer biomass")    
# hist(d$epifauna.total.plot.biomass.std.mg.g, col = "cyan", main = "Surveys by mobile epifauna biomass")    
# 


###################################################################################
# LOG TRANSFORMS                                                                  #
###################################################################################

# NOTE: For many variables I add a constant roughly equal to the smallest value recorded 

d$log10.Zostera.AG.mass <- log10(d$Zostera.aboveground.mean.mass + 1) 
d$log10.Zostera.BG.mass <- log10(d$Zostera.belowground.mean.mass + 1) 
d$log10.Zostera.shoots.core <- log10(d$Zostera.shoots.per.m2.core) 
d$log10.Zostera.sheath.width <- log10(d$Zostera.sheath.width) 
d$log10.Zostera.sheath.length <- log10(d$Zostera.sheath.length) 
d$log10.Zostera.longest.leaf.length <- log10(d$Zostera.longest.leaf.length) 

d$log10.epibiota.filter <- log10(d$epibiota.filter) 
d$log10.epibiota.zostera.marina <- log10(d$epibiota.zostera.marina) 
d$log10.periphyton.mass.per.g.zostera <- log10(d$periphyton.mass.per.g.zostera + 0.001) 
d$log10.periphyton.mass.per.area <- log10(d$periphyton.mass.per.area + 0.1) 

d$log10.mesograzer.abund.per.g.plant <- log10(d$mesograzer.total.plot.abund.std.g + 0.01) 
d$log10.crustacean.abund.per.g.plant <- log10(d$Malacostraca.mesograzer.plot.abund.std.g + 0.01) 
d$log10.gastropod.abund.per.g.plant <- log10(d$Gastropoda.mesograzer.plot.abund.std.g + 0.01) 

d$log10.mesograzer.mass.per.g.plant <- log10(d$mesograzer.total.plot.biomass.std.mg.g + 0.01) 
d$log10.crustacean.mass.per.g.plant <- log10(d$Malacostraca.mesograzer.plot.biomass.std.mg.g + 0.01) 
d$log10.gastropod.mass.per.g.plant <- log10(d$Gastropoda.mesograzer.plot.biomass.std.mg.g + 0.01) 

d$log10.mesograzer.abund.per.area <-  log10(d$mesograzer.abund.per.area + 1) 

d$log10.crustacean.mass.per.area <-  log10(d$crustacean.mass.per.area + 1) 
d$log10.gastropod.mass.per.area <-  log10(d$gastropod.mass.per.area + 1) 
d$log10.mesograzer.mass.per.area <-  log10(d$mesograzer.mass.per.area + 1) 

d$log10.grazer.richness.site <- log10(d$grazer.richness.site + 1) 

d$log10.day.length <- log10(d$day.length) 
d$log10.Leaf.PercN <- log10(d$Leaf.PercN) 
d$sqrt.nitrate <- sqrt(d$nitrate) 
d$log10.phosphate <- log10(d$phosphate) 
d$log10.chlomean <- log10(d$chlomean) 
d$log10.mean.fetch <- log10(d$mean.fetch) 



# hist(d$nitrate)
# hist(d$sqrt.nitrate)
# 
# hist(d$log10.Zostera.AG.mass)
# 

# Change values of NaN to NA:
d[d == "NaN"] = NA 


###################################################################################
# OBTAIN SITE MEANS                                                               #
###################################################################################

# CAN THIS GO AFTER IMPUTATION SECTION? SHOULD IT? 

# Obtain mean values per site
ZEN_2014_site_means <- d %>% 
  group_by(Site) %>% 
  dplyr::summarize( Zostera.AG.mass.site = mean(Zostera.aboveground.mean.mass, na.rm = T), 
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
             
             log10.mesograzer.abund.per.area.site = mean(log10.mesograzer.abund.per.area, na.rm = T),                  
             
             log10.mesograzer.mass.per.area.site = mean(log10.mesograzer.mass.per.area, na.rm = T),                      
             log10.crustacean.mass.per.area.site = mean(log10.crustacean.mass.per.area, na.rm = T),                      
             log10.gastropod.mass.per.area.site = mean(log10.gastropod.mass.per.area, na.rm = T),                      
             log10.periphyton.mass.per.area.site = mean(log10.periphyton.mass.per.area, na.rm = T),                      
             
             log10.Leaf.PercN.site = mean(log10.Leaf.PercN, na.rm = T)  )


ZEN_2014_site_means$grazer.richness.site <- d$grazer.richness.site[match(ZEN_2014_site_means$Site, d$Site)]

# Change values of NaN to NA:
ZEN_2014_site_means[ZEN_2014_site_means == "NaN"] = NA 

# Add site-level environmental (and other) variables back in
ZEN_2014_site_means$Ocean <- d$Ocean[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$Coast <- d$Coast[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$Latitude <- d$Latitude[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$Longitude <- d$Longitude[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$Temperature.C <- d$Temperature.C[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$Salinity.ppt <- d$Salinity.ppt[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$log10.mean.fetch <- d$log10.mean.fetch[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$day.length <- d$day.length[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$log10.day.length <- d$log10.day.length[match(ZEN_2014_site_means$Site, d$Site)]

ZEN_2014_site_means$sst.min <- d$sst.min[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$sst.mean <- d$sst.mean[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$sst.max <- d$sst.max[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$sst.range <- d$sst.range[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$salinity <- d$salinity[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$parmean <- d$parmean[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$cloudmean <- d$cloudmean[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$precipitation <- d$precipitation[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$nitrate <- d$nitrate[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$sqrt.nitrate <- d$sqrt.nitrate[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$ph <- d$ph[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$phosphate <- d$phosphate[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$log10.phosphate <- d$log10.phosphate[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$NP.ratio <- d$NP.ratio[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$chlomean <- d$chlomean[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$log10.chlomean <- d$log10.chlomean[match(ZEN_2014_site_means$Site, d$Site)]
ZEN_2014_site_means$pop.density.2015 <- d$pop.density.2015[match(ZEN_2014_site_means$Site, d$Site)]


# Add genetic data to site means data frame
ZEN_2014_site_means$FC1 <- d.gen_fca$FC1[match(ZEN_2014_site_means$Site, d.gen_fca$Site)]
ZEN_2014_site_means$FC2 <- d.gen_fca$FC2[match(ZEN_2014_site_means$Site, d.gen_fca$Site)]

# For boxplots, reorder variable 'Coast': WP to EA
ZEN_2014_site_means$Coast <- factor(ZEN_2014_site_means$Coast, levels = c("West Pacific", "East Pacific", "West Atlantic", "East Atlantic"))

# Create separate data sets by Ocean - SITE level
ZEN_2014_site_means_Atlantic <- droplevels(subset(ZEN_2014_site_means, Ocean == "Atlantic"))
ZEN_2014_site_means_Pacific <- droplevels(subset(ZEN_2014_site_means, Ocean == "Pacific"))
ZEN_2014_site_means_49_Atlantic <- droplevels(subset(ZEN_2014_site_means_Atlantic, Site != "SW.A"))



###################################################################################
# PCA - ENVIRONMENTAL VARIABLES (GLOBAL)                                          #
###################################################################################

# # Explore correlations among environmental drivers
# pairs.panels(ZEN_2014_site_means[,c("Latitude", "sst.mean", "sst.range", "sst.min", "sst.max", "Salinity.ppt", 
#   "parmean", "log10.day.length", "cloudmean", "precipitation", "sqrt.nitrate", "log10.phosphate", "log10.chlomean", 
#   "Leaf.PercN.site", "log10.mean.fetch")], 
#   smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 8)

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

# print(ZEN.env.pca)
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

# # plot cumulative proportion of variance explained by PC axes
# plot(ZEN.env.pca, type = "l")

# # Calculate proportion of variance explained by each PC
# summary(ZEN.env.pca)
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



###################################################################################
# PCA - ENVIRONMENTAL VARIABLES (ATLANTIC)                                        #
###################################################################################

# # Explore correlations among environmental drivers
# pairs.panels(ZEN_2014_site_means_Atlantic[,c("sst.mean", "sst.range", "Salinity.ppt", "parmean", 
#   "cloudmean", "log10.phosphate", "log10.chlomean", "Leaf.PercN.site"
#   # , "precipitation", "log10.day.length"
#   )], 
#   smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 8)

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

# print(ZEN.env.pca.atl)
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

# # plot cumulative proportion of variance explained by PC axes
# plot(ZEN.env.pca.atl, type = "l")

# # Calculate proportion of variance explained by each PC
# summary(ZEN.env.pca.atl)
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


###################################################################################
# PCA - ENVIRONMENTAL VARIABLES (PACIFIC)                                         #
###################################################################################

# # Explore correlations among environmental drivers
# pairs.panels(ZEN_2014_site_means_Pacific[,c("Latitude", "sst.mean", "sst.range", "sst.min", "sst.max", "Salinity.ppt", 
#   "parmean", "log10.day.length", "cloudmean", "precipitation", "sqrt.nitrate", "log10.phosphate", "log10.chlomean", 
#   "Leaf.PercN.site", "log10.mean.fetch")], 
#   smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 8)


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

# print(ZEN.env.pca.pac)
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

# # plot cumulative proportion of variance explained by PC axes
# plot(ZEN.env.pca.pac, type = "l")
# 
# # Calculate proportion of variance explained by each PC
# summary(ZEN.env.pca.pac)
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



###################################################################################
# EXPLORE DATA COMPLETENESS                                                       #
###################################################################################

# NOTE: AIC comparisons among models are invalid unless exactly the same number of plots
# are used in each comparison, because the DF influences calculation of the AIC score. 
# This means that we need data on all plots and need to impute missing data for 
# valid AIC model comparisons. 

# # How many observations are missing for each variable?
# sum(is.na(d$log10.Zostera.AG.mass)) # 24
# sum(is.na(d$log10.Zostera.shoots.core)) # 15
# sum(is.na(d$Zostera.longest.leaf.length)) # 0
# sum(is.na(d$Leaf.PercN)) # 14
# sum(is.na(d$Temperature.C)) # 0
# sum(is.na(d$Salinity.ppt)) # 0
# sum(is.na(d$pop.density.2015)) # 20 huh?
# sum(is.na(d$GenotypicRichness)) # 0
# sum(is.na(d$AllelicRichness)) # 0 
# sum(is.na(d$grazer.richness.site)) # 0
# sum(is.na(d$log10.periphyton.mass.per.g.zostera)) # 4
# sum(is.na(d$log10.mesograzer.abund.per.g.plant)) # 9
# sum(is.na(d$log10.crustacean.abund.per.g.plant)) # 9 
# sum(is.na(d$log10.gastropod.abund.per.g.plant)) # 9

# Look at percentage of values missing for each variable
# First create function to calculate % of missing values infor each variable in a data frame… 
pMiss <- function(x){sum(is.na(x))/length(x)*100}
# # Now apply it to the data frame: 
# apply(d,2,pMiss)




#










###################################################################################
# PCA - EELGRASS VARIABLES (GLOBAL)                                               #
###################################################################################

# NOTE: The PCA for eelgrass morphology uses imputed data (see impute_missing/R) 
d.imputed <- read.csv( "data/output/ZEN_2014_imputed.csv" )

# NOTE: This includes all available ZEN eelgrass morphological variables. We use the 
# first two axes, which together explain 83% of the variation in input variables, under 
# the (arbitrary) criterion of using those PC axes necessary to capture 75% of the variation. 

## PCA - EELGRASS VARIABLES (PLOT LEVEL)                                           

# Create data frame containing the ZEN 2014 eelgrass morphological variables
zos.morph.plot.2 <- d.imputed[c("log10.Zostera.AG.mass.imputed", "log10.Zostera.BG.mass.imputed", 
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
ZEN_2014_plot <- cbind(d.imputed, zos.morph.plot.2.pca.scores) 

# Rename PCA variables 1-2 and cull PC3-4
names(ZEN_2014_plot)[names(ZEN_2014_plot)=="PC1"] <- "PC1.zos"
names(ZEN_2014_plot)[names(ZEN_2014_plot)=="PC2"] <- "PC2.zos"
ZEN_2014_plot <- subset(ZEN_2014_plot, select = -c(PC3,PC4,PC5,PC6))




# NOTE: IS THIS WHERE THIS SHOULD BE?
# Obtain mean values per site: Eelgrass growth form PCz1 and PCz2 
add_means <- ddply(ZEN_2014_plot, c("Site"), summarize, 
                   PC1.zos.site = mean(PC1.zos, na.rm = T),
                   PC2.zos.site = mean(PC2.zos, na.rm = T)
)

# Add to site means data frame
ZEN_2014_site_means <- merge(ZEN_2014_site_means, add_means)

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
range01 <- function(x, ...){(x - min(x, na.rm = T, ...)) / (max(x, na.rm = T, ...) - min(x, na.rm = T, ...))}


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
ZEN_2014_site_means_Atlantic$zPC1.env.atl <- scale(ZEN_2014_site_means_Atlantic$PC1.env.atl)
ZEN_2014_site_means_Atlantic$zPC2.env.atl <- scale(ZEN_2014_site_means_Atlantic$PC2.env.atl)
ZEN_2014_site_means_Atlantic$zPC3.env.atl <- scale(ZEN_2014_site_means_Atlantic$PC3.env.atl)
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
ZEN_2014_site_means_Atlantic$rPC1.env.atl <- range01(ZEN_2014_site_means_Atlantic$PC1.env.atl)
ZEN_2014_site_means_Atlantic$rPC2.env.atl <- range01(ZEN_2014_site_means_Atlantic$PC2.env.atl)
ZEN_2014_site_means_Atlantic$rPC3.env.atl <- range01(ZEN_2014_site_means_Atlantic$PC3.env.atl)
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
ZEN_2014_site_means_Pacific$rperiphyton.area.pac <- range01(ZEN_2014_site_means_Pacific$log10.periphyton.mass.per.area.site)
ZEN_2014_site_means_Pacific$rperiphyton.perg.pac <- range01(ZEN_2014_site_means_Pacific$log10.periphyton.mass.per.g.zostera.site)
ZEN_2014_site_means_Pacific$rmesograzer.mass.area.pac <- range01(ZEN_2014_site_means_Pacific$log10.mesograzer.mass.per.area.site)
ZEN_2014_site_means_Pacific$rmesograzer.mass.perg.pac <- range01(ZEN_2014_site_means_Pacific$log10.mesograzer.mass.per.g.plant.site)




###################################################################################
# SUBSET DATA SETS BY GEOGRAPHY                                                   #
###################################################################################

# Create reduced data sets
# # Create separate data set excluding SW.A (no periphyton data)
ZEN_2014_site_means_49 <- droplevels(subset(ZEN_2014_site_means, Site != "SW.A"))




###################################################################################
# OUTPUT CURATED DATA SETS                                                        #
###################################################################################


# Export SITE-level data set
write.csv(ZEN_2014_site_means, "data/output/ZEN_2014_site_means.csv", row.names = F)

write.csv(ZEN_2014_site_means_Atlantic, "data/output/ZEN_2014_site_means_Atlantic.csv", row.names = F)
write.csv(ZEN_2014_site_means_49_Atlantic, "data/output/ZEN_2014_site_means_49_Atlantic.csv", row.names = F)
write.csv(ZEN_2014_site_means_Pacific, "data/output/ZEN_2014_site_means_Pacific.csv", row.names = F)






