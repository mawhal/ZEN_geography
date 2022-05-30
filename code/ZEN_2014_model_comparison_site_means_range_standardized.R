###################################################################################
#                                                                                ##
# ZEN 2014: Global eelgrass ecosystem structure: compare site-level models       ##
# Data are current as of 2017-04-24                                              ##
# Emmett Duffy (duffye@si.edu)                                                   ##  
# Last updated 2022-05-29                                                        ##
#                                                                                ##
###################################################################################

# NOTE: INTEGRATE THIS FILE WITH "DEFINITIVE" SCRIPT. RENAME THAT ONE TOO.
# PUT ALL SCRIPTS INTO A LOGICAL NAMING CONVENTION MAKING CLEAR THEIR TEMPORAL SEQUENCE

###################################################################################
# TABLE OF CONTENTS                                                               #
#                                                                                 #
# METADATA AND APPROACH                                                           #
# LOAD PACKAGES                                                                   #
# READ IN AND PREPARE DATA                                                        #
# EXPLORE DATA DISTRIBUTIONS: GLOBAL                                              #
# TEST FOR REDUNDANCY OF OCEAN AND FC1                                            #
# GLM: EELGRASS GROWTH FORM PC1 (GLOBAL) - MODEL USING SITE MEANS                 #
# GLM: EELGRASS GROWTH FORM PC2 (GLOBAL) - MODEL USING SITE MEANS                 #
# GLM: PERIPHYTON PER G EELGRASS (GLOBAL) - MODEL USING SITE MEANS                #
# GLM: PERIPHYTON PER BOTTOM AREA (GLOBAL) - MODEL USING SITE MEANS               #
# GLM: MESOGRAZERS PER G EELGRASS (GLOBAL) - MODEL USING SITE MEANS               $
# GLM: MESOGRAZERS PER BOTTOM AREA (GLOBAL) - MODEL USING SITE MEANS              #
# GLM: EELGRASS GROWTH FORM PC1 (ATLANTIC) - MODEL USING SITE MEANS               #
# GLM: EELGRASS GROWTH FORM PC2 (ATLANTIC) - MODEL USING SITE MEANS               #
# GLM: PERIPHYTON PER G EELGRASS (ATLANTIC) - MODEL USING SITE MEANS              #
# GLM: PERIPHYTON PER BOTTOM AREA (ATLANTIC) -  MODEL USING SITE MEANS            #
# GLM: MESOGRAZER MASS PER G EELGRASS (ATLANTIC) - MODEL USING SITE MEANS         #
# GLM: MESOGRAZER MASS PER BOTTOM AREA (ATLANTIC) - MODEL USING SITE MEANS        #
# GLM: EELGRASS GROWTH FORM PC1 (PACIFIC) - MODEL USING SITE MEANS                #
# GLM: EELGRASS GROWTH FORM PC2 (PACIFIC)  - MODEL USING SITE MEANS               #
# GLM: PERIPHYTON PER G EELGRASS (PACIFIC) - MODEL USING SITE MEANS               #
# GLM: PERIPHYTON PER BOTTOM AREA (PACIFIC) -  MODEL USING SITE MEANS             #
# GLM: MESOGRAZER MASS PER G EELGRASS (PACIFIC) - MODEL USING SITE MEANS          #
# GLM: MESOGRAZER MASS PER G BOTTOM AREA (PACIFIC) - MODEL USING SITE MEANS       #
# GLM: FIT BEST MODELS USING LOCALLY MEASURED PREDICTORS                          #
#                                                                                 #
# CALCULATION OF EELGRASS GROWTH (LEAF EXTENSION) FROM 2011 DATA                  #
# GLM: EELGRASS PRODUCTIVITY (GLOBAL) - MODEL USING SITE MEANS                    #
# GLM: EELGRASS PRODUCTIVITY (ATLANTIC) - MODEL USING SITE MEANS                  #
# GLM: EELGRASS PRODUCTIVITY (PACIFIC) - MODEL USING SITE MEANS                   #
# RANDOM FOREST ANALYSIS                                                          #
# MISCELLANEOUS STIATISCAL CHECKS                                                 #
# MISCELLANEOUS FIGURES                                                           #
#                                                                                 #
###################################################################################

###################################################################################
# METADATA AND APPROACH                                                           #
###################################################################################

# This script explores the data and fits models for ZEN 2014 global eelgrass survey. See also:
#   ZEN_2014_data_assembly_20210217.R series: for preparation of the data
#   ZEN_2014_figures_20210227.R series: for building figures for the MS


# NOTE: This script fits linear models to individual response variables (PCz1, PCz2, 
# periphyton, mesograzers) using SITE-level data, after a long period of attempting to 
# use plot-level data and getting very complicated results that are difficult to interpret. 
# Since all exogenous predictors are site-level variables (environmental PC axes, genetic
# FCA axes), conducting the analysis with site means seems the logical approach. 

# Beginning 5 April 2021 I begin using AICc for small sample sizes. Note that, when a null model 
# (estimating only mean of response variable) is included in candidate set is invariably chosen 
# as the best by AICc. Yet often the best NON-null model is highly significant and explains a 
# large proportion of the variance. Therefore AICc seems to be overly conservative, therefore 
# I do not include null models in these model comparisons. 

# Beginning 5 April 2021 I deleted all sections involving plot-level models. They can be found 
# in previous versions of this script. 

# Beginning 14 March 2021 I use range-standardized, centered variables as opposed to z scores
# (centered and standardized to +/- 1 standard deviation)

# This script uses ZEN 2014 data and linear models to find the best 
# explanation of the global variation in eelgrass ecosystem components as a function
# of interactions among environment, genetics, eelgrass growth form, periphyton, 
# and mesograzers. 


###################################################################################
# LOAD PACKAGES                                                                   #
###################################################################################

# Load packages:
library(psych) # for pairs panels
library(ggplot2)
library(piecewiseSEM) # for SEM fitting, and for partialResid
library(nlme) # needed to run mixed models with lme
library(MuMIn) # for model averaging and AICc


###################################################################################
# READ IN AND PREPARE DATA                                                        #
###################################################################################

# Read in zen2014 SITE-level data sets 
ZEN_2014_site_means <- read.csv("ZEN_2014_site_means_20220529.csv",  header = TRUE)
ZEN_2014_site_means_49 <- droplevels(subset(ZEN_2014_site_means, Site != "SW.A"))
# ZEN_2014_site_means_Atlantic <- read.csv("ZEN_2014_site_means_Atlantic_20210227.csv",  header = TRUE)
ZEN_2014_site_means_Pacific <- read.csv("ZEN_2014_site_means_Pacific_20210314.csv",  header = TRUE)
ZEN_2014_site_means_49_Atlantic <- read.csv("ZEN_2014_site_means_49_Atlantic_20210314.csv",  header = TRUE)

# Recode Ocean to 0/1 so SEM can handle it. SITE level
ZEN_2014_site_means_49$ocean.code <- ZEN_2014_site_means_49$Ocean
ZEN_2014_site_means_49[ZEN_2014_site_means_49$ocean.code == "Atlantic", "ocean.code"] = 0
ZEN_2014_site_means_49[ZEN_2014_site_means_49$ocean.code == "Pacific", "ocean.code"] = 1
ZEN_2014_site_means_49$ocean.code <- as.numeric(ZEN_2014_site_means_49$ocean.code)

ZEN_2014_site_means_49_Atlantic$ocean.code <- ZEN_2014_site_means_49_Atlantic$Ocean
ZEN_2014_site_means_49_Atlantic[ZEN_2014_site_means_49_Atlantic$ocean.code == "Atlantic", "ocean.code"] = 0
ZEN_2014_site_means_49_Atlantic[ZEN_2014_site_means_49_Atlantic$ocean.code == "Pacific", "ocean.code"] = 1
ZEN_2014_site_means_49_Atlantic$ocean.code <- as.numeric(ZEN_2014_site_means_49_Atlantic$ocean.code)

ZEN_2014_site_means_Pacific$ocean.code <- ZEN_2014_site_means_Pacific$Ocean
ZEN_2014_site_means_Pacific[ZEN_2014_site_means_Pacific$ocean.code == "Atlantic", "ocean.code"] = 0
ZEN_2014_site_means_Pacific[ZEN_2014_site_means_Pacific$ocean.code == "Pacific", "ocean.code"] = 1
ZEN_2014_site_means_Pacific$ocean.code <- as.numeric(ZEN_2014_site_means_Pacific$ocean.code)

# Read in data For estimating leaf growth rate from Ruesink et al. (2018. Oikos)
zmgrowth <- read.csv("ZEN_2011_ZRG_AllSites_Edit141102.csv",  header = TRUE)

# NOTE: Following is NOT the file with imputed values and it has all 50 sites. Figure this out ...
# # Read in zen2014 PLOT-level data sets for modeling, with missing data imputed:
# ZEN_2014_plot <- read.csv("ZEN_2014_plot_20210430.csv", header = TRUE)
# names(ZEN_2014_plot)


###################################################################################
# EXPLORE DATA DISTRIBUTIONS: GLOBAL                                              #
###################################################################################

ZEN_2014_site_means$Ocean <- as.factor(ZEN_2014_site_means$Ocean) # necessary to code pairs panels symbols by ocean ...

# Create new data frame with variable names friendlier to plotting
ZEN_2014_site_means_renamed <- ZEN_2014_site_means
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="sst.mean"] <- "SST mean"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="sst.range"] <- "SST range"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="Salinity.ppt"] <- "Salinity"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="parmean"] <- "PAR"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="day.length"] <- "day length"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="sqrt.nitrate"] <- "NO3 (sqrt)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.phosphate"] <- "PO4 (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.chlomean"] <- "chl (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Leaf.PercN.site"] <- "Leaf %N (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="PC1.env.global"] <- "Env PCe1"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="PC2.env.global"] <- "Env PCe2"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="PC3.env.global"] <- "Env PCe3"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="PC1.zos.site"] <- "Eelgrass form PCz1"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="PC2.zos.site"] <- "Eelgrass form PCz2"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="FC1"] <- "Genetics FCA1"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="FC2"] <- "Genetics FCA2"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Zostera.AG.mass.site"] <- "AG mass (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Zostera.BG.mass.site"] <- "BG mass (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Zostera.shoots.core.site"] <- "Shoot density"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Zostera.sheath.length.site"] <- "Sheath L (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Zostera.sheath.width.site"] <- "Sheath W (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.Zostera.longest.leaf.length.cm.site"] <- "Canopy ht (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.periphyton.mass.per.g.zostera.site"] <- "Periphyton mass (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.mesograzer.mass.per.g.plant.site"] <- "Mesograzer mass (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.periphyton.mass.per.area.site"] <- "Periphyton mass/area (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.mesograzer.mass.per.area.site"] <- "Mesograzer mass/area (log)"
names(ZEN_2014_site_means_renamed)[names(ZEN_2014_site_means_renamed)=="log10.mesograzer.abund.per.area.site"] <- "Mesograzer abund/area (log)"

# Correlates of environmental/oceanographic PC axes 
pairs.panels(ZEN_2014_site_means_renamed[,c("Ocean", "Latitude", "Env PCe1", "Env PCe2", 
  "Env PCe3", "SST mean", "SST range", "Salinity", "PAR", "day length", "NO3 (sqrt)", 
  "PO4 (log)", "Leaf %N (log)", "chl (log)")], 
  hist.col="gray", pch = 21, 
  smooth = T, ci = F, density = F, ellipses = F, lm = F, digits = 2, scale = F, cex = 12, 
  bg = c("blue","green")[ZEN_2014_site_means_renamed$Ocean])
# PCe1: latitude/climate: high = warmer, brighter, less cloudy (lower latitude)
# PCe2: nutrient status: high = high PO4, leaf N  
# PCe3: estuarine: low salinity, variable temp, high chl

# Correlations between environment and biology
pairs.panels(ZEN_2014_site_means_renamed[,c("Ocean", "Latitude", "Env PCe1", "Env PCe2", 
  "Env PCe3", "Genetics FCA1", "Genetics FCA2", "Eelgrass form PCz1", "Eelgrass form PCz2",
  "Periphyton mass (log)", "Mesograzer mass (log)", "Periphyton mass/area (log)", 
  "Mesograzer mass/area (log)")], hist.col="gray", pch = 21, 
  smooth = T, ci = F, density = F, ellipses = F, lm = F, digits = 2, scale = F, cex = 8, 
  bg = c("blue","green")[ZEN_2014_site_means_renamed$Ocean])

# Missing data? No.
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


###################################################################################
# TEST FOR REDUNDANCY OF OCEAN AND FC1                                            #
###################################################################################

# Genetic FC1 is bimodal and entirely non-overlapping between Atlantic and Pcific sites.
# Test for whether we need both or only one of these predictors. 

# Compare models with ocean vs FC1, # Main effects only

# EELGRASS FORM (PCz1)

# Both Ocean and FC1
pcz1.site.g.1 <- lm(rPC1.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2  
                    , data = ZEN_2014_site_means_49)
summary(pcz1.site.g.1)
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.80426    0.23255   3.458  0.00126 ** 
# OceanPacific    -0.42366    0.17728  -2.390  0.02142 *  
# rPC1.env.global -0.26854    0.08267  -3.248  0.00229 ** 
# rPC2.env.global -0.09434    0.10878  -0.867  0.39071    
# rPC3.env.global -0.07191    0.09673  -0.743  0.46137    
# rFC1            -0.17958    0.25131  -0.715  0.47882    
# rFC2             0.53142    0.10771   4.934 1.32e-05 ***
# Residual standard error: 0.1394 on 42 degrees of freedom
# Multiple R-squared:  0.7028,	Adjusted R-squared:  0.6604 
# F-statistic: 16.55 on 6 and 42 DF,  p-value: 1.114e-09 

# Omit FC1
pcz1.site.g.o <- lm(rPC1.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global # + rFC1 
                    + rFC2  
                    , data = ZEN_2014_site_means_49)
summary(pcz1.site.g.o)
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.64794    0.07845   8.259 2.04e-10 ***
# OceanPacific    -0.30220    0.05006  -6.037 3.24e-07 ***
# rPC1.env.global -0.25613    0.08036  -3.187  0.00268 ** 
# rPC2.env.global -0.07966    0.10621  -0.750  0.45732    
# rPC3.env.global -0.08965    0.09296  -0.964  0.34022    
# rFC2             0.53309    0.10707   4.979 1.08e-05 ***
# Residual standard error: 0.1386 on 43 degrees of freedom
# Multiple R-squared:  0.6992,	Adjusted R-squared:  0.6642 
# F-statistic: 19.99 on 5 and 43 DF,  p-value: 2.978e-10

# Omit Ocean
pcz1.site.g.fc1 <- lm(rPC1.zos.site ~ # Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2  
                    , data = ZEN_2014_site_means_49)
summary(pcz1.site.g.fc1)
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.31172    0.11346   2.747  0.00874 ** 
# rPC1.env.global -0.24338    0.08637  -2.818  0.00728 ** 
# rPC2.env.global -0.08894    0.11456  -0.776  0.44176    
# rPC3.env.global -0.12396    0.09928  -1.249  0.21854    
# rFC1             0.39626    0.07518   5.271 4.15e-06 ***
# rFC2             0.55302    0.11305   4.892 1.44e-05 ***
# Residual standard error: 0.1469 on 43 degrees of freedom
# Multiple R-squared:  0.6624,	Adjusted R-squared:  0.6232 
# F-statistic: 16.87 on 5 and 43 DF,  p-value: 3.301e-09


AICc(pcz1.site.g.1, pcz1.site.g.o, pcz1.site.g.fc1)
#                 df      AICc
# pcz1.site.g.1    8 -41.97843
# pcz1.site.g.o    7 -44.25458
# pcz1.site.g.fc1  7 -38.59952
# Best model excludes FC1. Model with both is better than Ocean only. 
# All models have strong effects of PCe1 and FC2.
# When ocean is excluded, FC1 is highly significant. 


# EELGRASS BIOMASS (PCz2)

# Both Ocean and FC1
pcz2.site.g.1 <- lm(rPC2.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2  
                    , data = ZEN_2014_site_means_49)
summary(pcz2.site.g.1)
#                 Estimate Std. Error t value Pr(>|t|)  
# (Intercept)      0.68250    0.29521   2.312   0.0258 *
# OceanPacific     0.03609    0.22505   0.160   0.8734  
# rPC1.env.global -0.04959    0.10495  -0.472   0.6390  
# rPC2.env.global -0.02545    0.13809  -0.184   0.8546  
# rPC3.env.global -0.31280    0.12279  -2.547   0.0146 *
# rFC1             0.06084    0.31902   0.191   0.8497  
# rFC2            -0.31134    0.13672  -2.277   0.0279 *
# Residual standard error: 0.177 on 42 degrees of freedom
# Multiple R-squared:   0.28,	Adjusted R-squared:  0.1771 
# F-statistic: 2.722 on 6 and 42 DF,  p-value: 0.02525


# Omit Ocean
pcz2.site.g.o <- lm(rPC2.zos.site ~ # Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2  
                    , data = ZEN_2014_site_means_49)
summary(pcz2.site.g.o)
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.72445    0.13517   5.359  3.1e-06 ***
# rPC1.env.global -0.05173    0.10290  -0.503   0.6178    
# rPC2.env.global -0.02591    0.13648  -0.190   0.8503    
# rPC3.env.global -0.30837    0.11828  -2.607   0.0125 *  
# rFC1             0.01179    0.08957   0.132   0.8959    
# rFC2            -0.31318    0.13469  -2.325   0.0249 *  
# Residual standard error: 0.175 on 43 degrees of freedom
# Multiple R-squared:  0.2795,	Adjusted R-squared:  0.1957 
# F-statistic: 3.336 on 5 and 43 DF,  p-value: 0.01237

# Omit FC1
pcz2.site.g.fc1 <- lm(rPC2.zos.site ~ Ocean
                        + rPC1.env.global + rPC2.env.global + rPC3.env.global # + rFC1 
                      + rFC2  
                      , data = ZEN_2014_site_means_49)
summary(pcz2.site.g.fc1)
#                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.735459   0.099027   7.427 3.12e-09 ***
# OceanPacific    -0.005066   0.063194  -0.080   0.9365    
# rPC1.env.global -0.053791   0.101445  -0.530   0.5987    
# rPC2.env.global -0.030428   0.134073  -0.227   0.8215    
# rPC3.env.global -0.306791   0.117344  -2.614   0.0123 *  
# rFC2            -0.311903   0.135152  -2.308   0.0259 *  
# Residual standard error: 0.175 on 43 degrees of freedom
# Multiple R-squared:  0.2793,	Adjusted R-squared:  0.1955 
# F-statistic: 3.333 on 5 and 43 DF,  p-value: 0.01243

AICc(pcz2.site.g.1, pcz2.site.g.o, pcz2.site.g.fc1)
#                 df      AICc
# pcz2.site.g.1    8 -18.59891
# pcz2.site.g.o    7 -21.43721
# pcz2.site.g.fc1  7 -21.42479
# Ocean and FC1 provide equivalent information. Either one is better than both together. 
# Otherwise all equivalent - all have significant effects of PCe3 and FC2.
# Neither ocean nor FC1 are significant in any model. 


# PERIPHYTON MASS PER UNIT AREA

# Both Ocean and FC1
peri.site.area.1 <- lm(rperiphyton ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2  
                    + rPC1.zos.site + rPC2.zos.site
                    , data = ZEN_2014_site_means_49)
summary(peri.site.area.1)
#                  Estimate Std. Error t value Pr(>|t|)  
# (Intercept)      0.826925   0.329250   2.512   0.0162 *
# OceanPacific    -0.127865   0.229689  -0.557   0.5808  
# rPC1.env.global -0.001243   0.112101  -0.011   0.9912  
# rPC2.env.global  0.052266   0.133046   0.393   0.6965  
# rPC3.env.global  0.004595   0.126177   0.036   0.9711  
# rFC1            -0.187937   0.306810  -0.613   0.5436  
# rFC2             0.388740   0.175099   2.220   0.0321 *
# rPC1.zos.site   -0.377849   0.189089  -1.998   0.0525 .
# rPC2.zos.site   -0.201705   0.148956  -1.354   0.1833  
# Residual standard error: 0.169 on 40 degrees of freedom
# Multiple R-squared:  0.4067,	Adjusted R-squared:  0.288 
# F-statistic: 3.427 on 8 and 40 DF,  p-value: 0.004309

peri.site.area.o <- lm(rperiphyton ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global # + rFC1 
                    + rFC2 + rPC1.zos.site + rPC2.zos.site
                    , data = ZEN_2014_site_means_49)
summary(peri.site.area.o)
#                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.657768   0.177926   3.697  0.00064 ***
# OceanPacific     0.003258   0.082644   0.039  0.96874    
# rPC1.env.global  0.014937   0.108111   0.138  0.89078    
# rPC2.env.global  0.068566   0.129361   0.530  0.59895    
# rPC3.env.global -0.014058   0.121511  -0.116  0.90846    
# rFC2             0.382087   0.173425   2.203  0.03326 *  
# rPC1.zos.site   -0.364529   0.186398  -1.956  0.05734 .  
# rPC2.zos.site   -0.205879   0.147662  -1.394  0.17075    
# Residual standard error: 0.1677 on 41 degrees of freedom
# Multiple R-squared:  0.4011,	Adjusted R-squared:  0.2989 
# F-statistic: 3.923 on 7 and 41 DF,  p-value: 0.002323

peri.site.area.fc1 <- lm(rperiphyton ~ # Ocean
                      + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2  
                      + rPC1.zos.site + rPC2.zos.site
                      , data = ZEN_2014_site_means_49)
summary(peri.site.area.fc1)
#                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.671190   0.172152   3.899 0.000351 ***
# rPC1.env.global  0.015055   0.107296   0.140 0.889103    
# rPC2.env.global  0.057033   0.131648   0.433 0.667122    
# rPC3.env.global -0.008423   0.122943  -0.069 0.945711    
# rFC1            -0.028763   0.110303  -0.261 0.795581    
# rFC2             0.372821   0.171288   2.177 0.035326 *  
# rPC1.zos.site   -0.340775   0.175477  -1.942 0.059030 .  
# rPC2.zos.site   -0.207880   0.147287  -1.411 0.165673    
# Residual standard error: 0.1676 on 41 degrees of freedom
# Multiple R-squared:  0.4021,	Adjusted R-squared:    0.3 
# F-statistic: 3.939 on 7 and 41 DF,  p-value: 0.002259


AICc(peri.site.area.1, peri.site.area.o, peri.site.area.fc1)
#                    df      AICc
# peri.site.area.1   10 -19.31829
# peri.site.area.o    9 -22.03487
# peri.site.area.fc1  9 -22.11421
# Ocean and FC1 provide equivalent information. Either one is better than both together. 
# Otherwise all equivalent - all have significant effects of FC2.
# Neither ocean nor FC1 are significant in any model. 


# SITE level test: Mesograzer mass per g eelgrass

# Main effects only
meso.site.area.1 <- lm(rmesograzer.mass ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2  
                    + rPC1.zos.site + rPC2.zos.site
                    , data = ZEN_2014_site_means_49)
summary(meso.site.area.1)
#                  Estimate Std. Error t value Pr(>|t|)   
# (Intercept)      0.209529   0.262366   0.799  0.42923   
# OceanPacific     0.135382   0.183030   0.740  0.46382   
# rPC1.env.global -0.018660   0.089329  -0.209  0.83560   
# rPC2.env.global  0.287859   0.106019   2.715  0.00973 **
# rPC3.env.global  0.340452   0.100545   3.386  0.00160 **
# rFC1             0.042117   0.244484   0.172  0.86409   
# rFC2             0.005721   0.139529   0.041  0.96750   
# rPC1.zos.site    0.098146   0.150677   0.651  0.51853   
# rPC2.zos.site   -0.171455   0.118697  -1.444  0.15639   
# Residual standard error: 0.1347 on 40 degrees of freedom
# Multiple R-squared:  0.5351,	Adjusted R-squared:  0.4421 
# F-statistic: 5.754 on 8 and 40 DF,  p-value: 6.855e-05

meso.site.area.o <- lm(rmesograzer.mass ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global # + rFC1 
                    + rFC2 + rPC1.zos.site + rPC2.zos.site
                    , data = ZEN_2014_site_means_49)
summary(meso.site.area.o)
#                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.247438   0.141174   1.753 0.087126 .  
# OceanPacific     0.105997   0.065573   1.616 0.113663    
# rPC1.env.global -0.022286   0.085780  -0.260 0.796318    
# rPC2.env.global  0.284207   0.102640   2.769 0.008407 ** 
# rPC3.env.global  0.344632   0.096412   3.575 0.000916 ***
# rFC2             0.007212   0.137603   0.052 0.958456    
# rPC1.zos.site    0.095161   0.147896   0.643 0.523525    
# rPC2.zos.site   -0.170519   0.117161  -1.455 0.153167    
# Residual standard error: 0.1331 on 41 degrees of freedom
# Multiple R-squared:  0.5347,	Adjusted R-squared:  0.4553 
# F-statistic: 6.731 on 7 and 41 DF,  p-value: 2.506e-05

meso.site.area.fc1 <- lm(rmesograzer.mass ~ # Ocean
                        + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2  
                      + rPC1.zos.site + rPC2.zos.site
                      , data = ZEN_2014_site_means_49)
summary(meso.site.area.fc1)
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.37442    0.13758   2.721 0.009498 ** 
# rPC1.env.global -0.03592    0.08575  -0.419 0.677523    
# rPC2.env.global  0.28281    0.10521   2.688 0.010341 *  
# rPC3.env.global  0.35423    0.09826   3.605 0.000837 ***
# rFC1            -0.12641    0.08815  -1.434 0.159151    
# rFC2             0.02258    0.13689   0.165 0.869826    
# rPC1.zos.site    0.05889    0.14024   0.420 0.676722    
# rPC2.zos.site   -0.16492    0.11771  -1.401 0.168732    
# Residual standard error: 0.1339 on 41 degrees of freedom
# Multiple R-squared:  0.5287,	Adjusted R-squared:  0.4482 
# F-statistic: 6.571 on 7 and 41 DF,  p-value: 3.177e-05

AICc(meso.site.area.1, meso.site.area.o, meso.site.area.fc1)
#                    df      AICc
# meso.site.area.1   10 -41.57172
# meso.site.area.o    9 -44.70947
# meso.site.area.fc1  9 -44.08014
# Ocean and FC1 provide equivalent information. Either one is better than both together. 
# Otherwise all equivalent - all have significant effects of PCe2, PCe3.
# Neither ocean nor FC1 are significant in any model. 


###################################################################################
# GLM: EELGRASS GROWTH FORM PC1 (GLOBAL) - MODEL USING SITE MEANS                 #
###################################################################################

# Using range-standardized variables.

# NOTE: This is an alternate approach to modeling eelgrass growth form axes. Because all 
# environmental and genetic predictors are site-level variables, it seems equally 
# appropriate to model site means of eelgrass growth form with simple linear models here. Also 
# there are few sites and many predictors so we build up from the no-interactions model

# Main effects only
pcz1.site.g.1 <- lm(rPC1.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2  
                    , data = ZEN_2014_site_means_49)
summary(pcz1.site.g.1)
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.80426    0.23255   3.458  0.00126 ** 
# OceanPacific    -0.42366    0.17728  -2.390  0.02142 *  
# rPC1.env.global -0.26854    0.08267  -3.248  0.00229 ** 
# rPC2.env.global -0.09434    0.10878  -0.867  0.39071    
# rPC3.env.global -0.07191    0.09673  -0.743  0.46137    
# rFC1            -0.17958    0.25131  -0.715  0.47882    
# rFC2             0.53142    0.10771   4.934 1.32e-05 ***
# Residual standard error: 0.1394 on 42 degrees of freedom
# Multiple R-squared:  0.7028,	Adjusted R-squared:  0.6604 
# F-statistic: 16.55 on 6 and 42 DF,  p-value: 1.114e-09

# Add interaction: ocean x PCe1
pcz1.site.g.2 <- lm(rPC1.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + Ocean*rPC1.env.global # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe2
pcz1.site.g.3 <- lm(rPC1.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + Ocean*rPC2.env.global # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe3
pcz1.site.g.4 <- lm(rPC1.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + Ocean*rPC3.env.global # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC1
pcz1.site.g.5 <- lm(rPC1.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + Ocean*rFC1 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC2
pcz1.site.g.6 <- lm(rPC1.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + Ocean*rFC2 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe1 x FC2
pcz1.site.g.7 <- lm(rPC1.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + rPC1.env.global*rFC2 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC2
pcz1.site.g.8 <- lm(rPC1.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + rPC2.env.global*rFC2 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC2
pcz1.site.g.9 <- lm(rPC1.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + rPC3.env.global*rFC2 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe1 x FC1
pcz1.site.g.10 <- lm(rPC1.zos.site ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                     + rPC1.env.global*rFC1 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC1
pcz1.site.g.11 <- lm(rPC1.zos.site ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                     + rPC2.env.global*rFC1 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC1
pcz1.site.g.12 <- lm(rPC1.zos.site ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                     + rPC3.env.global*rFC1 # added
                     , data = ZEN_2014_site_means_49)


AICc(pcz1.site.g.1, pcz1.site.g.2, pcz1.site.g.3, pcz1.site.g.4, pcz1.site.g.5, pcz1.site.g.6, pcz1.site.g.7, pcz1.site.g.8, pcz1.site.g.9, pcz1.site.g.10, pcz1.site.g.11, pcz1.site.g.12) 
#                df      AICc
# pcz1.site.g.1   8 -41.97843
# pcz1.site.g.2   9 -39.05397
# pcz1.site.g.3   9 -39.45692
# pcz1.site.g.4   9 -47.29079
# pcz1.site.g.5   9 -38.96368
# pcz1.site.g.6   9 -42.27455
# pcz1.site.g.7   9 -39.10851
# pcz1.site.g.8   9 -39.03086
# pcz1.site.g.9   9 -39.41174
# pcz1.site.g.10  9 -38.96803
# pcz1.site.g.11  9 -39.68549
# pcz1.site.g.12  9 -47.18834
# RESULTS: Best models are 4, 12. 

# Combine models 4 + 12
pcz1.site.g.13 <- lm(rPC1.zos.site ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                     + Ocean*rPC3.env.global # added
                     + rPC3.env.global*rFC1 # added
                     , data = ZEN_2014_site_means_49)

AICc(pcz1.site.g.4, pcz1.site.g.12, pcz1.site.g.13)
# No good. Average models 4 and 12

# Obtain model-averaged, standardized coefficients from models < 2 AIC from best 
summary(model.avg(pcz1.site.g.4, pcz1.site.g.12))
# Result: averaging simply adds two redundant interactions to model. They are both
# "significant" when conditionally averaged, and neither is significant in full
# average. Duh. Not helpful.

# Fit best model (4) without FC1 
pcz1.site.g.14 <- lm(rPC1.zos.site ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global # + rFC1 
                     + rFC2
                     + Ocean*rPC3.env.global # added
                     , data = ZEN_2014_site_means_49)


# Fit best model (12) without ocean 
pcz1.site.g.15 <- lm(rPC1.zos.site ~ # Ocean
                       + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC3.env.global*rFC1 # added
                     , data = ZEN_2014_site_means_49)

AICc(pcz1.site.g.4, pcz1.site.g.12, pcz1.site.g.14, pcz1.site.g.15)
#                df      AICc
# pcz1.site.g.4   9 -47.29079
# pcz1.site.g.12  9 -47.18834
# pcz1.site.g.14  8 -48.72557
# pcz1.site.g.15  8 -41.00106
# Excluding ocean is much worse than other models.

summary(pcz1.site.g.4)
#                              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   0.94962    0.22253   4.267 0.000114 ***
# OceanPacific                 -0.84453    0.22469  -3.759 0.000533 ***
# rPC1.env.global              -0.31300    0.07853  -3.986 0.000270 ***
# rPC2.env.global              -0.09104    0.10114  -0.900 0.373299    
# rPC3.env.global              -0.16235    0.09573  -1.696 0.097484 .  
# rFC2                          0.52789    0.10014   5.272 4.66e-06 ***
# rFC1                         -0.27374    0.23612  -1.159 0.253018    
# OceanPacific:rPC3.env.global  0.71071    0.25788   2.756 0.008693 ** 
# Residual standard error: 0.1296 on 41 degrees of freedom
# Multiple R-squared:  0.7493,	Adjusted R-squared:  0.7065 
# F-statistic:  17.5 on 7 and 41 DF,  p-value: 1.713e-10

summary(pcz1.site.g.12)
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           0.55045    0.23544   2.338 0.024350 *  
# OceanPacific         -0.48452    0.16648  -2.910 0.005810 ** 
# rPC1.env.global      -0.31294    0.07863  -3.980 0.000275 ***
# rPC2.env.global      -0.09154    0.10124  -0.904 0.371168    
# rPC3.env.global       0.63029    0.27185   2.319 0.025486 *  
# rFC2                  0.52752    0.10025   5.262 4.81e-06 ***
# rFC1                  0.18611    0.26934   0.691 0.493461    
# rPC3.env.global:rFC1 -0.90279    0.32979  -2.737 0.009115 ** 
# Residual standard error: 0.1298 on 41 degrees of freedom
# Multiple R-squared:  0.7487,	Adjusted R-squared:  0.7058 
# F-statistic: 17.45 on 7 and 41 DF,  p-value: 1.785e-10


# Examine residuals of best model
op <- par(mfrow = c(2,4))
ypred = predict(pcz1.site.g.4)
res = residuals(pcz1.site.g.4, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49$rPC1.env.global,res, xlab = "PCe1 (scaled)", ylab = "residuals",) # diagonal row - zeros?
plot(ZEN_2014_site_means_49$rPC2.env.global,res, xlab = "PCe2 (scaled)", ylab = "residuals",) # diagonal row - zeros?
plot(ZEN_2014_site_means_49$rPC3.env.global,res, xlab = "PCe3 (scaled)", ylab = "residuals",) # diagonal row - zeros?
plot(ZEN_2014_site_means_49$rFC2,res, xlab = "PCz1 (scaled)", ylab = "residuals",) # diagonal row - zeros?
par(op)
# RESULTS: Left-skewed BUT only one point well above the qq line 

# RESULTS and INTERPRETATION: Globally, eelgrass growth form (PCz1) is affected by ocean, 
# latitude/climate, and genetic FC2. Forest form is best developed in Pacific, at lower latitudes 
# (high PCe1, check in figure), and at high values of FC2. Details: Ocean and FC1 are 
# confounded. Best two models are equivalent and differ only in whether PCe3 interacts with 
# ocean vs FC1. Interaction of ocean and PCe3 is highly significant but doesn't affect value much 
# if I'm estimating right. 

# ACTION: Strong effect of ocean and significant interaction of PCe3 with ocean suggest
# we should proceed with separate analyses by ocean. Need partial plots of PCz1 vs main 
# predictors at site level. 


###################################################################################
# GLM: EELGRASS GROWTH FORM PC2 (GLOBAL) - MODEL USING SITE MEANS                 #
###################################################################################

# Main effects only
pcz2.site.g.1 <- lm(rPC2.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2  
                    , data = ZEN_2014_site_means_49)
summary(pcz2.site.g.1)
#                 Estimate Std. Error t value Pr(>|t|)  
# (Intercept)      0.68250    0.29521   2.312   0.0258 *
# OceanPacific     0.03609    0.22505   0.160   0.8734  
# rPC1.env.global -0.04959    0.10495  -0.472   0.6390  
# rPC2.env.global -0.02545    0.13809  -0.184   0.8546  
# rPC3.env.global -0.31280    0.12279  -2.547   0.0146 *
# rFC1             0.06084    0.31902   0.191   0.8497  
# rFC2            -0.31134    0.13672  -2.277   0.0279 *
# Residual standard error: 0.177 on 42 degrees of freedom
# Multiple R-squared:   0.28,	Adjusted R-squared:  0.1771 
# F-statistic: 2.722 on 6 and 42 DF,  p-value: 0.02525

# Add interaction: ocean x PCe1
pcz2.site.g.2 <- lm(rPC2.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + Ocean*rPC1.env.global # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe2
pcz2.site.g.3 <- lm(rPC2.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + Ocean*rPC2.env.global # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe3
pcz2.site.g.4 <- lm(rPC2.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + Ocean*rPC3.env.global # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC1
pcz2.site.g.5 <- lm(rPC2.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + Ocean*rFC1 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC2
pcz2.site.g.6 <- lm(rPC2.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + Ocean*rFC2 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe1 x FC2
pcz2.site.g.7 <- lm(rPC2.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + rPC1.env.global*rFC2 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC2
pcz2.site.g.8 <- lm(rPC2.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + rPC2.env.global*rFC2 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC2
pcz2.site.g.9 <- lm(rPC2.zos.site ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                    + rPC3.env.global*rFC2 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe1 x FC1
pcz2.site.g.10 <- lm(rPC2.zos.site ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                     + rPC1.env.global*rFC1 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC1
pcz2.site.g.11 <- lm(rPC2.zos.site ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                     + rPC2.env.global*rFC1 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC1
pcz2.site.g.12 <- lm(rPC2.zos.site ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                     + rPC3.env.global*rFC1 # added
                     , data = ZEN_2014_site_means_49)


AICc(pcz2.site.g.1, pcz2.site.g.2, pcz2.site.g.3, pcz2.site.g.4, pcz2.site.g.5, pcz2.site.g.6, pcz2.site.g.7, pcz2.site.g.8, pcz2.site.g.9, pcz2.site.g.10, pcz2.site.g.11, pcz2.site.g.12) 
#               df      AIC
# pcz2.site.g.1   8 -18.59891
# pcz2.site.g.2   9 -16.11037
# pcz2.site.g.3   9 -18.27818
# pcz2.site.g.4   9 -22.53802
# pcz2.site.g.5   9 -15.71361
# pcz2.site.g.6   9 -18.10942
# pcz2.site.g.7   9 -17.87445
# pcz2.site.g.8   9 -15.78955
# pcz2.site.g.9   9 -20.31331
# pcz2.site.g.10  9 -16.82209
# pcz2.site.g.11  9 -18.48339
# pcz2.site.g.12  9 -22.50389
# RESULTS: Best models are 4 and 12.These differ only in whether PCe3 interacts with ocean or FC1.

summary(pcz2.site.g.4)
#                              Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                   0.85227    0.28648   2.975  0.00489 **
# OceanPacific                 -0.45548    0.28925  -1.575  0.12301   
# rPC1.env.global              -0.10151    0.10110  -1.004  0.32121   
# rPC2.env.global              -0.02159    0.13020  -0.166  0.86910   
# rPC3.env.global              -0.41843    0.12324  -3.395  0.00153 **
# rFC2                         -0.31546    0.12891  -2.447  0.01877 * 
# rFC1                         -0.04914    0.30396  -0.162  0.87236   
# OceanPacific:rPC3.env.global  0.83011    0.33198   2.500  0.01649 * 
# Residual standard error: 0.1669 on 41 degrees of freedom
# Multiple R-squared:  0.3752,	Adjusted R-squared:  0.2686 
# F-statistic: 3.518 on 7 and 41 DF,  p-value: 0.004777

# Fit best model (4) without FC1 
pcz2.site.g.13 <- lm(rPC2.zos.site ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 # + rFC1
                     + Ocean*rPC3.env.global 
                     , data = ZEN_2014_site_means_49)


# Fit best model (12) without ocean 
pcz2.site.g.14 <- lm(rPC2.zos.site ~ # Ocean
                       + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                     + rPC3.env.global*rFC1 
                     , data = ZEN_2014_site_means_49)

AICc(pcz2.site.g.4, pcz2.site.g.12, pcz2.site.g.13, pcz2.site.g.14)
#                df       AIC
# pcz2.site.g.4   9 -22.53802
# pcz2.site.g.12  9 -22.50389
# pcz2.site.g.13  8 -25.52218
# pcz2.site.g.14  8 -25.48693
# Models 13 and 14 are equivalent, meaning ocean and FC1 are equivalent. Either one
# is better than having both predictors. Nevertheless, we keep both in because this helps 
# comparison with other response variables.

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(pcz2.site.g.4)
res = residuals(pcz2.site.g.4, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49$rPC1.env.global,res, xlab = "PCe1 (scaled)", ylab = "residuals",) # diagonal row - zeros?
plot(ZEN_2014_site_means_49$rPC2.env.global,res, xlab = "PCe2 (scaled)", ylab = "residuals",) # diagonal row - zeros?
plot(ZEN_2014_site_means_49$rPC3.env.global,res, xlab = "PCe3 (scaled)", ylab = "residuals",) # diagonal row - zeros?
plot(ZEN_2014_site_means_49$rFC2,res, xlab = "PCz1 (scaled)", ylab = "residuals",) # diagonal row - zeros?
par(op)
# RESULTS: 

# RESULTS and INTERPRETATION: Bottom line: Highest biomass (low PCz2) 
# occurs in productive estuaries (high PCe3), and at high values of FC2. Details: Ocean and FC1 are 
# confounded. Best two models are equivalent and differ only in whether PCe3 interacts with 
# ocean vs FC1, as was true of PCz1. Interaction of ocean and PCe3 is significant but doesn't 
# affect value much if I'm estimating right. 

# ACTION: Eelgrass biomass (PCz2): Significant interactive effect of PCe3 by ocean suggests we 
# should proceed with separate analyses by ocean. The partial plots below also show that relationships 
# are different in the two oceans.  


###################################################################################
# GLM: PERIPHYTON PER G EELGRASS (GLOBAL) - MODEL USING SITE MEANS                #
###################################################################################

# Main effects only
peri.site.g.1 <- lm(rperiphyton.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site
                    , data = ZEN_2014_site_means_49)
summary(peri.site.g.1)

# Add interaction: ocean x PCe1
peri.site.g.2 <- lm(rperiphyton.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site
                    + Ocean*rPC1.env.global # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe2
peri.site.g.3 <- lm(rperiphyton.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site
                    + Ocean*rPC2.env.global # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe3
peri.site.g.4 <- lm(rperiphyton.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site
                    + Ocean*rPC3.env.global # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC1
peri.site.g.5 <- lm(rperiphyton.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site
                    + Ocean*rFC1 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC2
peri.site.g.6 <- lm(rperiphyton.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site
                    + Ocean*rFC2 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe1 x FC1
peri.site.g.7 <- lm(rperiphyton.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site
                    + rPC1.env.global*rFC1 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC1
peri.site.g.8 <- lm(rperiphyton.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site
                    + rPC2.env.global*rFC1 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC1
peri.site.g.9 <- lm(rperiphyton.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site
                    + rPC3.env.global*rFC1 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe1 x FC2
peri.site.g.10 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site
                     + rPC1.env.global*rFC2 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC2
peri.site.g.11 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site
                     + rPC2.env.global*rFC2 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC2
peri.site.g.12 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site
                     + rPC3.env.global*rFC2 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x PCe1
peri.site.g.13 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.zos.site*rPC1.env.global # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x PCe2
peri.site.g.14 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.zos.site*rPC2.env.global # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x PCe3
peri.site.g.15 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.zos.site*rPC3.env.global # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x FC1
peri.site.g.16 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.zos.site*rFC1 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x FC2
peri.site.g.17 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.zos.site*rFC2 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x PCe1
peri.site.g.18 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC2.zos.site*rPC1.env.global # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x PCe2
peri.site.g.19 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC2.zos.site*rPC2.env.global # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x PCe3
peri.site.g.20 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC2.zos.site*rPC3.env.global # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x FC1
peri.site.g.21 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC2.zos.site*rFC1 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x FC2
peri.site.g.22 <- lm(rperiphyton.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC2.zos.site*rFC2 # added
                     , data = ZEN_2014_site_means_49)


AICc(peri.site.g.1, peri.site.g.2, peri.site.g.3, peri.site.g.4, peri.site.g.5, peri.site.g.6, 
    peri.site.g.7, peri.site.g.8, peri.site.g.9, peri.site.g.10, peri.site.g.11, peri.site.g.12,
    peri.site.g.13, peri.site.g.14, peri.site.g.15, peri.site.g.16, peri.site.g.17, peri.site.g.18, 
    peri.site.g.19, peri.site.g.20, peri.site.g.21, peri.site.g.22) 
#               df      AIC
# peri.site.g.1  10  -5.626465
# peri.site.g.2  11 -11.670728
# peri.site.g.3  11  -6.704004
# peri.site.g.4  11  -3.837470
# peri.site.g.5  11  -2.361773
# peri.site.g.6  11  -2.638429
# peri.site.g.7  11 -11.028504
# peri.site.g.8  11  -3.858982
# peri.site.g.9  11  -3.216732
# peri.site.g.10 11  -4.381094
# peri.site.g.11 11  -4.331242
# peri.site.g.12 11  -2.489553
# peri.site.g.13 11  -9.304805
# peri.site.g.14 11  -2.411154
# peri.site.g.15 11  -2.403725
# peri.site.g.16 11  -2.319564
# peri.site.g.17 11  -2.873401
# peri.site.g.18 11  -2.495797
# peri.site.g.19 11  -2.764250
# peri.site.g.20 11  -6.420954
# peri.site.g.21 11  -2.602693
# peri.site.g.22 11  -2.921126
# RESULTS: Best models are 2 and 7 (same except they switch ocean for FC1) 

summary(peri.site.g.2)
#                              Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                   0.34845    0.36364   0.958  0.34385   
# OceanPacific                  0.27353    0.28342   0.965  0.34043   
# rPC1.env.global               0.15459    0.13354   1.158  0.25403   
# rPC2.env.global               0.06980    0.14121   0.494  0.62388   
# rPC3.env.global               0.05341    0.13517   0.395  0.69493   
# rFC1                         -0.11079    0.32607  -0.340  0.73584   
# rFC2                          0.76730    0.22055   3.479  0.00125 **
# rPC1.zos.site                -0.15755    0.20017  -0.787  0.43598   
# rPC2.zos.site                 0.20488    0.15839   1.294  0.20345   
# OceanPacific:rPC1.env.global -0.73127    0.25478  -2.870  0.00660 **
# Residual standard error: 0.1789 on 39 degrees of freedom
# Multiple R-squared:  0.3427,	Adjusted R-squared:  0.191 
# F-statistic:  2.26 on 9 and 39 DF,  p-value: 0.03816
# NOTE: Barely significant. Not much explanatory power for global periphyton. 

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(peri.site.g.2)
res = residuals(peri.site.g.2, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49$rPC1.env.global,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC2.env.global,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC3.env.global,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rFC1,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rFC2,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC1.zos.site,res, xlab = "PCz1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC2.zos.site,res, xlab = "PCz2 (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: Tails depart from qq line.

# RESULTS and INTERPRETATION: On a global scale, no model explains site-level periphyton 
# mass well (R2 = 0.19). Periphyton increases with genetic FC2, and there is an interaction between 
# latitude/climate (PCe1) and ocean.


###################################################################################
# GLM: PERIPHYTON PER BOTTOM AREA (GLOBAL) - MODEL USING SITE MEANS               #
###################################################################################

# Main effects only
peri.area.site.g.1 <- lm(rperiphyton ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site
                         , data = ZEN_2014_site_means_49)
summary(peri.area.site.g.1)
#                  Estimate Std. Error t value Pr(>|t|)  
# (Intercept)      0.826925   0.329250   2.512   0.0162 *
# OceanPacific    -0.127865   0.229689  -0.557   0.5808  
# rPC1.env.global -0.001243   0.112101  -0.011   0.9912  
# rPC2.env.global  0.052266   0.133046   0.393   0.6965  
# rPC3.env.global  0.004595   0.126177   0.036   0.9711  
# rFC1            -0.187937   0.306810  -0.613   0.5436  
# rFC2             0.388740   0.175099   2.220   0.0321 *
# rPC1.zos.site   -0.377849   0.189089  -1.998   0.0525 .
# rPC2.zos.site   -0.201705   0.148956  -1.354   0.1833  
# Residual standard error: 0.169 on 40 degrees of freedom
# Multiple R-squared:  0.4067,	Adjusted R-squared:  0.288 
# F-statistic: 3.427 on 8 and 40 DF,  p-value: 0.004309

# Add interaction: ocean x PCe1
peri.area.site.g.2 <- lm(rperiphyton ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site
                         + Ocean*rPC1.env.global # added
                         , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe2
peri.area.site.g.3 <- lm(rperiphyton ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site
                         + Ocean*rPC2.env.global # added
                         , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe3
peri.area.site.g.4 <- lm(rperiphyton ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site
                         + Ocean*rPC3.env.global # added
                         , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC1
peri.area.site.g.5 <- lm(rperiphyton ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site
                         + Ocean*rFC1 # added
                         , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC2
peri.area.site.g.6 <- lm(rperiphyton ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site
                         + Ocean*rFC2 # added
                         , data = ZEN_2014_site_means_49)

# Add interaction:  PCe1 x FC1
peri.area.site.g.7 <- lm(rperiphyton ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site
                         + rPC1.env.global*rFC1 # added
                         , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC1
peri.area.site.g.8 <- lm(rperiphyton ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site
                         + rPC2.env.global*rFC1 # added
                         , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC1
peri.area.site.g.9 <- lm(rperiphyton ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site
                         + rPC3.env.global*rFC1 # added
                         , data = ZEN_2014_site_means_49)

# Add interaction:  PCe1 x FC2
peri.area.site.g.10 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site
                          + rPC1.env.global*rFC2 # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC2
peri.area.site.g.11 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site
                          + rPC2.env.global*rFC2 # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC2
peri.area.site.g.12 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site
                          + rPC3.env.global*rFC2 # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x PCe1
peri.area.site.g.13 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.zos.site*rPC1.env.global # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x PCe2
peri.area.site.g.14 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.zos.site*rPC2.env.global # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x PCe3
peri.area.site.g.15 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.zos.site*rPC3.env.global # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x FC1
peri.area.site.g.16 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.zos.site*rFC1 # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x FC2
peri.area.site.g.17 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.zos.site*rFC2 # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x PCe1
peri.area.site.g.18 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC2.zos.site*rPC1.env.global # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x PCe2
peri.area.site.g.19 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC2.zos.site*rPC2.env.global # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x PCe3
peri.area.site.g.20 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC2.zos.site*rPC3.env.global # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x FC1
peri.area.site.g.21 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC2.zos.site*rFC1 # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x FC2
peri.area.site.g.22 <- lm(rperiphyton ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC2.zos.site*rFC2 # added
                          , data = ZEN_2014_site_means_49)


AICc(peri.area.site.g.1, peri.area.site.g.2, peri.area.site.g.3, peri.area.site.g.4, peri.area.site.g.5, peri.area.site.g.6, 
    peri.area.site.g.7, peri.area.site.g.8, peri.area.site.g.9, peri.area.site.g.10, peri.area.site.g.11, peri.area.site.g.12,
    peri.area.site.g.13, peri.area.site.g.14, peri.area.site.g.15, peri.area.site.g.16, peri.area.site.g.17, peri.area.site.g.18, 
    peri.area.site.g.19, peri.area.site.g.20, peri.area.site.g.21, peri.area.site.g.22) 
#               df      AIC
# peri.area.site.g.1  10 -19.31829
# peri.area.site.g.2  11 -23.14907
# peri.area.site.g.3  11 -21.24154
# peri.area.site.g.4  11 -18.49845
# peri.area.site.g.5  11 -15.97459
# peri.area.site.g.6  11 -16.06597
# peri.area.site.g.7  11 -22.85434
# peri.area.site.g.8  11 -18.41580
# peri.area.site.g.9  11 -17.70419
# peri.area.site.g.10 11 -16.76792
# peri.area.site.g.11 11 -18.62932
# peri.area.site.g.12 11 -16.07364
# peri.area.site.g.13 11 -20.72647
# peri.area.site.g.14 11 -16.71695
# peri.area.site.g.15 11 -16.05726
# peri.area.site.g.16 11 -16.55457
# peri.area.site.g.17 11 -16.59543
# peri.area.site.g.18 11 -16.17146
# peri.area.site.g.19 11 -16.15917
# peri.area.site.g.20 11 -22.22782
# peri.area.site.g.21 11 -16.31189
# peri.area.site.g.22 11 -17.19929
# RESULTS: Best models are 2, 7 (same except they switch ocean for FC1), and 20 

summary(peri.area.site.g.2)
#                              Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                   0.59718    0.32345   1.846  0.07245 . 
# OceanPacific                  0.19371    0.25209   0.768  0.44686   
# rPC1.env.global               0.13403    0.11878   1.128  0.26602   
# rPC2.env.global               0.02811    0.12560   0.224  0.82406   
# rPC3.env.global               0.05110    0.12023   0.425  0.67315   
# rFC1                         -0.12122    0.29003  -0.418  0.67828   
# rFC2                          0.65266    0.19618   3.327  0.00192 **
# rPC1.zos.site                -0.36532    0.17805  -2.052  0.04695 * 
# rPC2.zos.site                -0.16737    0.14088  -1.188  0.24201   
# OceanPacific:rPC1.env.global -0.56207    0.22662  -2.480  0.01756 * 
# Residual standard error: 0.1591 on 39 degrees of freedom
# Multiple R-squared:  0.4875,	Adjusted R-squared:  0.3693 
# F-statistic: 4.122 on 9 and 39 DF,  p-value: 0.0008728
# NOTE: In all three best models, FC2 and PCz1 significantly affect periphyton.  

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(peri.area.site.g.2)
res = residuals(peri.area.site.g.2, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49$rPC1.env.global,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC2.env.global,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC3.env.global,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rFC1,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rFC2,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC1.zos.site,res, xlab = "PCz1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC2.zos.site,res, xlab = "PCz2 (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: Platykurtic. Tails depart from qq line.

# RESULTS and INTERPRETATION: Periphyton per area: Best model (2) explains 37% of variation. 
# All three best models are consistent in showing that periphyton per bottom area significantly 
# increases with genetic FC2, and is higher in forest-like stands (low PCz1).


###################################################################################
# GLM: MESOGRAZERS PER G EELGRASS (GLOBAL) - MODEL USING SITE MEANS               #
###################################################################################

# NOTE: Periphyton and its interactions have been removed. Rationale is that we are 
# near (arguably past) the limit of power with number of predictors and removing periphyton 
# thus simplifies the models, periphyton varies substantially within sites and should thus be 
# modeled at plot level, and inclusion of periphyton kills effects of local-regional 
# environmental predictors (especially PCe3) that are clearly important in bivariate plots. 

# Main effects only
meso.site.g.1 <- lm(rmesograzer.mass.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site 
                    , data = ZEN_2014_site_means_49)
summary(meso.site.g.1)

# Add interaction: ocean x PCe1
meso.site.g.2 <- lm(rmesograzer.mass.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site 
                    + Ocean*rPC1.env.global # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe2
meso.site.g.3 <- lm(rmesograzer.mass.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site 
                    + Ocean*rPC2.env.global # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe3
meso.site.g.4 <- lm(rmesograzer.mass.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site 
                    + Ocean*rPC3.env.global # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC1
meso.site.g.5 <- lm(rmesograzer.mass.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site 
                    + Ocean*rFC1 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC2
meso.site.g.6 <- lm(rmesograzer.mass.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site 
                    + Ocean*rFC2 # added
                    , data = ZEN_2014_site_means_49)


# Add interaction:  PCe1 x FC1
meso.site.g.7 <- lm(rmesograzer.mass.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site 
                    + rPC1.env.global*rFC1 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC1
meso.site.g.8 <- lm(rmesograzer.mass.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site 
                    + rPC2.env.global*rFC1 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC1
meso.site.g.9 <- lm(rmesograzer.mass.perg ~ Ocean
                    + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                    + rPC1.zos.site + rPC2.zos.site 
                    + rPC3.env.global*rFC1 # added
                    , data = ZEN_2014_site_means_49)

# Add interaction:  PCe1 x FC2
meso.site.g.10 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.env.global*rFC2 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC2
meso.site.g.11 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC2.env.global*rFC2 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC2
meso.site.g.12 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC3.env.global*rFC2 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x PCe1
meso.site.g.13 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.zos.site*rPC1.env.global # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x PCe2
meso.site.g.14 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.zos.site*rPC2.env.global # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x PCe3
meso.site.g.15 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.zos.site*rPC3.env.global # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x FC1
meso.site.g.16 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.zos.site*rFC1 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x FC2
meso.site.g.17 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC1.zos.site*rFC2 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x PCe1
meso.site.g.18 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC2.zos.site*rPC1.env.global # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x PCe2
meso.site.g.19 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC2.zos.site*rPC2.env.global # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x PCe3
meso.site.g.20 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC2.zos.site*rPC3.env.global # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x FC1
meso.site.g.21 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC2.zos.site*rFC1 # added
                     , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x FC2
meso.site.g.22 <- lm(rmesograzer.mass.perg ~ Ocean
                     + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                     + rPC1.zos.site + rPC2.zos.site 
                     + rPC2.zos.site*rFC2 # added
                     , data = ZEN_2014_site_means_49)


AICc(meso.site.g.1, meso.site.g.2, meso.site.g.3, meso.site.g.4, meso.site.g.5, 
    meso.site.g.6, meso.site.g.7, meso.site.g.8, meso.site.g.9, meso.site.g.10, meso.site.g.11, 
    meso.site.g.12, meso.site.g.13, meso.site.g.14, meso.site.g.15, meso.site.g.16, meso.site.g.17, 
    meso.site.g.18, meso.site.g.19, meso.site.g.20, meso.site.g.21, meso.site.g.22) 

#               df      AIC
# meso.site.g.1  10 -34.65561
# meso.site.g.2  11 -31.46953
# meso.site.g.3  11 -31.36904
# meso.site.g.4  11 -32.76033
# meso.site.g.5  11 -32.51865
# meso.site.g.6  11 -36.94808
# meso.site.g.7  11 -31.45683
# meso.site.g.8  11 -31.38716
# meso.site.g.9  11 -33.28747
# meso.site.g.10 11 -33.04097
# meso.site.g.11 11 -34.74199
# meso.site.g.12 11 -32.27862
# meso.site.g.13 11 -32.08211
# meso.site.g.14 11 -36.13550
# meso.site.g.15 11 -37.42215
# meso.site.g.16 11 -38.41982
# meso.site.g.17 11 -31.50373
# meso.site.g.18 11 -31.32481
# meso.site.g.19 11 -37.11937
# meso.site.g.20 11 -31.41673
# meso.site.g.21 11 -31.64239
# meso.site.g.22 11 -31.35495
# RESULTS: Best model is 16, but 6 and 19 are all close. 

# Obtain model-averaged, standardized coefficients from models with AICc < 2.0 from best 
summary(model.avg(meso.site.g.6, meso.site.g.16, meso.site.g.19)) # Strong effect of PCe3
# Model-averaged coefficients:  
# (full average) 
#                                 Estimate Std. Error Adjusted SE z value Pr(>|z|)    
# (Intercept)                    0.0187309  0.3880389   0.3951032   0.047 0.962188    
# OceanPacific                   0.3531320  0.3468472   0.3520422   1.003 0.315815    
# rPC1.env.global                0.0005085  0.1000915   0.1028800   0.005 0.996057    
# rPC2.env.global                0.1093142  0.3248382   0.3285235   0.333 0.739328    
# rPC3.env.global                0.3760452  0.1039627   0.1072161   3.507 0.000453 ***
# rFC1                          -0.2401905  0.3581213   0.3660171   0.656 0.511678    
# rFC2                           0.7313938  1.4480056   1.4587027   0.501 0.616089    
# rPC1.zos.site                  0.0798183  0.2745410   0.2789664   0.286 0.774785    
# rPC2.zos.site                 -0.1896438  0.3467644   0.3508893   0.540 0.588875    
# rFC1:rPC1.zos.site             0.4991415  0.5755490   0.5801554   0.860 0.389592    
# rPC2.env.global:rPC2.zos.site  0.3425036  0.6511957   0.6557544   0.522 0.601458    
# OceanPacific:rFC2             -0.7249399  1.4598035   1.4700481   0.493 0.621914 

summary(meso.site.g.16)
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.09283    0.26969   0.344 0.732528    
# OceanPacific        0.21653    0.18628   1.162 0.252133    
# rPC1.env.global     0.03672    0.09497   0.387 0.701091    
# rPC2.env.global     0.30920    0.10818   2.858 0.006804 ** 
# rPC3.env.global     0.39000    0.10240   3.809 0.000483 ***
# rFC1               -0.35152    0.30319  -1.159 0.253344    
# rFC2               -0.01962    0.14103  -0.139 0.890090    
# rPC1.zos.site      -0.11124    0.22539  -0.494 0.624395    
# rPC2.zos.site      -0.05137    0.13078  -0.393 0.696625    
# rFC1:rPC1.zos.site  0.99879    0.40473   2.468 0.018087 *  
# Residual standard error: 0.1361 on 39 degrees of freedom
# Multiple R-squared:   0.53,	Adjusted R-squared:  0.4216 
# F-statistic: 4.887 on 9 and 39 DF,  p-value: 0.0002098


# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(meso.site.g.16)
res = residuals(meso.site.g.16, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49$rPC1.env.global,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC2.env.global,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC3.env.global,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rFC1,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rFC2,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC1.zos.site,res, xlab = "PCz1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC2.zos.site,res, xlab = "PCz2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rperiphyton.perg,res, xlab = "periphyton (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: Looks pretty good.

# RESULTS AND INTERPRETATION: Global analysis shows that mesograzer biomass per g plant is 
# strongly higher in productive, nutrient rich estuaries (high PCe3, PCe2). Model averaging 
# confirms strong effect of PCe3 but not PCe2. Best model explains 42% of variation.


###################################################################################
# GLM: MESOGRAZERS PER BOTTOM AREA (GLOBAL) - MODEL USING SITE MEANS              #
###################################################################################

# The ecosytem-scale effects of eelgrass on associated organisms is best expressed
# on the basis of per unit area. 

# Main effects only
meso.area.site.g.1 <- lm(rmesograzer.mass ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site 
                         , data = ZEN_2014_site_means_49)
summary(meso.area.site.g.1)
#                  Estimate Std. Error t value Pr(>|t|)   
# (Intercept)      0.209529   0.262366   0.799  0.42923   
# OceanPacific     0.135382   0.183030   0.740  0.46382   
# rPC1.env.global -0.018660   0.089329  -0.209  0.83560   
# rPC2.env.global  0.287859   0.106019   2.715  0.00973 **
# rPC3.env.global  0.340452   0.100545   3.386  0.00160 **
# rFC1             0.042117   0.244484   0.172  0.86409   
# rFC2             0.005721   0.139529   0.041  0.96750   
# rPC1.zos.site    0.098146   0.150677   0.651  0.51853   
# rPC2.zos.site   -0.171455   0.118697  -1.444  0.15639   
# Residual standard error: 0.1347 on 40 degrees of freedom
# Multiple R-squared:  0.5351,	Adjusted R-squared:  0.4421 
# F-statistic: 5.754 on 8 and 40 DF,  p-value: 6.855e-05


# Add interaction: ocean x PCe1
meso.area.site.g.2 <- lm(rmesograzer.mass ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site 
                         + Ocean*rPC1.env.global # added
                         , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe2
meso.area.site.g.3 <- lm(rmesograzer.mass ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site 
                         + Ocean*rPC2.env.global # added
                         , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe3
meso.area.site.g.4 <- lm(rmesograzer.mass ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site 
                         + Ocean*rPC3.env.global # added
                         , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC1
meso.area.site.g.5 <- lm(rmesograzer.mass ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site 
                         + Ocean*rFC1 # added
                         , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC2
meso.area.site.g.6 <- lm(rmesograzer.mass ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site 
                         + Ocean*rFC2 # added
                         , data = ZEN_2014_site_means_49)


# Add interaction:  PCe1 x FC1
meso.area.site.g.7 <- lm(rmesograzer.mass ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site 
                         + rPC1.env.global*rFC1 # added
                         , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC1
meso.area.site.g.8 <- lm(rmesograzer.mass ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site 
                         + rPC2.env.global*rFC1 # added
                         , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC1
meso.area.site.g.9 <- lm(rmesograzer.mass ~ Ocean
                         + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                         + rPC1.zos.site + rPC2.zos.site 
                         + rPC3.env.global*rFC1 # added
                         , data = ZEN_2014_site_means_49)

# Add interaction:  PCe1 x FC2
meso.area.site.g.10 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.env.global*rFC2 # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC2
meso.area.site.g.11 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC2.env.global*rFC2 # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC2
meso.area.site.g.12 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC3.env.global*rFC2 # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x PCe1
meso.area.site.g.13 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.zos.site*rPC1.env.global # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x PCe2
meso.area.site.g.14 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.zos.site*rPC2.env.global # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x PCe3
meso.area.site.g.15 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.zos.site*rPC3.env.global # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x FC1
meso.area.site.g.16 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.zos.site*rFC1 # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz1 x FC2
meso.area.site.g.17 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.zos.site*rFC2 # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x PCe1
meso.area.site.g.18 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC2.zos.site*rPC1.env.global # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x PCe2
meso.area.site.g.19 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC2.zos.site*rPC2.env.global # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x PCe3
meso.area.site.g.20 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC2.zos.site*rPC3.env.global # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x FC1
meso.area.site.g.21 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC2.zos.site*rFC1 # added
                          , data = ZEN_2014_site_means_49)

# Add interaction:  PCz2 x FC2
meso.area.site.g.22 <- lm(rmesograzer.mass ~ Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC2.zos.site*rFC2 # added
                          , data = ZEN_2014_site_means_49)


AICc(meso.area.site.g.1, meso.area.site.g.2, meso.area.site.g.3, meso.area.site.g.4, meso.area.site.g.5, 
    meso.area.site.g.6, meso.area.site.g.7, meso.area.site.g.8, meso.area.site.g.9, meso.area.site.g.10, meso.area.site.g.11, 
    meso.area.site.g.12, meso.area.site.g.13, meso.area.site.g.14, meso.area.site.g.15, meso.area.site.g.16, meso.area.site.g.17, 
    meso.area.site.g.18, meso.area.site.g.19, meso.area.site.g.20, meso.area.site.g.21, meso.area.site.g.22) 

#               df      AIC
# meso.area.site.g.1  10 -41.57172
# meso.area.site.g.2  11 -38.60968
# meso.area.site.g.3  11 -38.22789
# meso.area.site.g.4  11 -39.19727
# meso.area.site.g.5  11 -39.05279
# meso.area.site.g.6  11 -44.44353
# meso.area.site.g.7  11 -38.53132
# meso.area.site.g.8  11 -38.44342
# meso.area.site.g.9  11 -39.60185
# meso.area.site.g.10 11 -39.18851
# meso.area.site.g.11 11 -42.16055
# meso.area.site.g.12 11 -39.52147
# meso.area.site.g.13 11 -39.33703
# meso.area.site.g.14 11 -44.54738
# meso.area.site.g.15 11 -45.06363
# meso.area.site.g.16 11 -47.23375
# meso.area.site.g.17 11 -38.38741
# meso.area.site.g.18 11 -38.22950
# meso.area.site.g.19 11 -44.11703
# meso.area.site.g.20 11 -39.06685
# meso.area.site.g.21 11 -38.73469
# meso.area.site.g.22 11 -38.33052
# RESULTS: Best model is 16. 

summary(meso.area.site.g.16)
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.335369   0.246491   1.361  0.18146    
# OceanPacific        0.191365   0.170257   1.124  0.26789    
# rPC1.env.global     0.056849   0.086801   0.655  0.51635    
# rPC2.env.global     0.249782   0.098877   2.526  0.01570 *  
# rPC3.env.global     0.372592   0.093588   3.981  0.00029 ***
# rFC1               -0.408358   0.277116  -1.474  0.14862    
# rFC2                0.002536   0.128903   0.020  0.98440    
# rPC1.zos.site      -0.327895   0.206003  -1.592  0.11953    
# rPC2.zos.site      -0.304920   0.119528  -2.551  0.01478 *  
# rFC1:rPC1.zos.site  1.037793   0.369918   2.805  0.00780 ** 
# Residual standard error: 0.1244 on 39 degrees of freedom
# Multiple R-squared:  0.6131,	Adjusted R-squared:  0.5239 
# F-statistic: 6.868 on 9 and 39 DF,  p-value: 7.506e-06
# RESULT: Strong increase in mesograzer biomass per area in productive estuaries,
# and where eelgrass biomass is high.

# Best model (16) had no main or interactive effect of ocean. Drop from model: 
meso.area.site.g.23 <- lm(rmesograzer.mass ~ # Ocean
                            + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                          + rPC1.zos.site + rPC2.zos.site 
                          + rPC1.zos.site*rFC1 # added
                          , data = ZEN_2014_site_means_49)
AICc(meso.area.site.g.16, meso.area.site.g.23) # Slightly better
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.55933    0.14557   3.842 0.000426 ***
# rPC1.env.global     0.02925    0.08353   0.350 0.728065    
# rPC2.env.global     0.24453    0.09909   2.468 0.017971 *  
# rPC3.env.global     0.39030    0.09256   4.217 0.000138 ***
# rFC1               -0.62216    0.20218  -3.077 0.003763 ** 
# rFC2                0.02618    0.12759   0.205 0.838454    
# rPC1.zos.site      -0.36261    0.20434  -1.775 0.083588 .  
# rPC2.zos.site      -0.28954    0.11913  -2.430 0.019661 *  
# rFC1:rPC1.zos.site  0.98906    0.36858   2.683 0.010546 *  
# Residual standard error: 0.1248 on 40 degrees of freedom
# Multiple R-squared:  0.6006,	Adjusted R-squared:  0.5207 
# F-statistic: 7.519 on 8 and 40 DF,  p-value: 4.506e-06

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(meso.area.site.g.16)
res = residuals(meso.area.site.g.16, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49$rPC1.env.global,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC2.env.global,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC3.env.global,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rFC1,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rFC2,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC1.zos.site,res, xlab = "PCz1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rPC2.zos.site,res, xlab = "PCz2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49$rperiphyton.perg,res, xlab = "periphyton (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: OK

# RESULTS: Mesograzer biomass per area increases strongly in productive estuaries (high PCe3 and to lesser
# extent high PCe2) and where eelgrass biomass is greater (low PCz2). There is a strong interaction between
# genetic FCA1 and eelgrass growth form (PCz1). Ocean has no main or interactive effect. If ocean 
# (non-significant) is removed, the model is slightly better and FCA1 becomes highly significant. Thus, the 
# effect of eelgrass genetics on mesograzer biomass is less than that of environment but still 
# significant. 


###################################################################################
# GLM: EELGRASS GROWTH FORM PC1 (ATLANTIC) - MODEL USING SITE MEANS               #
###################################################################################

# Main effects only
pcz1.site.a.1 <- lm(rPC1.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl  
                    , data = ZEN_2014_site_means_49_Atlantic)
summary(pcz1.site.a.1)

# Add interaction:  PCe1 x FC1
pcz1.site.a.2 <- lm(rPC1.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.env.global.atl*rFC1.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC1
pcz1.site.a.3 <- lm(rPC1.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC2.env.global.atl*rFC1.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC1
pcz1.site.a.4 <- lm(rPC1.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC3.env.global.atl*rFC1.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe1 x FC2
pcz1.site.a.5 <- lm(rPC1.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.env.global.atl*rFC2.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC2
pcz1.site.a.6 <- lm(rPC1.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC2.env.global.atl*rFC2.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC2
pcz1.site.a.7 <- lm(rPC1.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC3.env.global.atl*rFC2.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)


AICc(pcz1.site.a.1, pcz1.site.a.2, pcz1.site.a.3, pcz1.site.a.4, pcz1.site.a.5, pcz1.site.a.6, pcz1.site.a.7) 
#               df      AICc
# pcz1.site.a.1  7 -9.993391
# pcz1.site.a.2  8 -6.299913
# pcz1.site.a.3  8 -6.128481
# pcz1.site.a.4  8 -8.738088
# pcz1.site.a.5  8 -6.335890
# pcz1.site.a.6  8 -6.332210
# pcz1.site.a.7  8 -6.479100

# RESULTS: Best model is 1 (no interactions), but 4 is also close. Since
# the single interaction in model 4 is non-significant, we proceed with the main
# effects only model (1) for simplicity and ease of interpretation. 

summary(pcz1.site.a.1)
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           0.7208     0.1225   5.885 5.35e-06 ***
# rPC1.env.global.atl  -0.4515     0.1130  -3.994 0.000571 ***
# rPC2.env.global.atl  -0.3944     0.1943  -2.030 0.054079 .  
# rPC3.env.global.atl  -0.2709     0.1252  -2.164 0.041129 *  
# rFC1.global.atl      -0.2903     0.2301  -1.262 0.219760    
# rFC2.global.atl       0.6418     0.2087   3.076 0.005343 ** 
# Residual standard error: 0.1639 on 23 degrees of freedom
# Multiple R-squared:  0.609,	Adjusted R-squared:  0.524 
# F-statistic: 7.165 on 5 and 23 DF,  p-value: 0.0003583

# Test robustness of best model (1) by dropping correlated predictors (PCe2, FC2)
# Drop FC2
pcz1.site.a.8 <- lm(rPC1.zos.atl ~ 
  + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl # + rFC2.global.atl  
  , data = ZEN_2014_site_means_49_Atlantic)
summary(pcz1.site.a.8) # RESULT: PCe1 and PCe3 are significant
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           0.7976     0.1395   5.718 6.84e-06 ***
# rPC1.env.global.atl  -0.4853     0.1308  -3.709   0.0011 ** 
# rPC2.env.global.atl  -0.2422     0.2185  -1.109   0.2785    
# rPC3.env.global.atl  -0.3158     0.1446  -2.184   0.0390 *  
# rFC1.global.atl       0.1106     0.2206   0.501   0.6208    
# Residual standard error: 0.1906 on 24 degrees of freedom
# Multiple R-squared:  0.4482,	Adjusted R-squared:  0.3562 
# F-statistic: 4.873 on 4 and 24 DF,  p-value: 0.005086

# Drop PCe2
pcz1.site.a.9 <- lm(rPC1.zos.atl ~ 
  + rPC1.env.global.atl # + rPC2.env.global.atl 
  + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl  
  , data = ZEN_2014_site_means_49_Atlantic)
summary(pcz1.site.a.9) # RESULT: PCe1 and FC2 are significant
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           0.6420     0.1235   5.198 2.52e-05 ***
# rPC1.env.global.atl  -0.4555     0.1201  -3.791 0.000892 ***
# rPC3.env.global.atl  -0.2483     0.1326  -1.873 0.073331 .  
# rFC1.global.atl      -0.4631     0.2273  -2.037 0.052796 .  
# rFC2.global.atl       0.5340     0.2145   2.489 0.020120 *  
# Residual standard error: 0.1742 on 24 degrees of freedom
# Multiple R-squared:  0.539,	Adjusted R-squared:  0.4621 
# F-statistic: 7.014 on 4 and 24 DF,  p-value: 0.0006888

# Drop FC1
pcz1.site.a.10 <- lm(rPC1.zos.atl ~ 
  + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl # + rFC1.global.atl 
  + rFC2.global.atl  
  , data = ZEN_2014_site_means_49_Atlantic)
summary(pcz1.site.a.10) # 
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           0.7131     0.1238   5.758 6.19e-06 ***
# rPC1.env.global.atl  -0.4318     0.1133  -3.810  0.00085 ***
# rPC2.env.global.atl  -0.4850     0.1827  -2.654  0.01388 *  
# rPC3.env.global.atl  -0.3185     0.1208  -2.636  0.01448 *  
# rFC2.global.atl       0.4928     0.1741   2.831  0.00925 ** 
# Residual standard error: 0.1659 on 24 degrees of freedom
# Multiple R-squared:  0.582,	Adjusted R-squared:  0.5123 
# F-statistic: 8.353 on 4 and 24 DF,  p-value: 0.0002274

# Fit best model with unstandardized data to allow direct comparison betwen oceans
pcz1.site.a.1.raw <- lm(PC1.zos.site ~ 
                      + PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2
                    , data = ZEN_2014_site_means_49_Atlantic)

summary(pcz1.site.a.1.raw)
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     1.923493   1.032189   1.864 0.075209 .  
# PC1.env.global -0.334939   0.083857  -3.994 0.000571 ***
# PC2.env.global -0.335001   0.165025  -2.030 0.054079 .  
# PC3.env.global -0.222482   0.102830  -2.164 0.041129 *  
# FC1            -0.002196   0.001741  -1.262 0.219760    
# FC2             0.008304   0.002700   3.076 0.005343 ** 
# Residual standard error: 0.7251 on 23 degrees of freedom
# Multiple R-squared:  0.609,	Adjusted R-squared:  0.524 
# F-statistic: 7.165 on 5 and 23 DF,  p-value: 0.0003583

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(pcz1.site.a.1)
res = residuals(pcz1.site.a.1, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49_Atlantic$rPC1.env.global.atl,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC3.env.global.atl,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC1.global.atl,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC2.global.atl,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: Pretty good. 


# RESULTS and INTERPRETATION: In the Atlantic, meadow form (high PCz1) is best developed at 
# cooler, high-latitudes (low PCe1).There is also a strong positive effect of genetic FC2.

# PCe2, FC1, and FC2 are correlated. Tests of robustness show that PCe2 is not significant 
# in full model or when when FC2 is dropped (ad marginal when FC1 is dropped), so PCe2 effect
# is weak and unstable. PCe1 and FC2 are highly significant in all models so are strong and 
# little affected by collinear predictors.  


###################################################################################
# GLM: EELGRASS GROWTH FORM PC2 (ATLANTIC) - MODEL USING SITE MEANS               #
###################################################################################

# Main effects only
pcz2.test.a.1 <- lm(rPC2.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl  
                    , data = ZEN_2014_site_means_49_Atlantic)
summary(pcz2.test.a.1)

# Add interaction:  PCe1 x FC2
pcz2.test.a.2 <- lm(rPC2.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.env.global.atl*rFC2.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC2
pcz2.test.a.3 <- lm(rPC2.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC2.env.global.atl*rFC2.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC2
pcz2.test.a.4 <- lm(rPC2.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC3.env.global.atl*rFC2.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe1 x FC1
pcz2.test.a.5 <- lm(rPC2.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.env.global.atl*rFC1.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC1
pcz2.test.a.6 <- lm(rPC2.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC2.env.global.atl*rFC1.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC1
pcz2.test.a.7 <- lm(rPC2.zos.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC3.env.global.atl*rFC1.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

AICc(pcz2.test.a.1, pcz2.test.a.2, pcz2.test.a.3, pcz2.test.a.4, pcz2.test.a.5, pcz2.test.a.6, pcz2.test.a.7) 
#               df        AICc
# pcz2.test.a.1  7  -2.5426367
# pcz2.test.a.2  8 -12.3053003
# pcz2.test.a.3  8  -0.2554764
# pcz2.test.a.4  8   0.2837104
# pcz2.test.a.5  8  -3.7306188
# pcz2.test.a.6  8   1.2917508
# pcz2.test.a.7  8   1.3240189
# Model 2 is by far the best 

summary(pcz2.test.a.2)
#                                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           0.4049     0.1415   2.861 0.009075 ** 
# rPC1.env.global.atl                   0.7821     0.2575   3.037 0.006045 ** 
# rPC2.env.global.atl                  -0.2667     0.1786  -1.494 0.149471    
# rPC3.env.global.atl                  -0.4331     0.1153  -3.757 0.001089 ** 
# rFC1.global.atl                      -0.2617     0.2141  -1.222 0.234475    
# rFC2.global.atl                       1.1331     0.2849   3.977 0.000637 ***
# rPC1.env.global.atl:rFC2.global.atl  -1.8983     0.5225  -3.633 0.001470 ** 
# Residual standard error: 0.1506 on 22 degrees of freedom
# Multiple R-squared:  0.6048,	Adjusted R-squared:  0.497 
# F-statistic: 5.611 on 6 and 22 DF,  p-value: 0.001168

# Compare with main effects only model
summary(pcz2.test.a.1)
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.71627    0.13929   5.142 3.28e-05 ***
# rPC1.env.global.atl -0.07383    0.12853  -0.574  0.57128    
# rPC2.env.global.atl -0.27043    0.22089  -1.224  0.23324    
# rPC3.env.global.atl -0.40870    0.14238  -2.871  0.00864 ** 
# rFC1.global.atl     -0.14173    0.26167  -0.542  0.59327    
# rFC2.global.atl      0.36780    0.23726   1.550  0.13475    
# Residual standard error: 0.1863 on 23 degrees of freedom
# Multiple R-squared:  0.3677,	Adjusted R-squared:  0.2302 
# F-statistic: 2.675 on 5 and 23 DF,  p-value: 0.0478


# Test robustness of best model by dropping correlated predictors (PCe2, FC2)
# Drop FC2
pcz2.test.a.8 <- lm(rPC2.zos.atl ~ rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl 
  + rFC1.global.atl # + rFC2.global.atl
  # + rPC1.env.global.atl*rFC2.global.atl # added
  , data = ZEN_2014_site_means_49_Atlantic)
summary(pcz2.test.a.8) # RESULT: Strong effect of PCe3 only, but model P = 0.06
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.76022    0.14030   5.419 1.45e-05 ***
# rPC1.env.global.atl -0.09319    0.13161  -0.708  0.48570    
# rPC2.env.global.atl -0.18327    0.21977  -0.834  0.41255    
# rPC3.env.global.atl -0.43445    0.14548  -2.986  0.00641 ** 
# rFC1.global.atl      0.08798    0.22188   0.397  0.69522    
# Residual standard error: 0.1917 on 24 degrees of freedom
# Multiple R-squared:  0.3016,	Adjusted R-squared:  0.1852 
# F-statistic: 2.591 on 4 and 24 DF,  p-value: 0.0622

# Drop PCe2
pcz2.test.a.9 <- lm(rPC2.zos.atl ~ rPC1.env.global.atl # + rPC2.env.global.atl 
  + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl+ rPC1.env.global.atl*rFC2.global.atl 
  , data = ZEN_2014_site_means_49_Atlantic)
summary(pcz2.test.a.9) # RESULT: Strong effects of PCe1 and 3, FC2
#                                     Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                           0.3509     0.1404   2.499  0.02005 * 
# rPC1.env.global.atl                   0.7815     0.2643   2.957  0.00707 **
# rPC3.env.global.atl                  -0.4179     0.1179  -3.545  0.00173 **
# rFC1.global.atl                      -0.3788     0.2045  -1.853  0.07681 . 
# rFC2.global.atl                       1.0620     0.2883   3.684  0.00123 **
# rPC1.env.global.atl:rFC2.global.atl  -1.9028     0.5363  -3.548  0.00172 **
# Residual standard error: 0.1546 on 23 degrees of freedom
# Multiple R-squared:  0.5647,	Adjusted R-squared:  0.4701 
# F-statistic: 5.968 on 5 and 23 DF,  p-value: 0.001116

# Drop FC1
pcz2.test.a.10 <- lm(rPC2.zos.atl ~ rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl # + rFC1.global.atl 
  + rFC2.global.atl+ rPC1.env.global.atl*rFC2.global.atl 
  , data = ZEN_2014_site_means_49_Atlantic)
summary(pcz2.test.a.10) # RESULT: 
#                                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                           0.4143     0.1428   2.901 0.008055 ** 
# rPC1.env.global.atl                   0.7551     0.2593   2.912 0.007847 ** 
# rPC2.env.global.atl                  -0.3466     0.1679  -2.064 0.050454 .  
# rPC3.env.global.atl                  -0.4738     0.1116  -4.247 0.000305 ***
# rFC2.global.atl                       0.9622     0.2509   3.835 0.000846 ***
# rPC1.env.global.atl:rFC2.global.atl  -1.7998     0.5218  -3.449 0.002181 ** 
# Residual standard error: 0.1522 on 23 degrees of freedom
# Multiple R-squared:  0.5779,	Adjusted R-squared:  0.4862 
# F-statistic: 6.299 on 5 and 23 DF,  p-value: 0.0008062

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(pcz2.test.a.2)
res = residuals(pcz2.test.a.2, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49_Atlantic$rPC1.env.global.atl,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC3.env.global.atl,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC1.global.atl,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC2.global.atl,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: Model seems to systematically overpredict, i.e., many poiunts are below qq line. ? 

# RESULTS and INTERPRETATION: In the Atlantic, eelgrass biomass (inverse PCz2) is highest at 
# productive estuarine sites (high PCe3). Genetic FC2 also has strong effect but interaction 
# makes interpretation difficult. Best model explains 50% of variation. 

# PCe2, FC1, and FC2 are correlated. Tests of robustness show that FC1 and PCe2 are not 
# significant in any model. FC2 and PCe3 are significant in all models. 
# Thus, effects of PCe3 and FC2 are strong and unaffected by collinear predictors. 


###################################################################################
# GLM: PERIPHYTON PER G EELGRASS (ATLANTIC) - MODEL USING SITE MEANS              #
###################################################################################

# Main effects only
peri.test.a.1 <- lm(rperiphyton.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl 
                    + rPC1.zos.atl + rPC2.zos.atl
                    , data = ZEN_2014_site_means_49_Atlantic)
summary(peri.test.a.1)

# Add interaction:  PCe1 x FC1
peri.test.a.2 <- lm(rperiphyton.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl
                    + rPC1.env.global.atl*rFC1.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC1
peri.test.a.3 <- lm(rperiphyton.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl
                    + rPC2.env.global.atl*rFC1.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC1
peri.test.a.4 <- lm(rperiphyton.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl
                    + rPC3.env.global.atl*rFC1.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe1 x FC2
peri.test.a.5 <- lm(rperiphyton.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl
                    + rPC1.env.global.atl*rFC2.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC2
peri.test.a.6 <- lm(rperiphyton.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl
                    + rPC2.env.global.atl*rFC2.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC2
peri.test.a.7 <- lm(rperiphyton.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl
                    + rPC3.env.global.atl*rFC2.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)


# Add interaction:  PCz1 x PCe1
peri.test.a.8 <- lm(rperiphyton.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl
                    + rPC1.zos.atl*rPC1.env.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x PCe2
peri.test.a.9 <- lm(rperiphyton.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl
                    + rPC1.zos.atl*rPC2.env.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x PCe3
peri.test.a.10 <- lm(rperiphyton.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl
                     + rPC1.zos.atl*rPC3.env.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x FC1
peri.test.a.11 <- lm(rperiphyton.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl
                     + rPC1.zos.atl*rFC1.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x FC2
peri.test.a.12 <- lm(rperiphyton.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl
                     + rPC1.zos.atl*rFC2.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x PCe1
peri.test.a.13 <- lm(rperiphyton.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl
                     + rPC2.zos.atl*rPC1.env.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x PCe2
peri.test.a.14 <- lm(rperiphyton.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl
                     + rPC2.zos.atl*rPC2.env.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x PCe3
peri.test.a.15 <- lm(rperiphyton.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl
                     + rPC2.zos.atl*rPC3.env.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x FC1
peri.test.a.16 <- lm(rperiphyton.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl
                     + rPC2.zos.atl*rFC1.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x FC2
peri.test.a.17 <- lm(rperiphyton.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl
                     + rPC2.zos.atl*rFC2.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)


AICc(peri.test.a.1, peri.test.a.2, peri.test.a.3, peri.test.a.4, peri.test.a.5, 
    peri.test.a.6, peri.test.a.7, peri.test.a.8, peri.test.a.9, peri.test.a.10, peri.test.a.11, 
    peri.test.a.12, peri.test.a.13, peri.test.a.14, peri.test.a.15, peri.test.a.16, peri.test.a.17)

#                df      AICc
# peri.test.a.1   9 27.273077
# peri.test.a.2  10 30.280775
# peri.test.a.3  10 30.649207
# peri.test.a.4  10 31.593639
# peri.test.a.5  10 31.891958
# peri.test.a.6  10 28.342154
# peri.test.a.7  10 31.137947
# peri.test.a.8  10 23.698066
# peri.test.a.9  10 21.040834
# peri.test.a.10 10 30.168433
# peri.test.a.11 10 30.464230
# peri.test.a.12 10 31.992467
# peri.test.a.13 10 30.239225
# peri.test.a.14 10 30.655977
# peri.test.a.15 10 29.175472
# peri.test.a.16 10 29.239118
# peri.test.a.17 10 31.160401
# Model 9 is best.

summary(peri.test.a.9)
#                                  Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                        1.0620     0.4574   2.322  0.03092 * 
# rPC1.env.global.atl                0.3439     0.2208   1.558  0.13499   
# rPC2.env.global.atl               -1.4037     0.5764  -2.435  0.02438 * 
# rPC3.env.global.atl                0.2133     0.2212   0.965  0.34627   
# rFC1.global.atl                    0.4769     0.3546   1.345  0.19373   
# rFC2.global.atl                   -0.5073     0.3759  -1.350  0.19215   
# rPC1.zos.atl                      -1.5223     0.7021  -2.168  0.04237 * 
# rPC2.zos.atl                      -0.1184     0.3016  -0.393  0.69869   
# rPC2.env.global.atl:rPC1.zos.atl   3.6662     1.2083   3.034  0.00655 **
# Residual standard error: 0.2403 on 20 degrees of freedom
# Multiple R-squared:  0.4333,	Adjusted R-squared:  0.2067 
# F-statistic: 1.912 on 8 and 20 DF,  p-value: 0.1146


# Test robustness of best model by dropping correlated predictors (PCe2, FC2)
# Drop FC2
peri.test.a.18 <- lm(rperiphyton.perg.atl ~ 
  + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl # + rFC2.global.atl
  + rPC1.zos.atl + rPC2.zos.atl + rPC1.zos.atl*rPC2.env.global.atl 
  , data = ZEN_2014_site_means_49_Atlantic)
summary(peri.test.a.18) # RESULT: PCe2 abd PCz1 significant as in best model (but model not significant)
#                                  Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                        1.1032     0.4652   2.371   0.0274 *
# rPC1.env.global.atl                0.2804     0.2199   1.275   0.2162  
# rPC2.env.global.atl               -1.4319     0.5872  -2.438   0.0237 *
# rPC3.env.global.atl                0.1606     0.2219   0.724   0.4773  
# rFC1.global.atl                    0.1835     0.2856   0.642   0.5275  
# rPC1.zos.atl                      -1.5373     0.7156  -2.148   0.0435 *
# rPC2.zos.atl                      -0.1329     0.3072  -0.432   0.6698  
# rPC2.env.global.atl:rPC1.zos.atl   3.3079     1.2016   2.753   0.0119 *
# Residual standard error: 0.2449 on 21 degrees of freedom
# Multiple R-squared:  0.3817,	Adjusted R-squared:  0.1756 
# F-statistic: 1.852 on 7 and 21 DF,  p-value: 0.1296

# Drop PCe2
peri.test.a.19 <- lm(rperiphyton.perg.atl ~ 
  + rPC1.env.global.atl # + rPC2.env.global.atl 
  + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
  + rPC1.zos.atl + rPC2.zos.atl # + rPC1.zos.atl*rPC2.env.global.atl 
  , data = ZEN_2014_site_means_49_Atlantic)
summary(peri.test.a.19) # RESULT: Nothing
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)          0.09295    0.31756   0.293   0.7725  
# rPC1.env.global.atl  0.46625    0.24325   1.917   0.0683 .
# rPC3.env.global.atl  0.17896    0.24560   0.729   0.4739  
# rFC1.global.atl      0.30653    0.39264   0.781   0.4433  
# rFC2.global.atl     -0.22402    0.38404  -0.583   0.5656  
# rPC1.zos.atl         0.35522    0.34450   1.031   0.3137  
# rPC2.zos.atl         0.22093    0.31877   0.693   0.4955  
# Residual standard error: 0.2771 on 22 degrees of freedom
# Multiple R-squared:  0.1712,	Adjusted R-squared:  -0.05482 
# F-statistic: 0.7575 on 6 and 22 DF,  p-value: 0.6106

# Drop FC1
peri.test.a.20 <- lm(rperiphyton.perg.atl ~ 
  + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl # + rFC1.global.atl 
  + rFC2.global.atl + rPC1.zos.atl + rPC2.zos.atl + rPC1.zos.atl*rPC2.env.global.atl 
  , data = ZEN_2014_site_means_49_Atlantic)
summary(peri.test.a.20) # RESULT: 
#                                  Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                        1.0783     0.4660   2.314   0.0309 *
# rPC1.env.global.atl                0.2793     0.2196   1.272   0.2173  
# rPC2.env.global.atl               -1.1986     0.5665  -2.116   0.0465 *
# rPC3.env.global.atl                0.2495     0.2237   1.115   0.2774  
# rFC2.global.atl                   -0.1975     0.3026  -0.653   0.5210  
# rPC1.zos.atl                      -1.4740     0.7145  -2.063   0.0517 .
# rPC2.zos.atl                      -0.1075     0.3072  -0.350   0.7298  
# rPC2.env.global.atl:rPC1.zos.atl   3.3839     1.2126   2.791   0.0110 *
# Residual standard error: 0.2449 on 21 degrees of freedom
# Multiple R-squared:  0.3821,	Adjusted R-squared:  0.1761 
# F-statistic: 1.855 on 7 and 21 DF,  p-value: 0.129

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(peri.test.a.11)
res = residuals(peri.test.a.11, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49_Atlantic$rPC1.env.global.atl,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC3.env.global.atl,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC1.global.atl,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC2.global.atl,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC1.zos.atl,res, xlab = "PCz1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC2.zos.atl,res, xlab = "PCz2 (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: Generally OK. Possibly underpredicted. 

# RESULTS and INTERPRETATION: In the Atlantic, none of our predictors explains periphyton 
# well at the site level. "Best" model has P = 0.11, suggests slightly lower periphyton 
# load in meadow form of eelgrass and under low nutrient conditions (oddly).  

# PCe2, FC1, and FC2 are correlated. Tests of robustness show that none of these is very stable.  


###################################################################################
# GLM: PERIPHYTON PER BOTTOM AREA (ATLANTIC) -  MODEL USING SITE MEANS            #
###################################################################################

# Main effects only
peri.area.test.a.1 <- lm(rperiphyton.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl 
                         + rPC1.zos.atl + rPC2.zos.atl
                         , data = ZEN_2014_site_means_49_Atlantic)
summary(peri.area.test.a.1)

# Add interaction:  PCe1 x FC1
peri.area.test.a.2 <- lm(rperiphyton.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl
                         + rPC1.env.global.atl*rFC1.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC1
peri.area.test.a.3 <- lm(rperiphyton.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl
                         + rPC2.env.global.atl*rFC1.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC1
peri.area.test.a.4 <- lm(rperiphyton.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl
                         + rPC3.env.global.atl*rFC1.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe1 x FC2
peri.area.test.a.5 <- lm(rperiphyton.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl
                         + rPC1.env.global.atl*rFC2.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC2
peri.area.test.a.6 <- lm(rperiphyton.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl
                         + rPC2.env.global.atl*rFC2.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC2
peri.area.test.a.7 <- lm(rperiphyton.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl
                         + rPC3.env.global.atl*rFC2.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)


# Add interaction:  PCz1 x PCe1
peri.area.test.a.8 <- lm(rperiphyton.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl
                         + rPC1.zos.atl*rPC1.env.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x PCe2
peri.area.test.a.9 <- lm(rperiphyton.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl
                         + rPC1.zos.atl*rPC2.env.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x PCe3
peri.area.test.a.10 <- lm(rperiphyton.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl
                          + rPC1.zos.atl*rPC3.env.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x FC1
peri.area.test.a.11 <- lm(rperiphyton.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl
                          + rPC1.zos.atl*rFC1.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x FC2
peri.area.test.a.12 <- lm(rperiphyton.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl
                          + rPC1.zos.atl*rFC2.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x PCe1
peri.area.test.a.13 <- lm(rperiphyton.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl
                          + rPC2.zos.atl*rPC1.env.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x PCe2
peri.area.test.a.14 <- lm(rperiphyton.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl
                          + rPC2.zos.atl*rPC2.env.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x PCe3
peri.area.test.a.15 <- lm(rperiphyton.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl
                          + rPC2.zos.atl*rPC3.env.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x FC1
peri.area.test.a.16 <- lm(rperiphyton.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl
                          + rPC2.zos.atl*rFC1.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x FC2
peri.area.test.a.17 <- lm(rperiphyton.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl
                          + rPC2.zos.atl*rFC2.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)


AICc(peri.area.test.a.1, peri.area.test.a.2, peri.area.test.a.3, peri.area.test.a.4, peri.area.test.a.5, 
    peri.area.test.a.6, peri.area.test.a.7, peri.area.test.a.8, peri.area.test.a.9, peri.area.test.a.10, peri.area.test.a.11, 
    peri.area.test.a.12, peri.area.test.a.13, peri.area.test.a.14, peri.area.test.a.15, peri.area.test.a.16, peri.area.test.a.17)

#                     df       AICc
# peri.area.test.a.1   9  8.9072080
# peri.area.test.a.2  10 11.4603545
# peri.area.test.a.3  10 12.4702022
# peri.area.test.a.4  10 13.2777991
# peri.area.test.a.5  10 13.2917188
# peri.area.test.a.6  10  9.8691423
# peri.area.test.a.7  10 12.3083440
# peri.area.test.a.8  10  5.9795664
# peri.area.test.a.9  10  1.5734436
# peri.area.test.a.10 10 10.3105251
# peri.area.test.a.11 10 12.5447048
# peri.area.test.a.12 10 13.4449443
# peri.area.test.a.13 10 11.0074699
# peri.area.test.a.14 10 12.5915620
# peri.area.test.a.15 10  7.7587748
# peri.area.test.a.16 10 10.6095365
# peri.area.test.a.17 10 13.1886674
# Model 9 is best.

summary(peri.area.test.a.9)
#                                  Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                        1.1067     0.3270   3.384  0.00295 **
# rPC1.env.global.atl                0.3061     0.1578   1.939  0.06669 . 
# rPC2.env.global.atl               -1.1548     0.4121  -2.802  0.01100 * 
# rPC3.env.global.atl                0.1998     0.1581   1.263  0.22095   
# rFC1.global.atl                    0.3781     0.2535   1.491  0.15147   
# rFC2.global.atl                   -0.3548     0.2687  -1.321  0.20153   
# rPC1.zos.atl                      -1.2821     0.5019  -2.555  0.01889 * 
# rPC2.zos.atl                      -0.6188     0.2156  -2.870  0.00946 **
# rPC2.env.global.atl:rPC1.zos.atl   2.7772     0.8638   3.215  0.00434 **
# Residual standard error: 0.1718 on 20 degrees of freedom
# Multiple R-squared:  0.6013,	Adjusted R-squared:  0.4418 
# F-statistic: 3.771 on 8 and 20 DF,  p-value: 0.007583


# Test robustness of best model by dropping correlated predictors (PCe2, FC2)
# Drop FC2
peri.area.test.a.18 <- lm(rperiphyton.area.atl ~ 
  + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl # + rFC2.global.atl
  + rPC1.zos.atl + rPC2.zos.atl + rPC1.zos.atl*rPC2.env.global.atl 
  , data = ZEN_2014_site_means_49_Atlantic)
summary(peri.area.test.a.18)
#                                  Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                        1.1355     0.3320   3.420  0.00257 **
# rPC1.env.global.atl                0.2617     0.1569   1.667  0.11026   
# rPC2.env.global.atl               -1.1745     0.4191  -2.803  0.01066 * 
# rPC3.env.global.atl                0.1629     0.1584   1.028  0.31542   
# rFC1.global.atl                    0.1729     0.2038   0.848  0.40593   
# rPC1.zos.atl                      -1.2926     0.5106  -2.531  0.01941 * 
# rPC2.zos.atl                      -0.6289     0.2192  -2.868  0.00920 **
# rPC2.env.global.atl:rPC1.zos.atl   2.5265     0.8575   2.946  0.00771 **
# Residual standard error: 0.1748 on 21 degrees of freedom
# Multiple R-squared:  0.5665,	Adjusted R-squared:  0.4221 
# F-statistic: 3.921 on 7 and 21 DF,  p-value: 0.006932

# Drop PCe2
peri.area.test.a.19 <- lm(rperiphyton.area.atl ~ 
  + rPC1.env.global.atl # + rPC2.env.global.atl 
  + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
  + rPC1.zos.atl + rPC2.zos.atl # + rPC1.zos.atl*rPC2.env.global.atl 
  , data = ZEN_2014_site_means_49_Atlantic)
summary(peri.area.test.a.19)
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)           0.3270     0.2313   1.414   0.1715  
# rPC1.env.global.atl   0.4129     0.1772   2.331   0.0293 *
# rPC3.env.global.atl   0.1911     0.1789   1.068   0.2971  
# rFC1.global.atl       0.2263     0.2860   0.791   0.4373  
# rFC2.global.atl      -0.1852     0.2797  -0.662   0.5149  
# rPC1.zos.atl          0.1714     0.2509   0.683   0.5018  
# rPC2.zos.atl         -0.3507     0.2322  -1.510   0.1452  
# Residual standard error: 0.2018 on 22 degrees of freedom
# Multiple R-squared:  0.3946,	Adjusted R-squared:  0.2295 
# F-statistic:  2.39 on 6 and 22 DF,  p-value: 0.06255

# Drop FC1
peri.area.test.a.20 <- lm(rperiphyton.area.atl ~ 
  + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl # + rFC1.global.atl 
  + rFC2.global.atl + rPC1.zos.atl + rPC2.zos.atl + rPC1.zos.atl*rPC2.env.global.atl 
  , data = ZEN_2014_site_means_49_Atlantic)
summary(peri.area.test.a.20)
#                                  Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                        1.1196     0.3363   3.329  0.00318 **
# rPC1.env.global.atl                0.2549     0.1585   1.608  0.12269   
# rPC2.env.global.atl               -0.9922     0.4088  -2.427  0.02431 * 
# rPC3.env.global.atl                0.2284     0.1615   1.415  0.17182   
# rFC2.global.atl                   -0.1092     0.2184  -0.500  0.62215   
# rPC1.zos.atl                      -1.2439     0.5156  -2.412  0.02508 * 
# rPC2.zos.atl                      -0.6101     0.2217  -2.752  0.01195 * 
# rPC2.env.global.atl:rPC1.zos.atl   2.5533     0.8751   2.918  0.00823 **
# Residual standard error: 0.1767 on 21 degrees of freedom
# Multiple R-squared:  0.557,	Adjusted R-squared:  0.4093 
# F-statistic: 3.772 on 7 and 21 DF,  p-value: 0.00841

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(peri.area.test.a.9)
res = residuals(peri.area.test.a.9, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49_Atlantic$rPC1.env.global.atl,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC3.env.global.atl,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC1.global.atl,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC2.global.atl,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC1.zos.atl,res, xlab = "PCz1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC2.zos.atl,res, xlab = "PCz2 (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: Generally OK. 

# RESULTS and INTERPRETATION: In the Atlantic, best model for periphyton per area explains 44% of
# variation. Periphyton per area is higher in forest-like stands (low PCz1) of high biomass (low PCz2).
# Main and interactive effects of PCe2 but difficult to interpret. 

# PCe2 and its interaction are significant in all models where it is present.  
# PCz1 and PCz2 are significant in all models except where PCe2 is dropped. 


###################################################################################
# GLM: MESOGRAZER MASS PER G EELGRASS (ATLANTIC) - MODEL USING SITE MEANS         #
###################################################################################

# Periphyton removed (see explanation under global)

# Main effects only
meso.test.a.1 <- lm(rmesograzer.mass.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl 
                    + rPC1.zos.atl + rPC2.zos.atl  
                    , data = ZEN_2014_site_means_49_Atlantic)
summary(meso.test.a.1)

# Add interaction:  PCe1 x FC1
meso.test.a.2 <- lm(rmesograzer.mass.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl 
                    + rPC1.env.global.atl*rFC1.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC1
meso.test.a.3 <- lm(rmesograzer.mass.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl 
                    + rPC2.env.global.atl*rFC1.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC1
meso.test.a.4 <- lm(rmesograzer.mass.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl 
                    + rPC3.env.global.atl*rFC1.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe1 x FC2
meso.test.a.5 <- lm(rmesograzer.mass.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl 
                    + rPC1.env.global.atl*rFC2.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC2
meso.test.a.6 <- lm(rmesograzer.mass.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl 
                    + rPC2.env.global.atl*rFC2.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC2
meso.test.a.7 <- lm(rmesograzer.mass.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl 
                    + rPC3.env.global.atl*rFC2.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)


# Add interaction:  PCz1 x PCe1
meso.test.a.8 <- lm(rmesograzer.mass.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl 
                    + rPC1.zos.atl*rPC1.env.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x PCe2
meso.test.a.9 <- lm(rmesograzer.mass.perg.atl ~ 
                      + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                    + rPC1.zos.atl + rPC2.zos.atl 
                    + rPC1.zos.atl*rPC2.env.global.atl # added
                    , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x PCe3
meso.test.a.10 <- lm(rmesograzer.mass.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl 
                     + rPC1.zos.atl*rPC3.env.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x FC1
meso.test.a.11 <- lm(rmesograzer.mass.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl 
                     + rPC1.zos.atl*rFC1.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x FC2
meso.test.a.12 <- lm(rmesograzer.mass.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl 
                     + rPC1.zos.atl*rFC2.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x PCe1
meso.test.a.13 <- lm(rmesograzer.mass.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl 
                     + rPC2.zos.atl*rPC1.env.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x PCe2
meso.test.a.14 <- lm(rmesograzer.mass.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl 
                     + rPC2.zos.atl*rPC2.env.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x PCe3
meso.test.a.15 <- lm(rmesograzer.mass.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl 
                     + rPC2.zos.atl*rPC3.env.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x FC1
meso.test.a.16 <- lm(rmesograzer.mass.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl 
                     + rPC2.zos.atl*rFC1.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x FC2
meso.test.a.17 <- lm(rmesograzer.mass.perg.atl ~ 
                       + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                     + rPC1.zos.atl + rPC2.zos.atl 
                     + rPC2.zos.atl*rFC2.global.atl # added
                     , data = ZEN_2014_site_means_49_Atlantic)



AICc(meso.test.a.1, meso.test.a.2, meso.test.a.3, meso.test.a.4, meso.test.a.5, 
    meso.test.a.6, meso.test.a.7, meso.test.a.8, meso.test.a.9, meso.test.a.10, meso.test.a.11, 
    meso.test.a.12, meso.test.a.13, meso.test.a.14, meso.test.a.15, meso.test.a.16, meso.test.a.17)
#                df       AICc
# meso.test.a.1   9 -21.456048
# meso.test.a.2  10 -16.790762
# meso.test.a.3  10 -18.243300
# meso.test.a.4  10 -16.754072
# meso.test.a.5  10 -19.400349
# meso.test.a.6  10 -18.178491
# meso.test.a.7  10 -16.983702
# meso.test.a.8  10 -18.040114
# meso.test.a.9  10 -22.235008
# meso.test.a.10 10 -18.004158
# meso.test.a.11 10 -16.987964
# meso.test.a.12 10 -17.059189
# meso.test.a.13 10 -16.930336
# meso.test.a.14 10 -16.926627
# meso.test.a.15 10 -17.639560
# meso.test.a.16 10 -16.886885
# meso.test.a.17 10 -16.923887
# Model 9 is best but 1 (main effects only) is close.

summary(meso.test.a.9)
#                                  Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                      -0.11447    0.21691  -0.528  0.60350   
# rPC1.env.global.atl               0.08686    0.10470   0.830  0.41654   
# rPC2.env.global.atl               0.81013    0.27335   2.964  0.00768 **
# rPC3.env.global.atl               0.37256    0.10488   3.552  0.00200 **
# rFC1.global.atl                  -0.48491    0.16816  -2.884  0.00918 **
# rFC2.global.atl                   0.58413    0.17823   3.277  0.00377 **
# rPC1.zos.atl                      1.08850    0.33291   3.270  0.00383 **
# rPC2.zos.atl                     -0.14375    0.14300  -1.005  0.32680   
# rPC2.env.global.atl:rPC1.zos.atl -1.17419    0.57297  -2.049  0.05378 . 
# Residual standard error: 0.1139 on 20 degrees of freedom
# Multiple R-squared:  0.7969,	Adjusted R-squared:  0.7156 
# F-statistic: 9.807 on 8 and 20 DF,  p-value: 1.857e-05

summary(meso.test.a.1)
#                     Estimate Std. Error t value Pr(>|t|)   
# (Intercept)          0.20655    0.16106   1.282  0.21366   
# rPC1.env.global.atl  0.04437    0.11016   0.403  0.69117   
# rPC2.env.global.atl  0.33920    0.15890   2.135  0.04474 * 
# rPC3.env.global.atl  0.37952    0.11253   3.373  0.00288 **
# rFC1.global.atl     -0.42504    0.17777  -2.391  0.02625 * 
# rFC2.global.atl      0.50388    0.18665   2.700  0.01342 * 
# rPC1.zos.atl         0.47988    0.16149   2.972  0.00728 **
# rPC2.zos.atl        -0.25501    0.14202  -1.796  0.08696 . 
# Residual standard error: 0.1223 on 21 degrees of freedom
# Multiple R-squared:  0.7542,	Adjusted R-squared:  0.6723 
# F-statistic: 9.205 on 7 and 21 DF,  p-value: 3.332e-05

# Test robustness of best model by dropping correlated predictors (PCe2, FC2)
# Drop FC2
meso.test.a.18 <- lm(rmesograzer.mass.perg.atl ~ 
  + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl # + rFC2.global.atl
  + rPC1.zos.atl + rPC2.zos.atl + rPC1.zos.atl*rPC2.env.global.atl
  , data = ZEN_2014_site_means_49_Atlantic)
summary(meso.test.a.18) # RESULT: PCe2,3 and PCz1 significant as in best. FC1 no longer significant. 
#                                  Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                       -0.1620     0.2619  -0.619  0.54288   
# rPC1.env.global.atl                0.1600     0.1238   1.293  0.21000   
# rPC2.env.global.atl                0.8426     0.3305   2.549  0.01867 * 
# rPC3.env.global.atl                0.4333     0.1249   3.469  0.00229 **
# rFC1.global.atl                   -0.1471     0.1608  -0.915  0.37053   
# rPC1.zos.atl                       1.1058     0.4027   2.746  0.01211 * 
# rPC2.zos.atl                      -0.1271     0.1729  -0.735  0.47034   
# rPC2.env.global.atl:rPC1.zos.atl  -0.7616     0.6763  -1.126  0.27281   
# Residual standard error: 0.1379 on 21 degrees of freedom
# Multiple R-squared:  0.6878,	Adjusted R-squared:  0.5837 
# F-statistic: 6.608 on 7 and 21 DF,  p-value: 0.0003357

# Drop PCe2
meso.test.a.19 <- lm(rmesograzer.mass.perg.atl ~ 
  + rPC1.env.global.atl # + rPC2.env.global.atl 
  + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
  + rPC1.zos.atl + rPC2.zos.atl # + rPC1.zos.atl*rPC2.env.global.atl
  , data = ZEN_2014_site_means_49_Atlantic)
summary(meso.test.a.19) # RESULT: PCe3, PCz1, and FC2 significant as in best model, but not FC1
#                     Estimate Std. Error t value Pr(>|t|)   
# (Intercept)          0.37582    0.15109   2.487  0.02094 * 
# rPC1.env.global.atl -0.00814    0.11574  -0.070  0.94456   
# rPC3.env.global.atl  0.31522    0.11686   2.697  0.01315 * 
# rFC1.global.atl     -0.34077    0.18682  -1.824  0.08176 . 
# rFC2.global.atl      0.67054    0.18273   3.670  0.00135 **
# rPC1.zos.atl         0.36396    0.16392   2.220  0.03700 * 
# rPC2.zos.atl        -0.29588    0.15167  -1.951  0.06394 . 
# Residual standard error: 0.1318 on 22 degrees of freedom
# Multiple R-squared:  0.7009,	Adjusted R-squared:  0.6193 
# F-statistic: 8.591 on 6 and 22 DF,  p-value: 7.056e-05

# Drop FC1
meso.test.a.20 <- lm(rmesograzer.mass.perg.atl ~ 
  + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl # + rFC1.global.atl 
  + rFC2.global.atl + rPC1.zos.atl + rPC2.zos.atl + rPC1.zos.atl*rPC2.env.global.atl
  , data = ZEN_2014_site_means_49_Atlantic)
summary(meso.test.a.20) # RESULT: 
#                                  Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                       -0.1311     0.2518  -0.521   0.6081  
# rPC1.env.global.atl                0.1525     0.1187   1.286   0.2126  
# rPC2.env.global.atl                0.6016     0.3061   1.965   0.0627 .
# rPC3.env.global.atl                0.3358     0.1209   2.778   0.0113 *
# rFC2.global.atl                    0.2691     0.1635   1.646   0.1147  
# rPC1.zos.atl                       1.0395     0.3861   2.692   0.0136 *
# rPC2.zos.atl                      -0.1548     0.1660  -0.933   0.3615  
# rPC2.env.global.atl:rPC1.zos.atl  -0.8871     0.6552  -1.354   0.1901  
# Residual standard error: 0.1323 on 21 degrees of freedom
# Multiple R-squared:  0.7124,	Adjusted R-squared:  0.6165 
# F-statistic: 7.431 on 7 and 21 DF,  p-value: 0.0001529


AICc(meso.test.a.9, meso.test.a.18, meso.test.a.19, meso.test.a.20)
#                df      AICc
# meso.test.a.9  10 -22.23501
# meso.test.a.18  9 -14.51709
# meso.test.a.19  8 -20.03475
# meso.test.a.20  9 -16.90118
# Result: Model with both FC2 and PCe2 is best 

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(meso.test.a.9)
res = residuals(meso.test.a.9, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49_Atlantic$rPC1.env.global.atl,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC3.env.global.atl,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC1.global.atl,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC2.global.atl,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC1.zos.atl,res, xlab = "PCz1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC2.zos.atl,res, xlab = "PCz2 (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: OK. A bit right-skewed.


# RESULTS and INTERPRETATION: In the Atlantic, Mesograzers reach highest biomass in productive, 
# nutrient-rich estuaries (high PCe2 and 3), in the meadow form of eelgrass, and are strongly 
# associated with eelgrass genetic variation (both FC1 and 2). Best model explains 72% 
# of variation.

# Tests of robustness by dropping FC2, FC1, or PCe2: PCe3 and PCz1 are significant in 
# all models. PCe2 significant in all models except where FC1 is dropped (P = 0.06). 
# FC2 significant in all models except where FC1 is dropped. FC1 not stable to dropping
# other terms. 


###################################################################################
# GLM: MESOGRAZER MASS PER BOTTOM AREA (ATLANTIC) - MODEL USING SITE MEANS        #
###################################################################################

# Main effects only
meso.area.test.a.1 <- lm(rmesograzer.mass.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl 
                         + rPC1.zos.atl + rPC2.zos.atl  
                         , data = ZEN_2014_site_means_49_Atlantic)
summary(meso.area.test.a.1)

# Add interaction:  PCe1 x FC1
meso.area.test.a.2 <- lm(rmesograzer.mass.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl 
                         + rPC1.env.global.atl*rFC1.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC1
meso.area.test.a.3 <- lm(rmesograzer.mass.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl 
                         + rPC2.env.global.atl*rFC1.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC1
meso.area.test.a.4 <- lm(rmesograzer.mass.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl 
                         + rPC3.env.global.atl*rFC1.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe1 x FC2
meso.area.test.a.5 <- lm(rmesograzer.mass.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl 
                         + rPC1.env.global.atl*rFC2.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC2
meso.area.test.a.6 <- lm(rmesograzer.mass.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl 
                         + rPC2.env.global.atl*rFC2.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC2
meso.area.test.a.7 <- lm(rmesograzer.mass.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl 
                         + rPC3.env.global.atl*rFC2.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)


# Add interaction:  PCz1 x PCe1
meso.area.test.a.8 <- lm(rmesograzer.mass.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl 
                         + rPC1.zos.atl*rPC1.env.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x PCe2
meso.area.test.a.9 <- lm(rmesograzer.mass.area.atl ~ 
                           + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                         + rPC1.zos.atl + rPC2.zos.atl 
                         + rPC1.zos.atl*rPC2.env.global.atl # added
                         , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x PCe3
meso.area.test.a.10 <- lm(rmesograzer.mass.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl 
                          + rPC1.zos.atl*rPC3.env.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x FC1
meso.area.test.a.11 <- lm(rmesograzer.mass.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl 
                          + rPC1.zos.atl*rFC1.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz1 x FC2
meso.area.test.a.12 <- lm(rmesograzer.mass.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl 
                          + rPC1.zos.atl*rFC2.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x PCe1
meso.area.test.a.13 <- lm(rmesograzer.mass.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl 
                          + rPC2.zos.atl*rPC1.env.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x PCe2
meso.area.test.a.14 <- lm(rmesograzer.mass.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl 
                          + rPC2.zos.atl*rPC2.env.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x PCe3
meso.area.test.a.15 <- lm(rmesograzer.mass.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl 
                          + rPC2.zos.atl*rPC3.env.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x FC1
meso.area.test.a.16 <- lm(rmesograzer.mass.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl 
                          + rPC2.zos.atl*rFC1.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCz2 x FC2
meso.area.test.a.17 <- lm(rmesograzer.mass.area.atl ~ 
                            + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                          + rPC1.zos.atl + rPC2.zos.atl 
                          + rPC2.zos.atl*rFC2.global.atl # added
                          , data = ZEN_2014_site_means_49_Atlantic)



AICc(meso.area.test.a.1, meso.area.test.a.2, meso.area.test.a.3, meso.area.test.a.4, meso.area.test.a.5, 
    meso.area.test.a.6, meso.area.test.a.7, meso.area.test.a.8, meso.area.test.a.9, meso.area.test.a.10, meso.area.test.a.11, 
    meso.area.test.a.12, meso.area.test.a.13, meso.area.test.a.14, meso.area.test.a.15, meso.area.test.a.16, meso.area.test.a.17)
#                     df       AICc
# meso.area.test.a.1   9 -29.182613
# meso.area.test.a.2  10 -24.659125
# meso.area.test.a.3  10 -25.444133
# meso.area.test.a.4  10 -24.484656
# meso.area.test.a.5  10 -27.830932
# meso.area.test.a.6  10 -25.151876
# meso.area.test.a.7  10 -25.118061
# meso.area.test.a.8  10 -25.220946
# meso.area.test.a.9  10 -28.657720
# meso.area.test.a.10 10 -27.397875
# meso.area.test.a.11 10 -24.691606
# meso.area.test.a.12 10 -24.575077
# meso.area.test.a.13 10 -24.771305
# meso.area.test.a.14 10 -24.568900
# meso.area.test.a.15 10 -24.434130
# meso.area.test.a.16 10 -24.486408
# meso.area.test.a.17 10 -24.780227
# Model 1 is best but 9, 5, and 10 are close.

summary(meso.area.test.a.1)
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.34088    0.14097   2.418 0.024773 *  
# rPC1.env.global.atl  0.05528    0.09642   0.573 0.572502    
# rPC2.env.global.atl  0.23753    0.13908   1.708 0.102407    
# rPC3.env.global.atl  0.34178    0.09849   3.470 0.002288 ** 
# rFC1.global.atl     -0.36268    0.15560  -2.331 0.029807 *  
# rFC2.global.atl      0.44687    0.16337   2.735 0.012397 *  
# rPC1.zos.atl         0.31898    0.14135   2.257 0.034807 *  
# rPC2.zos.atl        -0.50811    0.12431  -4.088 0.000527 ***
# Residual standard error: 0.1071 on 21 degrees of freedom
# Multiple R-squared:  0.7839,	Adjusted R-squared:  0.7118 
# F-statistic: 10.88 on 7 and 21 DF,  p-value: 9.403e-06

# In all four top models, PCz1, PCz2, and PCe3 are significant.

# Test robustness of best model by dropping correlated predictors (PCe2, FC2)
# Drop FC2
meso.area.test.a.18 <- lm(rmesograzer.mass.area.atl ~ 
  + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl # + rFC2.global.atl
  + rPC1.zos.atl + rPC2.zos.atl + rPC1.zos.atl*rPC2.env.global.atl
  , data = ZEN_2014_site_means_49_Atlantic)
summary(meso.area.test.a.18)
#                                  Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                       0.05119    0.23224   0.220  0.82767   
# rPC1.env.global.atl               0.15191    0.10976   1.384  0.18090   
# rPC2.env.global.atl               0.63011    0.29313   2.150  0.04340 * 
# rPC3.env.global.atl               0.38932    0.11077   3.515  0.00206 **
# rFC1.global.atl                  -0.11467    0.14258  -0.804  0.43025   
# rPC1.zos.atl                      0.80485    0.35720   2.253  0.03506 * 
# rPC2.zos.atl                     -0.40757    0.15336  -2.658  0.01473 * 
# rPC2.env.global.atl:rPC1.zos.atl -0.54878    0.59982  -0.915  0.37063   
# Residual standard error: 0.1223 on 21 degrees of freedom
# Multiple R-squared:  0.7181,	Adjusted R-squared:  0.6241 
# F-statistic: 7.642 on 7 and 21 DF,  p-value: 0.0001261

# Drop PCe2
meso.area.test.a.19 <- lm(rmesograzer.mass.area.atl ~ 
  + rPC1.env.global.atl # + rPC2.env.global.atl 
  + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
  + rPC1.zos.atl + rPC2.zos.atl # + rPC1.zos.atl*rPC2.env.global.atl
  , data = ZEN_2014_site_means_49_Atlantic)
summary(meso.area.test.a.19)
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.45942    0.12794   3.591 0.001626 ** 
# rPC1.env.global.atl  0.01851    0.09800   0.189 0.851917    
# rPC3.env.global.atl  0.29675    0.09895   2.999 0.006609 ** 
# rFC1.global.atl     -0.30367    0.15819  -1.920 0.067957 .  
# rFC2.global.atl      0.56358    0.15472   3.643 0.001436 ** 
# rPC1.zos.atl         0.23780    0.13879   1.713 0.100703    
# rPC2.zos.atl        -0.53673    0.12843  -4.179 0.000389 ***
# Residual standard error: 0.1116 on 22 degrees of freedom
# Multiple R-squared:  0.7539,	Adjusted R-squared:  0.6867 
# F-statistic: 11.23 on 6 and 22 DF,  p-value: 9.404e-06

# Drop FC1
meso.area.test.a.20 <- lm(rmesograzer.mass.area.atl ~ 
  + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl # + rFC1.global.atl 
  + rFC2.global.atl + rPC1.zos.atl + rPC2.zos.atl + rPC1.zos.atl*rPC2.env.global.atl
  , data = ZEN_2014_site_means_49_Atlantic)
summary(meso.area.test.a.20)
#                                  Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                       0.07856    0.22164   0.354  0.72653   
# rPC1.env.global.atl               0.14355    0.10446   1.374  0.18385   
# rPC2.env.global.atl               0.42596    0.26946   1.581  0.12887   
# rPC3.env.global.atl               0.30541    0.10641   2.870  0.00916 **
# rFC2.global.atl                   0.24326    0.14395   1.690  0.10585   
# rPC1.zos.atl                      0.74839    0.33986   2.202  0.03897 * 
# rPC2.zos.atl                     -0.43141    0.14612  -2.952  0.00761 **
# rPC2.env.global.atl:rPC1.zos.atl -0.66613    0.57678  -1.155  0.26110   
# Residual standard error: 0.1165 on 21 degrees of freedom
# Multiple R-squared:  0.7442,	Adjusted R-squared:  0.6589 
# F-statistic: 8.728 on 7 and 21 DF,  p-value: 4.918e-05

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(meso.area.test.a.1)
res = residuals(meso.area.test.a.1, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49_Atlantic$rPC1.env.global.atl,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC3.env.global.atl,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC1.global.atl,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC2.global.atl,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC1.zos.atl,res, xlab = "PCz1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC2.zos.atl,res, xlab = "PCz2 (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: OK.  


# RESULTS and INTERPRETATION: In the Atlantic, mesograzers per area are greater in productive estuarine
# conditions (high PCe3, PCe2) as consistently seen in global analysis and when expressed per g 
# eelgrass. Mesograzer mass per area also responds to eelgrass form, being highest in meadows
# (high PCz1) with eelgrass biomass (low PCz2). Thus, there is clear influence of eelgrass form on
# system-level mesograzer biomass. 

# # Tests of robustness by dropping FC1, FC2, or PCe2: PCe3 and PCz2 are significant in 
# all models and PCz1 is significant in most. Genetic FC2 is significant in full model and when PCe2 
# is dropped but not when FC1 is dropped. PCe2 is significant only when FC2 is dropped from model. 
# FC1 is (weakly) significant only in full model. Thus, effects of FC1 and PCe2 seem weak and unstable.   
# Effects of FC2 and PCz2 are perhaps somewhat more stable.  


###################################################################################
# GLM: EELGRASS GROWTH FORM PC1 (PACIFIC) - MODEL USING SITE MEANS                #
###################################################################################

# Explore correlations among variables used in models
pairs.panels(ZEN_2014_site_means_Pacific[,c("Latitude", "PC1.env.global", "PC2.env.global", "PC3.env.global",
  "FC1", "FC2", "PC1.zos.site", "PC2.zos.site", 
  "log10.periphyton.mass.per.g.zostera.site", "log10.mesograzer.mass.per.g.plant.site")], 
  smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 8)

# Main effects only
pcz1.test.p.1 <- lm(rPC1.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac  
                    , data = ZEN_2014_site_means_Pacific)
summary(pcz1.test.p.1)

# Add interaction:  PCe1 x FC2
pcz1.test.p.2 <- lm(rPC1.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.env.global.pac*rFC2.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe2 x FC2
pcz1.test.p.3 <- lm(rPC1.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC2.env.global.pac*rFC2.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe3 x FC2
pcz1.test.p.4 <- lm(rPC1.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC3.env.global.pac*rFC2.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe1 x FC1
pcz1.test.p.5 <- lm(rPC1.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.env.global.pac*rFC1.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe2 x FC1
pcz1.test.p.6 <- lm(rPC1.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC2.env.global.pac*rFC1.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe3 x FC1
pcz1.test.p.7 <- lm(rPC1.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC3.env.global.pac*rFC1.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)


AICc(pcz1.test.p.1, pcz1.test.p.2, pcz1.test.p.3, pcz1.test.p.4, pcz1.test.p.5, pcz1.test.p.6, pcz1.test.p.7) 
#               df     AICc
# pcz1.test.p.1  7 15.73848
# pcz1.test.p.2  8 21.45282
# pcz1.test.p.3  8 21.17510
# pcz1.test.p.4  8 20.07942
# pcz1.test.p.5  8 17.90850
# pcz1.test.p.6  8 21.49303
# pcz1.test.p.7  8 21.49569
# RESULTS: Best model is 1 (no interactions). 


summary(pcz1.test.p.1)
#                     Estimate Std. Error t value Pr(>|t|)   
# (Intercept)           0.4057     0.4255   0.954  0.35647   
# rPC1.env.global.pac  -0.4481     0.2715  -1.650  0.12116   
# rPC2.env.global.pac  -0.2444     0.4382  -0.558  0.58581   
# rPC3.env.global.pac   0.4193     0.2376   1.765  0.09940 . 
# rFC1.global.pac      -0.3524     0.4776  -0.738  0.47281   
# rFC2.global.pac       0.7805     0.2557   3.052  0.00861 **  
# Residual standard error: 0.2392 on 14 degrees of freedom
# Multiple R-squared:  0.5048,	Adjusted R-squared:  0.3279 
# F-statistic: 2.854 on 5 and 14 DF,  p-value: 0.05568

# Test robustness of best model by dropping correlated predictors (PCe1, FC2)
# Drop PCe1
pcz1.test.p.8 <- lm(rPC1.zos.pac ~ # rPC1.env.global.pac 
  + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac  
  , data = ZEN_2014_site_means_Pacific)
summary(pcz1.test.p.8)
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)          0.12263    0.41113   0.298   0.7696  
# rPC2.env.global.pac  0.03859    0.42577   0.091   0.9290  
# rPC3.env.global.pac  0.23672    0.22202   1.066   0.3032  
# rFC1.global.pac     -0.08232    0.47372  -0.174   0.8644  
# rFC2.global.pac      0.46458    0.17899   2.596   0.0203 *
# Residual standard error: 0.2526 on 15 degrees of freedom
# Multiple R-squared:  0.4084,	Adjusted R-squared:  0.2507 
# F-statistic: 2.589 on 4 and 15 DF,  p-value: 0.07923

# Drop FC2
pcz1.test.p.9 <- lm(rPC1.zos.pac ~ rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac 
  + rFC1.global.pac # + rFC2.global.pac  
  , data = ZEN_2014_site_means_Pacific)
summary(pcz1.test.p.9)
#                     Estimate Std. Error t value Pr(>|t|)
# (Intercept)          0.13063    0.51845   0.252    0.804
# rPC1.env.global.pac  0.17244    0.22443   0.768    0.454
# rPC2.env.global.pac  0.14297    0.52291   0.273    0.788
# rPC3.env.global.pac  0.23305    0.28630   0.814    0.428
# rFC1.global.pac     -0.09868    0.58632  -0.168    0.869
# Residual standard error: 0.2982 on 15 degrees of freedom
# Multiple R-squared:  0.1752,	Adjusted R-squared:  -0.04473 
# F-statistic: 0.7966 on 4 and 15 DF,  p-value: 0.5457


# Fit best model with unstandardized data
pcz1.test.p.1.raw <- lm(PC1.zos.site ~ +PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2
  , data = ZEN_2014_site_means_Pacific)
summary(pcz1.test.p.1.raw)
#                  Estimate Std. Error t value Pr(>|t|)   
# (Intercept)    -3.1110782  2.3199478  -1.341  0.20127   
# PC1.env.global -0.5023905  0.3044489  -1.650  0.12116   
# PC2.env.global -0.3207520  0.5750617  -0.558  0.58581   
# PC3.env.global  0.9281119  0.5259200   1.765  0.09940 . 
# FC1            -0.0026872  0.0036422  -0.738  0.47281   
# FC2             0.0015080  0.0004941   3.052  0.00861 **
# Residual standard error: 1.281 on 14 degrees of freedom
# Multiple R-squared:  0.5048,	Adjusted R-squared:  0.3279 
# F-statistic: 2.854 on 5 and 14 DF,  p-value: 0.05568

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(pcz1.test.p.1)
res = residuals(pcz1.test.p.1, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_Pacific$rPC1.env.global.pac,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rPC2.env.global.pac,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rPC3.env.global.pac,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rFC1.global.pac,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rFC2.global.pac,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: Pretty Good. 

# RESULTS and INTERPRETATION: In the Pacific, meadow form (high PCz1) is not well explained
# by our models (best model P = 0.056). Only significant predictor is genetic FC2. 
# Best model explains only 33% of variation.

# Tests of robustness sshow that effect of FC2 remains significant (though overall model
# is not, P = 0.08) when correlated PCe1 is removed, so effect of FC2 is not an artifact
# of collinear predictors. 


###################################################################################
# GLM: EELGRASS GROWTH FORM PC2 (PACIFIC) - MODEL USING SITE MEANS                #
###################################################################################

# Main effects only
pcz2.test.p.1 <- lm(rPC2.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac  
                    , data = ZEN_2014_site_means_Pacific)
summary(pcz2.test.p.1)

# Add interaction:  PCe1 x FC2
pcz2.test.p.2 <- lm(rPC2.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.env.global.pac*rFC2.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe2 x FC2
pcz2.test.p.3 <- lm(rPC2.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC2.env.global.pac*rFC2.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe3 x FC2
pcz2.test.p.4 <- lm(rPC2.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC3.env.global.pac*rFC2.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe1 x FC1
pcz2.test.p.5 <- lm(rPC2.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.env.global.pac*rFC1.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe2 x FC1
pcz2.test.p.6 <- lm(rPC2.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC2.env.global.pac*rFC1.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe3 x FC1
pcz2.test.p.7 <- lm(rPC2.zos.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC3.env.global.pac*rFC1.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

AICc(pcz2.test.p.1, pcz2.test.p.2, pcz2.test.p.3, pcz2.test.p.4, pcz2.test.p.5, pcz2.test.p.6, pcz2.test.p.7) 
#               df     AICc
# pcz2.test.p.1  7 4.034721
# pcz2.test.p.2  8 7.701232
# pcz2.test.p.3  8 9.757615
# pcz2.test.p.4  8 5.037353
# pcz2.test.p.5  8 9.784571
# pcz2.test.p.6  8 9.738034
# pcz2.test.p.7  8 7.377881
# Model 1 is best but 4 (with no significant effects) is also close.

summary(pcz2.test.p.1)
#                     Estimate Std. Error t value Pr(>|t|)   
# (Intercept)          0.23005    0.31756   0.724  0.48073   
# rPC1.env.global.pac  0.03635    0.20266   0.179  0.86022   
# rPC2.env.global.pac  0.50320    0.32704   1.539  0.14618   
# rPC3.env.global.pac  0.16192    0.17733   0.913  0.37664   
# rFC1.global.pac      0.45236    0.35642   1.269  0.22507   
# rFC2.global.pac     -0.63352    0.19085  -3.319  0.00506 **
# Residual standard error: 0.1785 on 14 degrees of freedom
# Multiple R-squared:  0.6705,	Adjusted R-squared:  0.5528 
# F-statistic: 5.697 on 5 and 14 DF,  p-value: 0.004516


# Test robustness of best model by dropping correlated predictors (PCe1, FC2)
# Drop PCe1
pcz2.test.p.8 <- lm(rPC2.zos.pac ~ # rPC1.env.global.pac 
  + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac  
  , data = ZEN_2014_site_means_Pacific)
summary(pcz2.test.p.8)
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           0.2530     0.2811   0.900 0.382239    
# rPC2.env.global.pac   0.4802     0.2911   1.650 0.119751    
# rPC3.env.global.pac   0.1767     0.1518   1.164 0.262443    
# rFC1.global.pac       0.4305     0.3239   1.329 0.203665    
# rFC2.global.pac      -0.6079     0.1224  -4.968 0.000169 ***
# Residual standard error: 0.1727 on 15 degrees of freedom
# Multiple R-squared:  0.6697,	Adjusted R-squared:  0.5816 
# F-statistic: 7.603 on 4 and 15 DF,  p-value: 0.001485

# Drop FC2
pcz2.test.p.9 <- lm(rPC2.zos.pac ~ rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac 
  + rFC1.global.pac # + rFC2.global.pac  
  , data = ZEN_2014_site_means_Pacific)
summary(pcz2.test.p.9)
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)           0.4534     0.4008   1.131   0.2758  
# rPC1.env.global.pac  -0.4673     0.1735  -2.693   0.0167 *
# rPC2.env.global.pac   0.1888     0.4043   0.467   0.6472  
# rPC3.env.global.pac   0.3131     0.2213   1.415   0.1776  
# rFC1.global.pac       0.2465     0.4533   0.544   0.5946  
# Residual standard error: 0.2306 on 15 degrees of freedom
# Multiple R-squared:  0.4111,	Adjusted R-squared:  0.254 
# F-statistic: 2.618 on 4 and 15 DF,  p-value: 0.07699


# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(pcz2.test.p.1)
res = residuals(pcz2.test.p.1, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_Pacific$rPC1.env.global.pac,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rPC2.env.global.pac,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rPC3.env.global.pac,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rFC1.global.pac,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rFC2.global.pac,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: Residuals are strongly right-skewed. 

# RESULTS and INTERPRETATION: In the Pacific, eelgrass biomass is higher (low PCz2) at 
# sites with high genetic FC2. That's it. Best model explains 55% of variance (P = 0.005). 

# Tests of robustness show that effect of FC2 remains highly significant when correlated
# predictor PCe1 is removed from model. In contrast, PCe1 is only weakly significant
# (model P = 0.077) when FC2 is removed. 


###################################################################################
# GLM: PERIPHYTON PER G EELGRASS (PACIFIC) - MODEL USING SITE MEANS               #
###################################################################################

# Main effects only
peri.test.p.1 <- lm(rperiphyton.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
                    + rPC1.zos.pac + rPC2.zos.pac
                    , data = ZEN_2014_site_means_Pacific)
summary(peri.test.p.1)

# Add interaction:  PCe1 x FC1
peri.test.p.2 <- lm(rperiphyton.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac
                    + rPC1.env.global.pac*rFC1.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe2 x FC1
peri.test.p.3 <- lm(rperiphyton.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac
                    + rPC2.env.global.pac*rFC1.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe3 x FC1
peri.test.p.4 <- lm(rperiphyton.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac
                    + rPC3.env.global.pac*rFC1.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe1 x FC2
peri.test.p.5 <- lm(rperiphyton.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac
                    + rPC1.env.global.pac*rFC2.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe2 x FC2
peri.test.p.6 <- lm(rperiphyton.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac
                    + rPC2.env.global.pac*rFC2.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe3 x FC2
peri.test.p.7 <- lm(rperiphyton.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac
                    + rPC3.env.global.pac*rFC2.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)


# Add interaction:  PCz1 x PCe1
peri.test.p.8 <- lm(rperiphyton.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac
                    + rPC1.zos.pac*rPC1.env.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x PCe2
peri.test.p.9 <- lm(rperiphyton.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac
                    + rPC1.zos.pac*rPC2.env.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x PCe3
peri.test.p.10 <- lm(rperiphyton.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac
                     + rPC1.zos.pac*rPC3.env.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x FC1
peri.test.p.11 <- lm(rperiphyton.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac
                     + rPC1.zos.pac*rFC1.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x FC2
peri.test.p.12 <- lm(rperiphyton.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac
                     + rPC1.zos.pac*rFC2.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x PCe1
peri.test.p.13 <- lm(rperiphyton.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac
                     + rPC2.zos.pac*rPC1.env.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x PCe2
peri.test.p.14 <- lm(rperiphyton.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac
                     + rPC2.zos.pac*rPC2.env.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x PCe3
peri.test.p.15 <- lm(rperiphyton.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac
                     + rPC2.zos.pac*rPC3.env.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x FC1
peri.test.p.16 <- lm(rperiphyton.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac
                     + rPC2.zos.pac*rFC1.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x FC2
peri.test.p.17 <- lm(rperiphyton.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac
                     + rPC2.zos.pac*rFC2.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)


AICc(peri.test.p.1, peri.test.p.2, peri.test.p.3, peri.test.p.4, peri.test.p.5, 
    peri.test.p.6, peri.test.p.7, peri.test.p.8, peri.test.p.9, peri.test.p.10, peri.test.p.11, 
    peri.test.p.12, peri.test.p.13, peri.test.p.14, peri.test.p.15, peri.test.p.16, peri.test.p.17)

#                df      AICc
# peri.test.p.1   9 20.863276
# peri.test.p.2  10 29.166197
# peri.test.p.3  10 29.307720
# peri.test.p.4  10 29.278611
# peri.test.p.5  10 29.289996
# peri.test.p.6  10 28.232507
# peri.test.p.7  10 28.673314
# peri.test.p.8  10 27.719510
# peri.test.p.9  10 27.514226
# peri.test.p.10 10 29.249820
# peri.test.p.11 10 23.879051
# peri.test.p.12 10 28.426364
# peri.test.p.13 10 28.608565
# peri.test.p.14 10 29.161736
# peri.test.p.15 10 28.018926
# peri.test.p.16 10 28.511294
# peri.test.p.17 10 29.022788
# Best model by far is 1 (no interactions) 

summary(peri.test.p.1)
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)          0.20586    0.43652   0.472   0.6457  
# rPC1.env.global.pac -0.51410    0.27829  -1.847   0.0895 .
# rPC2.env.global.pac  0.89881    0.42953   2.093   0.0583 .
# rPC3.env.global.pac  0.04003    0.28640   0.140   0.8912  
# rFC1.global.pac      0.74178    0.45151   1.643   0.1263  
# rFC2.global.pac      0.77108    0.31289   2.464   0.0298 *
# rPC1.zos.pac        -0.51304    0.32606  -1.573   0.1416  
# rPC2.zos.pac        -0.34868    0.43688  -0.798   0.4403  
# Residual standard error: 0.214 on 12 degrees of freedom
# Multiple R-squared:  0.6362,	Adjusted R-squared:  0.424 
# F-statistic: 2.998 on 7 and 12 DF,  p-value: 0.04575


# Test robustness of best model by dropping correlated predictors (PCe1, FC2)
# Drop PCe1
peri.test.p.18 <- lm(rperiphyton.perg.pac ~ # rPC1.env.global.pac 
  + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
  + rPC1.zos.pac + rPC2.zos.pac
  , data = ZEN_2014_site_means_Pacific)
summary(peri.test.p.18)
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)          -0.2161     0.4050  -0.534   0.6026  
# rPC2.env.global.pac   1.0959     0.4530   2.419   0.0310 *
# rPC3.env.global.pac  -0.2816     0.2476  -1.137   0.2760  
# rFC1.global.pac       0.9715     0.4726   2.056   0.0605 .
# rFC2.global.pac       0.4197     0.2705   1.552   0.1448  
# rPC1.zos.pac         -0.2202     0.3102  -0.710   0.4904  
# rPC2.zos.pac         -0.1066     0.4538  -0.235   0.8179  
# Residual standard error: 0.233 on 13 degrees of freedom
# Multiple R-squared:  0.5328,	Adjusted R-squared:  0.3171 
# F-statistic: 2.471 on 6 and 13 DF,  p-value: 0.08093

# Drop PCe1
peri.test.p.19 <- lm(rperiphyton.perg.pac ~ rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac 
  + rFC1.global.pac # + rFC2.global.pac 
  + rPC1.zos.pac + rPC2.zos.pac
  , data = ZEN_2014_site_means_Pacific)
summary(peri.test.p.19)
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)          0.07363    0.51080   0.144   0.8876  
# rPC1.env.global.pac -0.09720    0.26054  -0.373   0.7151  
# rPC2.env.global.pac  1.32455    0.46369   2.857   0.0135 *
# rPC3.env.global.pac -0.07190    0.33342  -0.216   0.8326  
# rFC1.global.pac      1.09715    0.50449   2.175   0.0487 *
# rPC1.zos.pac        -0.34281    0.37572  -0.912   0.3782  
# rPC2.zos.pac        -0.70558    0.48600  -1.452   0.1703  
# Residual standard error: 0.2523 on 13 degrees of freedom
# Multiple R-squared:  0.4521,	Adjusted R-squared:  0.1993 
# F-statistic: 1.788 on 6 and 13 DF,  p-value: 0.1786


# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(peri.test.p.11)
res = residuals(peri.test.p.11, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_Pacific$rPC1.env.global.pac,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rPC2.env.global.pac,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rPC3.env.global.pac,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rFC1.global.pac,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rFC2.global.pac,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rPC1.zos.pac,res, xlab = "PCz1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_Pacific$rPC2.zos.pac,res, xlab = "PCz2 (scaled)", ylab = "residuals",) 
par(op)
# RESULTS: Generally OK.

# RESULTS and INTERPRETATION: In the Pacific, our models are not great for periphyton. 
# Best model (P = 0.045, R2  = 0.42) has significant association with genetic FC2, nearly
# significant increrase with PCe2. 

# Tests of robustness show that effect of PCe2 remains (weakly) significant in all models, 
# although overall model probabilities are poor. Thus, explanatory power is poor -- only 
# significant predictor in best model (FC2) is unstable to presence/absence of other predictors.  


###################################################################################
# GLM: PERIPHYTON PER BOTTOM AREA (PACIFIC) - MODEL USING SITE MEANS              #
###################################################################################

# Main effects only
peri.area.test.p.1 <- lm(rperiphyton.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
                         + rPC1.zos.pac + rPC2.zos.pac
                         , data = ZEN_2014_site_means_Pacific)
summary(peri.area.test.p.1)

# Add interaction:  PCe1 x FC1
peri.area.test.p.2 <- lm(rperiphyton.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac
                         + rPC1.env.global.pac*rFC1.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe2 x FC1
peri.area.test.p.3 <- lm(rperiphyton.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac
                         + rPC2.env.global.pac*rFC1.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe3 x FC1
peri.area.test.p.4 <- lm(rperiphyton.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac
                         + rPC3.env.global.pac*rFC1.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe1 x FC2
peri.area.test.p.5 <- lm(rperiphyton.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac
                         + rPC1.env.global.pac*rFC2.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe2 x FC2
peri.area.test.p.6 <- lm(rperiphyton.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac
                         + rPC2.env.global.pac*rFC2.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe3 x FC2
peri.area.test.p.7 <- lm(rperiphyton.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac
                         + rPC3.env.global.pac*rFC2.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)


# Add interaction:  PCz1 x PCe1
peri.area.test.p.8 <- lm(rperiphyton.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac
                         + rPC1.zos.pac*rPC1.env.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x PCe2
peri.area.test.p.9 <- lm(rperiphyton.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac
                         + rPC1.zos.pac*rPC2.env.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x PCe3
peri.area.test.p.10 <- lm(rperiphyton.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac
                          + rPC1.zos.pac*rPC3.env.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x FC1
peri.area.test.p.11 <- lm(rperiphyton.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac
                          + rPC1.zos.pac*rFC1.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x FC2
peri.area.test.p.12 <- lm(rperiphyton.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac
                          + rPC1.zos.pac*rFC2.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x PCe1
peri.area.test.p.13 <- lm(rperiphyton.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac
                          + rPC2.zos.pac*rPC1.env.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x PCe2
peri.area.test.p.14 <- lm(rperiphyton.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac
                          + rPC2.zos.pac*rPC2.env.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x PCe3
peri.area.test.p.15 <- lm(rperiphyton.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac
                          + rPC2.zos.pac*rPC3.env.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x FC1
peri.area.test.p.16 <- lm(rperiphyton.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac
                          + rPC2.zos.pac*rFC1.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x FC2
peri.area.test.p.17 <- lm(rperiphyton.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac
                          + rPC2.zos.pac*rFC2.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)


AICc(peri.area.test.p.1, peri.area.test.p.2, peri.area.test.p.3, peri.area.test.p.4, peri.area.test.p.5, 
    peri.area.test.p.6, peri.area.test.p.7, peri.area.test.p.8, peri.area.test.p.9, peri.area.test.p.10, peri.area.test.p.11, 
    peri.area.test.p.12, peri.area.test.p.13, peri.area.test.p.14, peri.area.test.p.15, peri.area.test.p.16, peri.area.test.p.17)

#                     df      AICc
# peri.area.test.p.1   9 17.719202
# peri.area.test.p.2  10 26.051381
# peri.area.test.p.3  10 26.162615
# peri.area.test.p.4  10 26.153675
# peri.area.test.p.5  10 26.147233
# peri.area.test.p.6  10 25.067536
# peri.area.test.p.7  10 25.189995
# peri.area.test.p.8  10 24.528569
# peri.area.test.p.9  10 24.067627
# peri.area.test.p.10 10 25.969967
# peri.area.test.p.11 10 20.380521
# peri.area.test.p.12 10 25.325278
# peri.area.test.p.13 10 25.353934
# peri.area.test.p.14 10 26.063247
# peri.area.test.p.15 10 25.261996
# peri.area.test.p.16 10 25.584376
# peri.area.test.p.17 10 25.881756
# Model 1 (no interactions) is best.

summary(peri.area.test.p.1)
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)          0.39889    0.40353   0.989   0.3424  
# rPC1.env.global.pac -0.31626    0.25725  -1.229   0.2425  
# rPC2.env.global.pac  0.80384    0.39706   2.024   0.0658 .
# rPC3.env.global.pac -0.04551    0.26475  -0.172   0.8664  
# rFC1.global.pac      0.68549    0.41738   1.642   0.1264  
# rFC2.global.pac      0.61102    0.28924   2.113   0.0563 .
# rPC1.zos.pac        -0.60995    0.30141  -2.024   0.0659 .
# rPC2.zos.pac        -0.52671    0.40386  -1.304   0.2166  
# Residual standard error: 0.1978 on 12 degrees of freedom
# Multiple R-squared:  0.6424,	Adjusted R-squared:  0.4338 
# F-statistic:  3.08 on 7 and 12 DF,  p-value: 0.04208

# Test robustness of best model by dropping correlated predictors (PCe1, FC2)
# Drop PCe1
peri.area.test.p.18 <- lm(rperiphyton.area.pac ~ # + rPC1.env.global.pac 
  + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
  + rPC1.zos.pac + rPC2.zos.pac
  , data = ZEN_2014_site_means_Pacific)
summary(peri.area.test.p.18)
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)           0.1393     0.3506   0.397   0.6976  
# rPC2.env.global.pac   0.9251     0.3921   2.359   0.0346 *
# rPC3.env.global.pac  -0.2433     0.2143  -1.135   0.2767  
# rFC1.global.pac       0.8268     0.4091   2.021   0.0643 .
# rFC2.global.pac       0.3949     0.2341   1.687   0.1155  
# rPC1.zos.pac         -0.4298     0.2685  -1.601   0.1335  
# rPC2.zos.pac         -0.3778     0.3928  -0.962   0.3537  
# Residual standard error: 0.2016 on 13 degrees of freedom
# Multiple R-squared:  0.5974,	Adjusted R-squared:  0.4115 
# F-statistic: 3.214 on 6 and 13 DF,  p-value: 0.03672

# Drop FC2
peri.area.test.p.19 <- lm(rperiphyton.area.pac ~ + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac 
  + rFC1.global.pac # + rFC2.global.pac 
  + rPC1.zos.pac + rPC2.zos.pac
  , data = ZEN_2014_site_means_Pacific)
summary(peri.area.test.p.19)
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)           0.2941     0.4507   0.653   0.5254  
# rPC1.env.global.pac   0.0141     0.2299   0.061   0.9520  
# rPC2.env.global.pac   1.1412     0.4091   2.790   0.0153 *
# rPC3.env.global.pac  -0.1342     0.2942  -0.456   0.6557  
# rFC1.global.pac       0.9671     0.4451   2.173   0.0489 *
# rPC1.zos.pac         -0.4751     0.3315  -1.433   0.1754  
# rPC2.zos.pac         -0.8095     0.4288  -1.888   0.0815 .
# Residual standard error: 0.2226 on 13 degrees of freedom
# Multiple R-squared:  0.5094,	Adjusted R-squared:  0.283 
# F-statistic:  2.25 on 6 and 13 DF,  p-value: 0.1039


# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(peri.area.test.p.1)
res = residuals(peri.area.test.p.1, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",)
plot(ypred,res, xlab = "predicted", ylab = "residuals",)
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "")
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_Pacific$rPC1.env.global.pac,res, xlab = "PCe1 (scaled)", ylab = "residuals",)
plot(ZEN_2014_site_means_Pacific$rPC2.env.global.pac,res, xlab = "PCe2 (scaled)", ylab = "residuals",)
plot(ZEN_2014_site_means_Pacific$rPC3.env.global.pac,res, xlab = "PCe3 (scaled)", ylab = "residuals",)
plot(ZEN_2014_site_means_Pacific$rFC1.global.pac,res, xlab = "FC1 (scaled)", ylab = "residuals",)
plot(ZEN_2014_site_means_Pacific$rFC2.global.pac,res, xlab = "FC2 (scaled)", ylab = "residuals",)
plot(ZEN_2014_site_means_Pacific$rPC1.zos.pac,res, xlab = "PCz1 (scaled)", ylab = "residuals",)
plot(ZEN_2014_site_means_Pacific$rPC2.zos.pac,res, xlab = "PCz2 (scaled)", ylab = "residuals",)
par(op)
# RESULTS:OK.

# RESULTS and INTERPRETATION: In Pacific, models have poor explanatopry power for 
# periphyton per area. Best model has P = 0.042 with marginally non-significant effects
# of PCe2, FC2, PCz1. 

# Tests of robustness show that effect of PCe2 becomes marginally significant after 
# dropping either PCe1 or FC2, althogh overall model P is weak or non-significant in both cases. 


###################################################################################
# GLM: MESOGRAZER MASS PER G EELGRASS (PACIFIC) - MODEL USING SITE MEANS          #
###################################################################################

# Main effects only
meso.test.p.1 <- lm(rmesograzer.mass.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
                    + rPC1.zos.pac + rPC2.zos.pac  
                    , data = ZEN_2014_site_means_Pacific)
summary(meso.test.p.1)

# Add interaction:  PCe1 x FC1
meso.test.p.2 <- lm(rmesograzer.mass.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac 
                    + rPC1.env.global.pac*rFC1.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe2 x FC1
meso.test.p.3 <- lm(rmesograzer.mass.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac 
                    + rPC2.env.global.pac*rFC1.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe3 x FC1
meso.test.p.4 <- lm(rmesograzer.mass.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac 
                    + rPC3.env.global.pac*rFC1.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe1 x FC2
meso.test.p.5 <- lm(rmesograzer.mass.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac 
                    + rPC1.env.global.pac*rFC2.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe2 x FC2
meso.test.p.6 <- lm(rmesograzer.mass.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac 
                    + rPC2.env.global.pac*rFC2.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe3 x FC2
meso.test.p.7 <- lm(rmesograzer.mass.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac 
                    + rPC3.env.global.pac*rFC2.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)


# Add interaction:  PCz1 x PCe1
meso.test.p.8 <- lm(rmesograzer.mass.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac 
                    + rPC1.zos.pac*rPC1.env.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x PCe2
meso.test.p.9 <- lm(rmesograzer.mass.perg.pac ~ 
                      + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                    + rPC1.zos.pac + rPC2.zos.pac 
                    + rPC1.zos.pac*rPC2.env.global.pac # added
                    , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x PCe3
meso.test.p.10 <- lm(rmesograzer.mass.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac 
                     + rPC1.zos.pac*rPC3.env.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x FC1
meso.test.p.11 <- lm(rmesograzer.mass.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac 
                     + rPC1.zos.pac*rFC1.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x FC2
meso.test.p.12 <- lm(rmesograzer.mass.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac 
                     + rPC1.zos.pac*rFC2.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x PCe1
meso.test.p.13 <- lm(rmesograzer.mass.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac 
                     + rPC2.zos.pac*rPC1.env.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x PCe2
meso.test.p.14 <- lm(rmesograzer.mass.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac 
                     + rPC2.zos.pac*rPC2.env.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x PCe3
meso.test.p.15 <- lm(rmesograzer.mass.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac 
                     + rPC2.zos.pac*rPC3.env.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x FC1
meso.test.p.16 <- lm(rmesograzer.mass.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac 
                     + rPC2.zos.pac*rFC1.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x FC2
meso.test.p.17 <- lm(rmesograzer.mass.perg.pac ~ 
                       + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                     + rPC1.zos.pac + rPC2.zos.pac 
                     + rPC2.zos.pac*rFC2.global.pac # added
                     , data = ZEN_2014_site_means_Pacific)

AICc(meso.test.p.1, meso.test.p.2, meso.test.p.3, meso.test.p.4, meso.test.p.5, 
    meso.test.p.6, meso.test.p.7, meso.test.p.8, meso.test.p.9, meso.test.p.10, meso.test.p.11, 
    meso.test.p.12, meso.test.p.13, meso.test.p.14, meso.test.p.15, meso.test.p.16, meso.test.p.17)
#                df      AICc
# meso.test.p.1   9 20.659551
# meso.test.p.2  10 27.441351
# meso.test.p.3  10 28.952643
# meso.test.p.4  10 29.096811
# meso.test.p.5  10 26.152577
# meso.test.p.6  10 28.215843
# meso.test.p.7  10 29.070943
# meso.test.p.8  10 28.521960
# meso.test.p.9  10 27.976852
# meso.test.p.10 10 29.055338
# meso.test.p.11 10 28.476004
# meso.test.p.12 10 28.375624
# meso.test.p.13 10 28.826613
# meso.test.p.14 10 28.462430
# meso.test.p.15 10 26.283175
# meso.test.p.16 10 28.616263
# meso.test.p.17 10 29.100026
# Model 1 (no interactions) is best.


summary(meso.test.p.1)
#                     Estimate Std. Error t value Pr(>|t|)
# (Intercept)           0.2221     0.4343   0.511   0.6183
# rPC1.env.global.pac  -0.2819     0.2769  -1.018   0.3286
# rPC2.env.global.pac   0.6111     0.4274   1.430   0.1782
# rPC3.env.global.pac   0.7274     0.2849   2.553   0.0253 *
# rFC1.global.pac       0.2292     0.4492   0.510   0.6192
# rFC2.global.pac       0.1080     0.3113   0.347   0.7348
# rPC1.zos.pac         -0.4312     0.3244  -1.329   0.2085
# rPC2.zos.pac         -0.5119     0.4347  -1.178   0.2618
# Residual standard error: 0.2129 on 12 degrees of freedom
# Multiple R-squared:  0.5434,	Adjusted R-squared:  0.2771
# F-statistic:  2.04 on 7 and 12 DF,  p-value: 0.1327
# Model has no significant explanatory power.

# Test robustness of best model by dropping correlated predictors (PCe1, FC2)
# Drop PCe1
meso.test.p.18 <- lm(rmesograzer.mass.perg.pac ~ # + rPC1.env.global.pac 
  + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
  + rPC1.zos.pac + rPC2.zos.pac  
  , data = ZEN_2014_site_means_Pacific)
summary(meso.test.p.18)
#                      Estimate Std. Error t value Pr(>|t|)  
# (Intercept)         -0.009272   0.370623  -0.025   0.9804  
# rPC2.env.global.pac  0.719204   0.414539   1.735   0.1064  
# rPC3.env.global.pac  0.551053   0.226583   2.432   0.0302 *
# rFC1.global.pac      0.355151   0.432461   0.821   0.4263  
# rFC2.global.pac     -0.084740   0.247530  -0.342   0.7376  
# rPC1.zos.pac        -0.270628   0.283878  -0.953   0.3578  
# rPC2.zos.pac        -0.379131   0.415240  -0.913   0.3778  
# Residual standard error: 0.2132 on 13 degrees of freedom
# Multiple R-squared:  0.504,	Adjusted R-squared:  0.275 
# F-statistic: 2.201 on 6 and 13 DF,  p-value: 0.1098

# Drop FC2
meso.test.p.19 <- lm(rmesograzer.mass.perg.pac ~ + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac 
  + rFC1.global.pac # + rFC2.global.pac 
  + rPC1.zos.pac + rPC2.zos.pac  
  , data = ZEN_2014_site_means_Pacific)
summary(meso.test.p.19)
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)           0.2036     0.4162   0.489   0.6328  
# rPC1.env.global.pac  -0.2236     0.2123  -1.053   0.3114  
# rPC2.env.global.pac   0.6707     0.3778   1.775   0.0992 .
# rPC3.env.global.pac   0.7117     0.2717   2.620   0.0212 *
# rFC1.global.pac       0.2789     0.4110   0.679   0.5093  
# rPC1.zos.pac         -0.4074     0.3061  -1.331   0.2061  
# rPC2.zos.pac         -0.5619     0.3960  -1.419   0.1794  
# Residual standard error: 0.2055 on 13 degrees of freedom
# Multiple R-squared:  0.5388,	Adjusted R-squared:  0.326 
# F-statistic: 2.532 on 6 and 13 DF,  p-value: 0.07563

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(meso.test.p.1)
res = residuals(meso.test.p.1, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",)
plot(ypred,res, xlab = "predicted", ylab = "residuals",)
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "")
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_Pacific$rPC1.env.global.pac,res, xlab = "PCe1 (scaled)", ylab = "residuals",)
plot(ZEN_2014_site_means_Pacific$rPC2.env.global.pac,res, xlab = "PCe2 (scaled)", ylab = "residuals",)
plot(ZEN_2014_site_means_Pacific$rPC3.env.global.pac,res, xlab = "PCe3 (scaled)", ylab = "residuals",)
plot(ZEN_2014_site_means_Pacific$rFC1.global.pac,res, xlab = "FC1 (scaled)", ylab = "residuals",)
plot(ZEN_2014_site_means_Pacific$rFC2.global.pac,res, xlab = "FC2 (scaled)", ylab = "residuals",)
plot(ZEN_2014_site_means_Pacific$rPC1.zos.pac,res, xlab = "PCz1 (scaled)", ylab = "residuals",)
plot(ZEN_2014_site_means_Pacific$rPC2.zos.pac,res, xlab = "PCz2 (scaled)", ylab = "residuals",)
par(op)
# RESULTS: Looks good generally. One point way off qq line.


# RESULTS and INTERPRETATION: In the Pacific, our models have little explanatory power 
# for mesograzers (best P = 0.13). In best model only PCe3 is significant predictor. 

# Tests of robustness show that PCe3 remains significant when potentially confounding
# variables are dropped from model. 


###################################################################################
# GLM: MESOGRAZER MASS PER BOTTOM AREA (PACIFIC) - MODEL USING SITE MEANS         #
###################################################################################

# Main effects only
meso.area.test.p.1 <- lm(rmesograzer.mass.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac 
                         + rPC1.zos.pac + rPC2.zos.pac  
                         , data = ZEN_2014_site_means_Pacific)
summary(meso.area.test.p.1)

# Add interaction:  PCe1 x FC1
meso.area.test.p.2 <- lm(rmesograzer.mass.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac 
                         + rPC1.env.global.pac*rFC1.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe2 x FC1
meso.area.test.p.3 <- lm(rmesograzer.mass.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac 
                         + rPC2.env.global.pac*rFC1.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe3 x FC1
meso.area.test.p.4 <- lm(rmesograzer.mass.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac 
                         + rPC3.env.global.pac*rFC1.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe1 x FC2
meso.area.test.p.5 <- lm(rmesograzer.mass.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac 
                         + rPC1.env.global.pac*rFC2.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe2 x FC2
meso.area.test.p.6 <- lm(rmesograzer.mass.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac 
                         + rPC2.env.global.pac*rFC2.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCe3 x FC2
meso.area.test.p.7 <- lm(rmesograzer.mass.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac 
                         + rPC3.env.global.pac*rFC2.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)


# Add interaction:  PCz1 x PCe1
meso.area.test.p.8 <- lm(rmesograzer.mass.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac 
                         + rPC1.zos.pac*rPC1.env.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x PCe2
meso.area.test.p.9 <- lm(rmesograzer.mass.area.pac ~ 
                           + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                         + rPC1.zos.pac + rPC2.zos.pac 
                         + rPC1.zos.pac*rPC2.env.global.pac # added
                         , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x PCe3
meso.area.test.p.10 <- lm(rmesograzer.mass.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac 
                          + rPC1.zos.pac*rPC3.env.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x FC1
meso.area.test.p.11 <- lm(rmesograzer.mass.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac 
                          + rPC1.zos.pac*rFC1.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz1 x FC2
meso.area.test.p.12 <- lm(rmesograzer.mass.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac 
                          + rPC1.zos.pac*rFC2.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x PCe1
meso.area.test.p.13 <- lm(rmesograzer.mass.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac 
                          + rPC2.zos.pac*rPC1.env.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x PCe2
meso.area.test.p.14 <- lm(rmesograzer.mass.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac 
                          + rPC2.zos.pac*rPC2.env.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x PCe3
meso.area.test.p.15 <- lm(rmesograzer.mass.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac 
                          + rPC2.zos.pac*rPC3.env.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x FC1
meso.area.test.p.16 <- lm(rmesograzer.mass.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac 
                          + rPC2.zos.pac*rFC1.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

# Add interaction:  PCz2 x FC2
meso.area.test.p.17 <- lm(rmesograzer.mass.area.pac ~ 
                            + rPC1.env.global.pac + rPC2.env.global.pac + rPC3.env.global.pac + rFC1.global.pac + rFC2.global.pac
                          + rPC1.zos.pac + rPC2.zos.pac 
                          + rPC2.zos.pac*rFC2.global.pac # added
                          , data = ZEN_2014_site_means_Pacific)

AICc(meso.area.test.p.1, meso.area.test.p.2, meso.area.test.p.3, meso.area.test.p.4, meso.area.test.p.5, 
    meso.area.test.p.6, meso.area.test.p.7, meso.area.test.p.8, meso.area.test.p.9, meso.area.test.p.10, meso.area.test.p.11, 
    meso.area.test.p.12, meso.area.test.p.13, meso.area.test.p.14, meso.area.test.p.15, meso.area.test.p.16, meso.area.test.p.17)
#                     df      AICc
# meso.area.test.p.1   9 21.938820
# meso.area.test.p.2  10 28.709725
# meso.area.test.p.3  10 30.193534
# meso.area.test.p.4  10 30.380800
# meso.area.test.p.5  10 27.700725
# meso.area.test.p.6  10 29.211367
# meso.area.test.p.7  10 30.245632
# meso.area.test.p.8  10 29.745203
# meso.area.test.p.9  10 28.710902
# meso.area.test.p.10 10 30.383251
# meso.area.test.p.11 10 29.216306
# meso.area.test.p.12 10 29.691991
# meso.area.test.p.13 10 30.056277
# meso.area.test.p.14 10 29.658989
# meso.area.test.p.15 10 27.290427
# meso.area.test.p.16 10 29.877321
# meso.area.test.p.17 10 30.382786
# Model 1 (no interactions) is best

summary(meso.area.test.p.1)
#                     Estimate Std. Error t value Pr(>|t|)  
# (Intercept)          0.33004    0.44842   0.736   0.4759  
# rPC1.env.global.pac -0.17858    0.28588  -0.625   0.5439  
# rPC2.env.global.pac  0.64683    0.44124   1.466   0.1684  
# rPC3.env.global.pac  0.66304    0.29420   2.254   0.0437 *
# rFC1.global.pac      0.27116    0.46382   0.585   0.5696  
# rFC2.global.pac      0.09865    0.32142   0.307   0.7642  
# rPC1.zos.pac        -0.60368    0.33494  -1.802   0.0967 .
# rPC2.zos.pac        -0.71243    0.44879  -1.587   0.1384  
# Residual standard error: 0.2198 on 12 degrees of freedom
# Multiple R-squared:  0.503,	Adjusted R-squared:  0.213 
# F-statistic: 1.735 on 7 and 12 DF,  p-value: 0.1917


# RESULTS and INTERPRETATION: In the Pacific, poor predictive power for mesograzers per area.
# Best model has P = 0.19. PCe3 is only significant predictor


###################################################################################
# GLM: FIT BEST MODELS USING LOCALLY MEASURED PREDICTORS                          #
###################################################################################

# Compare to models with main effects only 

# EELGRASS FORM (PCz1)

pcz1.site.g.1.raw <- lm(PC1.zos.site ~ Ocean
                    + day.length + Temperature.C + Salinity.ppt + log10.Leaf.PercN.site 
                    + FC1 + FC2  
                    , data = ZEN_2014_site_means_49)
summary(pcz1.site.g.1.raw)
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -1.9984582  2.0846681  -0.959  0.34336    
# OceanPacific          -1.8209600  1.3074463  -1.393  0.17120    
# Temperature.C          0.0105804  0.0334563   0.316  0.75342    
# Salinity.ppt          -0.0333630  0.0198576  -1.680  0.10054    
# log10.Leaf.PercN.site -0.6382489  1.3443790  -0.475  0.63748    
# day.length             0.2335368  0.0832223   2.806  0.00764 ** 
# FC1                    0.0001511  0.0009215   0.164  0.87053    
# FC2                    0.0011578  0.0002570   4.505 5.42e-05 ***
# Residual standard error: 1.026 on 41 degrees of freedom
# Multiple R-squared:  0.7254,	Adjusted R-squared:  0.6785 
# F-statistic: 15.47 on 7 and 41 DF,  p-value: 1.025e-09

# RESULT: Main effects model using PC axes (pcz1.site.g.1) had significant 
# effects of PCe1 (similar to day length) and FC2.  


# EELGRASS BIOMASS (PCz2)

pcz2.site.g.1.raw <- lm(PC2.zos.site ~ Ocean
                    + day.length + Temperature.C + Salinity.ppt + log10.Leaf.PercN.site 
                    + FC1 + FC2  
                    , data = ZEN_2014_site_means_49)
summary(pcz2.site.g.1.raw)
#                         Estimate Std. Error t value Pr(>|t|)   
# (Intercept)           -1.723e+00  2.063e+00  -0.835  0.40828   
# OceanPacific          -2.279e-03  1.294e+00  -0.002  0.99860   
# day.length             1.024e-01  8.235e-02   1.244  0.22072   
# Temperature.C         -4.957e-02  3.310e-02  -1.497  0.14195   
# Salinity.ppt           3.246e-02  1.965e-02   1.652  0.10619   
# log10.Leaf.PercN.site  8.126e-01  1.330e+00   0.611  0.54468   
# FC1                    6.497e-05  9.118e-04   0.071  0.94354   
# FC2                   -7.956e-04  2.543e-04  -3.129  0.00323 **
# Residual standard error: 1.016 on 41 degrees of freedom
# Multiple R-squared:  0.3106,	Adjusted R-squared:  0.1929 
# F-statistic: 2.639 on 7 and 41 DF,  p-value: 0.02393

# RESULT: Main effects model using PC axes (pcz2.site.g.1) had significant
# effects of PCe3 (correlated with salinity, chl) and FC2. Thus, FC2 
# significant with both PC axes and raw data.


# PERIPHYTON MASS PER UNIT AREA

peri.site.g.1.raw <- lm(log10.periphyton.mass.per.area.site ~ Ocean
                    + day.length + Temperature.C + Salinity.ppt + log10.Leaf.PercN.site 
                    + FC1 + FC2  
                    + PC1.zos.site + PC2.zos.site
                    , data = ZEN_2014_site_means_49)
summary(peri.site.g.1.raw)
#                         Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            1.0796526  1.2266968   0.880   0.3842  
# OceanPacific          -0.4091065  0.7745690  -0.528   0.6004  
# day.length            -0.0061691  0.0529413  -0.117   0.9078  
# Temperature.C         -0.0078469  0.0199489  -0.393   0.6962  
# Salinity.ppt           0.0004861  0.0123719   0.039   0.9689  
# log10.Leaf.PercN.site  0.9156965  0.7844805   1.167   0.2502  
# FC1                   -0.0002879  0.0005334  -0.540   0.5925  
# FC2                    0.0005006  0.0002027   2.470   0.0180 *
# PC1.zos.site          -0.1786137  0.0914108  -1.954   0.0579 .
# PC2.zos.site          -0.1471660  0.0923821  -1.593   0.1192  
# Residual standard error: 0.5939 on 39 degrees of freedom
# Multiple R-squared:  0.4297,	Adjusted R-squared:  0.2981 
# F-statistic: 3.265 on 9 and 39 DF,  p-value: 0.004743

# RESULT: Main effects model using PC axes (peri.area.site.g.1) has only
# FC2 significant, same as with the measured vaeriables here. 


# MESOGRAZER MASS PER UNIT AREA

meso.area.site.g.1.raw <- lm(log10.mesograzer.mass.per.area.site ~ Ocean
                             + day.length + Temperature.C + Salinity.ppt + log10.Leaf.PercN.site 
                             + FC1 + FC2  
                             + PC1.zos.site + PC2.zos.site
                             , data = ZEN_2014_site_means_49)
summary(meso.area.site.g.1.raw)
#                         Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            0.3651512  1.3972942   0.261   0.7952  
# OceanPacific           1.5058051  0.8822887   1.707   0.0958 .
# day.length             0.0431771  0.0603039   0.716   0.4783  
# Temperature.C          0.0549786  0.0227232   2.419   0.0203 *
# Salinity.ppt          -0.0228671  0.0140924  -1.623   0.1127  
# log10.Leaf.PercN.site  2.3481326  0.8935786   2.628   0.0122 *
# FC1                    0.0006146  0.0006075   1.012   0.3180  
# FC2                    0.0002070  0.0002308   0.897   0.3754  
# PC1.zos.site           0.0005910  0.1041234   0.006   0.9955  
# PC2.zos.site          -0.2039467  0.1052298  -1.938   0.0599 .
# Residual standard error: 0.6765 on 39 degrees of freedom
# Multiple R-squared:   0.48,	Adjusted R-squared:   0.36 
# F-statistic:     4 on 9 and 39 DF,  p-value: 0.001103

# RESULT: Main effects model using PC axes (meso.area.site.g.1) had significant 
# effects of PCe2 (nutrients) and PCe3 ( salinity, chl). Similar to the effect 
# of leaf N in thisthis model. 


###################################################################################
# CALCULATION OF EELGRASS GROWTH (LEAF EXTENSION) FROM 2011 DATA                  #
###################################################################################

# Use data from ZEN 2011 experiment to derive general equation predicting eelgrass 
# leaf extension from sheath length and sheath width (see Ruesink et al. 2018 Oikos). 

names(zmgrowth) 

# Streamline to only variables for use in modeling
zmgrowth_slim <- subset(zmgrowth, select = c(Site, Sum.extension.per.d, Sheath.total.length, Sheath.width, 
  Leaf.number, L1.old, L4.old))

# Remove ALL observations with NA for any variables. This is extreme -- come back and impute missing values later
zmgrowth_nona <- zmgrowth[complete.cases(zmgrowth_slim), ]
nrow(zmgrowth_nona) # [1] 832. Good. 

# Model total leaf extension per shoot (in cm) as function of sheath length and width, the two 
# variables in Jen's analysis of 2011 data for which we have data from all 50 ZEN 2014 sites. 
# Note that we need to fit as simple linear model (no random intercept) because we can't estimate 
# the intercept for sites not included in the 2011 data set. So this is a rough estimate ...

zen1_leaf_extension_model <- lm(Sum.extension.per.d ~ Sheath.total.length 
            + Sheath.width 
            # + Leaf.number + L1.old + L4.old, 
            , data = zmgrowth_nona)
summary(zen1_leaf_extension_model)
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         -5.24955    1.57637   -3.33 0.000906 ***
# Sheath.total.length  0.14391    0.01163   12.38  < 2e-16 ***
# Sheath.width         4.87641    0.47202   10.33  < 2e-16 ***
# Residual standard error: 19.92 on 829 degrees of freedom
# Multiple R-squared:  0.6328,	Adjusted R-squared:  0.6319 
# F-statistic: 714.4 on 2 and 829 DF,  p-value: < 2.2e-16

# Use zen1_leaf_extension_model to estimate leaf extension per (summer) day in each plot for ZEN 2014.
# NOTE: I have left out the (negative) intercept in the model since it results in the great majority of
# leaf extension rate in 20914 being negativbe numbers!
ZEN_2014_plot$leaf.extension.per.day <- 0.14391*ZEN_2014_plot$Zostera.sheath.length + 4.87641*ZEN_2014_plot$Zostera.sheath.width
ZEN_2014_plot$log10.leaf.extension.per.day <- log10(ZEN_2014_plot$leaf.extension.per.day)
hist(ZEN_2014_plot$leaf.extension.per.day) 
hist(ZEN_2014_plot$log10.leaf.extension.per.day) 

# Obtain estimated mean leaf extension rates by site
temp <- ddply(ZEN_2014_plot, c("Site"), summarize, 
              leaf.extension.per.day.site = mean(leaf.extension.per.day, na.rm = T),
              log10.leaf.extension.per.day.site = mean(log10.leaf.extension.per.day, na.rm = T)
)

# Add leaf extension rate back to site-level data ets for modeling
ZEN_2014_site_means_49$leaf.extension.per.day.site <- temp$leaf.extension.per.day.site[match(ZEN_2014_site_means_49$Site, temp$Site)]
ZEN_2014_site_means_49$log10.leaf.extension.per.day.site <- temp$log10.leaf.extension.per.day.site[match(ZEN_2014_site_means_49$Site, temp$Site)]

ZEN_2014_site_means_49_Atlantic$leaf.extension.per.day.site <- temp$leaf.extension.per.day.site[match(ZEN_2014_site_means_49_Atlantic$Site, temp$Site)]
ZEN_2014_site_means_49_Atlantic$log10.leaf.extension.per.day.site <- temp$log10.leaf.extension.per.day.site[match(ZEN_2014_site_means_49_Atlantic$Site, temp$Site)]

ZEN_2014_site_means_Pacific$leaf.extension.per.day.site <- temp$leaf.extension.per.day.site[match(ZEN_2014_site_means_Pacific$Site, temp$Site)]
ZEN_2014_site_means_Pacific$log10.leaf.extension.per.day.site <- temp$log10.leaf.extension.per.day.site[match(ZEN_2014_site_means_Pacific$Site, temp$Site)]

hist(ZEN_2014_site_means_49$leaf.extension.per.day.site) # highly non-normal. 
hist(ZEN_2014_site_means_49$log10.leaf.extension.per.day.site) # Better. Bimodal - Atl vs Pac?


###################################################################################
# GLM: EELGRASS PRODUCTIVITY (GLOBAL) - MODEL USING SITE MEANS                    #
###################################################################################

# productivity  = daily leaf extension rate (for now ...)

# standardize and center by the range of observed values. # The '...' allows it to work with NAs.
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

ZEN_2014_site_means_49$rlog10.leaf.extension.per.day.site <- range01(ZEN_2014_site_means_49$log10.leaf.extension.per.day.site)
ZEN_2014_site_means_49_Atlantic$rlog10.leaf.extension.per.day.site <- range01(ZEN_2014_site_means_49_Atlantic$log10.leaf.extension.per.day.site)
ZEN_2014_site_means_Pacific$rlog10.leaf.extension.per.day.site <- range01(ZEN_2014_site_means_Pacific$log10.leaf.extension.per.day.site)


# Main effects only
leafext.site.g.1 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                       + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1  
                       , data = ZEN_2014_site_means_49)
summary(leafext.site.g.1)

# Add interaction: ocean x PCe1
leafext.site.g.2 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                       + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                       + Ocean*rPC1.env.global # added
                       , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe2
leafext.site.g.3 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                       + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                       + Ocean*rPC2.env.global # added
                       , data = ZEN_2014_site_means_49)

# Add interaction: ocean x PCe3
leafext.site.g.4 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                       + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                       + Ocean*rPC3.env.global # added
                       , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC1
leafext.site.g.5 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                       + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                       + Ocean*rFC1 # added
                       , data = ZEN_2014_site_means_49)

# Add interaction: ocean x FC2
leafext.site.g.6 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                       + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                       + Ocean*rFC2 # added
                       , data = ZEN_2014_site_means_49)

# Add interaction:  PCe1 x FC2
leafext.site.g.7 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                       + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                       + rPC1.env.global*rFC2 # added
                       , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC2
leafext.site.g.8 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                       + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                       + rPC2.env.global*rFC2 # added
                       , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC2
leafext.site.g.9 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                       + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                       + rPC3.env.global*rFC2 # added
                       , data = ZEN_2014_site_means_49)

# Add interaction:  PCe1 x FC1
leafext.site.g.10 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                        + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                        + rPC1.env.global*rFC1 # added
                        , data = ZEN_2014_site_means_49)

# Add interaction:  PCe2 x FC1
leafext.site.g.11 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                        + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                        + rPC2.env.global*rFC1 # added
                        , data = ZEN_2014_site_means_49)

# Add interaction:  PCe3 x FC1
leafext.site.g.12 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                        + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                        + rPC3.env.global*rFC1 # added
                        , data = ZEN_2014_site_means_49)


AICc(leafext.site.g.1, leafext.site.g.2, leafext.site.g.3, leafext.site.g.4, leafext.site.g.5, leafext.site.g.6, leafext.site.g.7, leafext.site.g.8, leafext.site.g.9, leafext.site.g.10, leafext.site.g.11, leafext.site.g.12) 
#                   df      AICc
# leafext.site.g.1   8 -47.33987
# leafext.site.g.2   9 -46.21579
# leafext.site.g.3   9 -44.64493
# leafext.site.g.4   9 -55.98339
# leafext.site.g.5   9 -44.32775
# leafext.site.g.6   9 -47.61754
# leafext.site.g.7   9 -44.33958
# leafext.site.g.8   9 -44.45083
# leafext.site.g.9   9 -45.84038
# leafext.site.g.10  9 -45.03288
# leafext.site.g.11  9 -45.15334
# leafext.site.g.12  9 -56.25251
# RESULTS: Best models are 12, 4. 

# Combine models 4 + 12
leafext.site.g.13 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                        + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC2 + rFC1
                        + Ocean*rPC3.env.global # added
                        + rPC3.env.global*rFC1 # added
                        , data = ZEN_2014_site_means_49)

AICc(leafext.site.g.4, leafext.site.g.12, leafext.site.g.13)
# No good. Average models 4 and 12

# Obtain model-averaged, standardized coefficients from models < 2 AIC from best 
summary(model.avg(leafext.site.g.4, leafext.site.g.12))
# Result: averaging simply adds two redundant interactions to model. They are both
# "significant" when conditionally averaged, and neither is significant in full
# average. Duh. Not helpful.

# Fit best model (4) without FC1 
leafext.site.g.14 <- lm(rlog10.leaf.extension.per.day.site ~ Ocean
                        + rPC1.env.global + rPC2.env.global + rPC3.env.global # + rFC1 
                        + rFC2
                        + Ocean*rPC3.env.global # added
                        , data = ZEN_2014_site_means_49)


# Fit best model (12) without ocean 
leafext.site.g.15 <- lm(rlog10.leaf.extension.per.day.site ~ # Ocean
                          + rPC1.env.global + rPC2.env.global + rPC3.env.global + rFC1 + rFC2
                        + rPC3.env.global*rFC1 # added
                        , data = ZEN_2014_site_means_49)

AICc(leafext.site.g.4, leafext.site.g.12, leafext.site.g.14, leafext.site.g.15)
#                   df     AICc
# leafext.site.g.4   9 198.8717
# leafext.site.g.12  9 199.1198
# leafext.site.g.14  8 198.0293
# leafext.site.g.15  8 206.4216
# Excluding ocean is much worse than other models.

summary(leafext.site.g.12)
#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           0.42350    0.21464   1.973  0.05526 .  
# OceanPacific          0.48513    0.15177   3.196  0.00268 ** 
# rPC1.env.global       0.30907    0.07168   4.312 9.92e-05 ***
# rPC2.env.global       0.11030    0.09230   1.195  0.23893    
# rPC3.env.global      -0.75428    0.24784  -3.043  0.00407 ** 
# rFC2                 -0.51423    0.09139  -5.627 1.47e-06 ***
# rFC1                 -0.20199    0.24555  -0.823  0.41549    
# rPC3.env.global:rFC1  1.01068    0.30066   3.362  0.00169 ** 
# Residual standard error: 0.1183 on 41 degrees of freedom
# Multiple R-squared:  0.7666,	Adjusted R-squared:  0.7268 
# F-statistic: 19.24 on 7 and 41 DF,  p-value: 4.158e-11

summary(leafext.site.g.4)
#                              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  -0.02082    0.20365  -0.102 0.919076    
# OceanPacific                  0.88078    0.20561   4.284 0.000108 ***
# rPC1.env.global               0.30836    0.07187   4.291 0.000106 ***
# rPC2.env.global               0.10979    0.09255   1.186 0.242339    
# rPC3.env.global               0.13150    0.08760   1.501 0.141006    
# rFC2                         -0.51471    0.09164  -5.617 1.52e-06 ***
# rFC1                          0.31117    0.21608   1.440 0.157427    
# OceanPacific:rPC3.env.global -0.78319    0.23599  -3.319 0.001904 ** 
# Residual standard error: 0.1186 on 41 degrees of freedom
# Multiple R-squared:  0.7653,	Adjusted R-squared:  0.7253 
# F-statistic:  19.1 on 7 and 41 DF,  p-value: 4.635e-11


# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(leafext.site.g.12)
res = residuals(leafext.site.g.12, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49$rPC1.env.global,res, xlab = "PCe1 (scaled)", ylab = "residuals",) # diagonal row - zeros?
plot(ZEN_2014_site_means_49$rPC2.env.global,res, xlab = "PCe2 (scaled)", ylab = "residuals",) # diagonal row - zeros?
plot(ZEN_2014_site_means_49$rPC3.env.global,res, xlab = "PCe3 (scaled)", ylab = "residuals",) # diagonal row - zeros?
plot(ZEN_2014_site_means_49$rFC2,res, xlab = "PCz1 (scaled)", ylab = "residuals",) # diagonal row - zeros?
par(op)
# RESULTS: Good.

# RESULTS: Generally, results for leaf extension rate are almost identical to those 
# for eelgrass form (PCz1) because extension rate is strongly correlated with 
# sheath length, which is in turn highly correlated with canopy height, which is 
# the dominant component of PCz1. Strong effects of FC2 and differences by ocean. 

# Explore correlations 
pairs.panels(ZEN_2014_site_means_49[,c("log10.leaf.extension.per.day.site", "log10.Zostera.sheath.length.site", 
  "log10.Zostera.longest.leaf.length.cm.site", "log10.Zostera.AG.mass.site", "rPC1.zos.site" 
  )], smooth=T,density=F,ellipses=F,lm=F,digits=2,scale=F, cex.cor = 4)


###################################################################################
# GLM: EELGRASS PRODUCTIVITY (ATLANTIC) - MODEL USING SITE MEANS                  #
###################################################################################

# Main effects only
leafext.site.a.1 <- lm(log10.leaf.extension.per.day.site ~ 
                         + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl  
                       , data = ZEN_2014_site_means_49_Atlantic)
summary(leafext.site.a.1)

# Add interaction:  PCe1 x FC1
leafext.site.a.2 <- lm(log10.leaf.extension.per.day.site ~ 
                         + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                       + rPC1.env.global.atl*rFC1.global.atl # added
                       , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC1
leafext.site.a.3 <- lm(log10.leaf.extension.per.day.site ~ 
                         + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                       + rPC2.env.global.atl*rFC1.global.atl # added
                       , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC1
leafext.site.a.4 <- lm(log10.leaf.extension.per.day.site ~ 
                         + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                       + rPC3.env.global.atl*rFC1.global.atl # added
                       , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe1 x FC2
leafext.site.a.5 <- lm(log10.leaf.extension.per.day.site ~ 
                         + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                       + rPC1.env.global.atl*rFC2.global.atl # added
                       , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe2 x FC2
leafext.site.a.6 <- lm(log10.leaf.extension.per.day.site ~ 
                         + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                       + rPC2.env.global.atl*rFC2.global.atl # added
                       , data = ZEN_2014_site_means_49_Atlantic)

# Add interaction:  PCe3 x FC2
leafext.site.a.7 <- lm(log10.leaf.extension.per.day.site ~ 
                         + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl
                       + rPC3.env.global.atl*rFC2.global.atl # added
                       , data = ZEN_2014_site_means_49_Atlantic)


AICc(leafext.site.a.1, leafext.site.a.2, leafext.site.a.3, leafext.site.a.4, leafext.site.a.5, leafext.site.a.6, leafext.site.a.7) 
#                  df      AICc
# leafext.site.a.1  7 -45.14543
# leafext.site.a.2  8 -41.47838
# leafext.site.a.3  8 -41.96211
# leafext.site.a.4  8 -42.58317
# leafext.site.a.5  8 -41.87834
# leafext.site.a.6  8 -41.58105
# leafext.site.a.7  8 -42.34592
# RESULTS: Best model is 1 (no interactions),

summary(leafext.site.a.1)
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.21441    0.06682   3.209  0.00390 ** 
# rPC1.env.global.atl  0.31696    0.06166   5.140  3.3e-05 ***
# rPC2.env.global.atl  0.23862    0.10597   2.252  0.03418 *  
# rPC3.env.global.atl  0.13488    0.06830   1.975  0.06041 .  
# rFC1.global.atl      0.20255    0.12553   1.614  0.12026    
# rFC2.global.atl     -0.36888    0.11382  -3.241  0.00361 ** 
# Residual standard error: 0.08939 on 23 degrees of freedom
# Multiple R-squared:  0.6752,	Adjusted R-squared:  0.6046 
# F-statistic: 9.564 on 5 and 23 DF,  p-value: 4.858e-05

# Test robustness of best model (1) by dropping correlated predictors (PCe2, FC2)
# Drop FC2
leafext.site.a.8 <- lm(log10.leaf.extension.per.day.site ~ 
                         + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl + rFC1.global.atl # + rFC2.global.atl  
                       , data = ZEN_2014_site_means_49_Atlantic)
summary(leafext.site.a.8) # RESULT: PCe1 remainssignificant

# Drop PCe2
leafext.site.a.9 <- lm(log10.leaf.extension.per.day.site ~ 
                         + rPC1.env.global.atl # + rPC2.env.global.atl 
                       + rPC3.env.global.atl + rFC1.global.atl + rFC2.global.atl  
                       , data = ZEN_2014_site_means_49_Atlantic)
summary(leafext.site.a.9) # RESULT: PCe1 and FC2 remian significant

# Drop FC1
leafext.site.a.10 <- lm(log10.leaf.extension.per.day.site ~ 
                          + rPC1.env.global.atl + rPC2.env.global.atl + rPC3.env.global.atl # + rFC1.global.atl 
                        + rFC2.global.atl  
                        , data = ZEN_2014_site_means_49_Atlantic)
summary(leafext.site.a.10) # 

# Fit best model with unstandardized data to allow direct comparison betwen oceans
leafext.site.a.1.raw <- lm(log10.leaf.extension.per.day.site ~ 
                             + PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2
                           , data = ZEN_2014_site_means_49_Atlantic)
summary(leafext.site.a.1.raw)

# Examine residuals
op <- par(mfrow = c(2,4))
ypred = predict(leafext.site.a.1)
res = residuals(leafext.site.a.1, type = 'pearson')
hist(res, xlab = "residuals", ylab = "frequency",) 
plot(ypred,res, xlab = "predicted", ylab = "residuals",) 
qqnorm(res, xlab = "Model Quantiles", ylab = "Observation Quantiles", main = "") 
qqline(res, col = "blue", lwd = 2) # strong heavy tails
plot(ZEN_2014_site_means_49_Atlantic$rPC1.env.global.atl,res, xlab = "PCe1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC2.env.global.atl,res, xlab = "PCe2 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rPC3.env.global.atl,res, xlab = "PCe3 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC1.global.atl,res, xlab = "FC1 (scaled)", ylab = "residuals",) 
plot(ZEN_2014_site_means_49_Atlantic$rFC2.global.atl,res, xlab = "FC2 (scaled)", ylab = "residuals",) 
par(op)
# RESULTS:  Good. 

# RESULTS and INTERPRETATION: In the Atlantic, leaf extension rate decreases with
# latitude (increasing PCz1) and is strongly associated with genetic FC2. 


###################################################################################
# GLM: EELGRASS PRODUCTIVITY (PACIFIC) - MODEL USING SITE MEANS                   #
###################################################################################




###################################################################################
# RANDOM FOREST ANALYSIS                                                          #
###################################################################################

# Use random forest analysis as a check on robustness of the results from GLMs. Random forest
# analysis uses machine learning algorithms to rank the importance of predictor variables 
# in a regression problem but sidesteps for the overfitting common in high-dimensional 
# regression-type models, handles all types of variables, and is robust to different 
# forms of relationships among them. 

# First, test suggestion from co-authors that eelgrass canopy height may be influenced by 
# exposure (fetch) and by light levels, indexed here by C:N ratio. 

library(randomForest)

ZEN_2014_site_means$Ocean <- as.factor(ZEN_2014_site_means$Ocean)

# EELGRASS CANOPY HEIGHT as function of local variables:

canopy.rf = randomForest(Zostera.longest.leaf.length.site  ~  Temperature.C + Salinity.ppt 
  + day.length + leaf.CN.ratio.site + log10.mean.fetch + FC1 + FC2,
  na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = ZEN_2014_site_means)

# Examine summary output
canopy.rf # % Var explained: 45.12

# Plot error as a function of # of trees
plot(canopy.rf) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(canopy.rf) # NOTE: saving this as an object returns the numbers in the graph

# RESULT: FC2 and FC2 are substantially better repdictros of canopy height than fetch or leaf C:N. 



# EELGRASS GROWTH FORM PCz1 (FOREST-MEADOW AXIS) as function of local variables:

PC1.zos.local.rf = randomForest(PC1.zos.site  ~  Temperature.C + Salinity.ppt + day.length
  + leaf.CN.ratio.site + log10.mean.fetch + FC1 + FC2,
  na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = ZEN_2014_site_means)

# Examine summary output
PC1.zos.local.rf # % Var explained: 55.59

# Plot error as a function of # of trees
plot(PC1.zos.local.rf) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(PC1.zos.local.rf) # NOTE: saving this as an object returns the numbers in the graph

# RESULT: 
# Interpretation: 


# EELGRASS GROWTH FORM PCz1 (FOREST-MEADOW AXIS): GLOBAL

PC1.zos.rf = randomForest(PC1.zos.site  ~  Ocean
                          + PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2,
                          na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = ZEN_2014_site_means)

# Examine summary output
PC1.zos.rf # % Var explained: 59.43

# Plot error as a function of # of trees
plot(PC1.zos.rf) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(PC1.zos.rf) # NOTE: saving this as an object returns the numbers in the graph

# RESULT: FC2 > FC 1 > Ocean > PC1.env.all
# Interpretation: genetics is strongest predictor of eelgrass forest/meadow growth form GLOBALLY


# EELGRASS GROWTH FORM PCz1 (FOREST-MEADOW AXIS): ATLANTIC

ZEN_2014_site_means_49_Atlantic <- droplevels(subset(ZEN_2014_site_means_49, Ocean == "Atlantic"))
ZEN_2014_site_means_49_Atlantic$Ocean <- as.factor(ZEN_2014_site_means_49_Atlantic$Ocean)

PC1.zos.a.rf = randomForest(PC1.zos.site  ~  
                              + PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2 # global predictors
                            , 
                            na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = ZEN_2014_site_means_49_Atlantic)

# Examine summary output
PC1.zos.a.rf # % Var explained: 28.96

# Plot variable importance
varImpPlot(PC1.zos.a.rf) # 

# RESULT: PCe1 >> PCe3 = FC2 = PCe2 > FC1


# EELGRASS GROWTH FORM PCz1 (FOREST-MEADOW AXIS): PACIFIC

ZEN_2014_site_means_Pacific <- droplevels(subset(ZEN_2014_site_means_49, Ocean == "Pacific"))
ZEN_2014_site_means_Pacific$Ocean <- as.factor(ZEN_2014_site_means_Pacific$Ocean)

PC1.zos.p.rf = randomForest(PC1.zos.site  ~  
                              + PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2 # global predictors
                            ,
                            na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = ZEN_2014_site_means_Pacific)

# Examine summary output
PC1.zos.p.rf # % Var explained: 7.71. 

# Plot error as a function of # of trees
plot(PC1.zos.p.rf) # Good. 200 trees is plenty.

# Plot variable importance
varImpPlot(PC1.zos.p.rf) # 

# RESULT: Little resolution and litle explanatory power (7% of variance). This is attempting to 
# divide a set of 20 site means, so very low power. Probably should restrict the random forest
# analysis to global data set where we need to use site means.  


# EELGRASS GROWTH FORM PCz2 (BIOMASS): GLOBAL

ZEN_2014_site_means$Ocean <- as.factor(ZEN_2014_site_means$Ocean)
PC2.zos.rf = randomForest(PC2.zos.site  ~  Ocean
                          + PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2,
                          na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = ZEN_2014_site_means)

# Examine summary output
PC2.zos.rf # % Var explained: 13.99

# Plot error as a function of # of trees
plot(PC2.zos.rf) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(PC2.zos.rf) # 

# RESULT: PCe3 > FC2 > PCe1 > FC1 = PCe2 >> Ocean



# EELGRASS GROWTH FORM PCz2 (BIOMASS): ATLANTIC

PC2.zos.a.rf = randomForest(PC2.zos.site  ~  
                              + PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2,
                            , 
                            na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = ZEN_2014_site_means_Atlantic)

# Examine summary output
PC2.zos.a.rf # % Var explained: 10.92

# Plot error as a function of # of trees
plot(PC2.zos.a.rf) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(PC2.zos.a.rf) # 

# RESULT: PCe3 > FC2 > PCe1 = PCe2 = FC1 (not much resolution)


# EELGRASS GROWTH FORM PCz2 (BIOMASS): PACIFIC

PC2.zos.p.rf = randomForest(PC2.zos.site  ~  
                              + PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2
                            , 
                            na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = ZEN_2014_site_means_Pacific)

# Examine summary output
PC2.zos.p.rf # % Var explained: 20.27

# Plot error as a function of # of trees
plot(PC2.zos.p.rf) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(PC2.zos.p.rf) # 

# RESULT: FC2 = PCe1 >> PCe2 = PCe3 = FC1


# MESOGRAZER BIOMASS: GLOBAL

meso.rf = randomForest(log10.mesograzer.mass.per.area.site  ~  Ocean
                       + PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2
                       + PC1.zos.site + PC2.zos.site + log10.periphyton.mass.per.g.zostera.site
                       ,
                       na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = ZEN_2014_site_means)

# Examine summary output
meso.rf # % Var explained: 31.68

# Plot error as a function of # of trees
plot(meso.rf) # Good. 200 trees is plenty.

# Plot variable importance
varImpPlot(meso.rf) # 

# RESULT: PCe2 = PCe3 >> PC2.zos.site > FC1 = ...
# Interpretation: Nutrient status (PC2.env.all) is most important driver globally. This 
# flips in various models between PCe2 (nutrients per se) and PCe3 (productive estuarine 
# conditions). 


# MESOGRAZER BIOMASS: ATLANTIC

meso.a.rf = randomForest(log10.mesograzer.mass.per.area.site  ~  
                           + PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2
                         + PC1.zos.site + PC2.zos.site + log10.periphyton.mass.per.g.zostera.site
                         , 
                         na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = ZEN_2014_site_means_49_Atlantic)

# Examine summary output
meso.a.rf # % Var explained: 26.74

# Plot error as a function of # of trees
plot(meso.a.rf) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(meso.a.rf) # 

# RESULT: FC2 = PCe3 = PCz2 >> PCe2 > peri > FC1 > ...


# MESOGRAZER BIOMASS: PACIFIC

meso.p.rf = randomForest(log10.mesograzer.mass.per.area.site  ~  
                           + PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2
                         + PC1.zos.site + PC2.zos.site + log10.periphyton.mass.per.g.zostera.site
                         , 
                         na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = ZEN_2014_site_means_Pacific)

# Examine summary output
meso.p.rf # % Var explained: -14.87. Huh?

# Plot error as a function of # of trees
plot(meso.p.rf) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(meso.p.rf) # 

# RESULT: PCe3 > ...


# PERIPHYTON BIOMASS: GLOBAL

peri.rf = randomForest(log10.periphyton.mass.per.g.zostera.site  ~  Ocean
                       + PC1.env.global + PC2.env.global + PC3.env.global + FC1 + FC2
                       + PC1.zos.site + PC2.zos.site
                       ,
                       na.action = na.roughfix, corr.threshold = 0.7, ntree = 1000, data = ZEN_2014_site_means)

# Examine summary output
peri.rf # % Var explained: 8.5

# Plot error as a function of # of trees
plot(peri.rf) # Good. 100 trees is plenty.

# Plot variable importance
varImpPlot(peri.rf) # 

# RESULT: PCe2 > FC2 = ... (not much resolution)


###################################################################################
# MISCELLANEOUS STIATISCAL CHECKS                                                 #
###################################################################################

# Eelgrass leaf C:N ratio (low values) is considered an indicator of light limitation. 
# Any evidence that light limitation is stronger in Atlantic? Answer: No. 

#  leaf length vs C:N ratio
canopy.cn <- ggplot(ZEN_2014_site_means, aes(x = leaf.CN.ratio.site, y = log10.Zostera.longest.leaf.length.cm.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Leaf C:N ratio") +  
  ylab("Longest leaf length (log)") +  
  # scale_x_continuous(limits = c(-4,4.3), breaks = c(-4,-2,0,2,4)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) 
canopy.cn
pdf("canopy.cn.pdf",width = 6.5, height = 6); canopy.cn;  dev.off() 


# Does exposure (indexed by fetch) correlate with leaf length? Answer: No.

#  leaf length vs fetch
canopy.fetch <- ggplot(ZEN_2014_site_means, aes(x = log10.mean.fetch, y = log10.Zostera.longest.leaf.length.cm.site, group = Ocean, col = Ocean)) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue", "forestgreen")) +
  xlab("Fetch (log)") +  
  ylab("Longest leaf length (log)") +  
  # scale_x_continuous(limits = c(-4,4.3), breaks = c(-4,-2,0,2,4)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) 
canopy.fetch
pdf("canopy.fetch.pdf",width = 6.5, height = 6); canopy.fetch;  dev.off() 


###################################################################################
# MISCELLANEOUS FIGURES                                                           #
###################################################################################


# Eelgrass PCz1 predicted by genetics FC1  (Global) 
pcz1.fc1 <- ggplot(ZEN_2014_site_means_49, aes(x = FC1, y = PC1.zos.site, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Genetic FC1 (Global)") +  
  ylab("Eelgrass form PCz1 (Global)") +  
  # scale_x_continuous(limits = c(30,75)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz1.fc1

# Eelgrass PCz1 predicted by genetics FC2  (Global) 
pcz1.fc2 <- ggplot(ZEN_2014_site_means_49, aes(x = FC2, y = PC1.zos.site, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Genetic FC2 (Global)") +  
  ylab("Eelgrass form PCz1 (Global)") +  
  # scale_x_continuous(limits = c(30,75)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz1.fc2

# Eelgrass PCz1 predicted by environment PCe1  (Global) 
pcz1.pce1 <- ggplot(ZEN_2014_site_means_49, aes(x = PC1.env.global, y = PC1.zos.site, col = "blue")) +
  geom_point(size = 4) +
  geom_text(aes(label = unique(Site)), hjust = -0.25, vjust = 0, size = 3) +
  scale_color_manual(values = c("blue")) +
  xlab("Environment PCe1 (Global)") +  
  ylab("Eelgrass form PCz1 (Global)") +  
  # scale_x_continuous(limits = c(30,75)) +
  # scale_y_continuous(limits=c(-3.5, 2.5), breaks = c(2,1,0,-1,-2,-3), position = "right") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = lm, fullrange = F, se = T, lwd = 1.0, na.rm=T) 
pcz1.pce1 
