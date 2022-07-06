### ZEN 2014 geography
# code to merge datasets into the main data file
# by Matt Whalen
# started 5 July 2022

# packages
library(tidyverse)

# read data
d <- read.csv("data/input/archive/ZEN_2014_main_data_20220705.csv", header = TRUE)


# EPIBIOTA 
# Epibiota (periphyton) data were not calculated correctly in the summary 'main" data file, 
# specifically filtered material ("Epibiota filter") was not divided by the mass of Zostera scraped. 
# So we have to regenerate these numbers from scratch. 
epibiota <- read.csv("data/input/archive/ZEN_2014_epibiota_mass_recalc.csv", header = TRUE)
################################################################################

# First, Japan B separated "Chl large epiphytes filter" from "Epibiota filter", the latter of which are all zero. 
# We need to  recode those labeled "Chl large epiphytes filter" as "Epibiota filter", so we can add them together. 
epibiota$Species <- as.factor(epibiota$Species)
levels(epibiota$Species)[levels(epibiota$Species) == "Chl Large Epiphytes Filter"] <- "Epibiota filter"
# levels(epibiota$Species) #This appears to have worked ...

# Second, Oregon A misnamed the epibiota filters ...
levels(epibiota$Species)[levels(epibiota$Species) == "Epibiota Packet"] <- "Epibiota filter"

# Subset to only the variables of interest here
epibiota <- subset(epibiota, select = c(Site, Site.Code, Subsite, 
                                        Sampling.Time, Plot.ID, Unique.ID, Species, Taxa, Group, Type, Dry.Mass.g.))


# Recode second sampling time as 1 for LI
epibiota$Site.Code <- as.factor(epibiota$Site.Code)
epibiota[epibiota$Site.Code == "LI" & epibiota$Sampling.Time == "1", "Sampling.Time"] = 3
epibiota[epibiota$Site.Code == "LI" & epibiota$Sampling.Time == "2", "Sampling.Time"] = 1

# Regenerate Unique IDs
epibiota$Unique.ID = as.character(epibiota$Unique.ID)
epibiota[epibiota$Site.Code == "LI", "Unique.ID"] = 
  paste(epibiota[epibiota$Site.Code == "LI", "Site.Code"], 
        epibiota[epibiota$Site.Code == "LI", "Subsite"], 
        epibiota[epibiota$Site.Code == "LI", "Sampling.Time"], 
        epibiota[epibiota$Site.Code == "LI", "Plot.ID"],
        sep = ".")

epibiota$Unique.ID = as.factor(epibiota$Unique.ID)

# Include only first sampling time
epibiota = subset(epibiota, Sampling.Time == "1")


# # create new wide-form data frame containing only the dry mass data
# epibiota.temp <- epibiota %>%
#   # only take along columns that are unique, otherwise output is staggered in chunks
#   select(Unique.ID, Species, Dry.Mass.g.) %>%
#   # group the data by species and sample
#   group_by(Unique.ID, Species) %>%
#   # sum the dry mass for each species in each sample (i.e., sum measurements from the same unique.ID)
#   dplyr::summarize(Dry.mass.g = sum(Dry.Mass.g.)) %>%
#   # cast species as columns
#   spread(Species, Dry.mass.g, fill = 0)
################################################################################
# ----------  WHALEN UPGRADING THIS TO USE TIDYVERSE
# create new wide-form data frame containing only the dry mass data
epibiota.temp <- epibiota %>%
  # only take along columns that are unique, otherwise output is staggered in chunks
  select(Unique.ID, Species, Dry.Mass.g.) %>%
  # group the data by species and sample
  group_by(Unique.ID, Species) %>%
  # sum the dry mass for each species in each sample (i.e., sum measurements from the same unique.ID)
  dplyr::summarize(Dry.mass.g = sum(Dry.Mass.g.)) %>%
  # cast species as columns
  pivot_wider( names_from = Species, values_from = Dry.mass.g, values_fill =  0)
################################################################################

# rename some variables 
names(epibiota.temp)[names(epibiota.temp)=="Epibiota filter"] <- "epibiota.filter"
names(epibiota.temp)[names(epibiota.temp)=="Zostera marina"] <- "epibiota.zostera.marina"

# Subset to only the variables of interest
epibiota.temp <- subset(epibiota.temp, select = c(Unique.ID, epibiota.filter, epibiota.zostera.marina))

# Calculate periphyton mass per g Zostera marina:
epibiota.temp$periphyton.mass.per.g.zostera <- epibiota.temp$epibiota.filter / epibiota.temp$epibiota.zostera.marina

# Add raw and normalized periphyton data back into main data frame
d$epibiota.filter <- epibiota.temp$epibiota.filter[match(d$Unique.ID, epibiota.temp$Unique.ID)]
d$epibiota.zostera.marina <- epibiota.temp$epibiota.zostera.marina[match(d$Unique.ID, epibiota.temp$Unique.ID)]
d$periphyton.mass.per.g.zostera <- epibiota.temp$periphyton.mass.per.g.zostera[match(d$Unique.ID, epibiota.temp$Unique.ID)]

# Remove miscalculated periphyton variable from summary data set
d <- subset(d, select = -c(Epibiota.Periphyton))



# write to disk
write_csv( d, "data/input/ZEN_2014_main_data.csv")
