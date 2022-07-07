#########################################
# "ZEN GLOBAL PREDATION MAPPING 2014" #
########################################

# code by Matt Whalen
# last updated 2017.04.24

# This script prepares Bio-ORACLE and WorldClim data for ZEN sites,
# then merges these data with other relevant ZEN data
# http://www.oracle.ugent.be/download.html (Shift-click to follow link)
# downloaded .rar data were extracted in Ubuntu using unrar function
# Note that BioOracle offers many useful oceanographic predictors
# 

# Hat tip to Brian Cheng for getting me started
# https://bscheng.com/2016/10/22/spatial-data-using-raster-in-r/

# Run this code to get data for all 23 variables in the Bio-ORACLE dataset
# (or, at least how every many you have in your directory)

# load libraries
library(rgdal)
library(raster)

# function to extract raster data within a radius defined by buffer (in meters...10000m = 10km)
buff <- function(rast,buffer=10000) {
  unlist(lapply(extract( rast, zen.geo, buffer=buffer ),mean,na.rm=T))
}

########### Bio-ORACLE ##################
# get a list of files to work with, all stored in a folder named "BioOracle Data"
files <- list.files( "data/input/BioOracle/",pattern = ".asc")
# read in the raster data files and give appropriate projects (lat/lon WGS84)
r <- lapply( paste0("data/input/BioOracle/",files), raster ) 
# consider adding a projection (+init=epsg:4326 +proj=longlat +ellps=WGS84)
# crop all rasters to same extent, and get rid of southern hemisphere
e <- extent(-180,180,0,70)
r2 <- lapply( r, function(rast) crop(rast,e) )

########### WorldClim Precipitation ###############
precip <- raster( "data/output/WorldClim_precip_2-5.tif" )



########### ZEN Subsite Data ##################
# Read in ZEN site data to get Latitude and Longitude of each site
SiteData <- read.csv("data/input/Duffy_et_al_2022_site_metadata.csv")

# just isolate the Lat/Long
zen.geo  <- SiteData[ ,c("longitude","latitude") ] 



## Extract data for each site and add to the ZEN data we already have
# SiteData$sst.min <- extract( sst.min, zen.geo, method='bilinear' )
# No Bio-ORACLE data available for Croatia, France subsite B
# Use a buffer over which to look for surrounding cells that might have data and average them
# This will take a while for many variables
buffer <- 10000

# Apply the extract function to all rasters
zen.oracle <- lapply( r2, buff )
zen.precip <- buff( precip )

# combine all of these into a data.frame and give them names
Environmentals <- data.frame( do.call( cbind, zen.oracle ) )
names(Environmentals) <- unlist( strsplit( files, ".asc") )

# add the environmental data as additional columns on SiteData
SiteData <- cbind( SiteData, Environmentals, precip=zen.precip )

# # write the new data to disk
write.csv( SiteData, "data/output/Duffy_et_al_2022_environmental.csv" )

