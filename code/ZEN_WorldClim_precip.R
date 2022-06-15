library(rgdal)
library(raster)

# files are from WorldClim database (precipitation 2.5 arc-minute)
# download here: http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/prec_2-5m_bil.zip
# original data are from weather stations, and global, regional, national, local data sources from 1950s-2000 
# website says data range is 1960-1990

# ESRI format is global but split by month
r1 = raster("data/input/prec/prec_1")
r2 = raster("data/input/prec/prec_2")
r3 = raster("data/input/prec/prec_3")
r4 = raster("data/input/prec/prec_4")
r5 = raster("data/input/prec/prec_5")
r6 = raster("data/input/prec/prec_6")
r7 = raster("data/input/prec/prec_7")
r8 = raster("data/input/prec/prec_8")
r9 = raster("data/input/prec/prec_9")
r10= raster("data/input/prec/prec_10")
r11= raster("data/input/prec/prec_11")
r12= raster("data/input/prec/prec_12")
s <- stack( r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12)
r <- sum( s )
r

writeRaster( r, filename="data/output/WorldClim_precip_2-5.tif" )

plot( r )


# # You can also download 30 arc-seconds for 30 degree cells 
# # see http://www.worldclim.org/tiles.php
# # 11, 12, 22, 13, 15, 16, 6, 110
# w = getData('worldclim', var='prec', res=0.5, lon=0, lat=45)
# w
# #
