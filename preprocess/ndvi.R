### NDVI 
# last_updated: 2022/04/29

# Preparing NDVI predictors
# Aggregating derivatives to resolution of 1 km 
# NDVI time series 2002-2020 downloaded from 
# https://ivfl-info.boku.ac.at/satellite-data-processing (Vuolo et al., 2012)
# spatial resolution 250 m, temporal resolution 7 days

# Load libraries
library(raster)
library(rgdal)
library(rts)

# load downloaded data and combine all tifs in one stack
setwd("NDVI_time_series_2002_2021/filtered")
rlist=list.files("NDVI_BOKU/NDVI_time_series_2002_2021/filtered", pattern="tif$", full.names=FALSE)
for(i in rlist) { assign(unlist(strsplit(i, "[.]"))[1], raster(i)) } 
allrasters <- raster::stack(rlist)

# only select NDVI rasters (EVI rasters also present in the folder)
r_ndvi <- subset(allrasters,  grep("ndvi", names(allrasters)))
ndvi <- gsub("MODIS.ndvi.","",names(r_ndvi), fixed = TRUE) 
ndvi <- gsub(".yL6000.BOKU","",ndvi, fixed = TRUE) 
names(r_ndvi) <- ndvi

# scale dataset
ndvi_scaled <- r_ndvi/10000
# extract dates 
date_ndvi <- as.Date(ndvi,"%Y%j")
# create raster time series 
ndvi_rts <- rts(ndvi_scaled, date_ndvi)
# apply the mean function on each year
ndvi_yearly <- apply.yearly(ndvi_rts,mean) 

# calculate mean for each quarter
ndvi_quarterly <- apply.quarterly(ndvi_rts, mean)

#save the process
setwd("C:/Users/Admin/Documents/NDVI/preprocessed rasters")
writeRaster(ndvi.y@raster, filename="yearlymeans.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(ndvi.q@raster, filename="quarterlymeans.tif", options="INTERLEAVE=BAND", overwrite=TRUE)

# omit negative and >1 values 
ndvi_yearly[ndvi_yearly<0]=0
ndvi_yearly[ndvi_yearly>1]=0

# year 2021 is redundant - only 1 month available
# therefore, remove last year
ndvi_02_20 <- dropLayer(ndvi_yearly,20)

# calculate mean over whole time series 
mean_years <- stackApply(ndvi_02_20, indices =  rep(1,nlayers(ndvi_02_20)), fun = "mean", na.rm = T)

# range over whole time series 
min_years <- stackApply(ndvi_02_20, indices =  rep(1,nlayers(ndvi_02_20)), fun = "min", na.rm = T)
max_years <- stackApply(ndvi_02_20, indices =  rep(1,nlayers(ndvi_02_20)), fun = "max", na.rm = T)
years_range <- max_years - min_years

# extract quarter 2 and (=summer months, april-september)
indices =  c(seq(2,77,4),seq(3,77,4))
indices = sort(indices)
ndvi_q2q3 <- subset(ndvi_quarterly,indices)
# calculate mean and range for summer months for the whole time series
mean_summer <- stackApply(ndvi_q2q3, indices = rep(1,nlayers(ndvi_q2q3)), fun = "mean", na.rm = T)
min_summer <- stackApply(ndvi_q2q3, indices =  rep(1,nlayers(ndvi_q2q3)), fun = "min", na.rm = T)
max_summer <- stackApply(ndvi_q2q3, indices =  rep(1,nlayers(ndvi_q2q3)), fun = "max", na.rm = T)
summer_range <- max_summer-min_summer

# combine relevant ndvi predictors in a stack 
ndvi_predictors <- stack(mean_years,years_range,mean_summer, summer_range)
names(ndvi_predictors) <- c("mean_years","years_range","mean_summer", "summer_range")

# align raster stack with base raster (base raster was created in dem_derivatives.R script)
ndvi_proj <- projectRaster(ndvi_selected,base_raster)
ndvi_crop <- crop(ndvi_proj,base_raster)
ndvi_mask <- mask(ndvi_crop,base_raster)

# export ndvi predictors
writeRaster(x = ndvi_mask, filename = "ndvi_processed.tif", driver = "GeoTiff",options="INTERLEAVE=BAND", overwrite=TRUE)



