# Climate data 
# last_updated: 2022/04/29 HZ

# Preparing Climate predictors
# monthly mean temperature
# long term annual mean precipitation
# range of long term annual mean precipitation
# data basis = Spartacus (Hiebl and Frei, 2016; Hiebl and Frei, 2018) 



# Load libraries 
library(raster)
library(rgdal)
library(foreach)
library(rasterVis)
require(snow)
library(rts)

# set path to data 
setwd()
pathToNC_Backbone <- c("Spartacus/RR/","Spartacus/Tn/","Spartacus/Tx/") 

# generate names and path to all files
rastFiles <- matrix(0, nrow = length(list.files(path = pathToNC_Backbone[1], pattern = '\\.nc$')), ncol = 3)
rastFiles_short <- rastFiles

for (p in 1:length(pathToNC_Backbone)) {
  pathToNc <- pathToNC_Backbone[p]
  rastFiles[, p] <- list.files(path = pathToNc, pattern = '\\.nc$', full.names = TRUE)
  rastFiles_short[, p] <- list.files(path = pathToNc, pattern = '\\.nc$')
}

# define weights for calculating mean temperature
g1 = 0.5
g2 = 0.5

# loop over single year files and load data
beginCluster()
runtime_total <- proc.time()
for (f in 1:nrow(rastFiles)) {
  runtime <- proc.time()
  print(rastFiles_short[f, 1])
  
  
  # read single years as rasterbricks
  rr <- brick(rastFiles[f, 1]) # precip
  tn <- brick(rastFiles[f, 2]) # Tmin
  tx <- brick(rastFiles[f, 3]) # Tmax
  
  # calculate Tmean as a function of weights
  tm <- g1*tn + g2*tx
} # end loop over single year files

print(paste("Total Runtime:"))
proc.time() - runtime_total
endCluster()


beginCluster()
runtime_total <- proc.time()
# extract dates 
date <- gsub("X","",names(rr), fixed = TRUE) 
date <- gsub(".","-",date, fixed = TRUE) 
date <- as.Date(date)
# create time series raster and calculate yearly sum 
rr_rts <- rts(rr,date)
rr_rts_sum <- apply.yearly(rr_rts,sum) 

print(paste("Total Runtime:"))
proc.time() - runtime_total
endCluster()

# convert time series back to raster stack
rr_yearly_sum <- stack(rr_rts_sum@raster)
# name layers of raster stack after the corresponding year
years <- 1961:2020
names(rr_yearly_sum) <- years
# take out year 2020 because it's only available till September 
rr_yearly_sum <- dropLayer(rr_yearly_sum,60)

# calculate precipitation predictors
# mean of all years
rr_mean <- stackApply(rr_yearly_sum, indices =  rep(1,nlayers(rr_yearly_sum)), fun = "mean", na.rm = T)
# range (min and max) over whole time series 
rr_min <- stackApply(rr_yearly_sum, indices =  rep(1,nlayers(rr_yearly_sum)), fun = "min", na.rm = T)
rr_max <- stackApply(rr_yearly_sum, indices =  rep(1,nlayers(rr_yearly_sum)), fun = "max", na.rm = T)
rr_range <- rr_max - rr_min

# combine all perdictors in a stack and export them 
rr_predictors <- stack(rr_mean,rr_min, rr_max,rr_range)
names(rr_predictors) <- c("rr_mean", "rr_min", "rr_max","rr_range")
writeRaster(x = rr_predictors, filename = "rr_predictors.tif", driver = "GeoTiff",options="INTERLEAVE=BAND", overwrite=TRUE)


# Create time series for mean temperature
tm_rts <- rts(tm, date)

# Calculate monthly means
beginCluster()
runtime_total <- proc.time()

tm_rts_monthly_mean <- apply.monthly(tm_rts,mean) 


print(paste("Total Runtime:"))
proc.time() - runtime_total
endCluster()

# Convert raster time series back to raster stack
monthly_mean <- raster::stack(tm_rts_monthly_mean@raster)
# take out year 2020 because data is only available till 09/2020
tm_monthly_mean <- dropLayer(monthly_mean,c(709:717))

# name raster stack after months
months <- c("Jan", "Feb", "March", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec")
# create sequence for the 59 years of data
months_59 <- rep(months,59)
# name monthly mean raster stack
names(tm_monthly_mean) <- months_59

# separate  raster stack by month and calculate long- time monthly mean 
tm_01 <- subset(tm_monthly_59,  grep("Jan", names(tm_monthly_59)))
tm_02 <- subset(tm_monthly_59,  grep("Feb", names(tm_monthly_59)))
tm_03 <- subset(tm_monthly_59,  grep("March",names(tm_monthly_59)))
tm_04 <- subset(tm_monthly_59,  grep("Apr", names(tm_monthly_59)))
tm_05 <- subset(tm_monthly_59,  grep("May", names(tm_monthly_59)))
tm_06 <- subset(tm_monthly_59,  grep("Jun", names(tm_monthly_59)))
tm_07 <- subset(tm_monthly_59,  grep("Jul", names(tm_monthly_59)))
tm_08 <- subset(tm_monthly_59,  grep("Aug", names(tm_monthly_59)))
tm_09 <- subset(tm_monthly_59,  grep("Sept",names(tm_monthly_59)))
tm_10 <- subset(tm_monthly_59,  grep("Oct", names(tm_monthly_59)))
tm_11 <- subset(tm_monthly_59,  grep("Nov", names(tm_monthly_59)))
tm_12 <- subset(tm_monthly_59,  grep("Dec", names(tm_monthly_59)))

mean_01 <- stackApply(tm_01, indices =  rep(1,nlayers(tm_01)), fun = "mean", na.rm = T)
mean_02 <- stackApply(tm_02, indices =  rep(1,nlayers(tm_02)), fun = "mean", na.rm = T)
mean_03 <- stackApply(tm_03, indices =  rep(1,nlayers(tm_03)), fun = "mean", na.rm = T)
mean_04 <- stackApply(tm_04, indices =  rep(1,nlayers(tm_04)), fun = "mean", na.rm = T)
mean_05 <- stackApply(tm_05, indices =  rep(1,nlayers(tm_05)), fun = "mean", na.rm = T)
mean_06 <- stackApply(tm_06, indices =  rep(1,nlayers(tm_06)), fun = "mean", na.rm = T)
mean_07 <- stackApply(tm_07, indices =  rep(1,nlayers(tm_07)), fun = "mean", na.rm = T)
mean_08 <- stackApply(tm_08, indices =  rep(1,nlayers(tm_08)), fun = "mean", na.rm = T)
mean_09 <- stackApply(tm_09, indices =  rep(1,nlayers(tm_09)), fun = "mean", na.rm = T)
mean_10 <- stackApply(tm_10, indices =  rep(1,nlayers(tm_10)), fun = "mean", na.rm = T)
mean_11 <- stackApply(tm_11, indices =  rep(1,nlayers(tm_11)), fun = "mean", na.rm = T)
mean_12 <- stackApply(tm_12, indices =  rep(1,nlayers(tm_12)), fun = "mean", na.rm = T)

# stack monthly means together 
monthly_tm <- stack(mean_01,mean_02,mean_03,mean_04,mean_05,mean_06,mean_07,mean_08,mean_09,mean_10,mean_11,mean_12)

# export monthly meanss
writeRaster(x = monthly_tm, filename = "tm_monthly_mean.tif", driver = "GeoTiff",options="INTERLEAVE=BAND", overwrite=TRUE)
