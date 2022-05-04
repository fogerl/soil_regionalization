### DEM Derivatives
# last_updated: 2022/04/29

# Preparing DEM derivatives slope aspect and TWI 
# Aggregating derivatives to resolution of 1 km and
# calcuating statistical parameters mean, median, standard deviation Q5 and Q95

# Load libraries
library(raster)
library(rgdal)
library(envirem)
library(RSAGA)

# Load digital elevation model (resolution = 10x10 m)
setwd()
DEM_10 <- raster("dhm_at_lamb_10m_2018.tif")

# re-project raster to target CRS 
beginCluster()
dem10proj <- projectRaster(DEM_10, crs = 3416)
endCluster()

# load and preprocess base raster, reprojection will be based on cells of base raster
# dem of spartacus dataset
dem_spartacus <- raster("DEM_SPARTACUS.tif")
# shape of Austrian borders
at_shp <- readOGR("shp_file.shp")
# mask and crop base raster to Austrian borders
dem_mask <- mask(dem_spartacus,at_shp)
base_raster <- crop(dem_mask,extent(at_shp))


# calculate slope and aspect
slope_10 <- terrain(dem10proj)
aspect_10 <- terrain(dem10proj, opt = "aspect")

# aggregate to 1k resolution and calculate statistical parameters
beginCluster()

mean_slope <-raster::aggregate(slope_10,fact = 100, fun = mean)
sd_slope <-raster::aggregate(slope_10,fact = 100, fun = sd)
median_slope <-raster::aggregate(slope_10,fact = 100, fun = median)
q5_slope <- aggregate(x = slope_10, fact = 100, 
                      fun = function(i,...) quantile(i, probs=0.05, na.rm=T)) 
q95_slope <- aggregate(x = slope_10, fact = 100, 
                       fun = function(i,...) quantile(i, probs=0.95, na.rm=T))

mean_aspect <-raster::aggregate(aspect_10,fact = 100, fun = mean)
sd_aspect <-raster::aggregate(aspect_10,fact = 100, fun = sd)
median_aspect <-raster::aggregate(aspect_10,fact = 100, fun = median)
q5_aspect <- aggregate(x = aspect_10, fact = 100, 
                       fun = function(i,...) quantile(i, probs=0.05, na.rm=T)) 
q95_aspect <- aggregate(x = aspect_10, fact = 100, 
                        fun = function(i,...) quantile(i, probs=0.95, na.rm=T))

# stack all parameters together
aspect_predictors <- stack(mean_aspect, sd_aspect,median_aspect, q5_aspect,q95_aspect)
slope_predictors <- stack(mean_slope, sd_slope,median_slope, q5_slope,q95_slope)

endCluster()

# align generated parameters to cells of base raster 
slope_proj <- projectRaster(slope_predictors,base_raster)
aspect_proj <- projectRaster(aspect_predictors,base_raster)
# name layers of raster stack
names(slope_proj) <- c("mean_slope","sd_slope","median_slope","q5_slope","q95_slope")
names(aspect_proj) <- c("mean_aspect","sd_aspect","median_aspect","q5_aspect","q95_aspect")

# export slope and aspect derivatives 
writeRaster(x = aspect_proj, filename = "aspect_predictors.tif", driver = "GeoTiff",options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(x = slope_proj, filename = "slope_predictors.tif", driver = "GeoTiff",options="INTERLEAVE=BAND", overwrite=TRUE)


## calculate TWI 
# set-up sagaEnv
sagaEnv <- rsaga.env("C:/Program Files/QGIS 3.10/apps/saga-ltr", cores = 2, parallel = TRUE)

# aggregate to DEM to 50m resolution
# with finer resolution TWI-calculation was not possible for my computer
beginCluster()
dem50 <- aggregate(dem10proj, fact = 5, fun = mean)
endCluster()

# calculate TWI based on 50m resolution 
beginCluster()
runtime_total <- proc.time()
TWI_50 <- topoWetnessIndex(dem50, sagaEnv)
print(paste("Total Runtime:"))
proc.time() - runtime_total
endCluster() 

# aggregate to target resolution and calculate statistcal parameters
beginCluster()
mean_twi <-raster::aggregate(TWI_50,fact = 20, fun = mean)
sd_twi <-raster::aggregate(TWI_50,fact = 20, fun = sd)
median_twi <-raster::aggregate(TWI_50,fact = 20, fun = median)
q5_twi <- aggregate(x = TWI_50, fact = 20, 
                    fun = function(i,...) quantile(i, probs=0.05, na.rm=T)) 
q95_twi <- aggregate(x = TWI_50, fact = 20, 
                     fun = function(i,...) quantile(i, probs=0.95, na.rm=T))
# stack all parameters together
twi_predictors <- stack(mean_twi, sd_twi,median_twi, q5_twi,q95_twi)
endCluster()

# align generated parameters to cells of base raster 
twi_proj <- projectRaster(twi_predictors,base_raster, method = "ngb")
# name layers of raster stack
names(twi_proj) <- c("mean_twi","sd_twi","median_twi","q5_twi","q95_twi")

# export TWI derivatives
writeRaster(x = twi_proj, filename = "twi_predictors.tif", driver = "GeoTiff",options="INTERLEAVE=BAND", overwrite=TRUE)