### Preparing Input Dataframe for 3D Soil Mapping
# Data Preparation for indirect approach - XGBOOST algorithm 
# For FNN algorithm the same procedure is followed - only PCA Analysis is LEFT OUT
# Target Variables: Sand, Silt, Clay, Humus

## Loading Libraries 
library(raster)
library(fasterize)
library(sf)
library(rgdal)
library(dplyr)
library(tidyverse)
library(landmap)


### LOADING PREPROCESSED PREDICTOR VARIABLES
# set working directory (to be filled in)
setwd() 

# NDVI - preprocessed and resampled to 1k resolution 
ndvi <- stack("ndvi_processed_etrs.tif") 
names(ndvi) <- c("ndvi_mean_years","ndvi_years_range","ndvi_mean_summer", "ndvi_summer_range")

# DEM + derivatives slope, aspect and TWI - preprocessed and resampled to 1k resolution 
slope_predictors <- stack("slope_predictors.tif")
names(slope_predictors) <- c("mean_slope","sd_slope","median_slope","q5_slope","q95_slope")
aspect_predictors <-  stack("aspect_predictors.tif")
names(aspect_predictors) <- c("mean_aspect","sd_aspect","median_aspect","q5_aspect","q95_aspect")
dem_1k <- raster("dem1k.tif")
twi_predictors <- stack("twi_predictors.tif")
names(twi_predictors) <- c("mean_twi","sd_twi","median_twi","q5_twi","q95_twi")

# climate data - preprocessed to monthly mean temperature and long term annual mean precipitation
tm_monthly <- stack("tm_monthly_mean.tif")
months <- paste("tm","_",1:12,sep="")
names(tm_monthly) <- months
rr <- stack("rr_predictors.tif")
names(rr) <- c("rr_mean","rr_range")


# soil form polygons - eBOD database
layers_eBOD <- ogrListLayers("eBOD2_2018.gpkg")
bodenform <- st_read("eBOD2_2018.gpkg", layers_eBOD[1], fid_column_name = "fid")
# reproject to ETRS89 
bodenform <- st_transform(bodenform,crs =crs(tm_monthly))

# transform soil unit names (=character) to numeric code
bodenart_code <- seq(1:13)
bodenart <- unique(bodenform$bodenart)
bodenart_lookup <- cbind(bodenart_code,bodenart)
bodenart_df <- as.data.frame(bodenart_lookup)
bodenart_df <- bodenart_df %>% dplyr::mutate(across(c(where(is.character),-bodenart),function(x) as.numeric(as.character(x))))
merge_bodenart <- merge(bodenart_df,bodenform, by = "bodenart")
bodenart_shp <- st_sf(merge_bodenart)
# rasterize soil unit polygons
bodenart_r <- fasterize(st_collection_extract(bodenart_shp, "POLYGON"), tm_monthly[[1]], field = "bodenart_code")


# landcover map - preprocessed and resampled to 1k resolution  
lisa <- raster("prepL2_1k.tif")
lisa <- projectRaster(lisa, tm_monthly, method = "ngb")


# lithology
geology <- st_read("hydrogeology.shp")
geology_r <- fasterize(geology, tm_monthly[[1]], field = "litho_neu")


## Stacking together and renaming all Covariate layers 
cov <- stack(ndvi, slope_predictors, aspect_predictors, twi_predictors, dem_1k, tm_monthly, rr, bodenart_r, lisa, geology_r)
names(cov) <- c("ndvi_mean_years","ndvi_years_range","ndvi_mean_summer","ndvi_summer_range","mean_slope", "sd_slope", "median_slope", "q5_slope", 
                "q95_slope", "mean_aspect","sd_aspect",  "median_aspect","q5_aspect",  "q95_aspect", "mean_twi", "sd_twi","median_twi", "q5_twi",
                "q95_twi", "dem1k","tm_1","tm_2","tm_3","tm_4","tm_5","tm_6","tm_7", "tm_8", "tm_9","tm_10","tm_11","tm_12","rr_mean","rr_range", 
                "bodenart", "LISA", "lithology")

## PCA Analysis (only if XGBoost model is used) 
# Filtering out missing pixels and artefacts with PCA 
cov_px <- as(cov,"SpatialPixelsDataFrame")
# change data type of categorical variables
cov_px$bodenart <- as.factor(cov_px$bodenart)
cov_px$LISA <- as.factor(cov_px$LISA)
cov_px$lithology <- as.factor(cov_px$lithology)
# PCA transformation 
cov_spc <- landmap::spc(cov_px, ~ndvi_mean_years +ndvi_years_range + ndvi_mean_summer + ndvi_summer_range  +
                       mean_slope +sd_slope + median_slope + q5_slope+ q95_slope+ mean_aspect +sd_aspect+ median_aspect + q5_aspect+ q95_aspect+
                       mean_twi+sd_twi+median_twi+q5_twi+q95_twi+ dem1k+   tm_1+    tm_2+     tm_3+
                       tm_4+    tm_5+   tm_6+   tm_7+ tm_8+   tm_9+  tm_10+  tm_11+   tm_12+   
                       rr_mean+rr_range+bodenart +LISA + lithology )

# Filtering out missing pixels and artefacts with PCA 
plot(cumsum(cov_spc@pca$sdev^2/sum(cov_spc@pca$sdev^2)))
# use only info that explains >95% of variance
pc_var_expl <- cumsum(cov_spc@pca$sdev^2/sum(cov_spc@pca$sdev^2))
pc.use <- first(which(pc_var_expl>0.95)) 

# truncated data
trunc <- cov_spc@predicted[1:pc.use]
grid1k <- as(trunc, "SpatialGridDataFrame")


### LOADING AND PREPROCESSING TARGET VARIABLES (Soil data)
## Profile data from eBOD database - preprocessing and aggregating to depth levels 0-20,20-50,50-100 cm
# load spatial information of profiles
profil <- st_read("eBOD2_2018.gpkg",layers_eBOD[3], fid_column_name = "prof7_id")
profil_pre <- profil[,1]
profil_pre <- st_transform(profil_pre, crs(rr))
# load soil information of profiles
horizont_text <- st_read("eBOD2_2018.gpkg", layers_eBOD[22], fid_column_name = "hor8_id")

# select  attributes of interest and clean data
horizont_text <- horizont_text %>% dplyr::select(prof7_id, von_cm, bis_cm, entnahmetiefe_cm, sand_prozent, schluff_prozent, ton_prozent,humus)
clean_data <- horizont_text %>% filter(von_cm>bis_cm | von_cm<0 | bis_cm<0)
#drop rows with unplausible depth information 
horizont_text <- horizont_text[!(horizont_text$von_cm > horizont_text$bis_cm | horizont_text$von_cm<0 | horizont_text$bis_cm<0), ,drop=F]
# correct data type
horizont_text$humus <- as.numeric(horizont_text$humus)


#aggregate to depth 0-20cm
df0_20 <- horizont_text %>% filter(von_cm<20)
# new column "bis_20cm" and reorder
df0_20$bis_20cm <- ifelse(df0_20$bis_cm <=20, df0_20$bis_cm,20)
df0_20 <- df0_20[,c(1:3,9,4:8)]
# new column "horizon" = ("bis_20cm"-"von_cm") 
df0_20$horizon <- c((df0_20$bis_20cm-df0_20$von_cm))  
df0_20 <- df0_20[,c(1:4,10,5:9)]
#clean texture data from NA values (if target variable is humus drop_na(humus))
texture0_20 <- df0_20 %>% drop_na(sand_prozent)
# calculate weighted mean for top horizon
texture0_20 <- texture0_20 %>% group_by(prof7_id) %>% 
  dplyr::summarise(von_cm = min(von_cm),
                   bis_cm = max(bis_20cm),
                   sand = weighted.mean(sand_prozent,horizon), 
                   clay = weighted.mean(ton_prozent,horizon),
                   silt = weighted.mean(schluff_prozent,horizon),
                   total = sum(sand,clay,silt),
                   humus = weighted.mean(humus,horizon))
 
# add DEPTH information
texture0_20$DEPTH <- 10


# same procedure as above for depth aggregation 20-50cm and 50-100cm
df20_50 <- horizont_text %>% filter(von_cm<50 & bis_cm>20)
df20_50$von_20cm <- ifelse(df20_50$von_cm >=20, df20_50$von_cm,20)
df20_50$bis_50cm <- ifelse(df20_50$bis_cm <=20, df20_50$bis_cm,50)
df20_50$horizon <- c((df20_50$bis_50cm-df20_50$von_20cm))   
df20_50 <- df20_50[,c(1:3,9,10,11,4:8)]
texture20_50 <- df20_50 %>% drop_na(sand_prozent)
texture20_50 <- texture20_50 %>% group_by(prof7_id) %>% 
  dplyr::summarise(von_cm = min(von_20cm),
            bis_cm = max(bis_50cm),
            sand = weighted.mean(sand_prozent,horizon), 
            clay = weighted.mean(ton_prozent,horizon),
            silt = weighted.mean(schluff_prozent,horizon),
            total = sum(sand,clay,silt),
            humus = weighted.mean(humus,horizon))


texture20_50$DEPTH <- 35

df50_100 <- horizont_text %>% filter(von_cm<100 & bis_cm>50)
df50_100$von_50cm <- ifelse(df50_100$von_cm >=50, df50_100$von_cm,50)
df50_100$bis_100cm <- ifelse(df50_100$bis_cm <=100, df50_100$bis_cm,100)
df50_100$horizon <- c((df50_100$bis_100cm-df50_100$von_50cm))   
df50_100 <- df50_100[,c(1:3,9,10,11,4:8)]
texture50_100 <- df50_100 %>% drop_na(sand_prozent)
texture50_100 <- texture50_100 %>% group_by(prof7_id) %>% 
  dplyr::summarise(von_cm = min(von_50cm),
            bis_cm = max(bis_100cm),
            sand = weighted.mean(sand_prozent,horizon), 
            clay = weighted.mean(ton_prozent,horizon),
            silt = weighted.mean(schluff_prozent,horizon),
            total = sum(sand,clay,silt),
            humus = weighted.mean(humus,horizon))

texture50_100$DEPTH <- 75

# combine profile information for all three depth levels 
profile_depth <- rbind(texture0_20, texture20_50, texture50_100)
# combine profil information with spatial information 
profile_depth <- merge(profil_pre, profile_depth, by ="prof7_id")
# correct data type
profile_depth$prof7_id <- as.factor(profile_depth$prof7_id)


## Load preprocessed soil data from BORIS and aggregate to three depth levels
bzi <- read.csv("bzi.csv")

# Depth level 0-20cm 
bzi0_20 <- bzi %>% filter(von_cm<20)
bzi0_20$bis_20cm <- ifelse(bzi0_20$bis_cm <=20, bzi0_20$bis_cm,20)
bzi0_20$horizon <- c((bzi0_20$bis_20cm-bzi0_20$von_cm))  

bzi0_20_final <- bzi0_20 %>% group_by(prof7_id) %>% 
  dplyr::summarise(von_cm = min(von_cm),
            bis_cm = max(bis_20cm),
            sand = weighted.mean(sand,horizon), 
            clay = weighted.mean(clay,horizon),
            silt = weighted.mean(silt,horizon),
            total = sum(sand,clay,silt),
            humus = weighted.mean(toc,horizon),
            x = min(x), y = min(y))
bzi0_20_final$DEPTH <- 10

# Depth level 20-50cm 
bzi20_50 <- bzi %>% filter(von_cm<50 & bis_cm>20)
bzi20_50$von_20cm <- ifelse(bzi20_50$von_cm >=20, bzi20_50$von_cm,20)
bzi20_50$bis_50cm <- ifelse(bzi20_50$bis_cm <=50, bzi20_50$bis_cm,50)
bzi20_50$horizon <- c((bzi20_50$bis_50cm-bzi20_50$von_20cm))   

bzi20_50_final <- bzi20_50 %>% group_by(prof7_id) %>% 
  dplyr::summarise(von_cm = min(von_20cm),
            bis_cm = max(bis_50cm),
            sand = weighted.mean(sand,horizon), 
            clay = weighted.mean(clay,horizon),
            silt = weighted.mean(silt,horizon),
            total = sum(sand,clay,silt),
            humus = weighted.mean(toc,horizon),
            x = min(x), y = min(y))
bzi20_50_final$DEPTH <- 35

# Depth level 50-100cm 
bzi50_100 <- bzi %>% filter(von_cm<100 & bis_cm>50)
bzi50_100$von_50cm <- ifelse(bzi50_100$von_cm >=50, bzi50_100$von_cm,50)
bzi50_100$bis_100cm <- ifelse(bzi50_100$bis_cm <=100, bzi50_100$bis_cm,100)
bzi50_100$horizon <- c((bzi50_100$bis_100cm-bzi50_100$von_50cm))   
bzi50_100 <- bzi50_100 %>% group_by(prof7_id) %>% 
  dplyr::summarise(von_cm = min(von_50cm),
            bis_cm = max(bis_100cm),
            sand = weighted.mean(sand,horizon), 
            clay = weighted.mean(clay,horizon),
            silt = weighted.mean(silt,horizon),
            total = sum(sand,clay,silt),
            humus = weighted.mean(toc,horizon),
            x = min(x), y = min(y))
bzi50_100$DEPTH <- 75

# combine all BORIS soil information again and transform to spatial points dataframe
bzi_depth <- rbind(bzi0_20_final,bzi20_50_final, bzi50_100)
xy <- bzi_depth[,c(9,10)]
bzi_depth_sp <- SpatialPointsDataFrame(coords = xy, data = bzi_depth, 
                                       proj4string = CRS("+proj=lcc +lat_1=49 +lat_2=46 +lat_0=47.5 +lon_0=13.33333333333333 +x_0=400000 +y_0=400000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

# transform preprocessed eBOD information to spatial points data frame 
profile_depth_sp <- as_Spatial(profile_depth)

# check if duplicates in coordinates 
zerodist2(profile_depth_sp, bzi_depth_sp) # none

# merge relevant columns of soil data
bzi_depth_sp <- bzi_depth_sp[,-c(2,3,7:10)]
profile_depth_sp <- profile_depth_sp[,-c(2,3,7,8)]
merge_all <- rbind(bzi_depth_sp, profile_depth_sp)





## OVERLAY TARGET DATA WITH PREDICTOR DATA
ov <- data.frame(merge_all, raster::extract(stack(grid1k), merge_all), xy=TRUE)
# for FNN  model
# ov <- data.frame(merge_all, raster::extract(cov, merge_all), xy=TRUE)
merge_export <-rowid_to_column(merge,"rowid")

# export data
write_csv(merge_export,"data_2022-02-02_texture_xgboost.csv")


