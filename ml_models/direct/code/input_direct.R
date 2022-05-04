### Preparing Input Dataframe for 3D Soil Mapping
# Data Preparation for direct approach - XGBOOST algorithm 
# For FNN algorithm the same procedure is followed - only PCA Analysis is LEFT OUT
# Target Variables: Ksat
# Ksat Data Basis: Hydrobod Project II


# Load Libraries
library(raster)
library(sp)
library(dplyr)
library(terra)
library(sf)
library(rgdal)
library(tidyverse)

# Load target datasets
setwd()
ks_hydrobod0_20 <- raster("hybod_k1.tif")
ks_hydrobod20_50 <- raster("hybod_k2.tif")
ks_hydrobod50_100 <- raster("hybod_k3.tif")

# Load base raster
dem_1k <- raster("dem1k.tif")

# reproject data to target CRS and resolution 
hydrobod10_1k <- projectRaster(ks_hydrobod0_20, dem_1k)
hydrobod35_1k <- projectRaster(ks_hydrobod20_50, dem_1k)
hydrobod75_1k <- projectRaster(ks_hydrobod50_100, dem_1k)

## Predictor variables for XGBoost
# convert Spatial Points Dataframe to raster stack
# grid1k generated in input_indirect.R 
# (for FNN use cov stack instead of grid1k - also generated in input_indirect.R)
grid_stack <- stack(grid1k) 

# combine target variable with predictor variables
hydro_stack <- stack(hydrobod10_1k, grid_stack)
# convert to dataframe
hydro_df_at <- as.data.frame(hydro_stack, xy=TRUE, id = row.number())
# add depth information 
hydro_df_at$DEPTH <- 10
# rename columns
colnames(hydro_df_at)[3] <- "ksat"

# same procedure for other two depth levels
hydro_stack_35 <- stack(hydrobod35_1k, grid_stack)
hydro_df_35 <- as.data.frame(hydro_stack_35, xy=TRUE, id = row.number())
hydro_df_35$DEPTH <- 35
colnames(hydro_df_35)[3] <- "ksat"
hydro_stack_75 <- stack(hydrobod75_1k, grid_stack)
hydro_df_75 <- as.data.frame(hydro_stack_75, xy=TRUE, id = row.number())
hydro_df_75$DEPTH <- 75
colnames(hydro_df_75)[3] <- "ksat"

# combine input dataframes of all three depth levels to one dataframe
hydro_df <- rbind(hydro_df_at,hydro_df_35,hydro_df_75)

# clean data
hydro_df_depth <- drop_na(hydro_df)
hydro_df_depth <- rowid_to_column(hydro_df_depth)

# export dataframe
write_csv(hydro_df_depth, "data_2021-09-22_ks_PCA.csv")


