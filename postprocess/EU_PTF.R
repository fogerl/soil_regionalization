### Indirect Approach
## predict soil property maps for the whole study area
## use soil property maps as input data for EU-PTFs
## predict Ksat and uncertainty quantiles for the whole study area

# Load libraries
library(tidyvers)
library(dplyr)
library(raster)
library(terra)
library(sf)
library(sp)
library(euptf2)


# Load XGBOOST models 
setwd()
final_sand <- readRDS("indirect/result/run_2022-02-02_sand_PCA/model.rds") 
final_silt <- readRDS("indirect/result/run_2022-02-02_silt_PCA/model.rds") 
final_clay <- readRDS("indirect/result/run_2022-02-02_clay_PCA/model.rds") 
final_oc <- readRDS("indirect/result/run_2021-10-10_humus_PCA/model.rds")

# predict target variables for whole study area and top horizon 
# predictor variables from input_indirect.R
# define prediction depth
grid1k$DEPTH = 10  
# convert predictor variables to dataframe 
data_at <- as.data.frame(grid1k, xy =TRUE)
# clean and adapt data frame structure
data_at <- drop_na(data_at)
data_at$clay <- as.numeric(NA)
data_at$silt <- as.numeric(NA)
data_at$sand <- as.numeric(NA)

# predict soil properties
sand_pred <- predict(final_sand,data_at, na.action = na.pass)
silt_pred <- predict(final_silt,data_at, na.action = na.pass)
clay_pred <- predict(final_clay,data_at, na.action = na.pass)
oc_pred <- predict(final_oc,data_at, na.action = na.pass)

output <- data_at
output$sand_pred <- sand_pred
output$silt_pred <- silt_pred
output$clay_pred <- clay_pred
output$oc_pred <- oc_pred

# Correct unplausible data values
output$sand_pred <- ifelse(output$sand_pred > 100, 100, output$sand_pred)
output$sand_pred <- ifelse(output$sand_pred < 0, 0, output$sand_pred)
output$silt_pred <- ifelse(output$silt_pred > 100, 100, output$silt_pred)
output$silt_pred <- ifelse(output$silt_pred < 0, 0, output$silt_pred)
output$clay_pred <- ifelse(output$clay_pred > 100, 100, output$clay_pred)
output$clay_pred <- ifelse(output$clay_pred < 0, 0, output$clay_pred)
output$oc_pred <-  ifelse(output$oc_pred < 0, 0, output$oc_pred)
output$oc_pred <-  ifelse(output$oc_pred >100, 100, output$oc_pred)

# Scale data (texture has to sum up to 100 %)
output_cor <- output %>% mutate(sand_cor =sand_pred/(silt_pred+sand_pred+clay_pred)*100,
                                silt_cor = silt_pred/(silt_pred+sand_pred+clay_pred)*100,
                                clay_cor = clay_pred/(silt_pred+sand_pred+clay_pred)*100,
                                total_cor = sand_cor+silt_cor + clay_cor,
                                oc_pred  =  oc_pred)

# define position of coordinates
xy <- output_cor[,c(34,35)]
# convert to spatial dataframe
output_sp <- SpatialPointsDataFrame(coords = xy, data = output_cor,
                                    proj4string = CRS("+proj=lcc +lat_1=49 +lat_2=46 +lat_0=47.5 +lon_0=13.33333333333333 +x_0=400000 +y_0=400000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

# create raster maps (base raster = dem1k loaded in input_indirect.R)
sand_r <- raster::rasterize(output_sp, dem_1k, field = "sand_cor")
silt_r <- rasterize(output_sp, dem_1k, field = "silt_cor")
clay_r <- rasterize(output_sp, dem_1k, field = "clay_cor")
oc_r <-  rasterize(output_sp, dem_1k, field = "oc_pred")

# repeat same procedure for mid and bottom horizon (prediction depth = 35 and 75 cm)
#----------------------------------------------------------------------------------------


### Apply EU-PTFs

input_ptf <- output_cor
# rename columns according to required input names for EU-PTF
input_ptf <- input_ptf %>% dplyr::rename(DEPTH_M = DEPTH, OC = oc_pred, SAND = sand_cor, SILT = silt_cor, CLAY = clay_cor)
# only keep relevant columns
input_ptf <- input_ptf[,-c(1:32,36:42)]

# transform texture fractions to USDA texture fractions
input_ptf_transf <- TT.text.transf(
  tri.data = input_ptf,
  base.css.ps.lim = c(0,2,50,2000),
  dat.css.ps.lim = c(0,2,63,2000)
)

# transform structure for applying EU PTF
input_ptf_transf <- dplyr::rename(input_ptf_transf,
                                  USSAND = SAND, USSILT = SILT, USCLAY = CLAY)


# choose PTF
which_PTF(predictor = input_ptf_transf,target = c("KS")) # PTF02
# derive Ksat
pred_KS_AT <- euptfFun(ptf = "PTF02", predictor = input_ptf_transf, target = "KS", query = "predictions")

# transform unit (cm/d) to log10(cm/d)
pred_KS_AT$log_KS <- log10(pred_KS_AT$KS_PTF02) 

# define coordinates and transform to spatial data
xy <- pred_KS_AT[,c(4,5)]
pred_sp <- SpatialPointsDataFrame(coords = xy, data = pred_KS_AT, 
            proj4string = CRS("+proj=lcc +lat_1=49 +lat_2=46 +lat_0=47.5 +lon_0=13.33333333333333 +x_0=400000 +y_0=400000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
ks_pred <- terra::rasterize(vect(pred_sp), rast(dem_1k), field = "log_KS")


# derive uncertainty quantiles (Q25 and Q75)
quant_KS_top <- euptfFun(ptf = "PTF02", predictor = input_ptf_transf, target = "KS", query = "quantiles", quantiles = c(.25,.75))

# transform unit (cm/d) to log10(cm/d)
quant_KS_top$log25 <- log10(quant_KS_top$`KS_PTF02_quantile= 0.25`)
quant_KS_top$log75 <- log10(quant_KS_top$`KS_PTF02_quantile= 0.75`)

# define coordinates and transform to spatial data
xy <- quant_KS_top[,c(5,6)]
quant_sp <- SpatialPointsDataFrame(coords = xy, data = quant_KS_top, proj4string = CRS("+proj=lcc +lat_1=49 +lat_2=46 +lat_0=47.5 +lon_0=13.33333333333333 +x_0=400000 +y_0=400000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
quant25 <- rasterize(quant_sp, dem_1k, field = "log25")
quant75 <- rasterize(quant_sp, dem_1k, field = "log75")

# repeat same procedure for mid and bottom horizon 