### Direct Approach
# predict Ksat raster for study area
# based on trained XGBoost models

# load model 
setwd()
hydrobod_xgboost <- readRDS("direct/results/run_2021-09-22_xboost_hydrobod_PCA/model.rds")
input_data_hydrobod <- read.csv("direct/results/run_2021-09-22_xboost_hydrobod_PCA/test_data.csv")
test_pred_hydrobod <- read.csv("direct/results/run_2021-09-22_xboost_hydrobod_PCA/xgboost_test_prediction.csv")

# predict target variable for whole study area and top horizon 
# predictor variables from input_indirect.R
# define prediction depth
grid1k$DEPTH = 10 

# convert predictor variables to dataframe
data_at_top <- as.data.frame(grid1k, xy =TRUE)
# clean and adapt data frame structure
data_at_top <- drop_na(data_at_top)
data_at_top$ksat <- as.numeric(NA)
data_at_top <- data_at_top %>% 
  rename("x" = "s1",
         "y" = "s2",)
data_at_top$rowid <- as.numeric(NA)
pred_hydrobod_top <- predict(hydrobod_xgboost,data_at_top, na.action = na.pass)
data_at_top$pred <- pred_hydrobod_top

# correct prediction values (hydrobod10_1k loaded in input_direct.R)
data_at_top$pred <- ifelse(data_at_top$pred > max(values(hydrobod10_1k), na.rm = TRUE), 
                           max(values(hydrobod10_1k), na.rm = TRUE), data_at_top$pred)
data_at_top$pred <- ifelse(data_at_top$pred < 0, 0, data_at_top$pred)

# define coordinates and transform to spatial data (base raster = dem_1k loaded in input_indirect.R)
xy <- data_at_top[,c(34,35)]
hydrobod_top_sp <- SpatialPointsDataFrame(coords = xy, data = data_at_top, 
                      proj4string = CRS("+proj=lcc +lat_1=49 +lat_2=46 +lat_0=47.5 +lon_0=13.33333333333333 +x_0=400000 +y_0=400000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
hydrobod_top <- rasterize(hydrobod_top_sp, dem_1k, field = "pred")

# repeat procedure for mid and bottom horizon 


