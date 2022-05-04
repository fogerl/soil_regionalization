# Soil property prediction
# Hanna Zeitfogel, Moritz Feigl

setwd("ml_models/indirect")
# Load libraries and fnn model code
library(caret) 
library(magrittr)
setwd() # to be filled in 
source("code/utils_oc.R")

### Choose specific dataset or leave as "" to use the last one ###
dataset <- ""

# load data
if(dataset == ""){
  all_datasets <- list.files("data")
  chosen_dataset <- strsplit(all_datasets, "_") %>% 
    lapply(., tail, 1) %>% 
    lapply(., gsub, pattern = ".csv", replacement = "", fixed = TRUE) %>% 
    unlist %>% 
    as.Date %>% 
    which.max
  dataset <- all_datasets[chosen_dataset]
}
data <- read.csv(paste0("data/", dataset ))[,-c(1,2)]
if(colnames(data)[1] == "X") data <- data[, -1]

# pre-process categorical variables
data$soil_unit <- as.factor(data$soil_unit)
data$lithology <- as.factor(data$lithology)
data$landcover <- as.factor(data$landcover)
dummies <- dummyVars(~ ., data)
preprocessed_data <- predict(dummies, newdata = data) 

# train/test split
data_date <- gsub(".csv", "", dataset, fixed = TRUE) %>% 
  strsplit(., "_") %>% 
  unlist %>% 
  tail(., 1)
model_name <- paste0("run_", data_date, "_FNN")
set.seed(250690)
test_ids <- sample(nrow(preprocessed_data), floor(nrow(preprocessed_data)*0.1))
train_data <- as.data.frame(preprocessed_data[-test_ids, ])
test_data <- as.data.frame(preprocessed_data[test_ids, ])

if(!file.exists("results")) dir.create(file.path("results"))
if(!file.exists(paste0("results/", model_name))){
  dir.create(file.path(paste0("results/", model_name)), showWarnings = FALSE)
}
write.csv(train_data, paste0("results/", model_name, "/training_data.csv"))
write.csv(test_data, paste0("results/", model_name, "/test_data.csv"))

# Model run FNN
soil_fnn(train_data, 
         test_data,
         y_variables = c("humus"),
         model_name = model_name, 
         n_random_initial_points = 20, 
         n_iter = 40,
         seed = 220640,
         initial_grid_from_model_scores = FALSE)
