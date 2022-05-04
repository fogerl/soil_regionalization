# Soil property prediction
# Hanna Zeitfogel, Moritz Feigl

# Load libraries and xgboost model code
setwd() # to be filled in 
library(caret)
library(magrittr)
source("utils.R")

# Choose specific dataset or leave as "" to use the last one
dataset <- ""
variable <- "humus"
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
data <- read.csv(paste0("data/", dataset ))[, -c(1,2)]
if(colnames(data)[1] == "X") data <- data[, -1]

preprocessed_data <- data
# train/test split
data_date <- gsub(".csv", "", dataset, fixed = TRUE) %>% 
  strsplit(., "_") %>% 
  unlist %>% 
  tail(., 1)
model_name <- paste0("run_", data_date, "_", variable)
# train/test split
set.seed(250690)
test_ids <- sample(nrow(preprocessed_data), floor(nrow(preprocessed_data)*0.1))
train_data <- preprocessed_data[-test_ids, ]
test_data <- preprocessed_data[test_ids, ]
if(!file.exists("results")) dir.create(file.path("results"), showWarnings = FALSE)
if(!file.exists(paste0("results/", model_name))){
  dir.create(file.path(paste0("results/", model_name)), showWarnings = FALSE)
}
write.csv(train_data, paste0("results/", model_name, "/training_data.csv"))
write.csv(test_data, paste0("results/", model_name, "/test_data.csv"))

# Model run 
formular <- "humus ~."
xgboost_from_df(train_data, 
                formular, 
                test_data,
                model_name = model_name, 
                n_random_initial_points = 10, 
                n_iter = 40,
                seed = 220640,
                use_previous_grid = TRUE)

