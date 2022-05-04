# utility functions
# soil regionalization project
# Feb 2021
# Hanna Zeitfogel, Moritz Feigl

library(dplyr)
library(tensorflow)
library(keras)

# XGBoost with hyperparameter optimization
xgboost_from_df <- function(train_data,
                            formular,
                            test_data = NULL,
                            model_name = NULL,
                            no_cores = parallel::detectCores() - 1,
                            seed = NULL,
                            n_iter = 40,
                            n_random_initial_points = 20,
                            use_previous_grid = FALSE){
  
  # checks
  model_short <- "XGBoost"
  if(is.character(formular)) formular <- as.formula(formular)
  dependent_variable <- all.vars(formular)[1]
  cat("*** Starting XGBoost computation ***\n")
  # start time
  start_time <- Sys.time()
  # define model name with date if not selected
  if(is.null(model_name)) model_name <- paste0("xgboost_", format(Sys.time(), "%Y-%m-%d"))
  # create model_name folder and results folder
  if(!file.exists("results")) dir.create(file.path("results"))
  if(!file.exists(paste0("results/", model_name))){
    dir.create(file.path(paste0("results/", model_name)))
  }
  model_folder <- paste0("results/", model_name)
  # remove NA rows
  na_train <- which(is.na(train_data), arr.ind = TRUE)
  if(nrow(na_train) > 0) train_data <- train_data[-unique(na_train[,1]),]
  if(!is.null(test_data)){
    na_test <- which(is.na(test_data), arr.ind = TRUE)
    if(nrow(na_test) > 0) test_data <- test_data[-na_test[,1],]
  }
  # random seed
  if(is.null(seed)) seed <- sample(1000:100000, 1)
  
  # XGBoost
  # Train control settings
  seeds <- vector(mode = "list", length = 51)
  set.seed(seed)
  for(i in 1:51) seeds[[i]] <- sample(10000, 200)
  tc <- caret::trainControl(method = "repeatedcv",
                            number = 10,
                            repeats = 5,
                            allowParallel = TRUE,
                            verboseIter = TRUE,
                            seeds = seeds)
  
  # optimization parameter bounds
  parameter_bounds <- list(
    "nrounds" = c(300L, 3000L),#c(300L, 3000L),
    "eta" = c(0.001, 0.3),
    "max_depth" = c(3L, 20L),#c(3L, 12L),
    "min_child_weight" = c(1L, 15L),
    "subsample" = c(0.5, 1),#c(0.7, 1),
    "colsample_bytree" = c(0.5, 1),#c(0.7, 1),
    "gamma" = c(0L, 5L)
  )
  # initial sampled or previously computed grid
  if("hyperpar_opt_scores.csv" %in% list.files(model_folder) & use_previous_grid){
    # previous computed grid results
    hyperpar_opt_scores <- read.csv(paste0(model_folder, "/hyperpar_opt_scores.csv"))
    names(hyperpar_opt_scores)[8:9] <- c("RMSE", "MAE")
    initial_grid <- hyperpar_opt_scores[, 1:7]
    initial_grid$Value <- hyperpar_opt_scores$RMSE * -1
    if(nrow(initial_grid) > n_random_initial_points) {
      n_iter <- n_iter - (nrow(initial_grid) - n_random_initial_points)
    }
  } else {
    # new initial grid
    initial_grid <- parameter_bounds
    for(i in 1:length(parameter_bounds)){
      range <- parameter_bounds[[i]]
      if(is.integer(range)){
        initial_grid[[i]] <- sample(range[1]:range[2], n_random_initial_points, replace = TRUE)
      } else {
        initial_grid[[i]] <- round(sample(seq(range[1], range[2], length.out = 10),
                                          n_random_initial_points, replace = TRUE), 3)
      }
    }
    initial_grid <- as.data.frame(do.call(cbind, initial_grid))
    # run initial grid models
    nul_out <- ifelse(Sys.info()["sysname"] == "Windows", "NUL", "/dev/null")
    cl <- parallel::makeCluster(no_cores) # convention to leave 1 core for OS
    
   # parallel::clusterEvalQ(cl, .libPaths("/home/lv71468/hzeitfogel/Rlibs")) # hanna is testing this 
    
    doParallel::registerDoParallel(cl)
    model <- caret::train(formular,
                          data = train_data,
                          method = "xgbTree",
                          trControl = tc,
                          nthread = 1,
                          tuneGrid = initial_grid)
    parallel::stopCluster(cl)
    hyperpar_opt_scores <- model$results[, c(names(initial_grid), "RMSE", "MAE")]
    write.csv(hyperpar_opt_scores,
              file = paste0(model_folder, "/hyperpar_opt_scores.csv"),
              row.names = FALSE)
    initial_grid <- hyperpar_opt_scores[, 1:7]
    initial_grid$Value <- hyperpar_opt_scores$RMSE * -1
  }
  
  # define XGBoost function for optimization
  xgb_optim <- function(nrounds, eta, max_depth, min_child_weight,
                        subsample, colsample_bytree, gamma) {
    tune_grid <- expand.grid(nrounds = nrounds,
                             max_depth = max_depth,
                             eta = eta,
                             gamma = gamma,
                             colsample_bytree = colsample_bytree,
                             subsample = subsample,
                             min_child_weight = min_child_weight)
    cl <- parallel::makeCluster(no_cores) # convention to leave 1 core for OS
    doParallel::registerDoParallel(cl)
    model <- caret::train(formular,
                          data = train_data,
                          method = "xgbTree",
                          trControl = tc,
                          nthread = 1,
                          tuneGrid = tune_grid)
    parallel::stopCluster(cl)
    # save results of given hyperparameter set
    hyperpar_opt_scores <- read.csv(paste0(model_folder, "/hyperpar_opt_scores.csv"))
    hyperpar_opt_scores <- rbind(hyperpar_opt_scores,
                                 model$results[, names(hyperpar_opt_scores)])
    write.csv(hyperpar_opt_scores,
              file = paste0(model_folder, "/hyperpar_opt_scores.csv"),
              row.names = FALSE)
    return(list(Score = -caret::getTrainPerf(model)[, "TrainRMSE"], Pred = 0))
  }
  if(n_iter > 0){
    set.seed(seed)
    Bopt_xgboost <- rBayesianOptimization::BayesianOptimization(xgb_optim,
                                                                bounds = parameter_bounds,
                                                                init_grid_dt = initial_grid,
                                                                init_points = 0,
                                                                n_iter = n_iter,
                                                                acq = "ucb",
                                                                kappa = 2.576,
                                                                eps = 0.0,
                                                                verbose = TRUE)
  }
  initial_grid$Value <- NULL
  # update colnames of hyperpar_opt_scores
  hyperpar_opt_scores <- read.csv(paste0(model_folder, "/hyperpar_opt_scores.csv"))
  hyperpar_opt_scores[, c("RMSE", "MAE")] <- round(hyperpar_opt_scores[, c("RMSE", "MAE")], 3)
  best_par <- hyperpar_opt_scores[which.min(hyperpar_opt_scores$RMSE), ]
  # Optimized hyperparameters
  best_par <- data.frame(nrounds = best_par$nrounds,
                         eta = best_par$eta,
                         max_depth =  best_par$max_depth,
                         gamma = best_par$gamma,
                         colsample_bytree = best_par$colsample_bytree,
                         subsample = best_par$subsample,
                         min_child_weight = best_par$min_child_weight)
  cat("Finished hyperparameter optimization, optimized parameters:\n")
  for(i in 1:ncol(best_par)) cat(names(best_par)[i], ":", as.numeric(best_par[1, i]), "\n")
  # Run model with optimized hyperparameters
  cat("Running XGBoost with optimized hyperparameter set...")
  cl <- parallel::makeCluster(no_cores) # convention to leave 1 core for OS
  doParallel::registerDoParallel(cl)
  model <- caret::train(formular,
                        data = train_data,
                        method = "xgbTree",
                        trControl = tc,
                        nthread = 1,
                        tuneGrid = best_par)
  parallel::stopCluster(cl)
  cat("Done!\n")
  
  # save all hyperparameter set performances
  hyperpar_opt_scores2 <- model$results[, c(names(initial_grid), "RMSE", "MAE")]
  hyperpar_opt_scores <- rbind(hyperpar_opt_scores,
                               hyperpar_opt_scores2)
  hyperpar_opt_scores[, c("RMSE", "MAE")] <- round(hyperpar_opt_scores[, c("RMSE", "MAE")], 3)
  colnames(hyperpar_opt_scores) <- c(names(initial_grid), "cv_or_validation_RMSE", "cv_or_validation_MAE")
  cat("\nHyperparameter optimization results are saved in",
      paste0(model_folder, "/hyperpar_opt_scores.csv\n"))
  write.csv(hyperpar_opt_scores,
            paste0(model_folder, "/hyperpar_opt_scores.csv"),
            row.names = FALSE)
  
  cv_or_val_results <- model$results[which.min(model$results$RMSE), c("RMSE", "MAE")]
  
  # train prediction
  suppressWarnings({train_prediction <- predict(model, train_data)})
  write.csv(data.frame("observation" = train_data[, dependent_variable],
                       "prediction" = train_prediction),
            paste0(model_folder, "/xgboost_training_prediction.csv"), row.names = FALSE)
  # test prediction
  if(!is.null(test_data)) {
    suppressWarnings({test_prediction <- predict(model, test_data)})
    write.csv(data.frame("observation" = test_data[, dependent_variable],
                         "prediction" = test_prediction),
              paste0(model_folder, "/xgboost_test_prediction.csv"), row.names = FALSE)
  }
  # model diagnostics
  RMSE <- function(prediction, observation){
    round(sqrt(mean((prediction - observation)^2, na.rm = TRUE)), 3)
  }
  
  MAE <- function(prediction, observation){
    round(mean(abs(observation - prediction)), 3)
  }
  # training performance
  train_rmse <- RMSE(train_prediction, train_data[, dependent_variable])
  train_mae <- MAE(train_prediction, train_data[, dependent_variable])
  # CV or validation set performance
  cv_or_validation_rmse <- round(as.numeric(cv_or_val_results["RMSE"]), 3)
  cv_or_validation_mae <- round(as.numeric(cv_or_val_results["MAE"]), 3)
  # test performance (optional)
  if(!is.null(test_data)){
    test_rmse <- RMSE(test_prediction, test_data[, dependent_variable])
    test_mae <- MAE(test_prediction, test_data[, dependent_variable])
  } else {
    test_rmse <- NA
    test_mae <- NA
  }
  # run time
  run_time <-  paste0(
    round((as.numeric(Sys.time()) - as.numeric(start_time))/60, 2),
    " minutes")
  # output data frame and save as csv
  model_scores <- data.frame(start_time = as.character(start_time),
                             run_time = run_time,
                             model_name = model_name,
                             train_RMSE = train_rmse,
                             train_MAE = train_mae,
                             cv_or_validation_RMSE = cv_or_validation_rmse,
                             cv_or_validation_MAE = cv_or_validation_mae,
                             test_RMSE = test_rmse,
                             test_MAE = test_mae,
                             stringsAsFactors = FALSE)
  
  if("model_scores.csv" %in% list.files(model_folder)){
    model_scores_old <- read.csv(paste0(model_folder, "/model_scores.csv"),
                                 stringsAsFactors = FALSE)
    write.csv(rbind(model_scores_old, model_scores, stringsAsFactors = FALSE),
              paste0(model_folder, "/model_scores.csv"),
              row.names = FALSE)
  } else {
    write.csv(model_scores,
              paste0(model_folder, "/model_scores.csv"),
              row.names = FALSE)
  }
  cat("\nModel training performance:", paste0("RMSE = ", model_scores$train_RMSE, ","),
      paste0("MAE = ", model_scores$train_MAE))
  if(!is.null(test_data)){
    cat("\nModel testing performance:", paste0("RMSE = ", model_scores$test_RMSE, ","),
        paste0("MAE = ", model_scores$test_MAE), "\n")
  }
  cat("\nModel results are saved in",
      paste0(model_folder, "/model_scores.csv\n"))
  
  if(!is.null(model)){
    saveRDS(model, paste0(model_folder, "/model.rds"))
    cat("The trained model is saved in",
        paste0(model_folder, "/model.rds\n"))
  }
  
  # importance plot
  cat("\nSaving variable importance plot in",
      paste0(model_folder, "/"),"\n")
  importance <- caret::varImp(model)$importance
  v <- as.numeric(importance[,1])
  w <- rownames(importance)
  DF <- data.frame(w,v, stringsAsFactors = FALSE)
  ggplot(DF, aes(x = stats::reorder(w,v), y = v, fill = v))+
    geom_bar(stat="identity", position="dodge") + coord_flip() +
    ylab("Variable Importance") +
    xlab("") +
    ggtitle("Information Value Summary - XBoost") +
    guides(fill = FALSE) +
    scale_fill_gradient(low = "red", high = "blue") 
  ggsave(paste0(model_folder, "/importance_plot.png"), 
           height = 12, width = 8, units = "in")
  xgb_imp <- xgboost::xgb.importance(feature_names = model$finalModel$feature_names,
                                     model = model$finalModel)
  p <- xgboost::xgb.ggplot.importance(xgb_imp) 
  ggsave(paste0(model_folder, "/importance_plot_cluster.png"), 
           height = 12, width = 8, units = "in")
  write.csv(p[[1]], paste0(model_folder, "/importance_cluster_results.csv"), row.names = FALSE)
}

# Feedforward neural network with bayesian hyperaprameter optimization
soil_fnn <- function(train_data,
                     test_data = NULL,
                     y_variables,
                     model_name = NULL,
                     seed = NULL,
                     n_iter = 40,
                     n_random_initial_points = 20,
                     epochs = 100,
                     early_stopping_patience = 5,
                     ensemble_runs = 5,
                     bounds_layers = c(1, 5),
                     bounds_units = c(5, 200),
                     bounds_dropout = c(0, 0.2),
                     bounds_batch_size = c(5, 150),
                     initial_grid_from_model_scores = FALSE){
  
  # start time
  start_time <- Sys.time()
  # define model name with date if not selected
  if(is.null(model_name)) model_name <- paste0("fnn_", format(Sys.time(), "%Y-%m-%d"))
  # create model_name folder and results folder
  if(!file.exists("results")) dir.create(file.path("results"))
  if(!file.exists(paste0("results/", model_name))){
    dir.create(file.path(paste0("results/", model_name)))
  }
  # Check if data set are data.frames
  if(class(train_data) != "data.frame" | class(train_data) != "data.frame"){
    stop("Both train_data and test_data have to be data.frames!")
  }
  # model folder
  model_folder <- paste0("results/", model_name)
  # remove NA rows
  na_train <- which(is.na(train_data), arr.ind = TRUE)
  if(nrow(na_train) > 0) train_data <- train_data[-unique(na_train[,1]),]
  if(!is.null(test_data)) {
    na_test <- which(is.na(test_data), arr.ind = TRUE)
    if(nrow(na_test) > 0) test_data <- test_data[-na_test[,1],]
  }
  # train/val split
  part_training <- nrow(train_data)/4 * 3
  train_length <- floor(part_training)
  nn_val <- train_data[(train_length + 1):nrow(train_data), ]
  nn_train <- train_data[1:train_length, ]
  # define y and x
  y_train <- nn_train[, y_variables] / 100
  y_full_train <- train_data[, y_variables] / 100
  y_val <- nn_val[, y_variables] / 100
  if(!is.null(test_data)) y_test <- test_data[, y_variables] / 100
  x_train <- nn_train[, which(!(names(nn_train) %in% y_variables))]
  x_full_train <- train_data[, which(!(names(train_data) %in% y_variables))]
  x_val <- nn_val[, which(!(names(nn_val) %in% y_variables))]
  if(!is.null(test_data)) x_test <- test_data[, which(!(names(test_data) %in% y_variables))]
  # scale features with mean and sd from training data
  train_means <- apply(x_train, 2, mean)
  train_sd <- apply(x_train, 2, sd)
  cat("Mean and standard deviation used for feature scaling are saved under",
      paste0(model_folder, "/scaling_values.csv\n"))
  write.csv(rbind(train_means, train_sd),
            paste0(model_folder, "/scaling_values.csv"))
  # apply it on all data sets
  numeric_cols <- which(apply(x_train, 2, function(x) length(unique(x)) > 2))
  for (i in numeric_cols){
    x_train[, i] <- (x_train[, i] - train_means[i]) / train_sd[i]
    x_full_train[, i] <- (x_full_train[, i] - train_means[i]) / train_sd[i]
    x_val[, i] <- (x_val[, i] - train_means[i]) / train_sd[i]
    if(!is.null(test_data)) x_test[, i] <- (x_test[, i] - train_means[i]) / train_sd[i]
  }
  x_train <- x_train %>% as.matrix()
  y_train <- y_train %>% as.matrix()
  x_full_train <- x_full_train %>% as.matrix()
  x_val <- x_val %>% as.matrix()
  y_val <- y_val %>% as.matrix()
  if(!is.null(test_data)){
    x_test <- x_test %>% as.matrix()
  } else {
    x_test <- y_test <- test <- NULL
  }
  
  # initial value for flag -> if additional initial_grid points should be calculated
  ini_grid_cal_flag <- FALSE
  # if there are no model score available -> set initial_grid_from_model_scores = FALSE
  if(initial_grid_from_model_scores){
    if(!("hyperpar_opt_scores.csv" %in%
         list.files(paste0(model_folder)))) {
      initial_grid_from_model_scores <- FALSE
    }
  }
  if(initial_grid_from_model_scores){
    # get initial grid for optimization from the previous calculated model_scores
    cat("Using existing scores as initial grid for the Bayesian Optimization\n")
    hyperpar_opt_scores.csv <- read.csv(paste0(model_folder, "/hyperpar_opt_scores.csv"),
                                        stringsAsFactors = FALSE)
    initial_grid <- hyperpar_opt_scores.csv[c("layers", "units", "dropout",
                                              "batch_size", "cv_or_validation_RMSE")]
    if(nrow(initial_grid) < n_random_initial_points) ini_grid_cal_flag <- TRUE
    
  }
  # should a random grid be calculated first
  if(!initial_grid_from_model_scores | ini_grid_cal_flag) {
    set.seed(seed)
    n_random <- ifelse(ini_grid_cal_flag,
                       n_random_initial_points - nrow(initial_grid),
                       n_random_initial_points)
    grid <- data.frame(
      "layers" = replicate(n = n_random,
                           sample(
                             x = seq(bounds_layers[1], bounds_layers[2]),
                             size = 1)),
      "units" = replicate(n = n_random,
                          sample(
                            x = seq(bounds_units[1], bounds_units[2]), size = 1)),
      "dropout" = replicate(n = n_random,
                            sample(
                              x = seq(bounds_dropout[1], bounds_dropout[2], by = 0.025),
                              size = 1)),
      "batch_size" = replicate(n = n_random,
                               sample(
                                 x = seq(bounds_batch_size[1], bounds_batch_size[2]),
                                 size = 1))
    )
    cat(paste0("\nRandom hyperparameter sampling:\n"))
    # run initial grid models
    grid_results <- mapply(
      soil_nn,
      batch_size = grid$batch_size,
      layers = grid$layers,
      units = grid$units,
      dropout = grid$dropout,
      ensemble_runs = 1,
      MoreArgs = list(x_train = x_train, y_train = y_train,
                      x_val = x_val, y_val = y_val,
                      x_full_train = x_full_train,
                      y_full_train = y_full_train,
                      epochs = epochs,
                      early_stopping_patience = early_stopping_patience,
                      model_folder = model_folder,
                      model_name = model_name,
                      seed = seed, 
                      y_variables = y_variables))
    # Define initial_grid for bayesian hyperaparameter optimization
    if(ini_grid_cal_flag){
      additional_grid_points <- cbind(grid, grid_results)
      names(additional_grid_points) <- names(initial_grid)
      initial_grid <- rbind(initial_grid, additional_grid_points)
    } else {
      initial_grid <- cbind(grid, grid_results)
    }
  }
  # Bayesian Hyperparameter optimization
  cat("Bayesian Hyperparameter Optimization:\n")
  if(nrow(initial_grid) > n_random_initial_points){
    n_iter <- n_iter - (nrow(initial_grid) - n_random_initial_points)
    cat(nrow(initial_grid) - n_random_initial_points,
        "iterations were already computed\n")
  }
  colnames(initial_grid)[5] <- "Value"
  if(n_iter > 0){
    # function to be optimized
    Bopt_nnmodel <- function(layers, units, dropout,  batch_size) {
      results <- soil_nn(
        x_train = x_train,
        y_train = y_train,
        x_val = x_val,
        y_val = y_val,
        x_test = x_test,
        y_test = y_test,
        x_full_train = x_full_train,
        y_full_train = y_full_train,
        layers = layers,
        units = units,
        dropout = dropout,
        batch_size = batch_size,
        epochs = epochs,
        early_stopping_patience = early_stopping_patience,
        ensemble_runs = 1,
        model_name = model_name,
        model_folder = model_folder,
        y_variables = y_variables,
        seed = seed)
      return(list("Score" = results*-1, "Pred" = 0))
    }
    initial_grid$Value <- initial_grid$Value * -1
    set.seed(seed)
    Bopt_nn <- rBayesianOptimization::BayesianOptimization(
      Bopt_nnmodel,
      bounds = list(layers = as.integer(bounds_layers),
                    units = as.integer(bounds_units),
                    dropout = bounds_dropout,
                    batch_size = as.integer(bounds_batch_size)),
      n_iter = n_iter,
      init_grid_dt = initial_grid,
      acq = "ucb", kappa = 2.576, eps = 0.0,
      verbose = TRUE)
    all_model_results <- Bopt_nn$History
  } else {
    # if all iterations are done -> use initial grid
    all_model_results <- initial_grid
    all_model_results$Value <- all_model_results$Value * -1
  }
  
  # run best model as ensemble and save results
  cat("Run the best performing model as ensemble:\n")
  top_n_model_results <- all_model_results %>%
    top_n(n = 1, wt = Value) %>%
    mutate(dropout = round(dropout, 2)) %>%
    dplyr::select(-Value)
  if(nrow(top_n_model_results) != 1) top_n_model_results <- top_n_model_results[1, ]
  set.seed(seed)
  soil_nn(
    layers = top_n_model_results$layers,
    units = top_n_model_results$units,
    dropout = top_n_model_results$dropout,
    batch_size = top_n_model_results$batch_size,
    x_train = x_train,
    y_train = y_train,
    x_val = x_val,
    y_val = y_val,
    x_test = x_test,
    y_test = y_test,
    epochs = epochs,
    early_stopping_patience = early_stopping_patience,
    ensemble_runs = ensemble_runs,
    model_name = model_name,
    model_folder = model_folder,
    save_model_and_prediction = TRUE,
    x_full_train = x_full_train,
    y_full_train = y_full_train,
    start_time = start_time,
    seed = seed,
    y_variables = y_variables)
}

soil_nn <- function(x_train, x_val, x_test = NULL, x_full_train = NULL,
                    y_train, y_full_train = NULL, y_val, y_test = NULL, batch_size,
                    data_inputs, layers, units, dropout, ensemble_runs,
                    epochs, early_stopping_patience,
                    model_folder = model_folder, model_name,
                    save_model_and_prediction = FALSE,
                    start_time = NULL,
                    seed,
                    y_variables = y_variables){
  
  if(!save_model_and_prediction){
    cat("layers = ", layers,
        ", units = ", units,
        ", dropout = ", dropout,
        ", batch_size = ", batch_size,
        ", ensemble_runs = ", ensemble_runs, ", ", sep = "")
  }
  
  # Ensemble runs
  set.seed(seed)
  for(run in 1:ensemble_runs){
    # define model object
    model <- create_model(nn_type, x_train, layers, units, dropout, seed = seed + run,
                          y_variables = y_variables)
    # training
    tensorflow::tf$random$set_seed((seed + run))
    history <- model %>% keras::fit(x = x_train, y = y_train,
                                    epochs = epochs,
                                    batch_size = batch_size,
                                    callbacks = list(
                                      callback_early_stopping(patience = early_stopping_patience,
                                                              restore_best_weights = TRUE)),
                                    validation_data = list(x_val, y_val),
                                    verbose = 0)
    # store prediction results in lists
    if(run == 1){
      predict_val <- list(predict(model, x_val))
      if(save_model_and_prediction){
        predict_full_train <- list(predict(model, x_full_train))
        if(!is.null(x_test)) predict_test <- list(predict(model, x_test))
      }
    } else {
      predict_val[[run]] <- predict(model, x_val)
      if(save_model_and_prediction){
        predict_full_train[[run]] <- predict(model, x_full_train)
        if(!is.null(x_test)) predict_test[[run]] <- predict(model, x_test)
      }
    }
    # save model as rds
    if(save_model_and_prediction){
      saveRDS(model, paste0(model_folder, "/model", run, ".rds"))
    }
  }
  # get mean prediction
  if(ensemble_runs != 1){
    # reshape list to length = 1 with mean predictions
    mean_pred_results_val <- predict_val[[1]]
    for(col in 1:length(y_variables)){
      mean_pred_results_val[, col] <- predict_val %>%
        lapply(., function(x) x[, col]) %>% 
        do.call(rbind, .) %>%
        as.data.frame() %>%
        as.list() %>%
        lapply(matrix, ncol = ensemble_runs) %>%
        lapply(function(x) apply(x, 1, mean)) %>%
        unlist()
    }
    if(save_model_and_prediction){
      if(!is.null(x_test)){
        mean_pred_results_test <- predict_test[[1]]
        for(col in 1:length(y_variables)){
          mean_pred_results_test[, col] <- predict_test %>%
            lapply(., function(x) x[, col]) %>% 
            do.call(rbind, .) %>%
            as.data.frame() %>%
            as.list() %>%
            lapply(matrix, ncol = ensemble_runs) %>%
            lapply(function(x) apply(x, 1, mean)) %>%
            unlist()
        }
      }
      mean_pred_results_full_train <- predict_full_train[[1]]
      for(col in 1:length(y_variables)){
        mean_pred_results_full_train[, col] <- predict_full_train %>%
          lapply(., function(x) x[, col]) %>% 
          do.call(rbind, .) %>%
          as.data.frame() %>%
          as.list() %>%
          lapply(matrix, ncol = ensemble_runs) %>%
          lapply(function(x) apply(x, 1, mean)) %>%
          unlist()
      }
    }
  } else {
    mean_pred_results_val <- predict_val[[1]]
    if(save_model_and_prediction){
      mean_pred_results_full_train <- predict_full_train[[1]]
      if(!is.null(x_test)) mean_pred_results_test <- predict_test[[1]]
    }
  }
  # validation scores
  RMSE_val <- matrix(NA, ncol = length(y_variables) + 1)
  colnames(RMSE_val) <- c(paste0("val_RMSE_", y_variables), "val_RMSE_mean")
  MAE_val <- matrix(NA, ncol = length(y_variables) + 1)
  colnames(MAE_val) <- c(paste0("val_MAE_", y_variables), "val_MAE_mean")
  for(var in 1:length(y_variables)){
    RMSE_val[1,var] <- RMSE(prediction = mean_pred_results_val[, var]* 100,
                            observation = y_val[, var]*100) %>% round(3)
    MAE_val[1, var] <- MAE(prediction = mean_pred_results_val[, var] * 100,
                           observation = y_val[, var]*100) %>% round(3)
  }
  RMSE_val[1,length(y_variables) + 1] <- round(mean(RMSE_val[1, 1:length(y_variables)]), 3)
  MAE_val[1,length(y_variables) + 1] <- round(mean(MAE_val[1, 1:length(y_variables)]), 3)
  # save hyperparameter results
  model_scores <- data.frame("layers" = layers,
                             "units" = units,
                             "max_epochs" = epochs,
                             "early_stopping_patience" = early_stopping_patience,
                             "batch_size" = batch_size,
                             "dropout" = dropout,
                             "ensemble_runs" = ensemble_runs,
                             RMSE_val,
                             MAE_val,
                             stringsAsFactors = FALSE)
  
  # either combine with old results or write new
  if("hyperpar_opt_scores.csv" %in% list.files(model_folder)){
    existing_model_scores <- read.csv(paste0(model_folder, "/hyperpar_opt_scores.csv"))
    write.csv(rbind(existing_model_scores, model_scores),
              paste0(model_folder, "/hyperpar_opt_scores.csv"), row.names = FALSE)
  } else {
    write.csv(model_scores,
              paste0(model_folder, "/hyperpar_opt_scores.csv"), row.names = FALSE)
  }
  # save train and test prediction
  if(save_model_and_prediction){
    
    colnames(mean_pred_results_full_train) <- paste0("predicted_", y_variables)
    colnames(mean_pred_results_test)<- paste0("predicted_", y_variables)
    train_results <- data.frame(mean_pred_results_full_train, y_full_train) * 100
    test_results <- data.frame(mean_pred_results_test, y_test) * 100
    
    cat("Saving train_data and test_data prediction in", model_folder, "\n")
    write.csv(train_results, paste0(model_folder, "/train_prediction.csv"), 
              row.names = FALSE)
    write.csv(test_results, paste0(model_folder, "/test_prediction.csv"), 
              row.names = FALSE)
    
    
    # write model diagnostics
    # training performance
    RMSE_train <- matrix(NA, ncol = length(y_variables) + 1)
    colnames(RMSE_train) <- c(paste0("train_RMSE_", y_variables), "train_RMSE_mean")
    MAE_train <- matrix(NA, ncol = length(y_variables) + 1)
    colnames(MAE_train) <- c(paste0("train_MAE_", y_variables), "train_MAE_mean")
    for(var in 1:length(y_variables)){
      RMSE_train[1,var] <- RMSE(prediction = mean_pred_results_full_train[, var]* 100,
                                observation = y_full_train*100) %>% round(3)
      MAE_train[1, var] <- MAE(prediction = mean_pred_results_full_train[, var] * 100,
                               observation = y_full_train*100) %>% round(3)
    }
    RMSE_train[1,length(y_variables) + 1] <- round(mean(RMSE_train[1, 1:length(y_variables)]), 3)
    MAE_train[1,length(y_variables) + 1] <- round(mean(MAE_train[1, 1:length(y_variables)]), 3)
    # test performance (optional)
    if(!is.null(test_data)){
      # test scores
      RMSE_test <- matrix(NA, ncol = length(y_variables) + 1)
      colnames(RMSE_test) <- c(paste0("test_RMSE_", y_variables), "test_RMSE_mean")
      MAE_test <- matrix(NA, ncol = length(y_variables) + 1)
      colnames(MAE_test) <- c(paste0("test_MAE_", y_variables), "test_MAE_mean")
      for(var in 1:length(y_variables)){
        RMSE_test[1,var] <- RMSE(prediction = mean_pred_results_test[, var]* 100,
                                 observation = y_test*100) %>% round(3)
        MAE_test[1, var] <- MAE(prediction = mean_pred_results_test[, var] * 100,
                                observation = y_test*100) %>% round(3)
      }
      RMSE_test[1,length(y_variables) + 1] <- round(mean(RMSE_test[1, 1:length(y_variables)]), 3)
      MAE_test[1,length(y_variables) + 1] <- round(mean(MAE_test[1, 1:length(y_variables)]), 3)
    } else {
      RMSE_test <- matrix(NA, ncol = length(y_variables) + 1)
      MAE_test <- matrix(NA, ncol = length(y_variables) + 1)
    }
    # run time
    run_time <-  paste0(
      round((as.numeric(Sys.time()) - as.numeric(start_time))/60, 2),
      " minutes")
    # output data frame and save as csv
    model_scores <- data.frame(start_time = as.character(start_time),
                               run_time = run_time,
                               model_name = model_name,
                               RMSE_train,
                               MAE_train,
                               RMSE_val,
                               MAE_val,
                               RMSE_test,
                               MAE_test,
                               layers = layers,
                               units = units,
                               dropout = dropout,
                               batch_size = batch_size,
                               ensemble_runs = ensemble_runs,
                               stringsAsFactors = FALSE)
    
    if("model_scores.csv" %in% list.files(model_folder)){
      model_scores_old <- read.csv(paste0(model_folder, "/model_scores.csv"),
                                   stringsAsFactors = FALSE)
      write.csv(rbind(model_scores_old, model_scores, stringsAsFactors = FALSE),
                paste0(model_folder, "/model_scores.csv"),
                row.names = FALSE)
    } else {
      write.csv(model_scores, paste0(model_folder, "/model_scores.csv"), row.names = FALSE)
    }
    cat("\nModel training performance:", paste0("mean RMSE = ", model_scores$train_RMSE_mean, ","),
        paste0("mean MAE = ", model_scores$train_MAE_mean))
    cat("\nModel validation performance:", paste0("mean RMSE = ", model_scores$vak_RMSE_mean, ","),
        paste0("mean MAE = ", model_scores$val_MAE_mean))
    if(!is.null(test_data)){
      cat("\nModel testing performance:", paste0("mean RMSE = ", model_scores$test_RMSE_mean, ","),
          paste0("mean MAE = ", model_scores$test_MAE_mean), "\n")
    }
    cat("\nModel results are saved in", paste0(model_folder, "/model_scores.csv\n"))
    
  } else{
    cat("validation RMSE:", RMSE_val[, ncol(RMSE_val)], "\n")
    return(RMSE_val[, ncol(RMSE_val)])
  }
}


# create FNN model from hyperparameters
create_model <- function(nn_type, x_train, layers, units, dropout, seed, type = NULL,
                         timesteps = NULL, y_variables = y_variables){
  input <- keras::layer_input(shape = ncol(x_train))
  dropout_layers <- ifelse(dropout != 0, TRUE, FALSE)
  for(lay in 1:layers){
    # first time step
    if(lay == 1) {
      if(dropout_layers){
        output <- input %>%
          keras::layer_alpha_dropout(rate = dropout) %>%
          keras::layer_dense(units = units,
                             activation = 'selu',
                             kernel_initializer = "lecun_normal")
      } else {
        output <- input %>%
          keras::layer_dense(units = units,
                             activation = 'selu',
                             kernel_initializer = "lecun_normal")
      }
    } else {
      # all other time steps
      if(dropout_layers){
        output <- output %>% keras::layer_alpha_dropout(rate = dropout)
      }
      # add dense layers
      output <- output %>%
        keras::layer_dense(units = units,
                           activation = 'selu',
                           kernel_initializer = "lecun_normal")
    }
  }
  # add last layer
  output <- output %>% keras::layer_dense(length(y_variables))
  tensorflow::tf$random$set_seed(seed)
  model <- keras::keras_model(input, output)
  # compile
  model %>% compile(
    loss = "mse",
    optimizer = "adam"
  )
  return(model)
}

RMSE <- function(prediction, observation){
  round(sqrt(mean((prediction - observation)^2, na.rm = TRUE)), 3)
}

MAE <- function(prediction, observation){
  round(mean(abs(observation - prediction)), 3)
}
