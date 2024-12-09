
---
title: "sample script"
author: "Joanna Pranga"
date: "year 2019"
source: "master thesis"
---



# Define required packages
required_packages <- c("caret", "dplyr", "foreach", "doParallel", "ranger", "tidyverse", "sf", "blockCV", "rgdal","mlr", "terra", "raster")

# Loop to check, install if necessary, and load packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# to measure running time of R code the tictoc package will be used
# devtools::install_github("collectivemedia/tictoc")


##############################################################################
############################# DATA PREPARATION ###############################
##############################################################################

# Define the input file geodatabase
fgdb <- "E:/Joanna/other_data" 

# Read the spatial data layer from the file geodatabase specified in 'fgdb'
kantons <- readOGR(dsn= fgdb, layer= "kantons_merged_LV03")

# Create SpatialPointsDataFrame object
pts_comb <- SpatialPointsDataFrame(df_comb[,c("POINT_X", "POINT_Y")], df_comb, proj4string=crs(kantons))

# Plot data on the map
plot(kantons) # plot raster data
points(pts_comb[which(pts_comb$NDWI_0_1==1), ], col="red", pch= 20) # add presence points
points(pts_comb[which(pts_comb$NDWI_0_1==0), ], col="darkgreen", pch= 20) # add absence points
legend(x=630000, y=230000, legend=c("Drought-affected","Control"), col=c("red", "darkgreen"), pch=c(20,20), bty="n")

# Remove 20% of the data from the set for the final VALIDATION of spatial predictions
set.seed(12345)
data_VALIDATION <- sample_n(df_comb, 18000)
table(data_VALIDATION$NDWI_0_1)
write.table(data_VALIDATION, file = ".../data_VALIDATION_18000.csv", sep = ";", dec = ".", na = "NA", row.names = FALSE) # save dataset

data_CALIBRATION <- df_comb[!(df_comb$ID %in% data_VALIDATION$ID),]
table(data_CALIBRATION$NDWI_0_1)
write.table(data_CALIBRATION, file = ".../data_CALIBRATION_18000.csv", sep = ";", dec = ".", na = "NA", row.names = FALSE) # save dataset


##############################################################################
############################ VARIABLE SELECTION ##############################
##############################################################################

# This is what was used to genereate a matrix of predictor-quality:
quad.formula <- function(sp, env) {
  formula(paste(sp,'~',env,'+ I(',env,'^2)'))
}

## the function to be used
spitTSS <- function(spec, pred, nrep, split.prop, data) {
  # Determine the total number of observations in the dataset
  nobs <- nrow(data)
  # Calculate the number of samples for the training set based on the split proportion
  nsamp <- floor(split.prop * nobs)
  samps <- replicate(nrep, sample(seq_len(nobs), nsamp))
  # Apply the TSS calculation for each sample (column in 'samps')
  tsss <- sapply(1:ncol(samps), function(x) {
    # Extract the indices for the current sample
    samp.ind <- samps[, x]
    # Create the training dataset using the sampled indices
    train <- data.frame(occ = data[samp.ind, spec], data[samp.ind, ])
    # Create the testing dataset using the remaining indices
    test <- data.frame(occ = data[-samp.ind, spec], data[-samp.ind, ])
    # Train a logistic regression model using the training dataset
    mod_train <- glm(formula = quad.formula(sp = 'occ', env = pred), data = train, family = 'binomial')
    # Generate predictions on the testing dataset
    preds <- predict(mod_train, newdata = test, type = 'response')
    # Create a data frame to store observed and predicted values for the test set
    df <- data.frame(ID = 1:length(preds), obs = test$occ, pred = preds)
    # Determine the optimal threshold for TSS calculation
    thr_tss <- optimal.thresholds(df, threshold = 1001, opt.methods = 3)[, 2]
    # Calculate presence-absence accuracy metrics
    paa <- presence.absence.accuracy(df, 
                                     threshold = thr_tss, 
                                     find.auc = TRUE, 
                                     st.dev = FALSE)
    # Compute the True Skill Statistic (TSS)
    tss <- (paa$sensitivity + paa$specificity) - 1
    tss
  })
  
  # Compute the mean TSS across all replicates
  tss_mean <- mean(tsss)
  tss_mean
}

# define the data that will be used for tss matrix creation
target <- c('NDWI_0_1')
names(data_CALIBRATION)  # your df is your data frame that contains all predictors and the target variable
non_predictors <- c(1:7, 943:944)  # define columns that are not predictors!


## apply the function to create the matrix
# set up cluster
clust <- makeCluster(36)
registerDoParallel(clust)

set.seed(6789)
tss_matrix <- foreach(i.spec = target, .packages = c('dismo', 'PresenceAbsence'), .combine = 'data.frame') %:%
  foreach(i.pred = names(data_CALIBRATION[,-non_predictors]), .packages = c('dismo', 'PresenceAbsence'), .combine = 'c') %dopar% {
    spitTSS(spec = i.spec, pred = i.pred, nrep = 10, split.prop = 0.7, data = data_CALIBRATION)
  }

# clean up cluster environment
stopCluster(clust)

## create the matrix
tss_matrix <- data.frame(matrix(tss_matrix, ncol = 1))
names(tss_matrix) <- target
row.names(tss_matrix) <- names(data_CALIBRATION[,-non_predictors])
write.table(tss_matrix, file = ".../tss_matrix_dataset_reduced.csv", sep = ";", dec = ".", na = "NA", row.names = FALSE)

## create ranks matrix
ranking <- t(apply(-tss_matrix, 2, rank))
write.table(ranking, file = ".../ranking_tss_dataset_reduced.csv", sep = ";", dec = ".", na = "NA", row.names = FALSE)


# The function to select variables from the ranking created in previous step
## the required function
select.thr <- function(Rank = NULL, X = NULL, threshold = 0.7, method = NULL){
  cm <- cor(X, method = method, use = 'complete.obs')
  # get highly correlated pairs
  pairs <- which(abs(cm) >= threshold, arr.ind = TRUE)
  pairs <- pairs[apply(pairs, 1, function(row) row[1] != row[2]),]
  # get predictor lists sorted by importance
  sortings <- lapply(row.names(Rank), function(row){
    colnames(Rank)[order(Rank[row,])]
  })
  names(sortings) <- row.names(Rank)
  # function to select predictors
  selpred <- function(sort.imp, cm){
    exclude <- NULL
    for (i in 1:length(sort.imp)){
      tst1 <- sort.imp[i] %in% row.names(pairs)
      tst2 <- !sort.imp[i] %in% exclude
      tst <- tst1 & tst2
      if (tst) {
        cv<-cm[setdiff(row.names(cm),exclude),sort.imp[i]]
        cv<-cv[setdiff(names(cv),sort.imp[1:i])]
        exclude<-c(exclude,names(which((abs(cv)>=threshold))))
      }
    }
    sel <- sort.imp[!(sort.imp %in% unique(exclude))]
    sel
  }
  sels <- lapply(sortings, selpred, cm = cm)
  sels
}


## do the variable selection
data_predictors <- data_CALIBRATION[,-non_predictors]

tic("select variables with Pearson")
sel_vars_Pearson <- select.thr(Rank = ranking, X = data_predictors, threshold = 0.7, method = 'pearson')
toc()

write.table(sel_vars_Pearson, file = ".../selected_variables_dataset_reduced.csv", sep = ";", dec = ".", na = "NA", row.names = FALSE)


##############################################################################
######################### SPATIAL CROSS-VALIDATION ###########################
##############################################################################

# prepare data for applying in blockCV package
raster <- raster("..../masked_raster.tif")

pa_data <- data_CALIBRATION[,c(5, 943, 944)]
pa_data <- SpatialPointsDataFrame(pa_data[,c("POINT_X", "POINT_Y")], pa_data, proj4string=crs(raster))
head(pa_data)
crs(pa_data)

# creating a spatial block 
set.seed(12345)
sb <- spatialBlock(speciesData = pa_data,
                   species = "NDWI_0_1",
                   rasterLayer = raster,
                   theRange = 20000, # set the size of the blocks to 20 km
                   k = 5,
                   selection = "random",
                   iteration = 200, # find evenly dispersed folds
                   numLimit = NULL, # If numLimit = NULL, the most evenly dispersed number of records is chosen (given the number of iteration).
                   border = kantons,
                   biomod2Format = TRUE)

saveRDS(sb, ".../spatial_block_dataset_reduced.rds")


#############################################################
####################### APPLY INNER LOOP ####################
#############################################################

hyper_grid <- expand.grid(
  numtrees = c(500, 1000, 1500),
  mtry       = c(4, 8, 12, 16),
  minnodesize  = c(1, 5, 10)
)

# to check fold IDs
foldID <- sb$foldID
folds <- sb$folds

## apply the function to create the matrix
# set up cluster

tic()
clust <- makeCluster(36)

doParallel::registerDoParallel(clust)

# here apply the inner loop computations
set.seed(12345)
inner_cv <- foreach(i_folds = 1:5)%:%
  foreach(i_subfolds = 1:4, .combine = 'cbind')%:%
  foreach(i_grid = 1:nrow(hyper_grid), .combine = 'c')%dopar%{
    ext <- 1:5
    levels <- ext[-i_folds]
    subfold <- levels[i_subfolds]
    # define training and testing dataset
    testset <- data_REDUCED[which(sb$foldID == subfold),] # !!!!always change the name of dataset used!!!!
    trainset <- data_REDUCED[which(!sb$foldID %in% c(subfold, i_folds)),] # !!!!always change the name of dataset used!!!!
    # set hyperparameters
    par_numtrees <- hyper_grid[i_grid,'numtrees']
    par_mtry <- hyper_grid[i_grid,'mtry']
    par_minnodesize <- hyper_grid[i_grid, 'minnodesize']
    # train a number of models for each inner loop
    mod <- ranger::ranger(
      formula = NDWI_0_1~., 
      data = trainset, 
      num.trees = par_numtrees, 
      mtry = par_mtry, 
      min.node.size = par_minnodesize, 
      probability = TRUE,
      importance = 'impurity',
      replace = TRUE)
    observed <- testset$NDWI_0_1
    predicted <- predict(mod, testset, type="response")
    predicted <- predicted$predictions[,2]
    pred <- ROCR::prediction(predicted, observed)
    auc <- ROCR::performance(pred, "auc")
    auc <- unlist(slot(auc,"y.values"))
    auc
    # paste(i_folds, subfold, i_grid)
    # runif(1, 0, 1)
  }
toc()
# clean up cluster environment
stopCluster(clust)

# calculate averages for each inner loop and each hypeparameter combination 
averages <- lapply(inner_cv, function(x) apply(x, 1, mean))
averages


# extract hyperparameper combination with the highest TSS measure
max_tss <- lapply(averages, which.max)
max_tss

# check if the function works well and selects proper rows from hyper_grid
lapply(max_tss, function(x) hyper_grid[x,])


######################################################
################## APPLY OUTER LOOP ##################
######################################################

# create empty vectors for evaluation results
AUC <- vector()

tic()
set.seed(12345)
#outer_cv <- function(){
for(i in 1:length(folds)){
  best_numtrees <- hyper_grid[max_tss[[i]],'numtrees']
  best_mtry <- hyper_grid[max_tss[[i]],'mtry']
  best_minnodesize <- hyper_grid[max_tss[[i]],'minnodesize']
  calibset <- unlist(folds[[i]][1])
  validset <- unlist(folds[[i]][2])
  # model fitting on training set
  fit <- ranger(
    formula = NDWI_0_1~.,
    data = data_REDUCED[calibset, ],
    num.trees = best_numtrees,
    mtry = best_mtry,
    min.node.size = best_minnodesize,
    probability = TRUE,
    importance = 'impurity',
    replace = TRUE,
    num.threads = 36,
    verbose =  TRUE)
  # save model results
  path1 <- file.path(".../results_data_reduced", paste0("model_fold_", i, ".rds"))
  saveRDS(object = fit, file = path1)
  # save prediction tables
  observed <- data_REDUCED[validset, 1]
  predicted <- predict(fit, data_REDUCED[validset, ], type = "response", num.threads = 36)
  predicted <- predicted$predictions[,2]
  pred <- ROCR::prediction(predicted, observed)
  auc <- ROCR::performance(pred, "auc")
  auc <- unlist(slot(auc,"y.values"))
  AUC[i] <- as.numeric(auc)
  ROC <- ROCR::performance(pred,"tpr","fpr")
  plot(ROC, main = "ROC curve", col="blue", lwd=2)
  abline(a=0, b= 1, lwd=1, lty=2)
  ACC = performance(pred, measure = "acc")
  plot(ACC, main = "Accuracy", col="blue", lwd=2)
  }
toc()
print(AUC)
print(mean(AUC))


model_1 <- readRDS(".../model_fold_1.rds")


###################################################################
# read environmental variables (GeoTiff fies)

# List all .tif files in the directory specified by path and store their full paths 
filesClimate <- list.files(path = pathClimate, pattern =".tif$", full.names=TRUE)
filesLidarSoil <- list.files(path = pathLidarSoil, pattern =".tif$", full.names=TRUE)
filesTerrain <- list.files(path = pathTerrain1, pattern =".tif$", full.names=TRUE)

# Stack all the .tif files from into a single RasterStack object 
predictorsClimate <- stack(filesClimate)
predictorsLidarSoil <- stack(filesLidarSoil)
predictorsTerrain1 <- stack(filesTerrain)

# Combine all three RasterStack objects (climate, lidar soil, and terrain data) into one RasterStack named 'predictorsALL'
predictorsALL <- stack(predictorsClimate, predictorsLidarSoil, predictorsTerrain)

predictorNames <- colnames(data_REDUCED[,-1])

# Subset the combined RasterStack 'predictorsALL' to retain only the layers that match the names in 'predictorNames'
predictorsSELECTED <- raster::subset(predictorsALL, predictorNames)




####################################################################
# crop rasters within the prediction function
library(raster)
tic()
beginCluster(36)
sp_predictions <- raster::predict(predictorsSELECTED, model_1, ext = kantons, progress='window', fun = function(model, ...) predict(model, ...)$predictions[,2])
endCluster()
toc()

writeRaster(sp_predictions1, filename= ".../predictions_probabilities_dataset.tif", format = "GTiff", overwrite=TRUE)

plot(sp_predictions)

################################################################################



##############################################################################
############################ MODEL PRFORMANCE ################################
##############################################################################

coords = data_CALIBRATION[,943:944]

# create a classification task
task = makeClassifTask(data = data_REDUCED, target = "NDWI_0_1", coordinates = coords, positive = "1")

# create learner -  select a machine learning method and its properties
# If predict.type='prob' we set 'probability=TRUE' in ranger.
lrn = makeLearner(cl = "classif.ranger", predict.type = "response", replace = TRUE)

# performance estimation level - outer resampling strategy (loop) - five spatially disjoint partitions with three repetitions
outer = makeResampleDesc(method = "SpCV", iters = 5)

# tuning estimation level - inner resampling strategy (loop) - four spatially disjoint partitions with one repetition
inner = makeResampleDesc(method = "SpCV", iters = 4)

# specify the parameter settings - the lower and upper limit of tuning parameters
ps = makeParamSet(
  makeDiscreteParam("num.trees", values= c(500, 1000, 1500)),
  makeDiscreteParam("mtry", values= c(4, 8, 12, 16)),
  makeDiscreteParam("min.node.size", values= c(1, 5, 10)))

# In the case of discrete ps specified above, since we have manually specified the values, grid search will simply be the cross product. 
ctrl_grid = makeTuneControlGrid()


# fuse learner with tuning
wrapped_lrn = makeTuneWrapper(learner = lrn,
                              # inner loop (tunning level)
                              resampling = inner,
                              # hyperparameter search space
                              par.set = ps,
                              # random search
                              control = ctrl_grid,
                              # print output on the console
                              show.info = TRUE,
                              # performance measure
                              measures = list(tnr, tpr, tp, tn, fp, fn, mmce, kappa))

# to make sure that the modeling goes on even if one model fails
configureMlr(on.learner.error = "warn", on.error.dump = TRUE)


######################################################################################################
# Fit models according to a resampling strategy
tic("nested spatial cross-validation")

#install.packages("parallelMap")
library(parallelMap)

if (Sys.info()["sysname"] == "Windows") {
  parallelStartSocket(level = "mlr.tuneParams",
                      cpus =  round(parallel::detectCores() / 1.1))
}


set.seed(12345)
result_response = mlr::resample(learner = wrapped_lrn,
                                task = task,
                                resampling = outer,
                                extract = getTuneResult,
                                measures = list(kappa, mmce, tpr, tnr, fpr, fnr, ppv, npv, tp, tn, fp, fn)) # added a number of new measures

# stop parallelization
parallelStop()

toc() 

saveRDS(result_response, ".../result_dataset_reduced_MLR_response.rds")





