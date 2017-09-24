# Before start, make sure you have the following packages installed: dismo, dplyr, raster, rJava;
# To run MaxEnt in R, the file 'maxent.jar' (Phillips et al., 2006) is required, which can be downloaded from http://biodiversityinformatics.amnh.org/open_source/maxent/ 
# Make sure that the 'maxent.jar' file is put in the java folder of the dismo package, which is the folder returned by system.file("java", package="dismo");
# Download the 'case1.rar' (zipped file) from https://geoterry.github.io/GEOWAT-SDGinsights/downloads.html

# Import the surface water point sample
# Load the surface water point sample file, assuming the unzipped 'case1' folder is in the D: drive.
LBR_surface_raw <- read.table("D:/case1/sample/lbr_wpt_surface.csv", header=TRUE, sep=",")
# Extract the columns which contain the coordinates of the water points (columns 2 and 3)
LBR_surface <- LBR_surface_raw[,2:3]
# Select random background points incorporated with the bias layer
LBR_bias <- raster("D:/case1/bias/lbr_bias.asc")
LBR_background <- xyFromCell(LBR_bias, sample(which(!is.na(values(LBR_bias))), 1000, prob=values(LBR_bias)[!is.na(values(LBR_bias))]))
colnames(LBR_background) = c('Long', 'Lat')
# Load the predictive covariates in 'LBR_layers'
library(dismo)
LBR_layers <- list.files("D:/case1/covariates/", full.names=TRUE )
# Create a RasterStack of the predictive covariates
LBR_covariates <- stack(LBR_layers)
# Setup training/test samples and run it 50 times for aggregated results
library(dplyr)
LBR_model_list <- list()
LBR_e <- list()
LBR_Prediction_list <- list()
for (i in 1:50) {
  # Set up training and test samples
  # Randomly select 70% of the input surface water points (about 41 points) as training sample;
  # Set the remainder (30%, about 18 points) as test sample.
  LBR_surface_train <- sample_n(LBR_surface, 41, replace = FALSE)
  LBR_surface_test <- LBR_surface[!rownames(LBR_surface) %in% rownames(LBR_surface_train),]
  # Run MaxEnt if the file 'maxent.jar' is properly located, otherwise skip.
  jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
  if (file.exists(jar)) {
    LBR_MaxEnt_model <- maxent(LBR_covariates, LBR_surface_train, factors='lbr_cov10_land_cover')
  } else {
    cat('Error: maxent.jar is needed for running this model')
    plot(1)
  }
  LBR_model_list[[i]] <- LBR_MaxEnt_model
  LBR_e[[i]] <- evaluate(LBR_surface_test, LBR_background, LBR_MaxEnt_model, LBR_covariates)
  LBR_Prediction_list[[i]] <- predict(LBR_covariates, LBR_MaxEnt_model, progress='')
}
# Check out the results
# Check covariate contribution of the ith model (e.g. change 'i' to 3 if want to see the 3rd model)
plot(LBR_model_list[[i]])
# Check response curve of the ith model
response(LBR_model_list[[i]])
# Map the ith prediction
plot(LBR_Prediction_list[[i]])
# Map the aggregated (average) prediction
plot(mean(LBR_Prediction_list[[i]]))
# Extract AUC
LBR_auc <- sapply(LBR_e, function(x){slot(x, 'auc')})
# Check AUC of each single run
LBR_auc
# Calculate the mean AUC
mean(LBR_auc)
# Extract the maximum of the sum of the sensitivity and specificy, which could be employed as a threshold for defining potentially suitable or not
sapply(LBR_e, function(x){x@t[which.max(x@TPR + x@TNR)]})
# Export the ith prediction raster
library(raster)
writeRaster(LBR_Prediction_list[[i]], 'LBR_Prediction_i.tif')
# Export the average prediction raster
LBR_prediction_avg <- mean(LBR_Prediction_list[[i]])
writeRaster(LBR_prediction_avg, 'LBR_Prediction_avg.tif')



