# HydroMap

_HydroMap_ is a powerful R package for kriging interpolation of sparse groundwater level observations - that is, groundwater potentiometry. Key features include:
 * land surface elevation, topographic form (e.g. valleys and ridges) and the smoothness of the groundwater head (relative to the terrain) is account for in the mapping.
 * categorical land types such as geology can be included in the mapping.
 * fixed head boundary conditions, such as the ocean, can be included and importantly they only have an influence if there are no nearby observation points.
 * all mapping parameters can be calibrated using split-sample maximum likelihood estimation.
 * kriging variogram parameters can be calibrated - which trials have shown to significantly reduced the prediction error.

Importantly, this is a **beta** release. There is basic documentation (see [PDF MANUAL](https://github.com/peterson-tim-j/HydroMap/blob/master/hydroMap.pdf)) and the code has been tested. In the future, additional functions will be made accessible to the user and documented.

To illustrate the performance of _HydroMap_, Fig. 1 shows the interpolated watertable depth for Victoria, Australia, at April 2000. It was derived using the conventional kriging approach of the DEM elevation being the only predictor of the head. It shows the water table to be very smooth - which somewhat counter intuitively arises because the interpolated water level elevation is very noisy. That is, whenever the DEM elevation rises so does the heads. Fig. 2 shows that when the above _HydroMap_ features are included in the mapping then the heads become smooth and, consequently, the depth to the water table becomes "noisy" with significantly greater variability shown in the mountainous regions of the east.

![Conventional kriging approach](https://user-images.githubusercontent.com/8623994/44770420-57776580-abab-11e8-9b95-ff54604ba6e3.png)
Figure 1. Depth to water table for Victoria, Australia, derived only the DEM elevation as a predictor of groundwater head.

![HydroMap approach](https://user-images.githubusercontent.com/8623994/44770783-79bdb300-abac-11e8-9404-d0d7a4b4e9f4.png)
Figure 2. Depth to water table for Victoria, Australia, derived the key features of _HydroMap_.

# Installation

To install the package, following these steps:

1. Open R and Install the required packages using the following R command: `install.packages(c("sp", "grid", "gstat","raster", "parallel", "rgenoud", "snow","RSAGA"))`
1. Download the ![latest release](https://github.com/peterson-tim-j/HydroMap/releases).
1. Install _HydroMap_ using the R command: `install.packages('HydroMap.tar.gz',repos = NULL)`
1. Load _HydroMap_ using the R commands: `library("HydroMap")`
1. Open the help documentation using `?krige.head` and follow the example.  

# Examples: Data setup
The examples below all use the same gridded data and point data. Below the data is imported and cropped to central northern Victoria. Also, here the calibration training and prediction data are pre-defined. This is done to ensure that the results from each of the following examples are comparable.

```R
library(RSAGA)

# Setup RSAGA with the paths to the requuired modules. Note you will to do this yourself for your own
# installation of SAGA.
set.env(saga.path = 'C:/Program Files (x86)/saga-9.0.1_x64',saga.modules = 'C:/Program Files (x86)/saga-9.0.1_x64/tools')

# Load water table observations from  April 2000 for Victoria, Australia and a 250m state-wide DEM.
data('victoria.groundwater')

# Crop this stat-ewide DEM and data  points to a small in the centre north. 
DEM <- raster::crop(raster::raster(DEM), raster::extent(2400000, 2500000, 2550000, 2650000))
DEM = as(DEM,'SpatialGridDataFrame')
obs.data <- raster::crop(obs.data, DEM, inverse = F) 

# Load a model variogram and mapping parametyers found to be effective.
data('mapping.parameters')

# Calculate the depth to water table (DTWT)
obs.data$DTWT = obs.data$elev - obs.data$head

# Convert DTWT to catagories, to aid mapping
obs.data$DTWT.cats =cut(obs.data$DTWT,breaks=c(-Inf,0,2,5,10,25,50, Inf ), 
    labels=c('<0m','0-2m','2-5m','5-10m','10-25m','25-50m','>50m'),include.lowest=T)

# Enforce a minimum error variance of 5cm ^2 for the groundwater head elevation. 
obs.data$total_err_var = pmax(obs.data$total_err_var, 0.05^2)

# Define the prediction data by randomly sample 25% of the observed data points. 
# The remaining 75% of data points are used for the presictions.   
nObs = nrow(obs.data);
nObs.prediction = floor(0.25*nObs)
set.seed(123456, sample.kind='default')
predictionData.index = sample(1:nObs, nObs.prediction, replace=F)
predictionData = obs.data[predictionData.index,]
trainingData = obs.data[-predictionData.index,]

# Plot the depth to water table of the prediction and training data
png('Example1_pointData_A.png')
sp::spplot(predictionData, 'DTWT.cats',scales = list(draw = TRUE), main='Prediction data DTWT [m]')
dev.off()
png('Example1_pointData_B.png')
sp::spplot(trainingData, 'DTWT.cats',scales = list(draw = TRUE), main='Training data DTWT [m]')
dev.off()

# Define the variogram model. This can be either a full variogram object (see `gstat:vgm`),
# which is then calibrated, or simply a type of variogram, as used here. When the type is input then 
# the initial variogram parameters are estimated using the residuals of the observed head and the  
# covariates, using ordinary least squares. Latter calibration of these parameters uses a range of 0.1
# and 10 times the initial estimates.
variogram.model = 'Mat';
```
# Example 1: Kriging with only land surface elevation 

This is the simplest example. The point groundwater head observations are spatially interpolated using only the land surface elevation (from the digital elevation model, DEM). Only the variogram parameters are calibraed to minimise the prediction error. 
```R
# Define the covariates for the kriging.
f <- as.formula('head ~ elev')

# Calibrate the mapping parameters with 25% of the data randomly selected and using 2 cores.
# NOTE 1: The rigor of the calibration is best controlled using the pop.size.multiplier input. 
# Here the size of the population of random guesses equals four time the number of calibration 
# parameters.
calib.results.example1 <- krige.head.calib(formula=f, grid=DEM, data=trainingData, newdata=predictionData, 
                  nmin=0, nmax=Inf, maxdist=Inf, omax=0, data.errvar.colname='total_err_var', model =         
                  variogram.model,  fit.variogram.type=1,  smooth.std =NA,
                  pop.size.multiplier=4, debug.level=0, use.cluster = 2)

# Do the interpolation of the point data using the calibration results.  
# NOTE: All of the observed data is used for the calibration. All CPU cores are also used. 
head.grid.example1 <- krige.head(calibration.results = calib.results.example1, data=obs.data, use.cluster = T)

# Map the head elevation and kriging uncertainty.
png('Example1_head.png')
raster::plot(raster::raster(head.grid.example1),'head')
raster::contour(raster::raster(head.grid.example1,1), levels = seq(70,125,by=5), add=T)
dev.off()

# Categorise the DTWT to seven classes.
head.grid.example1$DTWT.cats =cut(head.grid.example1$DTWT,breaks=c(-Inf,0,2,5,10,25,50, Inf ), 
    labels=c('<0m','0-2m','2-5m','5-10m','10-25m','25-50m','>50m'),include.lowest=T)

# Map the categorised depth to water table.
png('Example1_DTWT.png')
sp::spplot(head.grid.example1,'DTWT.cats', scales = list(draw = TRUE))
dev.off()
```



# Example 2: Kriging with land surface elevation and smoothing

In this example point groundwater head observations are spatially interpolated using co-variates of the land surface elevation (from the digital elevation model, DEM) and a local smoothing of the DEM. The magnitude of the smoothing is calibrated along with the variogram parameters. Note, the only difference to example 1 is the formula and the parameters being calibrated (i.e. `smooth.std`)  

```R
# Define the covariates for the kriging.
f <- as.formula('head ~ elev + smoothing')

# Calibrate the mapping parameters with 25% of the data randomly selected and using 2 cores.
# NOTE: Here the smoothing parameter can between 0.5 and 2.5.   
calib.results.example2 <- krige.head.calib(formula=f, grid=DEM, data=trainingData, newdata=predictionData, 
                  nmin=0, nmax=Inf, maxdist=Inf, omax=0, data.errvar.colname='total_err_var', model =         
                  variogram.model,  fit.variogram.type=1, smooth.std = c(0.5, 5.0),
                  pop.size.multiplier=4, debug.level=0, use.cluster = 2)

# Do the interpolation of the point data using the calibration results.  
# NOTE: All of the observed data is used for the calibration. All CPU cores are also used. 
head.grid.example2 <- krige.head(calibration.results = calib.results.example2, data=obs.data, use.cluster = T)

# Map the head elevation and kriging uncertainty.
png('Example2_head.png')
raster::plot(raster::raster(head.grid.example2),'head')
raster::contour(raster::raster(head.grid.example2,1), levels = seq(70,125,by=5), add=T)
dev.off()

# Categorise the DTWT to seven classes.
head.grid.example2$DTWT.cats =cut(head.grid.example2$DTWT,breaks=c(-Inf,0,2,5,10,25,50, Inf ), 
    labels=c('<0m','0-2m','2-5m','5-10m','10-25m','25-50m','>50m'),include.lowest=T)

# Map the categorised depth to water table.
png('Example2_DTWT.png')
sp::spplot(head.grid.example2,'DTWT.cats', scales = list(draw = TRUE))
dev.off()
```

# Example 3: Kriging with land surface elevation, smoothing and MrVBF and MrRTF

In this example point groundwater head observations are spatially interpolated using co-variates of the land surface elevation (from the digital elevation model, DEM), MrVBF measure of the valley floor and a local smoothing of the DEM. The magnitude of the smoothing is calibrated along with the variogram parameters. 

```R
# Define the covariates for the kriging.
f <- as.formula('head ~ elev + smoothing + log(MrVBF) + log(MrRTF)')

# Calibrate the mapping parameters with 25% of the data randomly selected and using 2 cores.
calib.results.example3 <- krige.head.calib(formula=f, grid=DEM, data=trainingData, newdata=predictionData, 
                  nmin=0, nmax=Inf, maxdist=Inf, omax=0, data.errvar.colname='total_err_var', model =         
                  variogram.model,  fit.variogram.type=1, smooth.std = c(0.5, 5.0),
                  pop.size.multiplier=4, debug.level=0, use.cluster = 2)

# Do the interpolation of the point data using the calibration results.  
# NOTE: All of the observed data is used for the calibration. All CPU cores are also used. 
head.grid.example3 <- krige.head(calibration.results = calib.results.example3, data=obs.data, use.cluster = T)

# Map the head elevation and kriging uncertainty.
png('Example3_head.png')
raster::plot(raster::raster(head.grid.example3),'head')
raster::contour(raster::raster(head.grid.example3,1), levels = seq(70,125,by=5), add=T)
dev.off()

# Categorise the DTWT to seven classes.
head.grid.example3$DTWT.cats =cut(head.grid.example3$DTWT,breaks=c(-Inf,0,2,5,10,25,50, Inf ), 
    labels=c('<0m','0-2m','2-5m','5-10m','10-25m','25-50m','>50m'),include.lowest=T)

# Map the categorised depth to water table.
png('Example3_DTWT.png')
sp::spplot(head.grid.example3,'DTWT.cats', scales = list(draw = TRUE))
dev.off()
```