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

# Example 1: Kriging with smoothing

In this example point groundwater head observations are spatially interpolated using co-variates of the land surface elevation (from the digital elevation model, DEM) and a local smoothing of the DEM. The magnitude of the smoothing is calibrated along with the variogram parameters. The example is for northern central Victoria, Australia,   

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

# Define the covariates for the kriging.
f <- as.formula('head ~ elev + smoothing')

# Define the base variogram model. Note, only the structure is used - not
# the valuaes.
variogram.model = gstat::vgm(psill=10, model='Mat', range= 10000 , nugget=1, kappa=0.1);

# Enforce a minimum error variance of 5cm ^2 
obs.data$total_err_var = pmax(obs.data$total_err_var, 0.05^2)

# Calibrate the mapping parameters with 25% of the data randomly selected and using 2 cores.
# NOTE 1: The rigor of the calibration is best controlled using the pop.size.multiplier input. 
# Here the size of the population of random guesses equals four time the number of calibration 
# parameters.
# NOTE 2: The 25% of random data are used to calculate the prediction error - which is minimised. 
# The remaining 74% of data is used to make the predicts.   
# NOTE 3: Here the smoothing parameter can between 0.5 and 1.5.   
calib.results <- krige.head.calib(formula=f, grid=DEM, data=obs.data, newdata=0.25, 
                  nmin=0, nmax=Inf, maxdist=Inf, omax=0, data.errvar.colname='total_err_var', model =         
                  variogram.model,  fit.variogram.type=1, smooth.std = c(0.5, 1.5),
                  pop.size.multiplier=4, debug.level=0, use.cluster = 2)

# Do the interpolation of the point data using the calibration results.  
# NOTE: All of the observed data is used for the calibration. All CPU cores are also used. 
head.grid <- krige.head(calibration.results = calib.results, data=obs.data, use.cluster = T)

# Map the head elevation and kriging uncertainty.
sp::spplot(head.grid, scales = list(draw = TRUE))

# Calculate the depth to water table.
# NOTE, this requires getting the DEM elevation into the head grids - event if there are a 
# different number of finite values.
sp::gridded(head.grid)=F 
sp::gridded(DEM)=F
head.grid$DBNS = DEM$DEM - head.grid$head
sp::gridded(head.grid)=T

# Map the depth to water table.
sp::spplot(head.grid,3, scales = list(draw = TRUE))
```