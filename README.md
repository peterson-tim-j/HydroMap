# HydroMap

_HydroMap_ is a powerful R package for kriging interpolation of sparse groundwater level observatons - that is, groundwater potentiometry. Key features include:
 * land surface elevation, topographic form (e.g. valleys and ridges) and the smoothness of the groundwater head (relative to the terrain) is account for in the mapping.
 * catagorical land types such as geology can be included in the mapping.
 * fixed head boundary conditions, such as the ocean, can be included and importantly they only have an influence if there are no nearby observation points.
 * all mapping parameters can be calibrated using split-sample maximum likelihood estimation.
 * kriging variogram parameters can be calibrated - which trials have shown to significantly reduced the prediction error.

Importantly, this is a *beta* release. There is basic documentation (see pdf) and the code has been tested. In the future, additional functions will be made accessable to the user and documented.

To illustrate the perforamnce of _HydroMap_, Fig. 1 shows the interpolated watertable depth for Victoria, Australia, at April 2000. It was derived using the conventional kriging approach of the DEM elevation being the only predictor of the head. It shows the water table to be very soooth - which somewhat counter intuitively arises because the interpolated water level elevation is very noisy. That is, whenever the DEM elevation rises so does the heads. Fig. 2 shows that when the above _HydroMap_ features are included in the mapping then the heads become smooth and, consequently, the depth to the water table becomes "noisy" with significantly greater variability shown in the mountainous regions of the east.

![Convetional kriging approach](https://user-images.githubusercontent.com/8623994/44770420-57776580-abab-11e8-9b95-ff54604ba6e3.png)
Figure 1. Depth to water table for Victoria, Australia, derived only the DEM elevation as a predictor of groundwater head.

![HydroMap approach](https://user-images.githubusercontent.com/8623994/44770783-79bdb300-abac-11e8-9404-d0d7a4b4e9f4.png)
Figure 2. Depth to water table for Victoria, Australia, derived the key featues of _HydroMap_.

# Getting Started

To get started using the package, following these steps:

1. Download the zipped package.
1. Unzip the package on your local machine.
1. Open R
1. Install the required packages using the following R command: install(c("sp", "grid", "gstat","raster", "parallel", "rgenoud", "devtools","RSAGA"))
1. Load the required R packages using the following R commands: library(sp); library(grid); library(gstat); library(raster); library(parallel); library(rgenoud); library(RSAGA);library(devtools)
1. Within R, navigate to where you unzipped the HydroMap package.
1. Install the HydroMap package using the following command: install("HydroMap")
1. Open the help documentation for the main function: ?krige.head
