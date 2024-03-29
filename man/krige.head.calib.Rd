% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/krige.head.calib.R
\name{krige.head.calib}
\alias{krige.head.calib}
\title{Calibration of spatial interpolation parameters for mapping groundwater level observations.}
\usage{
krige.head.calib(
  formula = as.formula("head ~ elev + MrVBF + MrRTF + smoothing"),
  grid,
  grid.landtype.colname = NULL,
  data,
  newdata = 0.5,
  data.fixedHead = NULL,
  data.errvar.colname = NULL,
  model = c("Mat"),
  mrvbf.pslope = if (any(match(all.vars(as.formula(formula)), "MrVBF", nomatch = F) |
    match(all.vars(as.formula(formula)), "MrRTF", nomatch = F))) {
     seq(0.2, 2.5, by
    = 0.1)
 } else {
     NULL
 },
  mrvbf.ppctl = if (any(match(all.vars(as.formula(formula)), "MrVBF", nomatch = F) |
    match(all.vars(as.formula(formula)), "MrRTF", nomatch = F))) {
     seq(0.2, 2.5, by
    = 0.1)
 } else {
     NULL
 },
  smooth.std = seq(0.5, 2.5, length.out = 11),
  smooth.ncells = 11,
  nmax = if (is.character(data)) {
     -999
 } else {
     ceiling(seq(0.1, 0.2, 0.01) *
    length(data))
 },
  nmax.fixedHead = if (!is.null(data.fixedHead)) {
     seq(10, 110, length = 11)
 } else
    {
     NULL
 },
  maxdist = if (class(grid) == "SpatialPixelsDataFrame" || class(grid) ==
    "SpatialGridDataFrame") {
     ceiling(0.5 * sqrt((extent(grid)[2] -
    extent(grid)[1])^2 + (extent(grid)[4] - extent(grid)[3])^2) * seq(0.1, 1, 0.1))
 }
    else {
     -999
 },
  nmin = 0.2,
  omax = NULL,
  fit.variogram.type = 1,
  objFunc.type = 1,
  use.cluster = TRUE,
  pop.size.multiplier = 3,
  max.generations = 200,
  hard.generation.limit = FALSE,
  solution.tolerance = 1e-04,
  debug.level = 0
)
}
\arguments{
\item{\code{formula}}{see \code{\link{krige.head}}.}

\item{\code{MrVBF}}{see \code{\link{krige.head}}.}

\item{\code{MrRTF}}{see \code{\link{krige.head}}.}

\item{\code{smoothing}}{see \code{\link{krige.head}}.}

\item{\code{grid}}{see \code{\link{krige.head}}.}

\item{\code{grid.landtype.colname}}{see \code{\link{krige.head}}.}

\item{\code{data}}{see \code{\link{krige.head}}.#'}

\item{\code{data.fixedHead}}{see \code{\link{krige.head}}.}

\item{\code{newdata}}{is as for \code{data} but the points in a split-sample cross-validation scheme to estimate the interpolation error. Points listed within \code{newdata} should not be
listed within \code{data}. \code{newdata} can be (i) a real scalar >0 and <1  defining the fraction of data points within \code{data} to randomly remove and from \code{data} and use for cross-validation,
(i) an integer scalar >1 and <\code{length(data)} defining the number of data points within \code{data} to randomly remove and from \code{data} and use for cross-validation  (iii) a vector of indexes defining rows numbers within \code{data} to be extract and uses for the cross-validation; (iv) a vector of logicals
with \code{TRUE} defining rows numbers within \code{data} to be extract and used for the cross-validation, (v) a \code{SpatialPointsDataFrame} variable defining the complete data for cross validation - which must have identical columns to \code{data}.
The default is \code{newdata=0.5},}

\item{\code{data.errvar.colname}}{see \code{\link{krige.head}}.}

\item{\code{model}}{is either a character for the name of the variogram model type or a \code{gstat} variogram model object of type \code{variogramModel}. The available options are as per \code{gstat}, but
it is suggested to be \code{Mat}.}

\item{\code{mrvbf.pslope}}{defines the calibration type and range for the MrVFB shape parameter for the slope (see Gallant et al. 2003). A vector of two values defining the optimisation range when the parameter is
treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
not be optimised. If the \code{formula} includes either of the terms \code{MrVBF} or \code{MrRTF}, then the default is \code{seq(0.5, 1.5, length.out = 11)}. Else, the default is \code{NULL}.}

\item{\code{mrvbf.ppctl}}{defines the calibration type and range for the MrVFB shape parameter for elevation percentile (see Gallant et al. 2003). It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
not be optimised.  If the \code{formula} includes either of the terms \code{MrVBF} or \code{MrRTF}, then the default is \code{seq(0.5, 1.5, length.out = 11)}. Else, the default is \code{NULL}.}

\item{\code{smooth.std}}{defines the calibration type and range for the strength of the Gaussian kernel smoothing applied to the 5x5 grid cells surrounding each DEM grid cell. It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
not be optimised.  The default is \code{seq(0.5, 2.5, length.out = 11)}.}

\item{\code{smooth.ncells}}{defines the number of grid cells used (in each dimension) for Gaussian kernel smoothing. This input is current not able to be calibrated. Default is \code{11}.}

\item{\code{nmax}}{defines the calibration type and range for the maximum number of \code{data} observations to use when estimate each point using local kriging. It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
not be optimised.  The default is \code{ceiling(seq(0.1,0.20,0.01)*length(data))}.}

\item{\code{nmax.fixedHead}}{defines the calibration type and range for the maximum number of \code{data.fixedHead} observations to use when estimate each point using local kriging. It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
not be optimised.  The default is \code{seq(10,110,length=11)}.}

\item{\code{maxdist}}{defines the calibration type and range for the maximum search radius to use when estimate each point using local kriging. It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
not be optimised.  The default is from 10\% to 100\% of \code{grid} extend at increments of 10\%. If \code{grid} is \code{NULL}, then the user must input the search radius in one of the three accepted forms.}

\item{\code{nmin}}{see \code{\link{krige.head}}.}

\item{\code{omax}}{see \code{\link{krige.head}}.}

\item{\code{nsim}}{see \code{\link{krige.head}}.}

\item{\code{fit.variogram.type}}{defines the way the model variogram is to be derived/used.
For \code{fit.variogram.type==1} the input \code{model} then all of the variogram parameters will be calibrated. When \code{model} is a chacater for a variogram modle type, then
the variogram will be assumed isotropic. When \code{model} is a \code{gstat} variogram model object of type \code{variogramModel}, then all of the variogram parameters will be calibrated.
This approach allows significantly greater control of the variogram model form.
For \code{fit.variogram.type==2} the input \code{model} must be a \code{gstat} variogram model object of type \code{variogramModel}. The variogram model parameters will not be calibrated and instead
will be estimated by fitting the model variogram to an experimental variogram using multi-start local calibration. This approach is very similar to convential variogram model fitting.
For \code{fit.variogram.type==3} the input \code{model} must also be a variogram model object of type \code{variogramModel}. The variogram model parameters will not be calibrated or fit to an experimental variogram.
The default is \code{fit.variogram.type=1}.}

\item{\code{objFunc.type}}{defines the type of objective function to use in the calibration. For \code{objFunc.type==1}, the negative likelihood divided by the number of \code{newdata} observations is used
- which accounts for the expected error in the \code{newdata} estimates. That is, \code{newdata} points far from other \code{data} points will have a large kriging variance and hence the expected
error in the kriging estimate is larger. Conversely, \code{newdata} points near to \code{data}  points will have a low kriging variance and hence the kriging water level error should be low.
For details see Samper and Neumen (1989A,B) and Pardo-Iguzquiza & Dowd (2013). For \code{objFunc.type==2}, it is as for type 1 but with an added penalty when the estimate violates a physical constraint.
by first estimating the difference between the land surface elevation and the head. At those \code{newdata} points above the land surface, the absolute of this difference is added to the absolute error between
the predicated value and that within \code{newdata}. For \code{objFunc.type==3}, the root-mean-square error is used. For \code{objFunc.type==4}, the root-mean-square error is used but with the penalty from type 2.
The default is \code{objFunc.type=1}.}

\item{\code{pop.size.multiplier}}{defines the size of the calibration population, expressed as a multiplier for the number of parameters. Default =3.}

\item{\code{max.generations}}{defines the maximum number of calibration evolution loops. Default = 200.}

\item{\code{hard.generation.limit}}{defines if \code{max.generations} is a hard or soft limit (default). Default is \code{FALSE}.}

\item{\code{solution.tolerance}}{defines the calibration tolerance for convergence. Default is \code{1e-4}.}

\item{\code{use.cluster}}{see \code{\link{krige.head}}.}

\item{\code{debug.level}}{see \code{\link{krige.head}}.}
}
\value{
As per \code{genoud} plus the parameter values in the transformed calibration space "par.pretransformed" and "par" is replaced as a named list variable.
}
\description{
\code{krige.head.calib} calibrates parameters to minimise the interpolation error.
}
\details{
This function optimises the parameters for the spatial interpolation so as to minimise the
interpolation error, which by default is defined using a formal likelihood function. Once this function has return
the optimal parameter set, the mapping can be undertaken using \code{\link{krige.head}} (see examples below).

The mapping parameters are estimated using a mixed data-type (i.e. real and integer parameters) split-sample maximum likelihood global optimisation. The optimisation by default
includes the variogram parameters (e.g. range, sill and nuggest) and the search parameters for local kriging (e.g. radius, minimum and
maximum number of observations to use). Optimising these parameters is not common in kriging. It is done herein because trials for Victoria,
Australia, showed that calibrating these parameters produced significantly lower cross-validation errors (i.e. the error in predicting the observations
removed from the optimisation) compared to the standard approach using fitting to an experimental variogram. The optimisation is
numerically challenging and the following factors should be considered before use:

\itemize{
 \item{Optimisation of the parameters \code{mrvbf.pslope}, \code{mrvbf.ppctl} and \code{smooth.std} often required the creating of raster grids for every parameter combination. To ease the computation burden, these parameters should be treated as disrcete, not continuous, numbers.}
 \item{The optimisation package \code{rgeoud} is used herein. See the \code{rgeoud} documentation for details of the optimisation scheme.}
 \item{Trials have established default calibration parameters and settings that were effective for Victoria, Australia. There is no guarantee they will be effective for other regions.}
 \item{The available input point data must be split up into \code{data} and \code{newdata}. The point observations within \code{data} are used to estimate the water level at the locations defined
 within \code{newdata}. The difference between the observed and estimated values defines the interpolation prediction error. This process also gives the kriging variance at
 each \code{newdata} location. The points used for \code{data} and \code{newdata} should both cover all types of terrain, geology, landuse and the full mapping extend so as to avoid bias. Also, changing either may
 change the calibration solution.}
 \item{The \code{objFunc.type} allows for maximum likelihood estimation - that is optimisation that accounts for the expected kriging error. See \code{objFunc.type} below for details.}
}

In using this function, the primary user decisions are:
\itemize{
 \item{The kriging with external drift formula defining the independent gridded variables deemed to predict the groundwater elevation. See the input \code{formula}.}
 \item{The type of the objective function, defined by the input \code{objFunc.type}}
 \item{The type of variogram model, defined by the input \code{model}}
 }
}
\examples{

# Set enironment path for hydroMap
set.env()

# Load water table observations from  April 2000 for Victoria, Australia and a 250m state-wide DEM.
data('victoria.groundwater')

# Load a model variogram and mapping parameters found to be effective.
data('mapping.parameters')

# Define a simple kriging formula without MrVBF terms that does not require the package RSAGA.
f <- as.formula('head ~ elev + smoothing')

# Define an initial isotropic variogram model. All of the parameters will be calibrated.
varigram.model <- gstat::vgm(psill=25, model="Mat", range=5000, nugget=5, kappa = 0.5)

# Calibrate the mapping parameters with 25\% of the data randomly selected at 2 cores.
calib.results <- krige.head.calib(formula=f, grid=DEM, data=obs.data, newdata=0.25, nmin=0, nmax=Inf, maxdist=Inf, omax=0,
                      data.errvar.colname='total_err_var', model = variogram.model,  fit.variogram.type=1,
                      pop.size.multiplier=1, debug.level=0, use.cluster = 2)

# Grid the observed head using ONLY the training data from the calibration (and all cores) and then map..
head.grid.training <- krige.head(calibration.results = calib.results, use.cluster = T)
spplot(head.grid.training)

# Grid the observed head using ALL of the obsrved data from the calibration (and all cores).
head.grid.all <- krige.head(calibration.results = calib.results, data=obs.data, use.cluster = T)
spplot(head.grid.all)

}
\references{
Gallant, J.C., Dowling, T.I. (2003) A multiresolution index of valley bottom flatness for mapping depositional areas, Water Resources Research, 39/12:1347-1359

Rivoirard, J. & Romary, T. Math Geosci (2011) Continuity for Kriging with Moving Neighborhood, Mathematical Geosciences, 43: 469. DOI: 10.1007/s11004-011-9330-0

Pardo-Igúzquiza E. , Dowd P. A., (2013) Comparison of inference methods for estimating semivariogram model parameters and their uncertainty: The case of small data sets,
Computers & Geosciences, v50,pp 154-164, https://doi.org/10.1016/j.cageo.2012.06.002.

Samper, F. J., and S. P. Neuman (1989A), estimation of spatial covariance structures by adjoint state maximum likelihood cross validation: 1. Theory,
Water Resour. Res., 25(3), 351–362, doi: 10.1029/WR025i003p00351.

Samper, F. J., and S. P. Neuman (1989B), Estimation of spatial covariance structures by adjoint state maximum likelihood cross validation: 2. Synthetic experiments,
Water Resour. Res., 25(3), 363–371, doi: 10.1029/WR025i003p00363.
}
\seealso{
\code{\link{krige.head}} for undertaking the mapping.
}
