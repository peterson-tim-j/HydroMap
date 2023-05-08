#' Spatially interpolates sparse groundwater level observations, and if required, estimate mapping parameters.
#'
#' \code{krige.head} creates a groundwater level elevation map from sparse point observations for one time point.
#'
#' @details This function is the primary means for using the package. It interpolates sparse input groundwater elevation
#' observations for a single time point using a form of kriging with external drift that can account
#' for land surface elevation, topographic form (e.g. valleys and ridges), the smoothness of the groundwater
#' relative to the land surface and remote sensing gridded data. Also, a co-kriging features allows for
#' the inclusion of categorical land types and fixed head boundary conditions, such as the ocean. Each of
#' these features can be individually controlled by the user.
#'
#' Importantly, if the mapping parameters are not specified by the user, then this function estimates the parameters using
#' a mixed data-type (i.e. real and integer parameters) split-sample maximum likelihood global optimisation. The optimisation by default
#' includes the variogram parameters (e.g. range, sill and nugget) and the search parameters for local kriging (e.g. radius, minimum and
#' maximum number of observations to use). Optimising these parameters is not common in kriging. It is done herein because trials for Victoria,
#' Australia, showed that calibrating these parameters produced significantly lower cross-validation errors (i.e. the error in predicting the observations
#' removed from the optimisation) compared to the standard approach of graphical estimation from an experimental variogram. The optimisation is
#' numerically challenging and the following factors should be considered before use:
#'
#' \itemize{
#'  \item{Optimisation of the parameters \code{mrvbf.pslope}, \code{mrvbf.ppctl} and \code{smooth.std} often required the creating of raster grids for every parameter combination. To ease the computation burden, these parameters should be treated as discrete, not continuous, numbers.}
#'  \item{The optimisation package \code{rgeoud} is used herein. For control the optimisation process, consider directly using \code{\link{krige.head.calib}}.}
#'  \item{Trials have established default calibration parameters and settings that were effective for Victoria, Australia. There is no guarantee they will be effective for other regions.}
#' }
#'
#' In using this function, the primary user decisions are:
#' \itemize{
#'  \item{The kriging with external drift formula defining the independent gridded variables deemed to predict the groundwater elevation. See the input \code{formula}.}
#'  \item{The mapping extent and resolution, defined by the input \code{grid}, and the point observations of groundwater elevation, defined by the input \code{data}.}
#'  \item{The type of variogram model, defined by the input \code{model}}
#'  }
#'
#' @param \code{formula} defines the R formula (as a character or formula data type) to be used to interpolate the heads. The left hand side of the formula must be \code{head}.
#' The right hand side can contain any or all of the following terms: \code{elev} for the land surface elevation;
#' \code{MrVBF} for the Multiresolution Index of Valley Bottom Flatness as a measure of valley-ness at each DEM grd cell;
#' \code{MrRTF} for the Multiresolution Index of Ridge Top Flatness as a measure of ridge to plateaus at each DEM grd cell;
#' \code{smoothing} for a local smoothing factor derived from the DEM roughness. For any terms other than the prior, the data
#' for the variable must be listed within the inputs \code{grid} and \code{data}. The default is \code{as.formula("head ~ elev + MrVBF + MrRTF + smoothing")}.
#'
#' @param \code{grid} is either a character string to a ASCII grid digital elevation model (DEM) or a \code{SpatialPixelsDataFrame} or \code{SpatialGridDataFrame} containing the land
#' surface elevation, which must be named \code{elev}, and each formula variable other than \code{MrVBF}, \code{MrRTF} and \code{smoothing}; which are each derived from the DEM.
#'
#' @param \code{grid.landtype.colname} is a character of the column name within \code{grid} and \code{data} that define the land category. The land category data should be an integer and
#' be of a sufficiently small number of categories that multiple data points exists within each land category. If \code{NULL}, then land categories are not accounted for in the mapping. The default is \code{NULL}.
#'
#' @param \code{data} is either a character string to a .csv file of point data or a \code{SpatialPointsDataFrame} containing the columns \code{Easting}, \code{Northing} and \code{head} or \code{depth}.
#' If the formula includes the term \code{elev}, then the bore elevation should be included in \code{data} so as to account for difference between the bore and DEM elevation.
#' Each formula right hand side variable other than \code{MrVBF}, \code{MrRTF} and \code{smoothing} should also be provided within \code{data}.
#'
#' @param \code{data.fixedHead} is as for \code{data} but the points are treated as fixed head points within a cokriging approach. It can be used to guide the head estimates
#' toward zero along, say, the coastline. The fixed points will have a greater influence when no observation data is nearby. If \code{NULL}, then no fixed head points are used.
#' The default is \code{NULL}.
#'
#' @param \code{newdata} is as for \code{data} but the points are used in a split-sample cross-validation scheme to estimate the interpolation error. Points listed within \code{newdata} should not be
#' listed within \code{data}. If \code{is.null(newdata)==TRUE}, then \code{grid} should be \code{NULL}.
#'
#' @param \code{data.errvar.colname} is a character of the column name within \code{data} that define measurement error, as a variance. If \code{NULL}, then measurement error is not accounted for.
#' The default is \code{NULL}.
#'
#' @param \code{model} is either a character for the name of the variogram model type or a \code{gstat} variogram model object of type \code{variogramModel}. The available options are as per \code{gstat}, but
#' it is suggested to be \code{Mat}.
#'
#' @param  \code{mrvbf.pslope} defines the MrVFB shape parameter for the slope (see Gallant et al. 2003), a vector of two values defines the optimisation range when the parameter is
#' treated as a real number. A vector of length >2 values defines the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
#' not be optimised. If the \code{formula} includes either of the terms \code{MrVBF} or \code{MrRTF}, then the default is \code{seq(0.5, 1.5, length.out = 11)}. Else, the default is \code{NULL}.
#'
#' @param  \code{mrvbf.ppctl} defines the MrVFB shape parameter for elevation percentile (see Gallant et al. 2003). It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
#' treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
#' not be optimised.  If the \code{formula} includes either of the terms \code{MrVBF} or \code{MrRTF}, then the default is \code{seq(0.5, 1.5, length.out = 11)}. Else, the default is \code{NULL}.
#'
#' @param  \code{smooth.std} defines the strength of the Gaussian kernal smoothing applied to the 5x5 grid cells surrounding each DEM grid cell. It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
#' treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
#' not be optimised.  The default is \code{seq(0.5, 1.5, length.out = 11)}.
#'
#' @param \code{nmax} defines the maximum number of \code{data} observations to use when estimate each point using local kriging. It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
#' treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
#' not be optimised.  The default is \code{ceiling(seq(0.1,0.20,0.01)*length(data))}.
#'
#' @param \code{nmax.fixedHead} defines the maximum number of \code{data.fixedHead} observations to use when estimate each point using local kriging. It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
#' treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
#' not be optimised.  The default is \code{seq(10,110,length=11)}.
#'
#' @param \code{maxdist} defines the maximum search radius to use when estimate each point using local kriging. It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
#' treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
#' not be optimised.  The default is from 10\% to 100\% of \code{grid} extend at increments of 10\%. If \code{grid} is \code{NULL}, then the user must input the search radius in one of the three accepted forms.
#'
#' @param \code{nmin} defines the minimum number of \code{data} observations to use when estimating each point using local kriging. If the \code{nmin} observations cannot be located within the search radius \code{nmax}, then the search
#' radius is increased until \code{nmin} points are obtained. If \code{nmin} is between zero and one and \code{omax} is between zero and one, then \code{nmin} is treated as a fraction of \code{nmax}. Else, \code{nmax} is treated as an
#' integer number of data points and hence must be >0. This input cannot be optimised. The default value is 0.2.
#'
#' @param \code{omax} defines the maximum number of \code{data} observations to select per quadrant when estimating each point using local kriging. It must be either \code{NULL} (to not use quadrant search) or between zero and one,
#' which is treated as a fraction of \code{nmax}. This input cannot be optimised.  The default value is \code{NULL}.
#'
#' @param \code{nsim} defines the number of conditional simulations to undertake. If set to a non-zero value, conditional simulation is used instead of kriging interpolation. Importantly, this feature has not been tested.
#' The default value is 0.
#'
#' @param \code{fit.variogram.type} defines the way the model variogram is to be derived/used.
#' For \code{fit.variogram.type==1} the input \code{model} must be the variogram type (as character string) and optimisation must be undertaken. The variogram will be assumed isotropic. For more control of the variogram calibrate using \code{\link{krige.head.calib}}.
#' For \code{fit.variogram.type==2} the input \code{model} must be a \code{gstat} variogram model object of type \code{variogramModel}. The variogram model parameters will be estimated by fitting the model variogram to an experimental variogram using multi-start local calibration.
#' For \code{fit.variogram.type==3} the input \code{model} must also be a variogram model object of type \code{variogramModel}. If calibration is being undertaken, then the variogram model parameters will not be optimised or fit to an experimental variogram.
#'
#' @param \code{objFunc.type} defines the type of objective function to use in the optimisation. See \code{\link{krige.head.calib}} for details.
#'
#' @param \code{use.cluster} sets if the calibration and interpolation should be parallelised. If \code{TRUE}, then local all local cores will be used. An integer >0 sets the number of local cores to use.
#' An object from \code{makeCluster} allows for a user specified cluster.
#'
#' @param \code{debug.level} Control the user messages. A value >0 outputs \code{hydroMap} progress. See \code{gstat} for the influence of values >0 on the kriging.
#'
#' @seealso \code{\link{krige.head.calib}} for undertaking only the optimisation.
#'
#' @references
#' Gallant, J.C., Dowling, T.I. (2003): 'A multiresolution index of valley bottom flatness for mapping depositional areas', Water Resources Research, 39/12:1347-1359
#'
#' Rivoirard, J. & Romary, T. Math Geosci (2011) Continuity for Kriging with Moving Neighborhood, Mathematical Geosciences, 43: 469. DOI: 10.1007/s11004-011-9330-0
#'
#' @return
#' If \code{is.null(newdata)==TRUE}, then a Spatial object grid will be returned with "head" and "head.var" for the groundwater level and kriging variance respectively.
#'
#' Else, a point Spatial object will be returned with the estmates at the prediction locations and error estimates.
#'
#' @examples
#' # Load packages in case they have not loaded.
#' library(sp)
#' library(grid)
#' library(gstat)
#' library(raster)
#' library(RSAGA)
#' library(parallel)
#' library(rgenoud)
#'
#' # Set enironment path for hydroMap
#' set.env()
#'
#' # Load water table observations from  April 2000 for Victoria, Australia and a 250m state-wide DEM.
#' data('victoria.groundwater')
#'
#' # Load a model variogram and mapping parametyers found to be effective.
#' data('mapping.parameters')
#
#' # Define a simple kriging formula without MrVBF terms that does not require the package RSAGA.
#' f <- as.formula('head ~ elev + smoothing')
#'
#' # Interpolate the head data.
#' heads <- krige.head(formula=f, grid=DEM, data=obs.data, data.errvar.colname='total_err_var',
#' model=model, smooth.std=smooth.std, maxdist=maxdist, nmax=nmax, fit.variogram.type=3, debug.level=1)
#'
#' # Recalibrate the parameters and map using the default settings.
#' heads <- krige.head(formula=f, grid=DEM, data=obs.data, data.errvar.colname='total_err_var',
#' model = 'Mat',  fit.variogram.type=1, debug.level=1)
#'
#' @export
krige.head <- function(
  formula = as.formula("head ~ elev + MrVBF + MrRTF + smoothing"),
  grid=NULL,
  grid.landtype.colname = NULL,
  data=NULL,
  data.fixedHead=NULL,
  newdata=NULL,
  data.errvar.colname = NULL,
  model=NULL,
  mrvbf.pslope = if(any(match(all.vars(as.formula(formula)), 'MrVBF',nomatch=F) | match(all.vars(as.formula(formula)), 'MrRTF',nomatch=F))){seq(0.5, 1.5, length.out = 11)}else{NULL},
  mrvbf.ppctl  = if(any(match(all.vars(as.formula(formula)), 'MrVBF',nomatch=F) | match(all.vars(as.formula(formula)), 'MrRTF',nomatch=F))){seq(0.5, 1.5, length.out = 11)}else{NULL},
  smooth.std = seq(0.5,1.5,length.out=11),
  nmax = if(is.character(data)){-999}else{ceiling(seq(0.1,0.20,0.01)*length(data))},
  nmax.fixedHead = if(!is.null(data.fixedHead)) {seq(10,110,length=11)}else{NULL},
  maxdist = if(class(grid)=='SpatialPixelsDataFrame' || class(grid)=='SpatialGridDataFrame'){ceiling(0.5*sqrt((extent(grid)[2]-extent(grid)[1])^2 + (extent(grid)[4]-extent(grid)[3])^2)*seq(0.1,1,0.1))}else{-999},
  nmin = 0.2,
  omax = NULL,
  nsim = 0,
  fit.variogram.type=1,
  objFunc.type=1,
  use.cluster=TRUE,
  debug.level=0,
  ...) {

  # NOTE, to build pdf manual:
  # path <- find.package("hydroMap")
  # system(paste(shQuote(file.path(R.home("bin"), "R")),"CMD", "Rd2pdf", shQuote(path)))

  # Remove function details from warnings
  options(warning.expression=NULL)

  # Print initial message
  if (debug.level>0)
    message('Kriging groundwater level:');

  # Check enviro variables are setup
  if (!exists('pkg.env') || is.null(pkg.env))
    stop('The environment variable are not setup. Call set.env() with paths to SAGA.');

  # Check if point or grid kriging to be undertaken
  if (is.null(newdata)) {
    do.grid.est <- TRUE;
    if (debug.level>0)
      message('... Estimation type: grid cells')
  } else {
    do.grid.est <- FALSE;
    if (debug.level>0)
      message('... Estimation type: point cooridinates')
  }

  # Catch error if geostatistical simulations are requested for point estimation, rather than gridding.
  if (!do.grid.est && nsim!=0)
    stop('Point estimation cannot be undertaken for geostatistical simulations. To undertake gridded simulations, set newdata=NULL and nsim>0.')

  # Get variable names and check if using MrVBF or MrRTF or DEM
  if (debug.level>0)
    message('... Checking formula terms.');
  var.names = all.vars(formula);
  use.MrVBF = any(match(var.names, 'MrVBF'));
  use.MrRTF = any(match(var.names, 'MrRTF'));
  use.elev = any(match(var.names, 'elev')) || any(match(var.names, 'head'));
  use.DEMsmoothing = any(match(var.names, 'smoothing'));
  do.head.est = any(match(var.names, 'head'));
  do.depth.est = any(match(var.names, 'depth'));
  if (is.na(use.MrVBF))
    use.MrVBF = FALSE;
  if (is.na(use.MrRTF))
    use.MrRTF = FALSE;
  if (is.na(do.head.est))
    do.head.est = FALSE;
  if (is.na(do.depth.est))
    do.depth.est = FALSE;
  if (is.na(use.elev))
    use.elev = FALSE;
  if (is.null(smooth.std) || is.na(smooth.std) || is.na(use.DEMsmoothing))
    use.DEMsmoothing = FALSE;


  # Check variable names
  if (!use.elev && !use.MrVBF && !use.MrRTF && !use.DEMsmoothing && length(var.names)>1) {
    warning('Formula does not contain any of the expected terms: "elev", "MrVBF", "MrRTF", "smoothing"')
  }

  # Check if terms other tah than start options are listed in the formula
  use.extraTerms = FALSE;
  ind = match(var.names, c('head','depth','elev','smoothing','MrVBF','MrRTF'))
  if (any(is.na(ind))) {
    use.extraTerms = T;
    use.extraTerms.names = var.names[is.na(ind)]
  }

  # If there are extra terms to the formula, then check that the grid contauins the required extra data.
  if (use.extraTerms) {
    if (class(grid)!='SpatialPixelsDataFrame' && class(grid)!='SpatialGridDataFrame')
      stop('when the formula contains non-started extra terms, the input grid must be a SpatialPixelsDataFrame or SpatialGridDataFrame containing data for the extra terms.')

    # Check grid contains the data for the extra terms
    for (i in 1:length(use.extraTerms.names)) {
      if (all(is.na(match(names(grid), use.extraTerms.names[i]))))
        stop(paste('The input "grid" does not contain data for the extra formula term:',use.extraTerms.names[i]))
    }
  }

  # Check grid and data are provided.
  if (is.null(data) || (is.character(data) && nchar(data)==0) || ((is.null(grid) || (is.character(grid) && nchar(grid)==0)) && (!do.head.est || use.elev || use.MrVBF || use.MrRTF || use.DEMsmoothing || any(smooth.std>0))))
    stop('The inputs "data" and "grid" are both required when kriging is undertaken using the DEM.')

  # Check is MrVBF can be used.
  if ((use.MrVBF || use.MrRTF) && !pkg.env$saga.has.MrVFB)
    warning('MrVBF & MrRTF cannot calculated because RSAGA is not setup correctly. As a work around, copy the MRVBF and MrRTF grids for all combinations of pslope and ppctl to the working directory.');

  # Check the inputs 'data' contains columns called 'head'.
  if (!any(match(names(data),'head')))
    stop('Input "data" must contain a variable called "head".')


  # Check the inputs 'newdata' contains columns called 'head'.
  if (!is.null(newdata) && !any(match(names(data),'head')))
    stop('Input "newdata" must contain a variable called "head". If you do not require the residuals, it can be input as any value.')

  # Check variogram inputs.
  if (fit.variogram.type==1) {
    if (class(model)[1] =="variogramModel" || !is.character(model))
      stop(paste('When fit.variogram.type==1 the input to "model" must be the variogram type (as character string). The given input was of type',class(model)[1]))

  } else if (fit.variogram.type==2) {
    if (class(model)[1] !="variogramModel")
      stop(paste('When fit.variogram.type==2 the input to "model" must be a gstat variogram object. The given input was of type',class(model)[1]))
  } else if (fit.variogram.type==3) {
    if (class(model)[1] !="variogramModel")
      stop(paste('When fit.variogram.type==3 the input to "model" must be a gstat variogram object. The given input was of type',class(model)[1]))
  } else {
    stop('Unknown input fit.variogram.type.')
  }


  # Find   mid-point of MRVBF & smoothing parameters. This only does anything if
  # the parameters are vectors. It is required to that some initial data is available
  # for the calibration variogram estimation.
  if (use.MrVBF || use.MrRTF) {
    if (length(mrvbf.pslope)==2)
      ind_pslope = min(mrvbf.pslope) + range(mrvbf.pslope)/2
    else {
      ind_pslope = max(1,round(length(mrvbf.pslope)/2))
    }
    if (length(mrvbf.ppctl)==2)
      ind_ppctl = min(mrvbf.ppctl) + range(mrvbf.ppctl)/2
    else {
      ind_ppctl = max(1,round(length(mrvbf.ppctl)/2))
    }

    mrvbf.pslope.tmp = mrvbf.pslope[ind_pslope];
    mrvbf.ppctl.tmp = mrvbf.ppctl[ind_ppctl];

  } else {
    mrvbf.pslope.tmp = NULL;
    mrvbf.ppctl.tmp = NULL;
  }
  if (use.DEMsmoothing) {
    if (length(smooth.std)==2)
      ind_smooth = min(smooth.std) + range(smooth.std)/2
    else {
      ind_smooth = max(1,round(length(smooth.std)/2))
    }
    smooth.std.tmp = smooth.std[ind_smooth];
  } else {
    smooth.std.tmp = NULL;
  }

  # Delete unneccesary grid columns
  if (!is.character(grid)) {
    colnames = names(grid);
    filt = colnames !='DEM'
    if(!is.null(grid.landtype.colname))
      filt = filt & colnames !=grid.landtype.colname
    colnames = colnames[filt]
    if (length(colnames)>0)
      for (i in 1:length(colnames)) {
        if (!use.extraTerms || (use.extraTerms && all(is.na(match(use.extraTerms.names,colnames[i])))))
          grid[[colnames[i]]] = NULL;
      }
  }

  # Delete unneccesary data columns
  if (!is.character(data)) {
    colnames = names(data);
    filt = colnames !='head' & colnames!='depth' & colnames!='elev'
    if(!is.null(data.errvar.colname))
      filt = filt & colnames !=data.errvar.colname
    colnames = colnames[filt]
    if (length(colnames)>0)
      for (i in 1:length(colnames)) {
        #if (!use.extraTerms || (use.extraTerms && all(is.na(match(use.extraTerms.names,colnames[i])))))
          data[[colnames[i]]] = NULL;
      }
  }

  # Delete unneccesary data.fixedHead columns
  if (!is.null(data.fixedHead)) {
    if (!is.character(data.fixedHead)) {
      colnames = names(data.fixedHead);
      filt = colnames !='head' & colnames!='depth' & colnames!='elev'
      colnames = colnames[filt]
      if (length(colnames)>0)
        for (i in 1:length(colnames)) {
          #if (!use.extraTerms || (use.extraTerms && all(is.na(match(use.extraTerms.names,colnames[i])))))
            data.fixedHead[[colnames[i]]] = NULL;
        }
    }
  }

  # Delete unneccesary newdata columns
  if (!is.null(newdata) && !is.character(newdata)) {
    colnames = names(newdata);
    filt = colnames !='head' & colnames!='depth' & colnames!='elev'
    if(!is.null(data.errvar.colname))
      filt = filt & colnames !=data.errvar.colname
    colnames = colnames[filt]
    if (length(colnames)>0)
      for (i in 1:length(colnames)) {
        #if (!use.extraTerms || (use.extraTerms && all(is.na(match(use.extraTerms.names,colnames[i])))))
          newdata[[colnames[i]]] = NULL;
      }
  }

  # get grid data and remove grid columns not required.
  if (debug.level>0)
    message(' ... Getting grid data');
  grid = get.allData(formula = formula, NULL, grid, mrvbf.pslope = mrvbf.pslope.tmp, mrvbf.ppctl =mrvbf.ppctl.tmp, smooth.std = smooth.std.tmp, debug.level=debug.level)
  grid$smoothDEM = NULL;

  # Convert DEM to Raster object
  dem.asRaster = raster(grid,layer='DEM');

  # Setup newdata and get grid values ta points and then
  # check if any data points are outside the DEM area.
  if (!do.grid.est) {
    newdata <- import.pointData(newdata, debug.level=debug.level);
    row.names(newdata) = seq(1,length(newdata))
    tmp = extract(dem.asRaster, newdata, method='bilinear');
    filt = is.na(tmp);
    if (any(filt)) {
      warning(paste(sum(filt),' newdata points were outside the DEM extent and have been removed from the analysis.'),,immediate.=T)
      newdata = newdata[!filt,];
    }
  }

  # Setup data and get grid values at points anbd then
  # check if any data points are outside the DEM area.
  data  = import.pointData(data, debug.level=debug.level);
  row.names(data) = seq(length(newdata)+1,length(newdata)+length(data))
  tmp = extract(dem.asRaster, data, method='bilinear');
  filt = is.na(tmp);
  if (any(filt)) {
    warning(paste(sum(filt),' data points were outside the DEM extent and have been removed from the analysis.'),immediate.=T)
    data = data[!filt,];
  }

  # get point data for fixedHead points
  if (debug.level>0)
    message('... Getting point data');
  if (!is.null(data.fixedHead)) {
    # In order to store both obs points and fixed head points in 'pkg.env', merge the
    # data, get required point data and seperate.
    #---------------

    # Import data (does nothing if not a file name)
    data.fixedHead  = import.pointData(data.fixedHead, debug.level=debug.level);

    # Rewrite row names to avoid dubplicate row.names in latter rbind()
    row.names(data.fixedHead) = seq(length(newdata)+length(data)+1,length(newdata)+length(data)+length(data.fixedHead))

    # Check if any fixed-head points are outside the DEM area
    tmp = extract(dem.asRaster, data.fixedHead, method='bilinear');
    filt = is.na(tmp);
    if (any(filt)) {
      warning(paste(sum(filt),' fixed-head points were outside the DEM extent and have been removed from the analysis.'),immediate.=T)
      data.fixedHead = data.fixedHead[!filt,];
    }
    rm(tmp,filt);

    # Delete all columns from fixed heads data, except for 'head'
    colnames = names(data.fixedHead);
    filt = colnames !='head'
    colnames = colnames[filt]
    if (length(colnames)>0)
      for (i in 1:length(colnames)) data.fixedHead[[colnames[i]]] = NULL;

    # Add empty columns for fixed head
    colnames = names(data);
    filt = colnames !='head'
    colnames = colnames[filt]
    if (length(colnames)>0) {
      for (i in 1:length(colnames)) {
#        if (!use.extraTerms || (use.extraTerms && all(is.na(match(use.extraTerms.names,colnames[i])))))
          data.fixedHead[[colnames[i]]] = -9999.0;
      }
    }


    # Merge data
    if (do.grid.est) {
      data.all = rbind(data, data.fixedHead)
    } else {
      data.all = rbind(data, newdata, data.fixedHead)
    }
    data.all = get.allData(formula = formula, data.all, grid, mrvbf.pslope = mrvbf.pslope.tmp, mrvbf.ppctl =mrvbf.ppctl.tmp, smooth.std = smooth.std.tmp, debug.level=debug.level)

    # Disaggregate data
    data = data.all[seq(1,nrow(data)),]
    if (do.grid.est) {
      data.fixedHead = data.all[seq(nrow(data)+1,nrow(data.all)),]
    } else {
      newdata = data.all[seq(nrow(data)+1,nrow(data)+nrow(newdata)),]
      data.fixedHead = data.all[seq(nrow(data)+nrow(newdata)+1,nrow(data.all)),]
    }
    rm(data.all);

    # Check non-finite vales from grid
    filt = matrix(F, length(data.fixedHead),ncol(data.fixedHead))
    for (i in 1:ncol(data.fixedHead))
      filt[,i] = is.finite(data.fixedHead[[names(data.fixedHead)[i]]])
    filt = rowSums(filt)==ncol(data.fixedHead)
    if (any(!filt)) {
      warning(paste(sum(!filt),' non-finite fixed head values from the interpolated grid have been removed.'),immediate.=T)
      data.fixedHead = data.fixedHead[filt,]
    }
    #---------------
  } else {

    # Get point vakues from grid
    if (do.grid.est) {
      data.all = data;
    }else {
      data.all = rbind(data, newdata)
    }
    data.all = get.allData(formula = formula, data.all, grid, mrvbf.pslope = mrvbf.pslope.tmp, mrvbf.ppctl =mrvbf.ppctl.tmp, smooth.std = smooth.std.tmp, debug.level=debug.level)

    # Disaggregate data
    if (do.grid.est) {
      data = data.all;
    } else {
      data = data.all[seq(1,nrow(data)),]
      newdata = data.all[seq(nrow(data)+1,nrow(data.all)),]
    }
    rm(data.all);
  }

  # Free up memory
  rm(dem.asRaster);

  # Check non-finite vales from grid in newdata
  if (!do.grid.est) {
    filt = matrix(F, length(newdata),ncol(newdata))
    for (i in 1:ncol(newdata))
      filt[,i] = is.finite(newdata[[names(newdata)[i]]])
    filt = rowSums(filt)==ncol(newdata)
    if (any(!filt)) {
      warning(paste(sum(!filt),' non-finite newdata values from the interpolated grid have been removed.'),immediate.=T)
      newdata = newdata[filt,]
    }
  }

  # Check non-finite vales from grid in data
  filt = matrix(F, length(data),ncol(data))
  for (i in 1:ncol(data))
    filt[,i] = is.finite(data[[names(data)[i]]])
  filt = rowSums(filt)==ncol(data)
  if (any(!filt)) {
    warning(paste(sum(!filt),' non-finite data values from the interpolated grid have been removed.'),immediate.=T)
    data = data[filt,]
  }

  # Get calibration values for nmax and maxdist. This isshould only be undertaken if the user input strings for 'data' or 'grid'
  if (length(nmax)==1 && nmax==-999)
    nmax = ceiling(seq(0.04,0.20,0.02)*length(data));
  if (length(maxdist)==1 && maxdist==-999)
    maxdist = ceiling(0.5*sqrt((extent(grid)[2]-extent(grid)[1])^2 + (extent(grid)[4]-extent(grid)[3])^2)*seq(0.1,1,0.1));

  # Warn if there are major differences between the bore elevation and DEM elevation at the bore.
  nelev_diff = sum(abs(data$elev - data$DEM) > 25);
  if (nelev_diff>0) {
    warning(paste('the difference between the input bore elevation and DEM elevation is > 25 elevation units at ', nelev_diff, ' data points!'),immediate.=T);
  }

  # Check if any required mapping parameters are vectors.
  # If the parameters are vectors of two elemants, then continuous value
  # calibration is undertaken. If the parameters are discrete, the faster
  # discrete value calibration is undertaken at the provided parameter increments.
  if (!do.grid.est && (length(mrvbf.pslope)>1 || length(mrvbf.ppctl)>1 || length(smooth.std)>1 || length(nmax)>1 || length(maxdist)>1)) {
    do.calibration <- TRUE;
  } else {
    do.calibration <- FALSE;
  }

  if (!do.calibration && fit.variogram.type==1)
    stop('When fit.variogram.type==1 the variogram parameters must be calibrated. If not calibrating, then input a vgm() derived variogram model and set fit.variogram.type==3 to do no calibration.')

  if (do.calibration && nsim!=0) {
    stop('Parameter calibration cannot be undertaken for geostatistical simulations. Set nsim=0 or input the model parameter as scalars, not vectors.')
  }

  # Undertake calibration if required.
  if (do.calibration) {
    message('... Calling krige.head.calib() to calibrate parameters:')
    calib.solution <- krige.head.calib(formula = formula, grid=grid, grid.landtype.colname=grid.landtype.colname, data=data, newdata=newdata,  data.errvar.colname = data.errvar.colname, model=model,
                                       mrvbf.pslope = mrvbf.pslope,
                                       mrvbf.ppctl = mrvbf.ppctl,
                                       smooth.std = smooth.std,
                                       nmax = nmax,
                                       maxdist = maxdist ,
                                       nmin = nmin,
                                       omax = omax ,
                                       nsim = nsim, fit.variogram.type=fit.variogram.type, objFunc.type=objFunc.type, use.cluster=use.cluster, debug.level=debug.level, ...)

    # Exit if calibration failed.
    if (!is.list(calib.solution))
      stop('Calibration failed.')

    # Get model type and set default kappa and anisotropic ratio/angle in case not calibrated
    # Define/get initial variogram
    if (class(model)[1] =="variogramModel") {
      filt = model$model!='Nug';
      modelType = model$model[filt]
      range = model$range[filt]
      psill = model$psill[filt]
      nugget = model$psill[!filt]
      kappa= model$kappa[filt]
      ang1= model$ang1[filt]
      anis1= model$anis1[filt]
    } else {
      if (is.character(model)) {
        modelType = model;
      } else
        modelType ='Mat';

      # Set default kappa and anisotropic ratio/angle in case not calibrated
      anis1 = 1;
      ang1 = 0;
      kappa=0.5;
    }

    # Get list of calibrated parameter names
    param.Names = calib.solution$param.Names;
    for (i in 1:length(param.Names)) {

      if (param.Names[i] == 'mrvbf.pslope') {
        mrvbf.pslope = calib.solution$par[i];
      } else if (param.Names[i] == 'mrvbf.ppctl') {
        mrvbf.ppctl = calib.solution$par[i];
      } else if (param.Names[i] == 'smooth.std') {
        smooth.std = calib.solution$par[i];
      } else if (param.Names[i] == 'nmax') {
        nmax = calib.solution$par[i];
      } else if (param.Names[i] == 'maxdist') {
        maxdist = calib.solution$par[i];
      } else if (param.Names[i] == 'nug') {
        nug = calib.solution$par[i];
      } else if (param.Names[i] == 'psill_model_1') {
        psill = calib.solution$par[i];
      } else if (param.Names[i] == 'range_model_1') {
        range = calib.solution$par[i];
      } else if (param.Names[i] == 'kappa_model_1') {
        kappa = calib.solution$par[i];
      } else if (param.Names[i] == 'anis1_model_1') {
        anis1 = calib.solution$par[i];
      } else if (param.Names[i] == 'ang1_model_1') {
        ang1 = calib.solution$par[i];

      } else
        stop(paste('Unknown parameter names returned from calib: ',param.Names[i]))
    }

    if (fit.variogram.type==1){

      if (debug.level>0) {
        message('... Building the following variogram model:')
        message(paste('    Nugget = ',nug))
        message(paste('    Type = ',modelType))
        message(paste('    Partial sill = ',psill))
        message(paste('    Range = ',range))
        message(paste('    Kappa = ',kappa))
        message(paste('    Ang1 = ',ang1))
        message(paste('    Anis1 = ',anis1))
      }

      model = vgm(nugget=nug, psill=psill, range=range, kappa=kappa, anis=c(ang1,anis1),model=modelType)
    }
  }

  # Build variogram if not provided.
  if (fit.variogram.type==2) {
    if (debug.level>0)
      message('... Calling get.variogram():')

    model = get.variogram(formula=formula, data=data, model=model, nmin=nmin, nmax=nmax, maxdist = maxdist, debug.level=debug.level);

    if (debug.level>0)
      message('... Exiting get.variogram().')
  }

  # If 0 <= omax < 1, then nmin is treated as a fraction of nmax, else omax is just passed to gstat.
  omax.fixedHead = NULL
  if (!is.null(omax) && omax>=0 && omax<1) {
    if (!(is.null(nmax.fixedHead))) {
      omax.fixedHead = max(0,floor(omax * nmax.fixedHead))
    }
    omax = max(0,floor(omax * nmax))
    if (debug.level>0)
      message('... omax set to ', omax)
  }

  # If 0 <= nmin < 1, then nmin is treated as a fraction of nmax, else nmax is just passed to gstat.
  nmin.fixedHead = NULL;
  if (!is.null(nmin) && nmin>=0 && nmin<1) {
    if (!(is.null(nmax.fixedHead))) {
      nmin.fixedHead = max(0,floor(nmin * nmax.fixedHead))
    }
    nmin = max(0,floor(nmin * nmax))
    if (debug.level>0)
      message('... nmin set to ', nmin)
  }

  # Set-up observation errors, if required.
  data.weights = NULL;
  if (!is.null(data.errvar.colname)) {
    if (!is.character(data.errvar.colname)) {
      stop('The input "data.errvar.colname" must be NULL or a string for the name of a column within the input "data".')
    }

    # Extract and calc. weights for data.
    filt = data.errvar.colname %in% names(data)
    if (!any(filt)) {
      stop('The input string for "data.errvar.colname" is not a column within the input "data".')
    }

    data.weights = 1/data[[data.errvar.colname]];
  }

  # Remove columns not required.
  if (!do.depth.est)
    data$depth = NULL
  data$smoothDEM = NULL

  # Build gstat model if using land type of fixed head data
  if (!is.null(grid.landtype.colname) || !is.null(data.fixedHead)) {
    #g = gstat(NULL, 'trend', formula=formula, data=data, model=model, nmin=nmin, nmax=nmax, maxdist = maxdist, omax=omax, force=TRUE, weights=data.weights, set=list(trnd_threshdist=maxdist*trendMaxDistFrac))
    g = gstat(NULL, 'trend', formula=formula, data=data, model=model, nmin=nmin, nmax=nmax, maxdist = maxdist, omax=omax, force=TRUE, weights=data.weights)
  }

  # If land catagories, then extend the kriging to collocated univeral cokriging.
  # Doing so, requires estimation of the universal kriging residuals. This is done below.
  # The residuals are then used to find a threshold that, if not provided,
  # maximises the range of probabilities between the land types.
  # This threshold is then used with the Bayes-Markov theorm to derive a
  # variogram for the land catagory and then co-variogram with the primary variable.
  model.landtype = NULL;
  model.landtype.head = NULL;
  use.LandCatagory = F;
  if (!is.null(grid.landtype.colname)) {

    # Check inputs
    if (is.null(grid))
      stop('The grid must be input when including grid.landtype.colname');
    if (!is.character(grid.landtype.colname))
      stop('grid.landtype.colname must be a character string.');
    if (length(grid.landtype.colname)!=1)
      stop('grid.landtype.colname must be a character string containing one column name.');
    if (!(grid.landtype.colname %in% names(grid)))
      stop('The column name input to grid.landtype.colname is not a column within the input grid.');

    # Interpolate landtype grid column to residual coordinates
    tmp = extract(raster(grid,grid.landtype.colname), data, method='simple');
    if (any(is.na(tmp))) {
      # Apply a buffer to find the most frequent point in the region
      tmp = extract(raster(grid,grid.landtype.colname), data, method='simple');
    }
    data$tmpName =  tmp;
    ncols = length(names(data))
    names(data)[ncols] = grid.landtype.colname
    if (!is.null(newdata)) {
      tmp = extract(raster(grid,grid.landtype.colname), newdata, method='simple');
      newdata$tmpName =  tmp;
      ncols = length(names(newdata))
      names(newdata)[ncols] = grid.landtype.colname
    }

    # Find the unique land types.
    catagories =  unique(grid[[grid.landtype.colname]])
    catagories = catagories[!is.na(catagories)];

    # Randomly select observations from each land catagory.
    # This is, only a subset of observation points are used for
    # the Bayes-Markov variogram estimation. This is done to minimise
    # the computational load.
    if (is.null(pkg.env$landtype.indexes) || length(pkg.env$landtype.indexes) != length(data)) {
      ind.all = seq(1, length(data));
      ind.subset = c()
      landtype.sampleFrac = 0.1
      landtype.sampleMin = 100;
      pkg.env$landtype.indexes = rep(FALSE,  length(data));
      for (i in 1:length(catagories)) {
        filt = data[[grid.landtype.colname]] == catagories[i]
        filt = filt[!is.na(filt)]
        if (length(filt)==0)
          next
        sampleSize = max( min(sum(filt),landtype.sampleMin),ceiling(landtype.sampleFrac*sum(filt)))
        ind.tmp = sample(ind.all[filt], sampleSize,  replace=F)
        pkg.env$landtype.indexes[ind.tmp] = T;
      }
    }
    ind.subset = pkg.env$landtype.indexes;
    data.subset = data[ind.subset,];

    # estimate BLUE trend at points
    if (debug.level>0)
      message('... Starting BLUE trend estimation  of ', length(data.subset),' points for land catagories')

    trend <- predict(g, data.subset, BLUE=T, debug.level=debug.level)
    data.subset$trend.resid = data.subset$head - trend$trend.pred;

    # Find optimal residual thershold ie that which maximimises B.
    resid.percentiles = seq(0.025, 0.975, 0.025)
    B = vector(mode='numeric', length=length(resid.percentiles));
    for ( i in 1:length(resid.percentiles)) {
      residThreshold = quantile(data.subset$trend.resid, resid.percentiles[i], na.rm=T)
      data.landtype.prob = vector(mode='numeric', length=length(data.subset));
      for (j in 1:length(catagories)) {
        filt = data.subset[[grid.landtype.colname]] == catagories[j];
        filt = filt[!is.na(filt)]
        prob = sum(data.subset$trend.resid[filt] < residThreshold, na.rm=T)/sum(filt, na.rm=T)
        data.landtype.prob[filt]  = prob;
      }

      # Calculate Bayes-Markov parameter, B
      I = data.subset$trend.resid < residThreshold;
      m1 = sum(data.landtype.prob * I,na.rm=T)/sum(I,na.rm=T)
      m0 = sum(data.landtype.prob * (1-I),na.rm=T)/sum(1-I,na.rm=T)
      B[i] = m1 - m0;
    }

    # Report maximum B value and percentile.
    ind = seq(1,length(B),1);
    ind = ind[B == max(B)];
    if (length(ind)>1)
      ind = ind[1]
    B = B[ind];
    resid.percentiles = resid.percentiles[ind];
    if (debug.level>0)
      message(paste('... Maximum Markov-Bayes parameter, B, found =',B, ' at percentile ',resid.percentiles))

    if (B<1e-6) {
      # Turn off inclusion of land type
      use.LandCatagory = F;
      grid[[grid.landtype.colname]] = NULL;

    } else {

      # Add point data for the land type. The point data is probability of being below the thresholod residual vale,
      # but the mean is shifted to equal the mean of the residuals
      grid.landtype.colname.tmp = paste(grid.landtype.colname,'_tmp',sep='')
      grid[[grid.landtype.colname.tmp]] = NULL;
      grid[[grid.landtype.colname.tmp]] = NaN;
      data[[grid.landtype.colname.tmp]] = NULL;
      data[[grid.landtype.colname.tmp]] = NaN;
      residThreshold = quantile(data.subset$trend.resid, resid.percentiles, na.rm=T)
      for (j in 1:length(catagories)) {
        filt = data.subset[[grid.landtype.colname]] == catagories[j];
        filt = filt[!is.na(filt)]
        prob = sum(data.subset$trend.resid[filt] < residThreshold, na.rm=T)/sum(filt, na.rm=T)

        filt = data[[grid.landtype.colname]] == catagories[j];
        data[[grid.landtype.colname.tmp]][filt]  = prob;

        filt = grid[[grid.landtype.colname]] == catagories[j];
        filt = filt[!is.na(filt)]
        grid[[grid.landtype.colname.tmp]][filt]  = prob;

        if (!is.null(newdata)) {
          filt = newdata[[grid.landtype.colname]] == catagories[j];
          filt = filt[!is.na(filt)]
          newdata[[grid.landtype.colname.tmp]][filt]  = prob;
        }
      }

      # Remove 'data' column grid.landtype.colname to minimise RAM
      data[[grid.landtype.colname]] = NULL;
      grid[[grid.landtype.colname]] = NULL;
      data$trend.resid == NULL
      if (!is.null(newdata))
        newdata[[grid.landtype.colname]] = NULL;

      # Chnage name of column grid.landtype.colname.tmp to 'LandType'.
      names(data)[names(data) == grid.landtype.colname.tmp] = 'LandType'
      names(grid)[names(grid) == grid.landtype.colname.tmp] = 'LandType'
      if (!is.null(newdata))
        names(newdata)[names(newdata) == grid.landtype.colname.tmp] = 'LandType'
      data[[grid.landtype.colname.tmp]] = NULL;
      grid[[grid.landtype.colname.tmp]] = NULL;
      if (!is.null(newdata))
        newdata[[grid.landtype.colname.tmp]] = NULL;

      # Add variograms and cross variograms for the land type using Bayes Markov assumption
      model.landtype <- model
      vov <- var(grid$LandType, na.rm=T)/var(data.subset$trend.resid, na.rm=T)
      model.landtype$psill <- model$psill * vov
      model.landtype.head <- model
      model.landtype.head$psill <-  sqrt(model.landtype$psill*model$psill) * cor(data.subset$LandType,data.subset$trend.resid, use='na.or.complete');
      rm(data.subset);

      if (debug.level>0) {
        message('... Model Summary (compnent, nugget, sill, range, kappa):')
        message(paste('... Head,              ',model$psill[1], model$psill[2],model$range[2],model$kappa[2]))
        message(paste('... Land Type,         ',model.landtype$psill[1],model.landtype$psill[2],model.landtype$range[2],model.landtype$kappa[2]))
        message(paste('... Land Type Cross-Cor, ',model.landtype.head$psill[1],model.landtype.head$psill[2],model.landtype.head$range[2],model.landtype.head$kappa[2]))
      }
    }
  }

  # If using fixed heads, then build a variogram and cross-variogram for the fixed head.
  # The variograms are built using the input kriging formula and the Bayes-Markov theorm.
  # The deriving the variogram for the fixed heads, the formula and the residuals from the
  # primary variable are adopted and filtered to only those close to the fixed head points.
  # These points are then used to scale the primary variable variogram. In deriving the  cross-variogram,
  # the correlation is calculated from the primary variable residuals at fixed head points and the difference
  # between the fixed head value at the universal kriging trend.
  model.fixedHead = NULL;
  model.fixedHead.head = NULL;
  use.FixedHeads = F
  if (!is.null(data.fixedHead)) {

    use.FixedHeads = T

    # Get grid resolution for finding fixed points with same grid cell at obs points.
    # The for-loops combined find the minimum number UNIQUE of pairs of obs point-fxed head data.
    # That is, the first for-loop finds the closest fixed head point to each obs point. However,
    # multiple obs points may be closest to the same fixed head point. The second for-loop finds
    # only the obs points closest  each fixed head point.
    minPairs=25
    nPairs = 0
    nCells=1;
    maxCells=10;
    grid.DEM.params = gridparameters(grid)
    ind = matrix(nrow=0,ncol=2)
    data.fixedHead.easting = coordinates(data.fixedHead)[,1]
    data.fixedHead.northing = coordinates(data.fixedHead)[,2]
    data.easting = coordinates(data)[,1]
    data.northing = coordinates(data)[,2]
    while (nPairs<minPairs && nCells<=maxCells) {
      zero.dist = nCells *sqrt(grid.DEM.params$cellsize[1]^2+grid.DEM.params$cellsize[2]^2);

      # Find indexes to obs. points near to fixed head points.
      ind = zerodist2(data, data.fixedHead, zero = zero.dist)

      # Get list of unique obs points
      unique.bores = unique(ind[,1]);
      ind.unique = matrix(0, length(unique.bores),2)

      # Loop through each unique obs point to find the closest fixed head point.
      for (i in 1:length(unique.bores)) {
        # Add bore index indexes
        ind.unique[i,1] = unique.bores[i];

        # Find fixed head point ind values near to the current bore
        filt =  which(ind[,1] == unique.bores[i]);
        ind.fixedHead = ind[filt,2];

        if (length(ind.fixedHead)>1) {
          # Calc. distance from fixed points to obs point
          distSqr = (data.fixedHead.easting[ind.fixedHead] - data.easting[unique.bores[i]])^2 + (data.fixedHead.northing[ind.fixedHead] - data.northing[unique.bores[i]])^2;

          # Find the fixed point with the lowest distance to the obs point.
          distSqr.filt = which(distSqr == min(distSqr))[1]

          # Add index to the closest foxed point
          ind.unique[i,2] = ind.fixedHead[distSqr.filt];

        } else {
          ind.unique[i,2] = ind.fixedHead;
        }
      }

      # Replace ind and update while loop criteria etc
      ind = ind.unique

      # Get list of unique fixed head points
      unique.fixedHead = unique(ind[,2]);
      ind.unique = matrix(0, length(unique.fixedHead),2)

      # Loop through each unique fixed head point to find the closest obs point.
      for (i in 1:length(unique.fixedHead)) {
        # Add bore index indexes
        ind.unique[i,2] = unique.fixedHead[i];

        # Find obs point ind values near to the current bore
        filt =  which(ind[,2] == unique.fixedHead[i]);
        ind.obsHead = ind[filt,1];

        if (length(ind.obsHead)>1) {
          # Calc. distance from fixed points to obs point
          distSqr = (data.fixedHead.easting[unique.fixedHead[i]] - data.easting[ind.obsHead])^2 + (data.fixedHead.northing[unique.fixedHead[i]] - data.northing[ind.obsHead])^2;

          # Find the fixed point with the lowest distance to the obs point.
          distSqr.filt = which(distSqr == min(distSqr))[1]

          # Add index to the closest foxed point
          ind.unique[i,1] = ind.obsHead[distSqr.filt];

        } else {
          ind.unique[i,1] = ind.obsHead;
        }
      }

      # Replace ind and update while loop criteria etc
      ind = ind.unique
      nPairs = nrow(ind);
      nCells = nCells +1;
    }

    if (nrow(ind)==0) {
      warning(paste('No obs. points could be found near fixed points. Fixed head points will not be used.'),immediate.=T)
      use.FixedHeads = F;
    } else if (nrow(ind)<minPairs) {
      warning(paste('Only the following number of obs points could be found near fixed points: ',nrow(ind)),immediate.=T)
    }

    # Replace input elev with DEM. In development, the use of, say, zero elev for the coastline
    # at all fixed points caused singular kriging matrix.
    if (use.elev)
      data.fixedHead$elev = data.fixedHead$DEM

    # Build data sets fixed head and residuals near to each other.
    if (debug.level>0)
      message('... Kriging BLUE (trend) for obs. points near fixed head points.')
    data.nearFixedHead =  data[ind[,1],];
    trend.nearFixedHead = predict(g, data.nearFixedHead, BLUE=T, debug.level=debug.level);
    data.nearFixedHead$trend.resid = data.nearFixedHead$head - trend.nearFixedHead$trend.pred;

    if (debug.level>0)
      message('... Kriging BLUE (trend) for fixed head points near obs.')
    fixedHead.nearData =  data.fixedHead[ind[,2],];
    filt = unique(row.names(fixedHead.nearData));
    fixedHead.nearData =fixedHead.nearData[filt,];
    fixedHead.trend = predict(g, fixedHead.nearData, BLUE=T, debug.level=debug.level);
    fixedHead.nearData$trend.resid = fixedHead.nearData$head - fixedHead.trend$trend.pred;

    # Add variograms and cross variograms for the fixed head data using Bayes Markov assumption
    if (debug.level>0)
      message('... Calculating Bayes-Markov co-variograms.')
    model.fixedHead <- model
    vov <- var(fixedHead.nearData$trend.resid, na.rm=T)/var(data.nearFixedHead$trend.resid, na.rm=T)
    model.fixedHead$psill <- model$psill * vov
    model.fixedHead.head <- model
    crossCorr = cor(fixedHead.nearData$trend.resid,data.nearFixedHead$trend.resid, use='na.or.complete');
    model.fixedHead.head$psill <-  sqrt(model.fixedHead$psill*model$psill) * crossCorr;
    if (debug.level>0) {
      message('... Model Summary (compnent, nugget, sill, range, kappa):')
      message(paste('... Head,                ',model$psill[1], model$psill[2],model$range[2],model$kappa[2]))
      message(paste('... Fixed Head,          ',model.fixedHead$psill[1],model.fixedHead$psill[2],model.fixedHead$range[2],model.fixedHead$kappa[2]))
      message(paste('... Fixed Head Cross-Cor,',model.fixedHead.head$psill[1],model.fixedHead.head$psill[2],model.fixedHead.head$range[2],model.fixedHead.head$kappa[2]))
    }
  }

  # Remove columns not required.
  if (do.depth.est)
    data$depth = NULL
  data$smoothDEM = NULL
  data$BoreID = NULL

  # Do kriging. If doing grid est, then smooth after krigign head.
  # If doing, points, the est neighbouring 24 grid cells and smooth to get est est.
  if (do.grid.est) {
    message('... Setting up kriging of grid cells')

    colnames = names(grid);
    colnames[which(colnames =='DEM')]='elev';
    names(grid) = colnames
    gridded(grid) = FALSE;
    ncells = dim(coordinates(grid))[1];

    # Setup cluster
    is.cluster.setup=FALSE
    if (class(use.cluster)[1]=="SOCKcluster" | class(use.cluster)[1]=="PVMcluster" | class(use.cluster)[1]=="spawnedMPIcluster" | class(use.cluster)[1]=="MPIcluster") {
      cl = use.cluster;
      is.cluster.setup = TRUE
      nclus = length(use.cluster)
    } else {
      message('... Setting up cluster.');
      if (is.logical(use.cluster)) {
        nclus=1
        if (use.cluster) {nclus = detectCores(all.tests = TRUE, logical = FALSE)};
      } else if (is.numeric(use.cluster) && floor(use.cluster)>0) {
        nclus = min(detectCores(all.tests = TRUE, logical = FALSE),floor(use.cluster))
      } else {
        stop('The input use.cluster must be logical or an integer>0 or a cluster object from makeCLuster().')
      }
      if (nclus>1 && !is.cluster.setup) {
        clus <- c(rep("localhost", nclus));
        cl <- makeCluster(clus, type = "SOCK");
        clusterEvalQ(cl, library(gstat));
        clusterEvalQ(cl, library(raster));
        is.cluster.setup=TRUE
      }
    }

    # Creat gstat object.
    #gs = gstat(NULL, 'heads', formula=formula, data=data, model=model, nmin=nmin, nmax=nmax, maxdist = maxdist, omax=omax, force=TRUE, weights=data.weights, set=list(trnd_threshdist=maxdist*trendMaxDistFrac))
    gs = gstat(NULL, 'heads', formula=formula, data=data, model=model, nmin=nmin, nmax=nmax, maxdist = maxdist, omax=omax, force=TRUE, weights=data.weights)
    if (use.LandCatagory) {
      gs = gstat(gs,"land.type", formula = as.formula(paste(grid.landtype.colname,' ~ 1')), grid, nmax=1, nmin=1, model = model.landtype,merge=c("heads","land.type"))
      gs = gstat(gs, c("heads","land.type"), model = model.landtype.head)
    }
    if (use.FixedHeads) {
      nVars = length(all.vars(formula,unique=T))-1
      mergeVarIDs = c("heads", 2, "fixed.head", 2)
      for (i in 2:nVars) {
        mergeVarIDs = c(mergeVarIDs, "heads", i+1, "fixed.head", i+1)
      }

      gs = gstat(gs,"fixed.head", formula = formula, data=data.fixedHead, model = model.fixedHead, nmin=nmin.fixedHead, nmax=nmax.fixedHead, maxdist = maxdist, omax=omax.fixedHead, force=T,merge=c("heads","fixed.head"))
      gs = gstat(gs, c("heads","fixed.head"), model = model.fixedHead.head, merge=list(mergeVarIDs))
    }

    # Do prediction
    if (nclus>1) {
      clusterEvalQ(cl, library(gstat));
      clusterExport(cl, list("gs"), envir = environment());

      if (nsim==0) {
        splt = rep(1:nclus, each = ceiling(ncells/nclus), length.out = ncells)
        newgrid = lapply(as.list(1:nclus), function(w) grid[splt == w,])

        message('... Starting kriging of ', ncells,' grid cells using a cluster.')
        head <- do.call("rbind", parLapply(cl, newgrid, function(lst) predict(gs, lst, debug.level=debug.level)))

      } else {
        sims_perCore = rep(1:nclus, each = ceiling(nsim/nclus), length.out = nsim)
        sims_perCore = as.data.frame(table(sims_perCore))
        sims_perCore = sims_perCore$freq

        message('... Starting simulation of ', nsim,' grids using a cluster.')
        head <- do.call("cbind", parLapply(cl, sims_perCore, function(lst) predict(gs, newdata=grid, nsim=lst, debug.level=debug.level)))

      }


      # Close cluster of the input was not a cluster object
      if (is.logical(use.cluster) || is.numeric(use.cluster))
        stopCluster(cl)
    } else {
      if (nsim==0) {
        message('... Starting kriging of ', ncells,' grid cells using 1 core.')
        head <- predict(gs, grid, nsim=nsim,debug.level=debug.level)
      } else {
        message('... Starting ',nsim, 'simulations of ', ncells,' grid cells using 1 core.')
        head <- predict(gs, grid, nsim=nsim,debug.level=debug.level)
      }
    }

    # Rename simulation names
    if (nsim==0) {
      names(head) = c('head','head.var');
    } else {
      names(head) = paste('head.sim',seq(1,nsim),sep='')
    }

    gridded(head) = TRUE
    if (nsim==0) {
      message('... Finished kriging of ', ncells,' grid cells')
    } else {
      message('... Finished ',nsim, 'simulations of ', ncells,' grid cells.')
    }

    # Derive head
    if (do.depth.est) {
      if (nsim==0) {
        head$heads.pred = grid$elev - head$heads.pred
      } else {
        for (i in 1:nsim) {
          head[[i]] = grid$elev - head[[i]]
        }
      }
      gridded(grid) = TRUE;
    }

    # Smooth heads using Gaussian kernal
    if (is.na(smooth.std) || use.DEMsmoothing || smooth.std<=0) {
      return(head)
    } else {

      message('... Starting Gaussian smoothing of ', ncells,' grid cells')
      # Build Gaussian blur kernal
      cellDist = matrix(1,5,5);
      for (i in 1:5) {
        for(j in 1:5) {
          cellDist[i,j] = (i-3)^2 + (j-3)^2
        }
      }
      sigmaWeights = 1/(2*pi*smooth.std^2) * exp(-cellDist/(2*smooth.std^2) )
      sigmaWeights = sigmaWeights/sum(sigmaWeights);

      # Do smoothing
      if (nsim==0) {
        head.asRaster = focal(raster(head,1), sigmaWeights, na.rm=TRUE);
        var.asRaster = focal(raster(head,2), sigmaWeights, na.rm=TRUE);

        # Convert back from raster
        head = as(head.asRaster,'SpatialPixelsDataFrame')
        head.var = as(var.asRaster,'SpatialPixelsDataFrame')
        head$head.var = head.var[['layer']]

      } else {
        for (i in 1:nsim) {
          head.asRaster = focal(raster(head,i), sigmaWeights, na.rm=TRUE);
          head[[i]] = as(head.asRaster,'SpatialPixelsDataFrame')
        }
      }

      message('... Finished Gaussian smoothing of ', ncells,' grid cells')

      return(head)
    }
  }

  # Setup estimation points est at newdata locations
  nNewData = length(newdata);
  est = matrix(nrow=nNewData,ncol=6)
  easting  = coordinates(newdata)[,1]
  northing = coordinates(newdata)[,2]
  ind = seq(1,nNewData,1)

  # Is undertaking smoothing using predicted head, then convert to raster grids for efficient extraction of neighbouring cells
  grid.MrVBF = NULL;
  grid.MrRTF = NULL;
  grid.LandType = NULL;
  grid.params = NULL;
  grid.elev=NULL;
  if (!use.DEMsmoothing && !is.null(smooth.std)  && smooth.std>0) {
    if (debug.level>0)
      message('... Building raster grids.')
    grid.elev = raster(grid,layer='DEM')
    if (use.MrVBF)
      grid.MrVBF= raster(grid,layer='MrVBF')
    if(use.MrRTF)
      grid.MrRTF= raster(grid,layer='MrRTF')
    if (use.LandCatagory)
      grid.LandType= raster(grid,layer='LandType')
  }

  # Free up memory
  rm(grid);

  # Determining cluster size etc
  is.cluster.setup=FALSE
  if (class(use.cluster)[1]=="SOCKcluster" | class(use.cluster)[1]=="PVMcluster" | class(use.cluster)[1]=="spawnedMPIcluster" | class(use.cluster)[1]=="MPIcluster") {
    cl = use.cluster;
    is.cluster.setup = TRUE
    nclus = length(use.cluster)
  } else {
    if (is.logical(use.cluster)) {
      nclus=1
      if (use.cluster) {nclus = detectCores(all.tests = TRUE, logical = FALSE)};
    } else if (is.numeric(use.cluster) && floor(use.cluster)>0) {
      nclus = min(detectCores(all.tests = TRUE, logical = FALSE),floor(use.cluster))
    } else {
      stop('The input use.cluster must be logical or an integer>0 or a cluster object from makeCLuster().')
    }
  }

  # Do the kriging
  if (nclus>1) {
    if (!is.cluster.setup) {
      if (debug.level>0)
        message('... Setting up cluster.');
      clus <- c(rep("localhost", nclus));
      cl <- makeCluster(clus, type = "SOCK");
      is.cluster.setup=TRUE

      # Add libraries to nodes.
      clusterEvalQ(cl, library(gstat))
      clusterEvalQ(cl, library(raster))
    }

    # Split up obs point data for cluster
    njobs = min(nNewData,nclus*5);
    splitData = newdata;
    splt = rep(1:njobs, each = ceiling(nNewData/njobs), length.out = nNewData);
    splitData = lapply(as.list(1:njobs), function(w) splitData[splt == w,]);
    if (debug.level>0)
      message('... Starting kriging of ', nNewData,' points using a cluster.')
    est <- do.call("rbind", parLapplyLB(cl, splitData, krige.head.crossval,
                                        formula=formula,
                                        model=model,
                                        model.landtype=model.landtype,
                                        model.landtype.head=model.landtype.head,
                                        model.fixedHead=model.fixedHead,
                                        model.fixedHead.head=model.fixedHead.head,
                                        nmin=nmin,
                                        nmax=nmax,
                                        nmax.fixedHead=nmax.fixedHead,
                                        nmin.fixedHead=nmin.fixedHead,
                                        omax.fixedHead=omax.fixedHead,
                                        maxdist=maxdist,
                                        omax=omax,
                                        do.depth.est=do.depth.est,
                                        use.MrVBF=use.MrVBF,
                                        use.MrRTF=use.MrRTF,
                                        use.DEMsmoothing=use.DEMsmoothing,
                                        use.LandCatagory=use.LandCatagory,
                                        use.FixedHeads=use.FixedHeads,
                                        data=data,
                                        data.fixedHead=data.fixedHead,
                                        data.weights=data.weights,
                                        smooth.std=smooth.std,
                                        grid.elev=grid.elev,
                                        grid.MrVBF=grid.MrVBF,
                                        grid.MrRTF=grid.MrRTF,
                                        grid.LandType=grid.LandType,
                                        grid.params=grid.params,
                                        debug.level=debug.level))

    # Close cluster of the input was not a cluster object
    if (is.logical(use.cluster) || is.numeric(use.cluster))
      stopCluster(cl)

  } else {
    if (debug.level>0)
      message('... Starting kriging of ', nNewData,' points 1 core.')
    est <- krige.head.crossval(newdata=newdata,
                               formula=formula,
                               model=model,
                               model.landtype=model.landtype,
                               model.landtype.head=model.landtype.head,
                               model.fixedHead=model.fixedHead,
                               model.fixedHead.head=model.fixedHead.head,
                               nmin=nmin,
                               nmax=nmax,
                               nmax.fixedHead=nmax.fixedHead,
                               nmin.fixedHead=nmin.fixedHead,
                               omax.fixedHead=omax.fixedHead,
                               maxdist=maxdist,
                               omax=omax,
                               do.depth.est=do.depth.est,
                               use.MrVBF=use.MrVBF,
                               use.MrRTF=use.MrRTF,
                               use.DEMsmoothing=use.DEMsmoothing,
                               use.LandCatagory=use.LandCatagory,
                               use.FixedHeads=use.FixedHeads,
                               data=data,
                               data.fixedHead=data.fixedHead,
                               data.weights=data.weights,
                               smooth.std=smooth.std,
                               grid.elev=grid.elev,
                               grid.MrVBF=grid.MrVBF,
                               grid.MrRTF=grid.MrRTF,
                               grid.LandType=grid.LandType,
                               grid.params=grid.params,
                               debug.level=debug.level)
  }

  # Format outputs as data.frame and return
  est <- data.frame(Easting = est[,1], Northing= est[,2],head.obs= est[,3],head.pred= est[,4],head.var=est[,5],resid= est[,6],DEM= est[,7])
  coordinates(est) = ~Easting + Northing;
  if (debug.level>0)
    message('... Finished kriging. ', length(est),' points estimated.')
  return(est);
}
