#' Calibration of spatial interpolation parameters for mapping groundwater level observations.
#'
#' \code{krige.head.calib} calibrates parameters to minimise the interpolation error.
#'
#' @details This function optimises the parameters for the spatial interpolation so as to minimise the
#' interpolation error, which by default is defined using a formal likelihood function. Once this function has return
#' the optimal parameter set, the mapping can be undertaken using \code{\link{krige.head}} (see examples below).
#'
#' The mapping parameters are estimated using a mixed data-type (i.e. real and integer parameters) split-sample maximum likelihood global optimisation. The optimisation by default
#' includes the variogram parameters (e.g. range, sill and nuggest) and the search parameters for local kriging (e.g. radius, minimum and
#' maximum number of observations to use). Optimising these parameters is not common in kriging. It is done herein because trials for Victoria,
#' Australia, showed that calibrating these parameters produced significantly lower cross-validation errors (i.e. the error in predicting the observations
#' removed from the optimisation) compared to the standard approach using fitting to an experimental variogram. The optimisation is
#' numerically challenging and the following factors should be considered before use:
#'
#' \itemize{
#'  \item{Optimisation of the parameters \code{mrvbf.pslope}, \code{mrvbf.ppctl} and \code{smooth.std} often required the creating of raster grids for every parameter combination. To ease the computation burden, these parameters should be treated as disrcete, not continuous, numbers.}
#'  \item{The optimisation package \code{rgeoud} is used herein. See the \code{rgeoud} documentation for details of the optimisation scheme.}
#'  \item{Trials have established default calibration parameters and settings that were effective for Victoria, Australia. There is no guarantee they will be effective for other regions.}
#'  \item{The available input point data must be split up into \code{data} and \code{newdata}. The point observations within \code{data} are used to estimate the water level at the locations defined
#'  within \code{newdata}. The difference between the observed and estimated values defines the interpolation prediction error. This process also gives the kriging variance at
#'  each \code{newdata} location. The points used for \code{data} and \code{newdata} should both cover all types of terrain, geology, landuse and the full mapping extend so as to avoid bias. Also, changing either may
#'  change the calibration solution.}
#'  \item{The \code{objFunc.type} allows for maximum likelihood estimation - that is optimisation that accounts for the expected kriging error. See \code{objFunc.type} below for details.}
#' }
#'
#' In using this function, the primary user decisions are:
#' \itemize{
#'  \item{The kriging with external drift formula defining the independent gridded variables deemed to predict the groundwater elevation. See the input \code{formula}.}
#'  \item{The type of the objective function, defined by the input \code{objFunc.type}}
#'  \item{The type of variogram model, defined by the input \code{model}}
#'  }
#'
#' @param \code{formula} see \code{\link{krige.head}}.
#' @param \code{MrVBF} see \code{\link{krige.head}}.
#' @param \code{MrRTF} see \code{\link{krige.head}}.
#' @param \code{smoothing} see \code{\link{krige.head}}.
#' @param \code{grid} see \code{\link{krige.head}}.
#' @param \code{grid.landtype.colname} see \code{\link{krige.head}}.
#' @param \code{data} see \code{\link{krige.head}}.#'
#' @param \code{data.fixedHead} see \code{\link{krige.head}}.
#' @param \code{newdata} is as for \code{data} but the points in a split-sample cross-validation scheme to estimate the interpolation error. Points listed within \code{newdata} should not be
#' listed within \code{data}. \code{newdata} can be (i) a real scalar >0 and <1  defining the fraction of data points within \code{data} to randomly remove and from \code{data} and use for cross-validation,
#' (i) an integer scalar >1 and <\code{length(data)} defining the number of data points within \code{data} to randomly remove and from \code{data} and use for cross-validation  (iii) a vector of indexes defining rows numbers within \code{data} to be extract and uses for the cross-validation; (iv) a vector of logicals
#' with \code{TRUE} defining rows numbers within \code{data} to be extract and used for the cross-validation, (v) a \code{SpatialPointsDataFrame} variable defining the complete data for cross validation - which must have identical columns to \code{data}.
#' The default is \code{newdata=0.5},
#'
#' @param \code{data.errvar.colname} see \code{\link{krige.head}}.
#' @param \code{model} is either a character for the name of the variogram model type or a \code{gstat} variogram model object of type \code{variogramModel}. The available options are as per \code{gstat}, but
#' it is suggested to be \code{Mat}.
#' @param  \code{mrvbf.pslope} defines the calibration type and range for the MrVFB shape parameter for the slope (see Gallant et al. 2003). A vector of two values defining the optimisation range when the parameter is
#' treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
#' not be optimised. If the \code{formula} includes either of the terms \code{MrVBF} or \code{MrRTF}, then the default is \code{seq(0.5, 1.5, length.out = 11)}. Else, the default is \code{NULL}.
#' @param  \code{mrvbf.ppctl} defines the calibration type and range for the MrVFB shape parameter for elevation percentile (see Gallant et al. 2003). It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
#' treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
#' not be optimised.  If the \code{formula} includes either of the terms \code{MrVBF} or \code{MrRTF}, then the default is \code{seq(0.5, 1.5, length.out = 11)}. Else, the default is \code{NULL}.
#' @param  \code{smooth.std} defines the calibration type and range for the strength of the Gaussian kernal smoothing applied to the 5x5 grid cells surrounding each DEM grid cell. It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
#' treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
#' not be optimised.  The default is \code{seq(0.5, 1.5, length.out = 11)}.
#' @param \code{nmax} defines the calibration type and range for the maximum number of \code{data} observations to use when estimate each point using local kriging. It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
#' treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
#' not be optimised.  The default is \code{ceiling(seq(0.1,0.20,0.01)*length(data))}.
#' @param \code{nmax.fixedHead} defines the calibration type and range for the maximum number of \code{data.fixedHead} observations to use when estimate each point using local kriging. It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
#' treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
#' not be optimised.  The default is \code{seq(10,110,length=11)}.
#' @param \code{maxdist} defines the calibration type and range for the maximum search radius to use when estimate each point using local kriging. It can be a scalar number, a vector of two values defining the optimisation range when the parameter is
#' treated as a real number or a vector of length >2 values defining the optimisation increments when the parameter is treated as not continuous but discrete. If a single number is input, then the parameter will
#' not be optimised.  The default is from 10\% to 100\% of \code{grid} extend at increments of 10\%. If \code{grid} is \code{NULL}, then the user must input the search radius in one of the three accepted forms.
#' @param \code{nmin} see \code{\link{krige.head}}.
#' @param \code{omax} see \code{\link{krige.head}}.
#' @param \code{nsim} see \code{\link{krige.head}}.
#' @param \code{fit.variogram.type} defines the way the model variogram is to be derived/used.
#' For \code{fit.variogram.type==1} the input \code{model} then all of the variogram parameters will be calibrated. When \code{model} is a chacater for a variogram modle type, then
#' the variogram will be assumed isotropic. When \code{model} is a \code{gstat} variogram model object of type \code{variogramModel}, then all of the variogram parameters will be calibrated.
#' This approach allows significantly greater control of the variogram model form.
#' For \code{fit.variogram.type==2} the input \code{model} must be a \code{gstat} variogram model object of type \code{variogramModel}. The variogram model parameters will not be calibrated and instead
#' will be estimated by fitting the model variogram to an experimental variogram using multi-start local calibration. This approach is very similar to convential variogram model fitting.
#' For \code{fit.variogram.type==3} the input \code{model} must also be a variogram model object of type \code{variogramModel}. The variogram model parameters will not be calibrated or fit to an experimental variogram.
#' The default is \code{fit.variogram.type=1}.
#'
#' @param \code{objFunc.type} defines the type of objective function to use in the calibration. For \code{objFunc.type==1}, the negative likelihood divided by the number of \code{newdata} observations is used
#' - which accounts for the expected error in the \code{newdata} estimates. That is, \code{newdata} points far from other \code{data} points will have a large kriging variance and hence the expected
#' error in the kriging estimate is larger. Conversely, \code{newdata} points near to \code{data}  points will have a low kriging variance and hence the kriging water level error should be low.
#' For details see Samper and Neumen (1989A,B) and Pardo-Iguzquiza & Dowd (2013). For \code{objFunc.type==2}, it is as for type 1 but with an added penalty when the estimate violates a physical constraint.
#  For simplicity, a penalty is applied only when the groundwater level is above the land surface elevation; and hence is only appropriate when the unconfined water level is being estimated. The penalty is applied
#' by first estimating the difference between the land surface elevation and the head. At those \code{newdata} points above the land surface, the absolute of this difference is added to the absolute error between
#' the predicated value and that within \code{newdata}. For \code{objFunc.type==3}, the root-mean-square error is used. For \code{objFunc.type==4}, the root-mean-square error is used but with the penalty from type 2.
#' The default is \code{objFunc.type=1}.
#'
#' @param \code{use.cluster} see \code{\link{krige.head}}.
#'
#' @param \code{debug.level} see \code{\link{krige.head}}.
#'
#' @seealso \code{\link{krige.head}} for undertaking the mapping.
#'
#' @return
#' As per \code{genoud} plus the parameter values in the transformed calibration space "par.pretransformed" and "par" is replaced as a named list variable.
#'
#' @references
#' Gallant, J.C., Dowling, T.I. (2003) A multiresolution index of valley bottom flatness for mapping depositional areas, Water Resources Research, 39/12:1347-1359
#'
#' Rivoirard, J. & Romary, T. Math Geosci (2011) Continuity for Kriging with Moving Neighborhood, Mathematical Geosciences, 43: 469. DOI: 10.1007/s11004-011-9330-0
#'
#' Pardo-Igúzquiza E. , Dowd P. A., (2013) Comparison of inference methods for estimating semivariogram model parameters and their uncertainty: The case of small data sets,
#' Computers & Geosciences, v50,pp 154-164, https://doi.org/10.1016/j.cageo.2012.06.002.
#'
#' Samper, F. J., and S. P. Neuman (1989A), estimation of spatial covariance structures by adjoint state maximum likelihood cross validation: 1. Theory,
#' Water Resour. Res., 25(3), 351–362, doi: 10.1029/WR025i003p00351.
#'
#' Samper, F. J., and S. P. Neuman (1989B), Estimation of spatial covariance structures by adjoint state maximum likelihood cross validation: 2. Synthetic experiments,
#' Water Resour. Res., 25(3), 363–371, doi: 10.1029/WR025i003p00363.
#'
#' @examples
#'
#' # Set enironment path for hydroMap
#' set.env()
#'
#' # Load water table observations from  April 2000 for Victoria, Australia and a 250m state-wide DEM.
#' data('victoria.groundwater')
#'
#' # Load a model variogram and mapping parameters found to be effective.
#' data('mapping.parameters')
#'
#' # Define a simple kriging formula without MrVBF terms that does not require the package RSAGA.
#' f <- as.formula('head ~ elev + smoothing')
#'
#' # Define an initial isotropic variogram model. All of the parameters will be calibrated.
#' varigram.model <- gstat::vgm(psill=25, model="Mat", range=5000, nugget=5, kappa = 0.5)
#'
#' # Calibrate the mapping parameters with 25% of the data randomly selected at 2 cores.
#' calib.results <- krige.head.calib(formula=f, grid=DEM, data=obs.data, newdata=0.25, nmin=0, nmax=Inf, maxdist=Inf, omax=0,
#'                       data.errvar.colname='total_err_var', model = variogram.model,  fit.variogram.type=1,
#'                       pop.size.multiplier=1, debug.level=0, use.cluster = 2)
#'
#' # Grid the observed head using ONLY the training data from the calibration (and all cores) and then map..
#' head.grid.training <- krige.head(calibration.results = calib.results, use.cluster = T)
#' spplot(head.grid.training)
#' 
#' # Grid the observed head using ALL of the obsrved data from the calibration (and all cores).
#' head.grid.all <- krige.head(calibration.results = calib.results, data=obs.data, use.cluster = T)
#' spplot(head.grid.all)
#'
#' @export
krige.head.calib <-
  function(formula = as.formula("head ~ elev + MrVBF + MrRTF + smoothing"),
           grid,
           grid.landtype.colname = NULL,
           data,
           newdata = 0.5,
           data.fixedHead=NULL,
           data.errvar.colname = NULL,
           model = c('Mat'),
           mrvbf.pslope = if(any(match(all.vars(as.formula(formula)), 'MrVBF',nomatch=F) | match(all.vars(as.formula(formula)), 'MrRTF',nomatch=F))){seq(0.20, 2.5, by = 0.1)}else{NULL},
           mrvbf.ppctl  = if(any(match(all.vars(as.formula(formula)), 'MrVBF',nomatch=F) | match(all.vars(as.formula(formula)), 'MrRTF',nomatch=F))){seq(0.20, 2.5, by = 0.1)}else{NULL},
           smooth.std = seq(0.5, 1.5, length.out = 11),
           nmax = if(is.character(data)){-999}else{ceiling(seq(0.1,0.20,0.01)*length(data))},
           nmax.fixedHead = if(!is.null(data.fixedHead)) {seq(10,110,length=11)}else{NULL},
           maxdist = if(class(grid)=='SpatialPixelsDataFrame' || class(grid)=='SpatialGridDataFrame'){ceiling(0.5*sqrt((extent(grid)[2]-extent(grid)[1])^2 + (extent(grid)[4]-extent(grid)[3])^2)*seq(0.1,1,0.1))}else{-999},
           nmin = 0.2,
           omax = NULL,
           fit.variogram.type = 1,
           objFunc.type = 1,
           use.cluster = TRUE,
           pop.size.multiplier = 3,
           max.generations = 200,
           hard.generation.limit=FALSE,
           solution.tolerance=1e-4,
           debug.level=0) {

    # Check gstst is loaded
    if (!require('gstat',quietly = T))
      stop('Please install the gstat package. It is required by HydroMap.')
    
    if (debug.level>0)
      message('Calibrating mapping parameters.')

    # Store input parameter values (exluding variogram params). These are updated if they're calibrated.
    solution = list()
    solution$parameter.Values$mrvbf.pslope = mrvbf.pslope
    solution$parameter.Values$mrvbf.ppctl = mrvbf.ppctl
    solution$parameter.Values$smooth.std = smooth.std
    solution$parameter.Values$nmax = nmax
    solution$parameter.Values$nmax.fixedHead = nmax.fixedHead
    solution$parameter.Values$maxdist = maxdist
    solution$parameter.Values$nmin = nmin
    solution$parameter.Values$omax = omax
    
    # Assess formula componants
    var.names = all.vars(as.formula(formula));
    use.MrVBF = any(match(var.names, 'MrVBF'));
    use.MrRTF = any(match(var.names, 'MrRTF'));
    use.DEMsmoothing = any(match(var.names, 'smoothing'));
    if (is.na(use.MrVBF))
      use.MrVBF = FALSE;
    if (is.na(use.MrRTF))
      use.MrRTF = FALSE;
    if (is.na(use.DEMsmoothing))
      use.DEMsmoothing = FALSE;

    # Check RSAGA is loaded
    if (use.MrVBF || use.MrRTF)
      if (!require('RSAGA',quietly = T))
        stop('Please install the RSAGA package. It is required by HydroMap.')
    
    # Reset the package environment.
    clear.env()
    
    # Add a NA buffer if observations are on the boundary AND smoothing or MrVBF are used.
    # This step was added because smoothing and MrVBF can have NAs at boundaries,
    # resulting in points and grid cells at the boundary not being able to be estimated. 
    if (use.MrVBF || use.MrRTF || use.DEMsmoothing) {
      x.colbuffer = 0;
      y.rowbuffer = 0;
      kernal.maxdim = 7
      grid.asRaster = raster::raster(grid);
      if (any(!is.na(grid.asRaster[,1])) || any(!is.na(grid.asRaster[,ncol(grid.asRaster)]))) {
        x.colbuffer = kernal.maxdim
      }  
      if (any(!is.na(grid.asRaster[1,])) || any(!is.na(grid.asRaster[,nrow(grid.asRaster)]))) {
        y.rowbuffer = kernal.maxdim
      }     
      if (x.colbuffer>0 || y.rowbuffer>0) {
        warning('Buffer added around the input-grid of 1-gridcell. This is required to allow estimation of points at the end of the DEM.',immediate.=T);
        grid.asRaster = raster::extend(grid.asRaster, c(y.rowbuffer, x.colbuffer),NA)
        
        # Infill NA DEM values of grid by taking the local average. This was essential to ensure
        # DEM values at fixed head points beyond the mappng area (eg coastal points
        # with a fixed head of zero just beyond the DEM extent)
        dem.asRaster = raster::raster(grid,layer='DEM');
        for (i in 1:kernal.maxdim)
          grid.asRaster = raster::focal(grid.asRaster, w=matrix(1,kernal.maxdim,kernal.maxdim), fun=mean, na.rm=TRUE, NAonly=TRUE)
        
        grid.input = grid
        grid = as(grid.asRaster, class(grid))
        names(grid)='DEM'
      }
      rm('grid.asRaster')
    }
    
    # Import data if string names are input.
    data  = import.pointData(data);
    grid  = import.DEM(grid);
    if (is.character(newdata))
      newdata  = import.pointData(newdata);

    # Get calibration values for nmax and maxdist. This should only be undertaken if the user input strings for 'data' or 'grid'
    if (length(nmax)==1 && nmax==-999)
      nmax = ceiling(seq(0.04,0.20,0.02)*length(data));
    if (length(maxdist)==1 && maxdist==-999)
      maxdist = ceiling(0.5*sqrt((extent(grid)[2]-extent(grid)[1])^2 + (extent(grid)[4]-extent(grid)[3])^2)*seq(0.1,1,0.1));


    # Set up cluster for parallel calculations.
    use.local.cluster = FALSE
    if (is.logical(use.cluster)) {
      if (use.cluster)
        nclus = parallel::detectCores(all.tests = TRUE, logical = FALSE)
      else {
        nclus = 1

      }
    } else if (is.numeric(use.cluster)) {
      if (use.cluster > 0) {
        nclus = min(parallel::detectCores(all.tests = TRUE, logical = FALSE),floor(use.cluster));
      } else {
        stop('When use.cluster is a numerical value is must be >0')
      }
    }
    if (class(use.cluster)[1] == "SOCKcluster" |
        class(use.cluster)[1] == "PVMcluster" |
        class(use.cluster)[1] == "spawnedMPIcluster" |
        class(use.cluster)[1] == "MPIcluster") {
      use.local.cluster = use.cluster

    } else if (nclus>1) {
      if (debug.level>0)
        message('... Setting up local cores for parallel calibration.')
      
      clus <- c(rep("localhost", nclus))
      use.local.cluster <- parallel::makeCluster(clus, type = "SOCK")
    } else {
      use.local.cluster = F;
    }

    # Add libraries to nodes.
    if (nclus>1) {
      parallel::clusterEvalQ(use.local.cluster, library(gstat))
      parallel::clusterEvalQ(use.local.cluster, library(raster))
    }

    # Get parameter bounds and assess if discrete
    param.Names = vector(mode = "character", length = 0)
    param.bounds = matrix(Inf, 0, 2)
    param.IntegerLookupTable = list()
    doIntegerCalib = vector(mode = "logical", length = 0)

    message('... Assessing parameters to calibrate and, if discrete, define lookup table.')
    nParams = 0

    if (!is.null(mrvbf.pslope)) {
      if (length(mrvbf.pslope)==1) {
        # Add fixed parameter value to the lookup table This is requied so that only a subset of parameters can be calibrated.
        param.IntegerLookupTable[['mrvbf.pslope']] = cbind(1, mrvbf.pslope)
      } else if (length(mrvbf.pslope) == 2) {
        nParams = nParams + 1

        param.Names[nParams] = 'mrvbf.pslope'
        param.bounds = rbind(param.bounds, c(mrvbf.pslope[1], mrvbf.pslope[2]))
        doIntegerCalib[nParams] = FALSE
      } else {
        nParams = nParams + 1

        param.Names[nParams] = 'mrvbf.pslope'
        param.bounds = rbind(param.bounds, c(1, length(mrvbf.pslope)))
        doIntegerCalib[nParams] = TRUE


        param.IntegerLookupTable[[param.Names[nParams]]] = cbind(seq(1, length(mrvbf.pslope), 1), mrvbf.pslope)
      }
    }
    if (!is.null(mrvbf.ppctl)) {
      if (length(mrvbf.ppctl)==1) {
        # Add fixed parameter value to the lookup table This is requied so that only a subset of parameters can be calibrated.
        param.IntegerLookupTable[['mrvbf.ppctl']] = cbind(1, mrvbf.ppctl)
      } else if (length(mrvbf.ppctl) == 2) {
        nParams = nParams + 1

        param.Names[nParams] = 'mrvbf.ppctl'
        param.bounds = rbind(param.bounds, c(mrvbf.ppctl[1], mrvbf.ppctl[2]))
        doIntegerCalib[nParams] = FALSE
      } else {
        nParams = nParams + 1

        param.Names[nParams] = 'mrvbf.ppctl'
        param.bounds = rbind(param.bounds, c(1, length(mrvbf.ppctl)))
        doIntegerCalib[nParams] = TRUE

        param.IntegerLookupTable[[param.Names[nParams]]] = cbind(seq(1, length(mrvbf.ppctl), 1), mrvbf.ppctl)
      }
    }
    if (!is.null(smooth.std)) {
      if (length(smooth.std)==1) {
        # Add fixed parameter value to the lookup table This is requied so that only a subset of parameters can be calibrated.
        param.IntegerLookupTable[['smooth.std']] = cbind(1, smooth.std)
      } else if (length(smooth.std) == 2) {
        nParams = nParams + 1

        param.Names[nParams] = 'smooth.std'
        param.bounds = rbind(param.bounds, c(smooth.std[1], smooth.std[2]))
        doIntegerCalib[nParams] = FALSE
      } else {
        nParams = nParams + 1

        param.Names[nParams] = 'smooth.std'
        param.bounds = rbind(param.bounds, c(1, length(smooth.std)))
        doIntegerCalib[nParams] = TRUE

        param.IntegerLookupTable[[param.Names[nParams]]] = cbind(seq(1, length(smooth.std), 1), smooth.std)
      }
    }
    if (!is.null(nmax)) {
      if (length(nmax)==1) {
        # Add fixed parameter value to the lookup table This is requied so that only a subset of parameters can be calibrated.
        param.IntegerLookupTable[['nmax']] = cbind(1, nmax)
      } else if (length(nmax) == 2) {
        nParams = nParams + 1

        param.Names[nParams] = 'nmax'
        param.bounds = rbind(param.bounds, c(nmax[1], nmax[2]))
        doIntegerCalib[nParams] = FALSE
      } else {
        nParams = nParams + 1

        param.Names[nParams] = 'nmax'
        param.bounds = rbind(param.bounds, c(1, length(nmax)))
        doIntegerCalib[nParams] = TRUE

        param.IntegerLookupTable[[param.Names[nParams]]] = cbind(seq(1, length(nmax), 1), nmax)
      }
    }
    if (!is.null(nmax.fixedHead)) {
      if (length(nmax.fixedHead)==1) {
        # Add fixed parameter value to the lookup table This is requied so that only a subset of parameters can be calibrated.
        param.IntegerLookupTable[['nmax.fixedHead']] = cbind(1, nmax.fixedHead)
      } else if (length(nmax.fixedHead) == 2) {
        nParams = nParams + 1

        param.Names[nParams] = 'nmax.fixedHead'
        param.bounds = rbind(param.bounds, c(nmax.fixedHead[1], nmax.fixedHead[2]))
        doIntegerCalib[nParams] = FALSE
      } else {
        nParams = nParams + 1

        param.Names[nParams] = 'nmax.fixedHead'
        param.bounds = rbind(param.bounds, c(1, length(nmax.fixedHead)))
        doIntegerCalib[nParams] = TRUE

        param.IntegerLookupTable[[param.Names[nParams]]] = cbind(seq(1, length(nmax.fixedHead), 1), nmax.fixedHead)
      }
    }

    # Split the data create a prediction data set
    #-------------------
    if (is.null(newdata)) {
      stop('newdata cannot be NULL. It is required to define a prediction set of the calibration.')

    } else if (length(newdata)==1) {
      if (as.integer(newdata)==newdata && newdata>=1) {
        newdata = sample(1:length(data), size=as.integer(newdata), replace=F)
      } else if (newdata>0 && newdata<1) {
        newdata = sample(1:length(data), size=floor(newdata*length(data)), replace=F)
      } else {
        stop('newdata input is invalid. When a scalar, it must be >0 and <length(data) or >0 and <1')
      }

      data.IDs.logical = logical(length(data))
      data.IDs.logical[newdata] = TRUE
      newdata = data[data.IDs.logical,]
      data = data[!data.IDs.logical,]

    } else if (is.numeric(newdata)) {
      # Split up the input 'data' using the row indexes within 'newdata'
      # - which define the rows within 'data' to remove from 'data ' and input to 'newdata'

      if (min(newdata)<1 || max(newdata)>length(data))
        stop('newdata input vector of indexes contain values <1 and/or > length(data)')

      if (unique(newdata)!=length(newdata))
        stop('newdata input vector of indexes contain duplicate indexes.')

      data.IDs.logical = logical(length(data))
      data.IDs.logical[newdata] = TRUE
      newdata = data[data.IDs.logical,]
      data = data[!data.IDs.logical,]

    } else if (is.logical(newdata)) {
      # Split up the input 'data' using the a vector of logicals within 'newdata'
      # - which define the rows within 'data' to remove from 'data ' and input to 'newdata'

      if (length(newdata)!=length(data))
        stop('When newdata is a logical vector its length must be equal to length(data)')

      data.IDs.logical = newdata
      newdata = data[data.IDs.logical,]
      data = data[!data.IDs.logical,]

    } else if (class(newdata)!='SpatialPointsDataFrame') {
      stop('newdata input type is invalid. It must be in integer scalar >1, a real scalar >0 and <1, a vector of data row IDs or T/F or or a SpatialPointsDataFrame')
    } else {

      # Check that newdata and data have the same colum names.
      if (any(!(names(data) == names(newdata))))
        stop('The inputs data and newdata must have the same column names when both are input as SpatialPointsDataFrames.')
      
      # Reset package enviro data      
      pkg.env$DEM.data <<- NULL; pkg.env$MrVBF.data <<- NULL; pkg.env$MrRTF.data <<- NULL; pkg.env$smoothDEM.data <<- NULL
    
    }
    message(paste('... Number of obs points used for kriging:',length(data)))
    message(paste('... Number of obs points used for obj. function:',length(newdata)))
    #-------------------


    if (!is.null(maxdist)) {

      # Determine the minimum search distance to ensure the search criteria is met.
      # This is undertaken to minimise the probability of discontinuities.
      if (length(maxdist) >1) {
        message('... Calculating the mapping extent.');
        if (is.null(grid)) {
          maxMapRadius = sqrt((extent(data)[2]-extent(data)[1])^2 + (extent(data)[4]-extent(data)[3])^2);
        } else {
          maxMapRadius = sqrt((extent(grid)[2]-extent(grid)[1])^2 + (extent(grid)[4]-extent(grid)[3])^2);
        }
        maxMapRadius = maxMapRadius/2;

        # Set search criteria nmin and omax if to be set as a function of nmax.
        message('... Calculating the nmin and omax (if required).');
        nmax_max = max(nmax, na.rm=T)
        if (nmax_max<Inf && !is.null(nmax_max)) {

          # If 0 <= omax < 1, then omax is treated as a fraction of nmax, else omax is just passed to gstat.
          omax.tmp = omax;
          if (!is.null(omax) && omax>=0 && omax<1) {
            omax.tmp = max(0,floor(omax * nmax_max))
          }

          # If 0 <= nmin < 1, then nmin is treated as a fraction of nmax, else nmax is just passed to gstat.
          nmin.tmp = nmin;
          if (!is.null(nmin) && nmin>=0 && nmin<1) {
            nmin.tmp = max(0,floor(nmin * nmax_max))
          }

        } else {
          omax.tmp = omax;
          nmin.tmp = nmin;
        }

        # use gstat to find the maximum distance to the minimum number of data points.
        message('... Checking the maxdist to nmin observations from grid cells to assess if all grid cells can be estimated.');
        g = gstat(NULL, 'mindist.est', formula=as.formula('head~1'), data=rbind(data,newdata), nmin=nmin.tmp, nmax=nmin.tmp, maxdist = maxMapRadius, omax=omax.tmp, set=list(method="distance"))
        if (!is.logical(use.local.cluster)) {
          grid.tmp = grid;
          sp::gridded(grid.tmp) = FALSE;
          nData = length(grid.tmp);
          splt = rep(1:nclus, each = ceiling(nData/nclus), length.out = nData)
          grid.tmp = lapply(as.list(1:nclus), function(w) grid.tmp[splt == w,])
          clusterExport(use.local.cluster, list("g"), envir = environment());
          grid.dist <- do.call("rbind", parallel::parLapply(use.local.cluster, grid.tmp, function(lst) predict(g, lst)))
        } else {
          grid.dist <- predict(g, grid)
        }
        maxdist.min = ceiling(max(grid.dist$mindist.est.var,na.rm=T));
        message(paste('    Minimum search distance = ',maxdist.min));
        remove(list=c('g','grid.tmp','grid.dist','splt'))

        # Check input search distance values.
        if (maxdist.min >= max(maxdist))
          stop('The input search distance range "maxdist" is less than the minimum search distance required for estimating all grid cells.')
      }

      # Assign updated maxdist
      if (length(maxdist)==1) {
        # Add fixed parameter value to the lookup table This is requied so that only a subset of parameters can be calibrated.
        param.IntegerLookupTable[['maxdist']] = cbind(1, maxdist)
      } else if (length(maxdist) == 2) {
        nParams = nParams + 1
        param.Names[nParams] = 'maxdist'

        # Update max dist values
        if (maxdist.min > maxdist[1] && maxdist.min < maxdist[2]) {
          message('    The input lower search distance range "maxdist" has been increased to the minimum search distance required for estimating all grid cells.');
          maxdist[1] = maxdist.min;
        }

        # Add parameter bounds to calib. data.
        param.bounds = rbind(param.bounds, c(maxdist[1], maxdist[2]))
        doIntegerCalib[nParams] = FALSE
      } else {
        nParams = nParams + 1
        param.Names[nParams] = 'maxdist'

        if (any(maxdist<maxdist.min)) {
          message('    Scaling the input search distance ranges to be greater than or equal to  the minimum search distance...')
          maxdist = seq(maxdist.min, max(maxdist), length.out=length(maxdist))
        }

        # Add parameter bounds to calib. data.
        param.bounds = rbind(param.bounds, c(1, length(maxdist)))
        doIntegerCalib[nParams] = TRUE
        param.IntegerLookupTable[[param.Names[nParams]]] = cbind(seq(1, length(maxdist), 1), maxdist)
      }
    }

    # If all of the input data is not provided (ie if the function was not called from krige.heads()) then extract point data
    #-------------------

    # Find   mid-point of MRVBF & smoothing parameters. This only does anything if
    # the parameters are vectors. It is required so that some initial data is available
    # for the calibration variogram estimation.
    if (use.MrVBF || use.MrRTF) {
      if (length(mrvbf.pslope)==2)
        mrvbf.pslope.tmp = min(mrvbf.pslope) + (max(mrvbf.pslope)-min(mrvbf.pslope))/2
      else {
        ind_pslope = max(1,round(length(mrvbf.pslope)/2))
        mrvbf.pslope.tmp = mrvbf.pslope[ind_pslope];
      }
      if (length(mrvbf.ppctl)==2)
        mrvbf.ppctl.tmp = min(mrvbf.ppctl) + range(mrvbf.ppctl)/2
      else {
        ind_ppctl = max(1,round(length(mrvbf.ppctl)/2))
        mrvbf.ppctl.tmp = mrvbf.ppctl[ind_ppctl];
      }
    } else {
      mrvbf.pslope.tmp = NULL;
      mrvbf.ppctl.tmp = NULL;
    }
    if (use.DEMsmoothing) {
      if (length(smooth.std)==2)
        smooth.std.tmp = min(smooth.std) + (max(smooth.std)-min(smooth.std))/2
      else {
        ind_smooth = max(1,round(length(smooth.std)/2))
        smooth.std.tmp = smooth.std[ind_smooth];
      }
    } else {
      smooth.std.tmp = NULL;
    }

    # Get point data (if required)
    if ((use.MrVBF && !('MrVBF' %in% names(data))) ||
        (use.MrRTF && !('MrRTF' %in% names(data))) ||
        (use.DEMsmoothing && !('smoothing' %in% names(data))) ||
        (('DEM' %in% names(data)))) {

      if (debug.level>0)
        message(' ... Getting point data');

      data = get.allData(formula = formula, data, grid, mrvbf.pslope = mrvbf.pslope.tmp, mrvbf.ppctl =mrvbf.ppctl.tmp, smooth.std = smooth.std.tmp)

    }
    #-------------------

    # Setup the variogram parameters if they are to be calibrated (ie fit.variogram.type=1).
    if (fit.variogram.type == 1) {
      if (all(doIntegerCalib)) {
        useParamTransform = TRUE

      } else
        useParamTransform = FALSE


      doIntegerCalib = FALSE

      use.gradient.check = FALSE

      if (debug.level>0)
        message('... Assessing variogram parameters to calibrate routine.')

      # If the input variogram model is empty, then get an initial estimate.
      if (class(model)[1] != "variogramModel") {
        message('... No input variogram. Deriving OLS initial estimate.')
        model = get.variogram(formula, data, model = model, max.Its = 10)
      }

      # Extract the model parameters.
      nModels = dim(model)[1]

      variogram.type = as.character(model$model)

      variogram.range = model$range

      variogram.psill = model$psill

      variogram.kappa = model$kappa

      variogram.ang1 = model$ang1

      variogram.anis1 = model$anis1


      #message('DBG: nModels=',nModels)

      # Add variogram parameters to the list of parameters.
      for (i in 1:nModels) {
        if (any(model$model[i] == 'Nug')) {
          #if (nchar(grid.landtype.colname) == 0) {
          nParams = nParams + 1

          param.Names[nParams] = 'nug'
          param.bounds = rbind(param.bounds,
                               c(0.1 * variogram.psill[i], 10 * variogram.psill[i]));
          #}
        } else {
          #if (nchar(grid.landtype.colname) == 0) {
          nParams = nParams + 1

          param.Names[nParams] = paste('psill_model_', i - 1, sep = '')

          param.bounds = rbind(param.bounds,
                               c(0.1 * variogram.psill[i], 10 * variogram.psill[i]))
          #}

          nParams = nParams + 1

          param.Names[nParams] = paste('range_model_', i - 1, sep = '')

          param.bounds = rbind(param.bounds,
                               c(0.1 * variogram.range[i], 10 * variogram.range[i]))

          if (variogram.ang1[i] != 0) {
            nParams = nParams + 1

            param.Names[nParams] = paste('ang1_model_', i - 1, sep = '')

            param.bounds = rbind(param.bounds, c(0, 180))
          }

          if (variogram.anis1[i] != 1) {
            nParams = nParams + 1

            param.Names[nParams] = paste('anis1_model_', i - 1, sep = '')

            param.bounds = rbind(param.bounds, c(0, 1))
          }

          if (model$model[i] == 'Mat' || model$model[i] == 'Ste') {
            nParams = nParams + 1

            param.Names[nParams] = paste('kappa_model_', i - 1, sep = '')

            param.bounds = rbind(param.bounds,
                                 c(0.1 * variogram.kappa[i], 10 * variogram.kappa[i]))
          }
        }
      }

      # Remove 'nug' from list of model types and add to 'model' variable
      filt = variogram.type != 'Nug'

      model = variogram.type[filt]


    } else {
      # Check that if all parameters have discrete values then all are discrete.
      if (all(doIntegerCalib)) {
        doIntegerCalib = TRUE

        useParamTransform = FALSE

        use.gradient.check = FALSE

      } else if (any(doIntegerCalib)) {
        doIntegerCalib = FALSE

        useParamTransform = TRUE

        use.gradient.check = FALSE

      } else {
        doIntegerCalib = FALSE

        useParamTransform = FALSE

        use.gradient.check = TRUE

      }
    }

    # Check and set population size
    if (is.null(pop.size.multiplier)) {
      pop.size = 3 * nParams
    } else if (!is.numeric(pop.size.multiplier) || pop.size.multiplier < 1) {
      stop('The input calibration population size multiplier (see pop.size.multiplier) must be numerical and >=1')
    } else {
      pop.size = ceiling(pop.size.multiplier) * nParams
    }

    # Assess if smoothing is to be undertaken.
    var.names = all.vars(formula)
    use.DEMsmoothing = any(match(var.names, 'smoothing'))

    # Set as minimisation problem
    doObjMaximisation = FALSE

    # Turn off gradient estimation
    use.BFGS = FALSE

    # Do calibration
    if (doIntegerCalib) {
      message('... Starting integer value calibration of the followng parameters:',paste(param.Names,sep="",collapse=', '))
    } else if (useParamTransform) {
      message('... Starting mixed data-type calibration of the following parameters:',paste( param.Names,sep="",collapse=', '))
    } else {
      message('... Starting continuous value calibration of the following parameters:',paste( param.Names,sep="",collapse=', '))
    }
    genoud.solution = rgenoud::genoud(
      HydroMap:::get.objFunc,
      nvars = nParams,
      max = doObjMaximisation,
      pop.size = pop.size,
      max.generations = max.generations,
      wait.generations = 4,
      hard.generation.limit = hard.generation.limit,
      MemoryMatrix = TRUE,
      Domains = param.bounds,
      solution.tolerance = solution.tolerance,
      boundary.enforcement = 2,
      data.type.int = doIntegerCalib,
      unif.seed = unclass(Sys.time()),
      int.seed = floor(sqrt(unclass(Sys.time()))),
      transform = useParamTransform,
      BFGS = use.BFGS,
      gradient.check = use.gradient.check,
      cluster = F,
      param.names = param.Names,
      formula = formula,
      grid = grid,
      grid.landtype.colname = grid.landtype.colname,
      data = data,
      data.fixedHead = data.fixedHead,
      newdata = newdata,
      data.errvar.colname = data.errvar.colname,
      model = model,
      nmin = nmin,
      omax = omax,
      fit.variogram.type = fit.variogram.type,
      objFunc.type = objFunc.type,
      lookup.table = param.IntegerLookupTable,
      return.transformed.params = useParamTransform,
      use.cluster = use.local.cluster,
      debug.level=debug.level
    )


    message('... Finished calibration.')

    if (!is.logical(use.local.cluster)) {
      if (debug.level>0)
        message('... Closing cluster.')
      parallel::stopCluster(use.local.cluster)
    }

    # Store the returned value for use when calling get.objFunc directly.
    # This is required if some of the parameters have discrete parameters (which are used by the lookup table to transform to meaningful values)
    genoud.solution$par.pretransformed = genoud.solution$par;

    # Back transform parameters
    lookup.table.names = names(param.IntegerLookupTable)

    for (i in 1:length(param.Names)) {
      # Assess if the parameter name is listed in the lookup table. If so, get the non-transformed value.
      isIntegerParam =  param.Names[i] %in% lookup.table.names

      # Get non-transformed value from table.
      if (isIntegerParam) {
        genoud.solution$param.Values[i] = param.IntegerLookupTable[[param.Names[i]]][genoud.solution$par[i], 2]
      } else {
        genoud.solution$param.Values[i] = genoud.solution$par[i];
      }
    }

    # Print summary
    message(paste('    Best objective function value:', genoud.solution$value))
    for (i in 1:nParams)
      message(paste('    Parameter ', param.Names[i], ' = ', genoud.solution$param.Values[i]))
    
    # Add formula
    solution$inputs$formula = formula;     
    
    # Add data and new data
    solution$inputs$data = data
    solution$inputs$newdata = newdata
    if (exists('grid.input')) {
      solution$inputs$grid = grid.input 
    } else {
      solution$inputs$grid = grid
    }
    
    # Add remaining inputs 
    solution$inputs$grid.landtype.colname = grid.landtype.colname
    solution$inputs$data.fixedHead = data.fixedHead
    solution$inputs$data.errvar.colname = data.errvar.colname  

    ind = match(param.Names,'mrvbf.pslope')
    ind = na.omit(ind)
    if (length(ind)==1 && is.finite(ind))
      solution$parameter.Values$mrvbf.pslope = genoud.solution$param.Values[ind]
    
    ind = match(param.Names,'mrvbf.ppctl')
    ind = na.omit(ind)
    if (length(ind)==1 && is.finite(ind))
      solution$parameter.Values$mrvbf.ppctl = genoud.solution$param.Values[ind]
    
    ind = match(param.Names,'smooth.std')
    ind = na.omit(ind)
    if (length(ind)==1 && is.finite(ind))
      solution$parameter.Values$smooth.std = genoud.solution$param.Values[ind]
    
    ind = match(param.Names,'nmax')
    ind = na.omit(ind)
    if (length(ind)==1 && is.finite(ind))
      solution$parameter.Values$nmax = genoud.solution$param.Values[ind]
    
    ind = match(param.Names,'nmax.fixedHead')
    ind = na.omit(ind)
    if (length(ind)==1 && is.finite(ind))
      solution$parameter.Values$nmax.fixedHead = genoud.solution$param.Values[ind]
    
    ind = match(param.Names,'maxdist')
    ind = na.omit(ind)
    if (length(ind)==1 && is.finite(ind))
      solution$parameter.Values$maxdist = genoud.solution$param.Values[ind]
    
    ind = match(param.Names,'nmin')
    ind = na.omit(ind)
    if (length(ind)==1 && is.finite(ind))
      solution$parameter.Values$nmin = genoud.solution$param.Values[ind]
    
    ind = match(param.Names,'omax')
    ind = na.omit(ind)
    if (length(ind)==1 && is.finite(ind))
      solution$parameter.Values$omax = genoud.solution$param.Values[ind]
    
    # Add variogram model
    if (fit.variogram.type == 1)
      solution$variogramModel = get.variogramModel(genoud.solution$par, param.Names, model)
    
    
    # Add the lookup table to the solution. This is required if the the objectiveget.objFunc() is to be directly called.
    solution$lookup.table = param.IntegerLookupTable;
    
    # Get point data for derived covariates
    solution$inputs$newdata <- get.allData(formula =solution$inputs$formula, 
                              data = solution$inputs$newdata, 
                              grid = grid, 
                              mrvbf.pslope = solution$parameter.Values$mrvbf.pslope, 
                              mrvbf.ppctl = solution$parameter.Values$mrvbf.ppctl,
                              smooth.std = solution$parameter.Values$smooth.std,
                              debug.level=debug.level )
      
    # Get prediction errors
    xval.est = krige.head.crossval(solution)
    xval.est = data.frame(xval.est)
    sp::coordinates(xval.est) = ~Easting + Northing
    solution$predictions = xval.est
    
    # Add calibration settings
    solution$genoud.output = genoud.solution
    
    # Reset the package environment.
    clear.env()
    
    return(solution)

  }
