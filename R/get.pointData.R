get.pointData <- function(formula = as.formula("head ~ elev + MrVBF + MrRTF"), data, grid, mrvbf.pslope = 1.0, mrvbf.ppctl = 1.0, smooth.std = 1.0, smooth.ncells=11, smoothingKernal = NULL, debug.level=0 ) {

    # Check enviro variables are setup
    if (!exists('pkg.env') || is.null(pkg.env))
      stop('The environment variable are not setup. Call set.env() with paths to SAGA.');	

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
    if (is.na(use.DEMsmoothing))
      use.DEMsmoothing = FALSE;
      
    # Check is MrVBF can be used.
    if ((use.MrVBF || use.MrRTF) && !pkg.env$saga.has.MrVFB)
      stop('MrVBF & MrRTF cannot be used because RSAGA is not setup correctly. Call set.env() with paths to SAGA.');
    
    # Check is MrVBF can be used.
    if ((use.MrVBF || use.MrRTF) && !pkg.env$saga.has.MrVFB)
      stop('MrVBF & MrRTF cannot be used because RSAGA is not setup correctly. Call set.env() with paths to SAGA.');	  
      
    # Get DEM
    data = get.elev(data, grid)
    
    # Get MRVBF grid and point data.
    if (use.MrVBF || use.MrRTF) {
      if (debug.level>0) 
        message('... Calling get.MrVBF() to get grid values:')
      
      data = get.MrVBF(data, grid, mrvbf.pslope, mrvbf.ppctl, use.MrVBF, use.MrRTF)	    
      
      if (debug.level>0) 
        message('... Exiting get.MrVBF().')
    }
	    
    # Get smoothed DEM.
    if (use.DEMsmoothing) {
    
      # Calc. smoothed DEM
      data = get.smoothedDEM(data, grid, smooth.std = smooth.std, smooth.ncells=smooth.ncells, smoothingKernal=smoothingKernal)
    
      # Calculate difference between DEM and smoothed DEM.
      data$smoothing = data$smoothDEM - data$DEM;      
      # 
      # deltaDEMsmoothed = as(as.raster(grid$smoothing) - as.raster(grid$elev), 'SpatialPixelsDataFrame');
      # sp::gridded(deltaDEMsmoothed) = TRUE;
      # sp::fullgrid(deltaDEMsmoothed) = TRUE;
      # grid$smoothing = deltaDEMsmoothed[['layer']]	        
    }
}