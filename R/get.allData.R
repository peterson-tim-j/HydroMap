get.allData <- function(formula = as.formula("head ~ elev + MrVBF + MrRTF"), data, grid, grid_buffer=NULL, mrvbf.pslope = 1.0, mrvbf.ppctl = 1.0, smooth.std = 1.0, smoothingKernal = NULL,  smooth.ncells=11, min.sepDist = 5, maxStoredGrids=5, use.cluster=TRUE, debug.level=0) {

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

	use.extraTerms = FALSE;
	ind = match(var.names, c('head','depth','elev','smoothing','MrVBF','MrRTF'))
	if (any(is.na(ind))) {
	  use.extraTerms = T;
	  use.extraTerms.names = var.names[is.na(ind)]
	}
	
	# Get input data
	if (!is.null(data))	  
	  data = import.pointData(data, min.sepDist = min.sepDist, orderby.head = do.head.est, debug.level=debug.level)
	  
	# Get grids
	grid = import.DEM(grid, debug.level=debug.level);
	  
	# Get DEM point values
	if (!is.null(data))
	  newdata= get.elev(data, grid, debug.level=debug.level)
	else {
	  newdata = grid;
	}

	# If there are extra terms to the formula, then check that the grid contains the required extra data.
	if (use.extraTerms && !is.null(data)) {
	  if (debug.level>0) {
	    message('... Examining existence of point data for the extra formula terms.');
	    message(paste('    Grid column names are:',paste(names(grid),collapse=", ")));
	    message(paste('    Point data column names are:',paste(names(newdata),collapse=", ")));
	    message(paste('    User defined formula terms:',paste(use.extraTerms.names,collapse=", ")));    
	  }
	  for (i in 1:length(use.extraTerms.names)) {
	    if (all(is.na(match(names(newdata), use.extraTerms.names[i])))) {
	      
	      # data column nor in data. Interpolate grid.
	      grid.asRaster = raster::raster(grid,layer=use.extraTerms.names[i]);
	      #newdata[[use.extraTerms.names[i]]] = NULL;
	      if (debug.level>0) 
	        message(paste('... Interpolating grid to point locations for formula term:',use.extraTerms.names[i]));    
	      
	      tmp = raster::extract(grid.asRaster, newdata, method='bilinear');    

	      if (debug.level>0) 
	        message(paste('... Add interpolated point data locations for formula term',use.extraTerms.names[i], 'to data set.'));    
	      
	      newdata$tmpName =  tmp;
	      ncols = length(names(newdata))
	      names(newdata)[ncols] = use.extraTerms.names[i];

	    }
	  }
	}
	
	
	# Get MRVBF grid and point data.
	if (use.MrVBF || use.MrRTF) {
	  if (is.null(data))
	    newdata = get.MrVBF(NULL, newdata, mrvbf.pslope, mrvbf.ppctl, use.MrVBF, use.MrRTF, maxStoredGrids, debug.level=debug.level)
	  else {
	    newdata = get.MrVBF(newdata, grid, mrvbf.pslope, mrvbf.ppctl, use.MrVBF, use.MrRTF, maxStoredGrids, debug.level=debug.level)	  
	  }
	}

	# Get smoothed DEM.
	if (use.DEMsmoothing) {
	
	  # Calc. smoothed DEM
	  if (is.null(data)) {
	    newdata = get.smoothedDEM(NULL, newdata, grid_buffer, smooth.std = smooth.std, smoothingKernal=smoothingKernal, smooth.ncells=smooth.ncells, maxStoredGrids, use.cluster=use.cluster ,debug.level=debug.level)
	  } else {
	    newdata = get.smoothedDEM(newdata, grid, grid_buffer, smooth.std = smooth.std, smoothingKernal=smoothingKernal, smooth.ncells=smooth.ncells, maxStoredGrids, use.cluster=use.cluster, debug.level=debug.level)	  
	  }
	
	  # Calculate difference between DEM and smoothed DEM.
	  if (is.null(data)) {
	    if (debug.level>0) 
	      message('... Calculate difference between DEM and smoothed DEM grids.')
	    deltaDEMsmoothed = as(raster::raster(newdata,layer='smoothDEM') - raster::raster(newdata,layer='DEM'), 'SpatialPixelsDataFrame');
	    sp::gridded(deltaDEMsmoothed) = TRUE;
	    sp::fullgrid(deltaDEMsmoothed) = TRUE;
	    newdata$smoothing = deltaDEMsmoothed[['layer']]	        
	  } else {
	    if (debug.level>0) 
	      message('... Calculate difference between DEM and smoothed DEM points.')
	    newdata$smoothing = newdata$smoothDEM - newdata$DEM;
	  }
	}
	
	# Return
	if (debug.level>0) 
	  message('... Finished extracting data from grids.')
	return(newdata)
}
