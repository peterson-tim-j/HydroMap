get.MrVBF <- function(data = NULL, grid, pslope=NULL, ppctl=NULL, return.MrVBF=TRUE, return.MrRTF=TRUE, maxStoredGrids=10, debug.level=0) {

  if (debug.level>0)
    message('Getting MrVBF data:');

  # Check enviro variables are setup
  if (!exists('pkg.env') || is.null(pkg.env))
    stop('    The environment variable are not setup. Call set.env() with paths to SAGA.');

  # Check if 'data', 'DEM' or 'newdata' are strings. If so, import the data.
  if (is.character(data)) {
    if (debug.level>0)
      message('... Reading in point file');

    fname.data <- paste(pkg.env$working.path,data,sep=.Platform$file.sep);
    data <- read.table(fname.data, sep=',',header=T);
    sp::coordinates(data) = ~Easting + Northing;
  }
  if (is.character(grid)) {
    if (debug.level>0)
      message('... Reading in grid and writing DEM.asc');

    # Build path to DEM asc file
    fname.dem <- paste(pkg.env$working.path,grid,sep=.Platform$file.sep);

    # Read in grid
    grid = sp::read.asciigrid(fname.dem,colname='elev');
    sp::coordinates(grid) = ~Easting + Northing;
    sp::gridded(grid)=TRUE;
    sp::fullgrid(grid)=TRUE;
  }

  # Build name of MrVBF grid
  var.name.MrVBF <- paste('mrvbf_slope',pslope,'_pctl',ppctl,sep='');
  var.name.MrRTF <- paste('mrrtf_slope',pslope,'_pctl',ppctl,sep='');

  # Build path to DEM asc file
  fname.MrVBF<- paste(pkg.env$working.path,var.name.MrVBF,sep=.Platform$file.sep);
  fname.MrRTF<- paste(pkg.env$working.path,var.name.MrRTF,sep=.Platform$file.sep);

  # Check if the required MrVBF point data has already been extracted.
  if (return.MrVBF) {
    do.MrVBF <- TRUE;
  } else {
    do.MrVBF <- FALSE;
  }

  if (return.MrRTF) {
    do.MrRTF <- TRUE;
  } else {
    do.MrRTF <- FALSE;
  }
  if (!is.null(data) && ((!is.null(pkg.env$MrVBF.data) && return.MrVBF) || (!is.null(pkg.env$MrRTF.data) && return.MrRTF))) {
    if (return.MrVBF) {
	    var.name.all = names(pkg.env$MrVBF.data)
	    if (any(!is.na(match(var.name.all, var.name.MrVBF))) && length(data) == length(pkg.env$MrVBF.data))
	    do.MrVBF <- FALSE;
    }
    if (return.MrRTF) {
    	var.name.all = names(pkg.env$MrRTF.data)
    	if (any(!is.na(match(var.name.all, var.name.MrRTF))) && length(data) == length(pkg.env$MrRTF.data))
  	  do.MrRTF <- FALSE;
    }
  } else if (is.null(data)) {
	  var.name.all = names(pkg.env$MrVBF.grid)
  	if (any(!is.na(match(var.name.all, var.name.MrVBF))) && return.MrVBF && length(grid) == length(pkg.env$MrVBF.grid))
	    do.MrVBF <- FALSE;

	  var.name.all = names(pkg.env$MrRTF.grid)
  	if (any(!is.na(match(var.name.all, var.name.MrRTF))) && return.MrRTF && length(grid) == length(pkg.env$MrRTF.grid))
	    do.MrRTF <- FALSE;
  }


  # Calculate MrVBF and MrRTF if required and if the file does not already exist
  grid.MrVBF <- NULL
  grid.MrRTF <- NULL
  if (do.MrVBF || do.MrRTF) {

    # Check if fname.dem exists before re-running SAGA.
    if (!file.exists(paste(fname.MrVBF,'.asc',sep='')) || !file.exists(paste(fname.MrRTF,'.asc',sep=''))) {

      if (!pkg.env$saga.has.MrVFB)
        stop('MrVBF & MrRTF cannot be calculated because RSAGA is not setup correctly. Install the RSAGA package and call set.env() with the paths to SAGA. Alternativelly, copy the required MrVBF and MrRTF grids to the working directory.');


      if (pkg.env$saga.has.MrVFB) {
        # Define SAGA DEM file names
        fname.dem <- paste(pkg.env$working.path,'DEM.asc',sep=.Platform$file.sep);
        fname.dem.sgrd <- paste(pkg.env$working.path,'DEM.sgrd',sep=.Platform$file.sep);
        fname.dem.sdat <- paste(pkg.env$working.path,'DEM.sdat',sep=.Platform$file.sep);

        # Write DEM grid to .asc if not already done and then convert ASCII DEM to SAGA format.
        if (!file.exists(fname.dem)) {
          if (debug.level>0)
            message('... Writing DEM.asc');

          write.asciigrid(grid,fname.dem, attr=1);
        }

        # Write DEM grid to .asc if not already done and then convert ASCII DEM to SAGA format.
        if (!file.exists(fname.dem.sgrd) || !file.exists(fname.dem.sdat)) {
          if (debug.level>0)
            message('... Writing DEM in SAGA format.');

          RSAGA::rsaga.import.gdal(fname.dem,env = pkg.env$saga.settings);
        }
      }

      if (debug.level>0)
        message('... Calling SAGA to calculate MrVBF.');

    	# Calc. MrVBF
      RSAGA::rsaga.geoprocessor(lib = 'ta_morphometry', module='Multiresolution Index of Valley Bottom Flatness (MRVBF)',
    	              param = list(DEM = 'DEM.sgrd', MRVBF=paste(fname.MrVBF,'.sgrd',sep=''),MRRTF=paste(fname.MrRTF,'.sgrd',sep=''),
    	                           P_SLOPE=pslope,P_PCTL=ppctl),env = pkg.env$saga.settings, show.output.on.console=F, warn=F, flags='s');

    	# Convert SAGA MrVBF format to ASCII DEM.
      if (debug.level>0)
    	  message('... Converting SAGA grids to .ASC.');

    	if (do.MrVBF) {
    	  RSAGA::rsaga.sgrd.to.esri(paste(fname.MrVBF,'.sgrd',sep=''),env = pkg.env$saga.settings,  show.output.on.console=F, warn=F, flags='s');
    	  file.remove(paste(fname.MrVBF,'.sgrd',sep=''));
    	  file.remove(paste(fname.MrVBF,'.mgrd',sep=''));
    	  file.remove(paste(fname.MrVBF,'.sdat',sep=''));
    	}
    	if (do.MrRTF) {
    	  RSAGA::rsaga.sgrd.to.esri(paste(fname.MrRTF,'.sgrd',sep=''),env = pkg.env$saga.settings,  show.output.on.console=F, warn=F, flags='s');
    	  file.remove(paste(fname.MrRTF,'.sgrd',sep=''));
    	  file.remove(paste(fname.MrRTF,'.mgrd',sep=''));
    	  file.remove(paste(fname.MrRTF,'.sdat',sep=''));
    	}
    }

    # Check if the MrVBF / MrRTF data already exists in the enviro.
    # If not then read in Mr VBF grid and get grid attributes
    store.grid.MrVBF <- FALSE;
    store.grid.MrRTF <- FALSE;
    if (do.MrVBF) {
  	  var.name.all = names(pkg.env$MrVBF.grid)
  	  ind <- which(names(pkg.env$MrVBF.grid)==var.name.MrVBF);
	    if (!is.na(ind) && sum(ind)>0) {
    	  # Extract grid data
    	  grid.MrVBF = pkg.env$MrVBF.grid[ind];
  	  } else if (!is.null(grid.MrVBF)) {
	      store.grid.MrRTF <- TRUE;
	    } else {
	      if (debug.level>0)
    	    message('... Reading in MrVBF .asc file.');

    	  # read in file with grid data
    	  grid.MrVBF = sp::read.asciigrid(paste(fname.MrVBF,'.asc',sep=''),colname=var.name.MrVBF);
    	  sp::gridded(grid.MrVBF) = TRUE;
    	  sp::fullgrid(grid.MrVBF) = TRUE;

    	  store.grid.MrVBF <- TRUE;
	    }

  	  # Assess if the maximum number of stored grids is exceeded. If so, then delete the first grid.
  	  if (store.grid.MrVBF) {
  	    if (maxStoredGrids<=0)
  	      store.grid.MrVBF <- FALSE;

  	    if (length(var.name.all)>= maxStoredGrids) {
  	      if (debug.level>0)
  	        message('... Removing least recently created MrVBF grid from enviroment memory.');

  	      filt = seq(length(var.name.all),1,-1) < maxStoredGrids
  	      pkg.env$MrVBF.grid = pkg.env$MrVBF.grid[var.name.all[filt]];
  	    }
  	  }

    	# Get grid geometry
    	grid.MrVBF.params = sp::gridparameters(grid.MrVBF)
    }
    if (do.MrRTF) {
    	var.name.all = names(pkg.env$MrRTF.grid)
    	ind <- which(names(pkg.env$MrRTF.grid)==var.name.MrRTF);
    	if (!is.na(ind) && sum(ind)>0) {
    	  # Extract grid data
    	  grid.MrRTF = pkg.env$MrRTF.grid[ind];
    	} else if (!is.null(grid.MrRTF)) {
    	  store.grid.MrRTF <- TRUE;
    	} else {
    	  if (debug.level>0)
    	    message('... Reading in MrRTF .asc file.');

    	  # read in file with grid data
    	  grid.MrRTF = sp::read.asciigrid(paste(fname.MrRTF,'.asc',sep=''),colname=var.name.MrRTF);
    	  sp::gridded(grid.MrRTF) = TRUE;
    	  sp::fullgrid(grid.MrRTF) = TRUE;

    	  store.grid.MrRTF <- TRUE;
  	  }

    	# Assess if the maximum number of stored grids is exceeded. If so, then delete the first grid.
    	if (store.grid.MrRTF) {
    	  if (maxStoredGrids<=0)
    	    store.grid.MrRTF <- FALSE;

    	  if (length(var.name.all)>= maxStoredGrids) {
    	    if (debug.level>0)
    	      message('... Removing least recently created MrRTF grid from enviroment memory.');

    	    filt = seq(length(var.name.all),1,-1) < maxStoredGrids
    	    pkg.env$MrRTF.grid = pkg.env$MrRTF.grid[var.name.all[filt]];
    	  }
    	}

    	# Get grid geometry
    	grid.MrRTF.params = sp::gridparameters(grid.MrRTF)
    }

    # Get DEM grid attributes
    grid.DEM.params = sp::gridparameters(grid)

  # Check if MrVBF grid is of equal res and extent to DEM
  if (do.MrVBF) {
  	if (grid.DEM.params$cellcentre.offset[1] != grid.MrVBF.params$cellcentre.offset[1] ||
  	grid.DEM.params$cellcentre.offset[2] != grid.MrVBF.params$cellcentre.offset[2] ||
  	grid.DEM.params$cellsize[1] != grid.MrVBF.params$cellsize[1] ||
  	grid.DEM.params$cellsize[2] != grid.MrVBF.params$cellsize[2] ||
  	grid.DEM.params$cells.dim[1] != grid.MrVBF.params$cells.dim[1] ||
  	grid.DEM.params$cells.dim[2] != grid.MrVBF.params$cells.dim[2])
  	  stop('    MrVBF pre-calculated grid(s) are of a different resolution and or extent to the input DEM. Reset using set.env().')
  }

      # Check if MrVBF grid is of equal res and extent to DEM
  if (do.MrRTF) {
  	if (grid.DEM.params$cellcentre.offset[1] != grid.MrRTF.params$cellcentre.offset[1] ||
  	grid.DEM.params$cellcentre.offset[2] != grid.MrRTF.params$cellcentre.offset[2] ||
  	grid.DEM.params$cellsize[1] != grid.MrRTF.params$cellsize[1] ||
  	grid.DEM.params$cellsize[2] != grid.MrRTF.params$cellsize[2] ||
  	grid.DEM.params$cells.dim[1] != grid.MrRTF.params$cells.dim[1] ||
  	grid.DEM.params$cells.dim[2] != grid.MrRTF.params$cells.dim[2])
  	  stop('    MrRTF pre-calculated grid(s) are of a different resolution and or extent to the input DEM. Reset using set.env().')
  }

  # Store the grids in the eniro space
  if (store.grid.MrVBF && do.MrVBF) {
    if (debug.level>0)
  	  message('... Storing MrVBF grid into envir. variable.');

  	if (is.null(pkg.env$MrVBF.grid)) {
  	  pkg.env$MrVBF.grid = grid.MrVBF;
  	} else {
  	  pkg.env$MrVBF.grid[[var.name.MrVBF]] = grid.MrVBF[[var.name.MrVBF]];
  	}
  }
  if (store.grid.MrRTF && do.MrRTF) {
  	if (is.null(pkg.env$MrRTF.grid)) {
  	  pkg.env$MrRTF.grid = grid.MrRTF;
  	} else {
  	  pkg.env$MrRTF.grid[[var.name.MrRTF]] = grid.MrRTF[[var.name.MrRTF]];
  	}
  }

  # Interpolate MrVBF to the required points.
  if (!is.null(data)) {
  	if (do.MrVBF) {

      # Interpole grid. If any NA values then take the local average. This was essential to ensure
      # DEM values at fixed head points beyond the mappng area (eg coastal points
      # with a fixed head of zero just beyond the DEM extent)
      MrVBF.asRaster = raster::raster(grid.MrVBF,layer=var.name.MrVBF);
      data$MrVBF = NULL;
      if (debug.level>0)
        message('... Interpolating grid to point locations.');

      tmp = raster::extract(MrVBF.asRaster, data, method='bilinear');
      data$tmpName =  tmp;
      ncols = length(names(data))
      names(data)[ncols] = 'MrVBF'

  	  # Append interpolated points to enviro variable
  	  if (is.null(pkg.env$MrVBF.data) || length(data) != length(pkg.env$MrVBF.data)) {
  	    Easting = sp::coordinates(data)[1];
  	    Northing = sp::coordinates(data)[2];
  	    tmp = data.frame(Easting, Northing, tmp);
  	    names(tmp) = c('Easting','Northing',var.name.MrVBF)
  	    sp::coordinates(tmp) = ~Easting + Northing
  	    pkg.env$MrVBF.data = tmp;
  	  } else {
  	    pkg.env$MrVBF.data$tmp = tmp;
  	  }
  	  ncols = length(names(pkg.env$MrVBF.data))
  	  names(pkg.env$MrVBF.data)[ncols] = var.name.MrVBF
  	  rm(tmp);
  	}

  	# Interpolate MrRTF to the required points.
  	if (do.MrRTF) {
      # Interpole grid. If any NA values then take the local average. This was essential to ensure
      # DEM values at fixed head points beyond the mappng area (eg coastal points
      # with a fixed head of zero just beyond the DEM extent)
      MrRTF.asRaster = raster::raster(grid.MrRTF,layer=var.name.MrRTF);
      data$MrRTF = NULL;
      if (debug.level>0)
        message('... Interpolating grid to point locations.');

      tmp = raster::extract(MrRTF.asRaster, data, method='bilinear');
      data$tmpName =  tmp;
      ncols = length(names(data))
      names(data)[ncols] = 'MrRTF'

  	  # Append interpolated points to enviro variable
  	  if (is.null(pkg.env$MrRTF.data) || length(data) != length(pkg.env$MrRTF.data)) {
  	    Easting = sp::coordinates(data)[1];
  	    Northing = sp::coordinates(data)[2];
  	    tmp = data.frame(Easting, Northing, tmp);
  	    names(tmp) = c('Easting','Northing',var.name.MrRTF)
  	    sp::coordinates(tmp) = ~Easting + Northing
  	    pkg.env$MrRTF.data = tmp;
  	  } else {
  	    pkg.env$MrRTF.data$tmp = tmp;
  	  }
  	  ncols = length(names(pkg.env$MrRTF.data))
  	  names(pkg.env$MrRTF.data)[ncols] = var.name.MrRTF
  	  rm(tmp);
    }

	  return(data);

    } else {
      if (debug.level>0)
        message('... Retrieving previously estimated MrVBF / MrRTF grid vals.');

    	if (return.MrVBF)
    	  grid$MrVBF <- grid.MrVBF[[var.name.MrVBF]];

    	if (return.MrRTF)
    	  grid$MrRTF <- grid.MrRTF[[var.name.MrRTF]];

    	return(grid)
    }

  } else {

    if (!is.null(data)) {
      if (debug.level>0)
        message('... Retrieving previously estimated MrVBF / MrRTF point vals.');

      if (return.MrVBF)
        data$MrVBF <- pkg.env$MrVBF.data[[var.name.MrVBF]];

      if (return.MrRTF)
        data$MrRTF <- pkg.env$MrRTF.data[[var.name.MrRTF]];

      return(data);
    } else {
      if (debug.level>0)
        message('... Retrieving previously estimated MrVBF / MrRTF grid.');

      if (return.MrVBF)
      	grid$MrVBF <- pkg.env$MrVBF.grid[[var.name.MrVBF]];

      if (return.MrRTF)
      	grid$MrRTF <- pkg.env$MrRTF.grid[[var.name.MrRTF]];

      return(grid);

      }
  }
}
