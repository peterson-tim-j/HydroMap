#' Reads in an ARCMAP ASCII grid file containing the DEM and prepares it for use in the head mapping.
#'
#' \code{import.DEM} reading in a an ARCMAP ASCII grid file.
#'
#' @details This function reads in an ARCMAP ASCII grid file containing the DEM
#' and prepares it for use in the mapping of groundwater level data. The preparations
#' include removing NA grid cells and converting it to a fully gridded Spatial object.
#' Also, the coordinates must be projected.
#'
#' @param \code{grid} character string for the ARCMAP ASCII grid file name.
#'
#' @param \code{debug.level} Control the user messages. A value >0 outputs progress.
#'
#' @return
#' A gridded Spatial object will be returned with one field named "DEM".
#'
#' @export
import.DEM <- function(grid, debug.level=0 ) {

	if (is.character(grid)) {
	  if (debug.level>0)
	    message('... Reading in grid and writing DEM.asc');

	  # Build path to DEM asc file
	  if (dirname(grid) == '.') {
	    fname.dem <- paste(pkg.env$working.path,grid,sep=.Platform$file.sep);
	  } else {
	    fname.dem <- grid;
	  }

	  # Read in grid
	  grid = read.asciigrid(fname.dem,colname='DEM');

	  # Remove NULLS
	  grid = data.frame(grid)
	  names(grid) = c('DEM',"Easting","Northing");
	  grid[['NA']] <- NULL
	  coordinates(grid) = ~Easting + Northing;

	  gridded(grid)=TRUE;
	  fullgrid(grid)=TRUE;
	} else if (class(grid)!='SpatialPixelsDataFrame' && class(grid)!='SpatialGridDataFrame') {
	  stop('The input grid must be a file name or a SpatialPixelsDataFrame or SpatialGridDataFrame.')
	}

	return(grid)

}
